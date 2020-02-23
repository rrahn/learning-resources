#include <atomic>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/algorithm/count_if.hpp>

template <typename value_t, typename init_t>
concept bool accumulable_with = requires (value_t && value, init_t && init)
{
    { std::move(init) + std::forward<value_t>(value)} -> init_t;
};

template <std::ranges::input_range range_t>
    requires accumulable_with<std::ranges::range_value_t<range_t>, int32_t>
int32_t sum(range_t && data)
{
    auto common_range = data | std::views::common;
    return std::accumulate(std::ranges::begin(common_range), std::ranges::end(common_range), int32_t{});
}

seqan3::phred42 read_user_base_quality()
{
    std::cout << "Please specify the minimum base quality as a probability between 0 and 41\n";
    int32_t user_base_quality{};
    std::cin >> user_base_quality;

    if (user_base_quality < 0 || user_base_quality > 41)
        throw std::invalid_argument{"Only values in the interval [0, 41] can be used."};

    return seqan3::phred42{}.assign_rank(user_base_quality);
}

using character_string = const char *;

int main(int const , character_string argv[])
{
    std::string_view fastq_input_path1{argv[1]};
    std::string_view fastq_input_path2{argv[2]};

    seqan3::phred42 minimum_phred_quality = read_user_base_quality();
    std::cout << "The corresponding phred score is: " << static_cast<int32_t>(minimum_phred_quality.to_rank()) << "\n";
    std::cout << "The corresponding ASCII mapping is: " << minimum_phred_quality.to_char() << "\n";

    seqan3::sequence_file_input seq_file_in1{fastq_input_path1};
    seqan3::sequence_file_input seq_file_in2{fastq_input_path2};
    // seqan3::sequence_file_output seq_file_out{std::cout, seqan3::format_fasta{}};

    // Find average quality.
    auto by_average_quality = [minimum_phred_quality] (auto && fastq_record) // ["ACGTGAATGAC...", "id", "==+!456AA..."]
    {
        auto && quality_sequence = seqan3::get<seqan3::field::qual>(fastq_record); // ["==+!456AA..."]
        auto quality_rank_sequence = quality_sequence | seqan3::views::to_rank; // ["56,56,34, 33, ..."]
        int32_t rank_sum = sum(quality_rank_sequence);
        int32_t average_rank = rank_sum / std::ranges::distance(quality_sequence);

        assert(average_rank < seqan3::alphabet_size<seqan3::phred42>);
        return seqan3::phred42{}.assign_rank(average_rank) >= minimum_phred_quality;
    };

    std::atomic<size_t> total_sequence_count{};
    // Only print sequences with average quality filter.
    auto filtered_count = ranges::count_if(seqan3::views::zip(seq_file_in1, seq_file_in2), [=, total_sequence = std::ref(total_sequence_count)] (auto && fastq_record_pair) {
        ++total_sequence.get();
        using std::get;
        return by_average_quality(get<0>(fastq_record_pair)) && by_average_quality(get<1>(fastq_record_pair));
    });

    std::cout << "Filtered " << std::setprecision(3) << static_cast<double>(total_sequence_count - filtered_count) / static_cast<double>(total_sequence_count) * 100 << "%\n";
    std::cout << "Total sequences read " << total_sequence_count << "\n";

    return EXIT_SUCCESS;
}
