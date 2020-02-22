#include <cmath>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

template <typename range_t>
int32_t sum(range_t const & data)
{
    return std::accumulate(std::ranges::begin(data), std::ranges::end(data), 0);  // Use algorithms
}

int32_t convert_to_sanger_quality(double probability)
{
    return std::lround(-10 * std::log10(1.0 - probability));
}

using character_string = const char *;

int main(int const , character_string argv[])
{
    std::string_view fastq_input_path{argv[1]};

    std::cout << "Please specify the minimum base quality as a value between 0.00 and 1.00\n";
    double minimum_base_quality{};
    std::cin >> minimum_base_quality;

    seqan3::phred42 minimum_phred_quality{};
    minimum_phred_quality.assign_rank(convert_to_sanger_quality(minimum_base_quality));

    std::cout << "The corresponding phred score is: " << static_cast<int32_t>(minimum_phred_quality.to_rank()) << "\n";
    std::cout << "The corresponding ASCII mapping is: " << minimum_phred_quality.to_char() << "\n";

    seqan3::sequence_file_input seq_file_in{fastq_input_path};
    seqan3::sequence_file_output seq_file_out{std::cout, seqan3::format_fasta{}};

    // Find average quality.
    auto by_average_quality = [minimum_phred_quality] (auto const & fastq_record) // ["ACGTGAATGAC...", "id", "==+!456AA..."]
    {
        auto quality_sequence = seqan3::get<seqan3::field::qual>(fastq_record); // ["==+!456AA..."]
        auto quality_rank_sequence = quality_sequence | seqan3::views::to_rank; // ["56,56,34, 33, ..."]
        auto rank_sum = sum(quality_rank_sequence);
        auto average_rank = rank_sum / std::ranges::distance(quality_sequence);
        return seqan3::phred42{}.assign_rank(average_rank) >= minimum_phred_quality;
    };

    // Only print sequences with average quality filter.
    seq_file_out = seq_file_in | std::views::filter(by_average_quality); // untie elements of a tuple

    return EXIT_SUCCESS;
}
