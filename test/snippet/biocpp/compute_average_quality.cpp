#include <cmath>
#include <iostream>
#include <numeric>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/algorithm>

int32_t sum(std::vector<int32_t> const & data)
{
    return std::accumulate(data.begin(), data.end(), 0);  // Use algorithms
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
    // Only print files
    for (auto && [seq, id, qual] : seq_file_in) // untie elements of a tuple
    {
        if (std::ranges::any_of(qual, [&](auto current_quality){ return current_quality < minimum_phred_quality; }))
            seq_file_out.emplace_back(seq, id);
    }

    return EXIT_SUCCESS;
}
