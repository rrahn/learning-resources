#include <iostream>
#include <numeric>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

template <std::ranges::input_range range_t>
    requires seqan3::arithmetic<std::ranges::range_value_t<range_t>>
// Concept: forwarding reference.
// Concept: auto return type
// Concept: type_traits
// Concept: range algorithms
// Concept: trailing return type
auto average(range_t && rng) -> std::ranges::range_value_t<range_t>
{
    using arithmetic_t = std::ranges::range_value_t<range_t>;
    using index_sum_t = std::pair<arithmetic_t, arithmetic_t>;

    auto indexed_range = seqan3::views::zip(std::views::iota(1), rng) | std::views::common;
    index_sum_t result = std::accumulate(std::ranges::begin(indexed_range),
                                         std::ranges::end(indexed_range),
                                         index_sum_t{},
                                         [] (auto && initial, auto && current_value) -> index_sum_t
                                         {
                                             return {current_value.first, initial.second + current_value.second};
                                         });
    return result.second / result.first;
}

// Concept: range concept
// Concept: requires clause
// Concept: arithmetic
template <std::ranges::forward_range range_t>
    requires seqan3::arithmetic<std::ranges::range_value_t<range_t>>
// Concept: forwarding reference.
// Concept: auto return type
// Concept: type_traits
// Concept: range algorithms
// Concept: trailing return type
auto average(range_t && rng) -> std::ranges::range_value_t<range_t>
{

    using arithmetic_t = std::ranges::range_value_t<range_t>;
    auto common_range = rng | std::views::common;
    return std::accumulate(std::ranges::begin(common_range),
                           std::ranges::end(common_range),
                           arithmetic_t{}) /
           std::ranges::distance(rng);
}

// Concept: Program entry
int main(int const argc, char const * argv[])
{
    // Concept: exception
    if (argc < 3)
        throw std::invalid_argument{"Not enough arguments specified."};

    // Concept: strong type
    int32_t min_average{};
    // Concept: string_view
    std::string_view user_average{argv[2]};
    std::cout << "user value: " << user_average << "\n";

    // Concept: Char conversion from C++17
    // Concept: error codes from C++11
    if (auto [p, ec] = std::from_chars(user_average.begin(), user_average.end(), min_average); ec != std::errc{})
        throw std::make_error_code(ec);

    // Concept: Sequence file.
    seqan3::sequence_file_input input{argv[1]};

    // Concept: lambda-function -> invocable types C++11
    auto by_average_quality = [min_average] (auto const & fastq_record)
    {
        return average(seqan3::get<seqan3::field::qual>(fastq_record) | seqan3::views::to_rank | seqan3::views::single_pass_input) > min_average;
    };

    // Concept: range based for-loop
    // Concept: filter view
    for (auto const & [seq, id, qual] : input | std::views::filter(by_average_quality))
    {
        // Concept debug stream
        seqan3::debug_stream << "\n@" << id << "\n";
        seqan3::debug_stream << seq << "\n";
        seqan3::debug_stream << "+" << "\n";
        seqan3::debug_stream << qual;
    }
}
