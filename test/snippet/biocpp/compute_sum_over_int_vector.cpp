#include <iostream>
#include <numeric>
#include <vector>

int32_t sum(std::vector<int32_t> const & data)
{
    int32_t intermediate_sum{};
    for (int32_t element : data) // Calls in the background data.begin() and data.end().
        intermediate_sum += element;

    return intermediate_sum;
}

int32_t sum_short(std::vector<int32_t> const & data)
{
    return std::accumulate(data.begin(), data.end(), 0);  // Use algorithms
}

int main()
{
    std::vector<int32_t> int_data{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    std::cout << "The sum of my vector is " << sum(int_data) << "\n";
    std::cout << "The sum of my vector is " << sum_short(int_data) << "\n";

    return EXIT_SUCCESS;
}
