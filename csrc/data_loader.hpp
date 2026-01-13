// Note to self: Generated using Gemini

#ifndef DATA_LOADER_HPP
#define DATA_LOADER_HPP

#include <vector>
#include <string>

// Use type aliases for clarity.
using Point = std::vector<float>;
using Dataset = std::vector<Point>;

/**
 * @brief Reads a dataset of variable-dimension points from a specified file.
 * * Each line in the file is expected to contain space-separated floating-point numbers.
 *
 * @param filename The path to the input text file.
 * @param dataset The vector to store the loaded points. Cleared before loading.
 * @return true if the file was successfully opened and processing finished (even if lines were skipped), false otherwise.
 */
bool readPointsFromFile(const std::string& filename, Dataset& dataset);

/**
 * @brief Prints the contents of the loaded dataset to standard output.
 *
 * @param points The dataset to print.
 */
void print_points(const Dataset& points);

#endif // DATA_LOADER_HPP
