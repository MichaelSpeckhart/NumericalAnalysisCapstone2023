#pragma once

#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include <cstddef>

//#include "../src/functions.cc"

enum RETURN {
    SUCCESS = 0,
    FAILURE = 1
};

struct result_t {
    bool succeeded;
    std::string msg;
    std::vector<uint8_t> data;
};

struct Received {
    size_t func_id;
    std::string args;
    std::vector<std::byte> data;
};

/**
 * @brief Serialize the 
 * 
 * @tparam T 
 * @param matrix 
 * @return std::vector<byte> 
 */
// template <typename T>
// std::vector<std::byte> serialize_matrix(const std::vector<std::vector<T>> matrix) {
//     std::vector<std::byte> byte_data;

//     for (const auto& row : matrix) {
//         for (const auto& element : row) {
//             byte* byte_ptr = reinterpret_cast<byte*>(&element);
//             for (size_t i = 0; i < sizeof(T); ++i) {
//                 byte_data.push_back(byte_ptr[i]);
//             }
//         }   
//     }

//     return byte_data;
// }

// /**
//  * @brief Deserialize the matrix from a vector of bytes 
//  * 
//  * @tparam T 
//  * @param byte_data 
//  * @param rows 
//  * @param cols 
//  * @return std::vector<std::vector<T>> 
//  */
template <typename T>
std::vector<std::vector<T>> deserialize_matrix(const std::vector<std::byte> byte_data, size_t rows, size_t cols) {
    std::vector<std::vector<T>> matrix(rows, std::vector<T>(cols));
    size_t byte_index = 0;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            T element = 0;
            for (size_t k = 0; k < sizeof(T); ++k) {
                element |= static_cast<T>(byte_data[byte_index++]) << (k * 8);
            }
            matrix[i][j] = element;
        }
    }

    return matrix;
}