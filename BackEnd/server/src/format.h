#pragma once

#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>



enum RETURN {
    SUCCESS = 0,
    FAILURE = 1
};

struct result_t {
    bool succeeded;
    std::string msg;
    std::vector<double> vecData;
    std::vector<std::vector<double>> matData;
};

struct Received {
    size_t func_id;
    std::string args;
    std::vector<std::byte> data;
};

/**
* @brief Storing all the function data in a struct to easily pass around rather than individual variables
* 
*/
struct FunctionData {
   size_t mFuncId;
   std::vector<std::vector<double>> mFirstMatrix;
   std::vector<std::vector<double>> mSecondMatrix;
   size_t mScalar;
   std::vector<double> mVector;
};

/**
 * @brief Serialize the 
 * 
 * @tparam T 
 * @param matrix 
 * @return std::vector<byte> 
 */

std::vector<std::vector<double>> deserialize_matrix(const std::vector<double> matrixVals, size_t rows, size_t cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix[i][j] = matrixVals[i * cols + j];
        }
    }

    return matrix;
}