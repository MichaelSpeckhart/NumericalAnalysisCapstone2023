#pragma once

#ifndef PARSE_REQUESTS_H
#define PARSE_REQUESTS_H


#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <string>



#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

    struct result_t {
        bool succeeded;
        std::string msg;
        std::vector<double> vecData;
        std::vector<std::vector<double>> matData;
        std::string clientMatrix;
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

namespace Capstone {

    

    result_t parse_request(std::string receivedData, std::size_t bytes);

    std::string extract_json(std::string receivedData);

    std::vector<std::vector<double>> parse_matrix(std::string matrixData);

    std::vector<std::vector<double>> deserialize_matrix(const std::vector<double> matrixVals, size_t rows, size_t cols);
    
    std::string serialize_matrix(const std::vector<std::vector<double>> matrix);

    template <typename T> T mapIdToFunction(FunctionData& data);
}

#endif
