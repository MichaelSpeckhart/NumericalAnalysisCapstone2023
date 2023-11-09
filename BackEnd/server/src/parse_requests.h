#pragma once

#ifndef PARSE_REQUESTS_H
#define PARSE_REQUESTS_H


#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <string>
#include <type_traits>
#include <tuple>
#include <list>
#include <sstream>



#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

inline const std::string MAGIC_NUMBER = "XXXXXXXX";

typedef struct result_t {
    bool succeeded;
    std::string client_response;
}result_t;

typedef struct received_t {
    std::string func_id;
    std::string exp_resp;
    std::string data;
}received_t;

typedef std::vector<std::vector<double>> matrix;

namespace Capstone {

    result_t parse_request(std::string receivedData, std::size_t bytes);

    std::string extract_json(std::string receivedData);

    std::vector<double> extract_vector(std::string vec_str);

    matrix extract_matrix(std::string mat_str);

    std::tuple<std::vector<double>, std::vector<matrix>> parse_data(received_t msg); 

    std::string vector_to_matrix(std::vector<double> vec);

    std::string serialize_matrix(matrix mat);

    std::string serialize_vector(std::vector<double> vec);

    void map_func(uint32_t id, std::tuple<std::vector<double>, std::vector<matrix>> data, result_t *resp);
}

#endif
