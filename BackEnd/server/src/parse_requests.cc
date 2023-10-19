#include <typeinfo>

#include "parse_requests.h"
#include "../test/test_functions.h"

const std::string EMPTY = "";
const size_t MAX_SIZE = 1024;

namespace pt = boost::property_tree;


/**
 * @brief Once the service handles the client and receives the data, parse the json to extract
 * the function id, the matrix data, and etc.
 * 
 * @param receivedData 
 * @param bytes 
 * @return true 
 * @return false 
 */
result_t Capstone::parse_request(std::string receivedData, std::size_t bytes) {
    // try {
        result_t result;
        std::string jsonData = extract_json(receivedData);

        if (jsonData == EMPTY) {
            result.succeeded = false;
            result.client_response = "Error: Ill formed JSON, cannot properly extract\n";
            return result;
        }

        std::cout << jsonData << std::endl;

        std::istringstream jsonStream(jsonData);
        pt::ptree jsonTree;
        pt::read_json(jsonStream, jsonTree);
    
        received_t msg;
        msg.func_id = jsonTree.get<std::string>("operation"); 
        /*  EXPECTED RESPONSE
         _______________
        |  Scalar --> 0 |
        |  Vector --> 1 |
        |  Matrix --> 2 |
        |_______________|
        */
        msg.exp_resp = jsonTree.get<std::string>("exp_resp");  
        msg.data = jsonTree.get<std::string>("data");
        
        int id = std::stoi(msg.func_id); /* convert from string to int */
        std::tuple<std::vector<double>,std::vector<std::vector<double>>, std::vector<matrix>> data = Capstone::parse_data(msg); /* data tuple */
        int exp_resp = std::stoi(msg.exp_resp);

        switch (exp_resp) {
            case 0:{
                Scalar res_scalar(Capstone::mapIdToFunction(id, data));
                result.succeeded = true;
                result.client_response = std::to_string(res_scalar);
                break;
            }
            case 1:{
                Vector res_vector(Capstone::mapIdToFunction(id, data));
                result.succeeded = true;
                result.client_response = Capstone::serialize_vector(res_vector.values);
                break;
            }
            case 2:{
                Matrix res_matrix(Capstone::mapIdToFunction(id, data));
                result.succeeded = true;
                result.client_response = Capstone::serialize_matrix(res_matrix.values);
                break;
            }
            default:{
                result.succeeded = false;
                result.client_response = "Error in mapping backend function for case " + std::to_string(exp_resp) + "\n";
                break;
            }
        }
    return result;
}

/**
 * @brief Since the client sends over more than just JSON data, need to extract the JSON and return it so that Boost throws
 * no errors
 * 
 * @param receivedData 
 * @return std::string 
 */
std::string Capstone::extract_json(std::string receivedData) {
    size_t startPos = receivedData.find('{');
    
    if (startPos != std::string::npos) { 
        size_t endPos = receivedData.find('}', startPos);

        if (endPos != std::string::npos) { 
            size_t length = endPos - startPos - 1;
            std::string result = receivedData.substr(startPos, length + 2);

            return result;
        } else {
            std::cout << "Matching '}' not found." << std::endl;
        }
    } else {
        std::cout << "'{' not found." << std::endl;
    }

    return EMPTY;
}

std::vector<double> Capstone::extract_vector(std::string vec_str){
    size_t pos = 0;
    std::vector<double> result;
    while((pos = vec_str.find(',')) != std::string::npos){
        result.push_back(std::stod(vec_str.substr(0, pos + 1)));
        vec_str.erase(0, pos + 1);
    }
    result.push_back(std::stod(vec_str)); //get the last element
    return result;
}

matrix Capstone::extract_matrix(std::string mat_str){
    size_t newline = mat_str.find('\n');
    std::string dim = mat_str.substr(0, newline);
    std::string data = mat_str.substr(newline + 1, mat_str.length() - 1);
    size_t curr = 0;

    curr = dim.find(',');
    int num_rows = std::stoi(dim.substr(0, curr + 1));
    dim.erase(0, curr + 1);
    int num_cols = std::stoi(dim);

    vector<double> data_vec;

    while((curr = data.find(',')) != std::string::npos){
        data_vec.push_back(std::stod(data.substr(0, curr + 1)));
        data.erase(0, curr + 1);
    }


    matrix result(num_rows, std::vector<double>(num_cols));

    for(int row = 0; row < num_rows; row++){
        for(int col = 0; col < num_cols; col++){
            result[row][col] = data_vec[row * num_cols + col];
        }
    }
    return result;
}


std::tuple<std::vector<double>,std::vector<std::vector<double>>, std::vector<matrix>> Capstone::parse_data(received_t msg) {
    size_t d_pos = 0;
    std::string data = msg.data;

    std::vector<double> scalars;
    std::vector<std::vector<double>> vectors; 
    std::vector<matrix> matrices;

    while((d_pos = data.find(MAGIC_NUMBER)) != std::string::npos){
        std::string obj = data.substr(0, d_pos + 1);
        if(obj.length() == 1){ /* it is a scalar */
            scalars.push_back(std::stod(obj));
        }else if(obj.find('\n') == std::string::npos){ /* it is a vector */
            vectors.push_back(Capstone::extract_vector(obj));
        }else{ /* it is a matrix */
            matrices.push_back(Capstone::extract_matrix(obj));
        }
    }
        
    return std::make_tuple(scalars, vectors, matrices);
}

std::string Capstone::serialize_matrix(matrix mat){
    std::string result = "";
    int num_rows = mat.size();
    int num_cols = mat[0].size();
    result += std::to_string(num_rows) + ", " + std::to_string(num_cols) + '\n';
    /* converting the matrix type to a string */
    for(int row = 0; row < (int) num_rows; row++){
        for(int col = 0; col < (int) num_cols; col++){
            result += std::to_string(mat[row][col]);
            if(row != num_rows - 1 || col != num_cols - 1){ /* only add comma if not the end */
                result += ", ";
            }
        }
    }
    return result + '\n'; /* add newline to the end */
}

std::string Capstone::serialize_vector(std::vector<double> vec){
    std::string result = "";
    for(int i = 0; i < vec.size(); i ++){
        result += std::to_string(vec[i]);
        if(i != vec.size() - 1){ /* only add comma if not the end */
            result += ", ";
        }
    }
    return result;
}


/**
    * @brief Match the function ID given in the JSON to one of the matrix functions
    * 
    * @param funcID 
    */

   Result Capstone::mapIdToFunction(int id, std::tuple<std::vector<double>,std::vector<std::vector<double>>, std::vector<matrix>> data) {
       switch (id) {
            case 0x10:{
                std::cout << "Adding two matrices\n" << std::endl;
                std::vector<matrix> m_list = std::get<2>(data);
                matrix m1 = m_list[0];
                matrix m2 = m_list[1];
                return Matrix(sum_matrix(m1, m2));
            }
            case 0x11:{ 
                std::cout << "";
                break;
            }
            case 0x12:{
                break;
            }
            case 0x13:{
                break;
            }
            case 0x14:{ 
                break;
            }
            case 0x15:{
                break;
            }
            default:{ 
                break;
            }
       }
       return Result;
    }   

