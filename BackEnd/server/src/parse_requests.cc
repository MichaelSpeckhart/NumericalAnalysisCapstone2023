#include <typeinfo>

#include "parse_requests.h"
#include "function_map.h"

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
    std::cout << " parse_request called\n" << std::endl;
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
    /* msg.exp_resp = jsonTree.get<std::string>("exp_resp");  */
    msg.data = jsonTree.get<std::string>("data");
    
    uint32_t id = std::stoi(msg.func_id); /* convert from string to uint32_t */
    std::tuple<std::vector<double>, std::vector<matrix>> data = Capstone::parse_data(msg); /* data tuple */
    /* uint32_t exp_resp = std::stoul(msg.exp_resp); */
    Capstone::map_func(id, data, &result);
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
    std::cout << " extract_json called\n" << std::endl;
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
    std::cout << " extract vector called\n" << std::endl;
    std::cout << vec_str << std::endl;
    std::vector<double> result;
    result.reserve(vec_str.length() / 2);
    std::istringstream ss(vec_str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        std::cout << "Token: " << token << std::endl;
        double num = std::stod(token);
        result.push_back(num);
    }
    return result;
}

matrix Capstone::extract_matrix(std::string mat_str){
    std::cout << " extract matrix called\n" << std::endl;
    std::cout << mat_str << std::endl;
    size_t newline = mat_str.find('\n');
    /* extract dimensions without newline */
    std::string dim = mat_str.substr(0, newline); 
    /* extract matrix data without the newline */
    std::string data = mat_str.substr(newline + 1, mat_str.length() - (newline + 1));

    int num_rows, num_cols = 0;
    /* spliting the dimensions of variable length */
    std::vector<double> dim_vec = Capstone::extract_vector(dim);
    num_rows = (int) dim_vec[0];
    num_cols = (int) dim_vec[1];

    std::vector<double> data_vec = Capstone::extract_vector(data);

    matrix result(num_rows, std::vector<double>(num_cols));

    for(int row = 0; row < num_rows; row++){
        for(int col = 0; col < num_cols; col++){
            result[row][col] = data_vec[row * num_cols + col];
        }
    }
    return result;
}


std::tuple<std::vector<double>, std::vector<matrix>> Capstone::parse_data(received_t msg) {
    std::cout << " parse_data called\n" << std::endl;
    size_t d_pos = 0;
    std::string data = msg.data;

    std::vector<double> scalars;
    std::vector<std::vector<double>> vectors; 
    std::vector<matrix> matrices;

    /* TODO: modifiy to account for the fact that vectors are treated matrices */

    while((d_pos = data.find(MAGIC_NUMBER)) != std::string::npos){
        std::cout << "PASS: " + data << std::endl;
        std::string obj = data.substr(0, d_pos);
        if(obj.length() == 1){ /* it is a scalar */
            scalars.push_back(std::stod(obj));
            data.erase(0, d_pos + MAGIC_NUMBER.length());
        }else{ /* it is a matrix */
            matrices.push_back(Capstone::extract_matrix(obj));
            data.erase(0, d_pos + MAGIC_NUMBER.length());
        }
    }
    return std::make_tuple(scalars, matrices);
}

/**
 * Converts a matrix to a column vector.
 *
 * @param mat The input matrix.
 *
 * @returns A column vector representation of the matrix.
 */
std::vector<double> Capstone::matrix_to_colvector(matrix mat){
    /* assumes matrix is in the right format */
    std::vector<double> result;
    for(int row = 0; row < (int) mat.size(); row++){
        result.push_back(mat[row][0]);
    }
    return result;
}

std::string Capstone::serialize_matrix(matrix mat){
    std::cout << " serialize_matrix called\n" << std::endl;
    std::string result = "";
    std::stringstream token;
    int num_rows = mat.size();
    int num_cols = mat[0].size();
    result += std::to_string(num_rows) + "," + std::to_string(num_cols) + '\n';
    /* converting the matrix type to a string */
    for(int row = 0; row < (int) num_rows; row++){
        for(int col = 0; col < (int) num_cols; col++){
            token << std::fixed << std::setprecision(2) << mat[row][col];
            result += token.str();
            token.str(std::string()); /* clear the stringstream */
            if(row != num_rows - 1 || col != num_cols - 1){ /* only add comma if not the end */
                result += ",";
            }
        }
    }
    return result; /* add newline to the end */
}

std::string Capstone::serialize_vector(std::vector<double> vec){
    std::cout << " serialize_vector called\n" << std::endl;
    std::string result = "";
    for(int i = 0; i < (int) vec.size(); i ++){
        result += std::to_string(vec[i]);
        if(i != (int) vec.size() - 1){ /* only add comma if not the end */
            result += ",";
        }
    }
    return result;
}

void Capstone::map_func(uint32_t id, std::tuple<std::vector<double>, std::vector<matrix>> data, result_t *resp){
    std::cout << "map_func called\n" << std::endl;
    switch(id){
        case 0x10:{ 
            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            matrix m1 = mat_list[0];
            matrix m2 = mat_list[1];
            std::cout << "Matrix 1: " << Capstone::serialize_matrix(m1) << std::endl;
            std::cout << "Matrix 2: " << Capstone::serialize_matrix(m2) << std::endl;
            matrix sum = sum_matrix(m1, m2);
            std::string result = Capstone::serialize_matrix(sum);
            resp->client_response = result; /* attach the sum to the result struct in other scope */
            resp->succeeded = true; /* TODO: add a checking mechanism */
            break;
        }
        case 0x11:{ /* multiply */
            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            std::vector<double> scalars = std::get<0>(data); /* access the list of scalars from tuple */
            matrix m1 = mat_list[0];
            double s1 = scalars[0];
            std::cout << "Matrix 1: " << Capstone::serialize_matrix(m1) << std::endl;
            std::cout << "Scalar 1: " << s1 << std::endl;
            matrix mat = scalar_multiply(m1, s1);
            std::string result = Capstone::serialize_matrix(mat);
            resp->client_response = result; /* attach the sum to the result struct in other scope */
            resp->succeeded = true; /* TODO: add a checking mechanism */
            break;
        }
        case 0x12:{ /* transpose */
            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            matrix result = transpose(mat_list[0]);
            std::string result_str = Capstone::serialize_matrix(result);
            resp->client_response = result_str;
            resp->succeeded = true;
            break;
        }
        case 0x13:{ /* inverse */
            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            matrix result = mat_list[0];
            if(matrix_inverse(result)){
                std::string result_str = Capstone::serialize_matrix(result); 
                std::cout << "Inverse: " << result_str << std::endl;
                resp->client_response = result_str;
                resp->succeeded = true;
            }else{
                resp->client_response = "Error: Matrix is not invertible\n";
                resp->succeeded = false;
            }
            break;
        }
        case 0x20:{ /* gauss elimination */
            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            std::cout << "Im working 1" << std::endl;
            matrix m1 = mat_list[0];
            std::vector<double> v1 = matrix_to_colvector(mat_list[1]);;

            if(gaussian_elimination(m1, v1)){ /* result is stored in v1 */
                std::cout << "Result: " << Capstone::serialize_vector(v1) << std::endl;
                std::string result_str = Capstone::serialize_vector(v1);
                resp->client_response = result_str;
                resp->succeeded = true;
            }else{
                resp->client_response = "Error completing gauss elimination\n";
                resp->succeeded = false;
            }
            break;
        }
        case 0x21:{ /* lu factorization */
            std::vector<matrix> mat_list = std::get<1>(data);
            matrix m1 = mat_list[0];
            std::vector<int> res = lu_factorization_inplace(m1);
            std::vector<double> convert(res.begin(), res.end());
            resp->client_response = serialize_vector(convert);
            resp->succeeded = true;
            break;
        }
        case 0x30:{ /* jacobi */
            const double TOLERANCE = 1e-6;
            const int MAX_ITER = 100;

            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            matrix m1 = mat_list[0];
            /* grab the vector */
            std::vector<double> v1 = matrix_to_colvector(mat_list[1]);

            std::vector<double> res = jacobi_iteration(m1, v1, TOLERANCE, MAX_ITER);
            resp->client_response = serialize_vector(res);
            resp->succeeded = true;
            break;
        }
        case 0x31:{ /* gauss sidel */
            const int MAX_ITER = 100;

            std::vector<matrix> mat_list = std::get<1>(data); /* access the list of matrices from tuple */
            matrix m1 = mat_list[0];
            std::vector<double> v1 = matrix_to_colvector(mat_list[1]);
            std::vector<double> x(v1.size(), 0.0);

            if(gauss_seidel(m1, v1, x, MAX_ITER)){
                resp->client_response = serialize_vector(x);
                resp->succeeded = true;
            }else{
                resp->client_response = "Error completing gauss sidel \n";
                resp->succeeded = false;
            }
            break;
        }
        default:
            resp->client_response = "Error generating function mapping\n";
            break;
    }
}


