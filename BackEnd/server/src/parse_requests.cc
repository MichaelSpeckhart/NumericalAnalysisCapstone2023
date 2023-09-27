#include "parse_requests.h"
#include "format.h"
#include "http_server.h"
#include "function_map.h"
#include <typeinfo>

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
bool Capstone::parse_request(std::string receivedData, std::size_t bytes) {
    // try {
        std::string jsonData = extract_json(receivedData);

        if (jsonData == EMPTY) {
            std::cerr << "Error: Ill formed JSON, cannot properly extract: " << std::endl;
            return false;
        }

        std::cout << jsonData << std::endl;

        std::istringstream jsonStream(jsonData);
        pt::ptree jsonTree;
        pt::read_json(jsonStream, jsonTree);
    
        FunctionData data;

        data.mFuncId = jsonTree.get<size_t>("operation"); 
        size_t numArguments = jsonTree.get<size_t>("args"); 

        // If the number of arguments is only 1 (inverse, transpose, etc), just need to get one matrix
        if (numArguments == 1) {
            data.mFirstMatrix = parse_matrix(jsonTree.get<std::string>("matrixData"));
        } else {
            // Check if (matrix, matrix), (matrix, scalar), (matrix, vector), etc
            if (jsonTree.find("secondMatrixData") != jsonTree.not_found()) {
                data.mSecondMatrix = parse_matrix(jsonTree.get<std::string>("secondMatrixData"));
            } else if (jsonTree.find("scalar") != jsonTree.not_found()) {
                data.mScalar = jsonTree.get<size_t>("scalar");
            } else if (jsonTree.find("vectorData") != jsonTree.not_found()) {
                // Handle a vector
            }

        }

        Capstone::mapIdToFunction(data);

    // } catch (boost::system::system_error& bException) {
    //     std::cerr << "Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
    //     return false;
    // }
    
    return true;
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

/**
 * @brief Extract the matrix into a standard vector<vector<>> format from the string
 * 2,2 \r\n 2,1,2,3
 * @param matrixData 
 * @return std::vector<double> 
 */
std::vector<std::vector<double>> Capstone::parse_matrix(std::string matrixData) {
    size_t rows = static_cast<size_t>(matrixData[0]) - '0';
    size_t cols = static_cast<size_t>(matrixData[3]) - '0';

    size_t delimiter = matrixData.find("\r\n");
    std::vector<double> matrix;
    if (delimiter != std::string::npos) {
        std::string matrixVals = matrixData.substr(delimiter + 2);
        std::istringstream matrixStream(matrixVals);

        size_t value;
        while (matrixStream  >> value) {
            matrix.push_back(value);
            char comma;
            matrixStream >> comma;

            if (comma != ',') {
                break;
            }
        }
    }
    
    auto result = deserialize_matrix(matrix, rows, cols);

    return result;
}