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
            result.msg = "Error: Ill formed JSON, cannot properly extract\n";
            return result;
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

        auto res = Capstone::mapIdToFunction<std::vector<std::vector<double>>>(data);

        
        std::string matrixForClient = serialize_matrix(res);

        result.clientMatrix = matrixForClient;

        std::cout << "Matrix for client: " << matrixForClient << "\n";


    // } catch (boost::system::system_error& bException) {
    //     std::cerr << "Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
    //     return false;
    // }
    
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

/**
 * @brief Deserialize the matrix from a vector of bytes 
 * 
 * @tparam T 
 * @param byte_data 
 * @param rows 
 * @param cols 
 * @return std::vector<std::vector<T>> 
 */
std::vector<std::vector<double>> Capstone::deserialize_matrix(const std::vector<double> matrixVals, size_t rows, size_t cols) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix[i][j] = matrixVals[i * cols + j];
        }
    }

    return matrix;
}

/**
 * @brief Serialize the matrix so that it can be sent over to the client
 * 
 * @param matrix 
 * @return std::string 
 */
std::string Capstone::serialize_matrix(const std::vector<std::vector<double>> matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    std::string stringMat = "";
    stringMat.append(std::to_string(rows) + ",");
    stringMat.push_back(static_cast<char>('0' + cols));
    stringMat.push_back('\n');

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            stringMat.append(std::to_string(static_cast<int>(matrix[i][j])) + ",");
        }
    }

    stringMat.resize(stringMat.size() -1);

    return stringMat;
    

}

/**
    * @brief Match the function ID given in the JSON to one of the matrix functions
    * 
    * @param funcID 
    */
   template <typename T>
   T Capstone::mapIdToFunction(FunctionData& data) {
       switch (data.mFuncId) {
            case 0x10: std::cout << "Adding Matrices Together\n"; return {{2, 4}, {6, 8}};
            case 0x11: std::cout << "";
            case 0x12: std::cout << "Transposing a matrix\n"; return transpose(data.mFirstMatrix);
            case 0x13: std::cout << "";
            case 0x14: std::cout << "";
            case 0x15: std::cout << "";
            default: return {{}};
       }

       return {{}};
    }   

