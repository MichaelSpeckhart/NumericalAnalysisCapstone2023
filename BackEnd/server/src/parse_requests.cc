#include "parse_requests.h"
#include "format.h"
#include "http_server.h"

const std::string EMPTY = "";


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
    try {
        std::string jsonData = extract_json(receivedData);

        if (jsonData == EMPTY) {
            std::cerr << "Error: ill formed JSON, cannot properly extract: " << std::endl;
        }

        std::istringstream jsonStream(jsonData);
        pt::ptree jsonTree;
        pt::read_json(jsonStream, jsonTree);






    } catch (boost::system::system_error& bException) {
        std::cerr << "Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
        return false;
    }
    
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
 * @brief 
 * 
 * @param matrixData 
 * @return std::vector<double> 
 */
std::vector<double> Capstone::parse_matrix(std::string matrixData) {

}