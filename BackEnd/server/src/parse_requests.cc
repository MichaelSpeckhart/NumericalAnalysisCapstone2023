#include "parse_requests.h"
#include "format.h"

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
    std::istringstream jsonStream(receivedData);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(jsonStream, pt);

    std::string name = pt.get<std::string>("name");

    std::cout << name << std::endl;

    return true;
}