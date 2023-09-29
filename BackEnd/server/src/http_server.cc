#include "http_server.h"
#include "parse_requests.h"

#include <iostream>
#include <cstddef>

/**
 * @brief 
 * 
 */
namespace Capstone {

    /**
     * @brief Construct a new HTTPServer::HTTPServer object with the defined port number and IP Address
     * 
     * @param port 
     * @param ipAddress 
     */
    HTTPServer::HTTPServer(int port, std::string ipAddress) 
        : mPort(port), mIpAddress(ipAddress){ }


    /**
     * @brief Destroy the HTTPServer::HTTPServer object
     * 
     */
    HTTPServer::~HTTPServer() {
        
    }

    std::string HTTPServer::constructResponse(std::string responseBody) {
        std::string response = "HTTP/1.1 200 OK\r\n";
        response += "Content-Type: text/plain\r\n";
        response += "Content-Length: " + std::to_string(responseBody.length()) + "\r\n";
        response += "\r\n";
        response += responseBody;

        return response;
    }

    /**
     * @brief Handle the client once the connection is accepted and asynchronously read the data being inputted into a buffer
     * and then that buffer will be parsed and processed.
     * 
     * @param bSocket socket that will be used to connect to the client
     */
    void HTTPServer::handleClients(std::shared_ptr<boost::asio::ip::tcp::socket> bSocket) {
        //Create a vector of characters to load the received data into, shared pointer needs to be passed to async functions
        auto bReceivedData = std::make_shared<std::vector<char>>(1024);
        //Directly from the socket, asynchronously read the data coming in, this lambda is non-blocking and allows the main
        //thread to continue doing work
        // Asynchronously read data from the client into the buffer
        bSocket->async_read_some(boost::asio::buffer(*bReceivedData),
            [this, bSocket, bReceivedData](const boost::system::error_code& bErrorCode, std::size_t bytesRead) {
                if (!bErrorCode) {
                    //Converting received data to a string (decoding)
                    std::string clientData(bReceivedData->begin(), bReceivedData->begin() + bytesRead);
                    //Parse the incoming json
                    //TODO: validate the json to make sure fields match up
                    parse_request(clientData, bytesRead);
                    //Recursively call handleClients() 
                    std::string writeData = "2,2\n2,4,6,8";
                    std::string response = constructResponse(writeData);
                    boost::asio::async_write(*bSocket, boost::asio::buffer(response), [bSocket] (const boost::system::error_code& bError, 
                        std::size_t bytesWritten) {
                            if (!bError) {
                                std::cout << "Written Succesfully\n";
                            } else {
                                std::cerr << "Nothing was written\n";
                            }
                            bSocket->close();
                    });
                    handleClients(bSocket);
                } else if (bErrorCode == boost::asio::error::eof) {
                    std::cout << "Client disconnected." << std::endl;
                } else {
                    std::cout << "In handle clients: Boost Error: " << " (" << bErrorCode.value() << ") -> " << bErrorCode.message() << std::endl;
                }
            });
            bReceivedData.reset();
    }

    /**
     * @brief Start asynchronously accepting connections, once a connection is accepted, we want to handle it. If their is an error binding this
     * bSocket to the accept, then an error is thrown. Keep calling start accepts
     * 
     * @param acceptor 
     * @param context 
     */
    void HTTPServer::startAccepts(asio::ip::tcp::acceptor &acceptor, asio::io_context &context) {
        try {
            //Shared pointers work very well with async functions, wrapping the socket with a shared pointer
            std::shared_ptr<boost::asio::ip::tcp::socket> bSocket(new boost::asio::ip::tcp::socket(context));
            std::cout << "Waiting for client" << std::endl;
            acceptor.async_accept(*bSocket, [this, bSocket, &acceptor, &context](const boost::system::error_code& bErrorCode) {
                if (!bErrorCode) {
                    std::cout << "Client connected" << std::endl;
                    handleClients(bSocket);
                } else {
                    std::cout << "Acceptor error: Boost Error: " << " (" << bErrorCode.value() << ") -> " << bErrorCode.message() << std::endl;
                }
                startAccepts(acceptor, context);
            });
            bSocket.reset();
            
        } catch (const boost::system::system_error& bException) {
            std::cerr << "In accept clients: Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
        }
    
    }

    /**
     * @brief Initialize the endpoint with the ip address and port passed to the server
     * 
     */
    void HTTPServer::init() {
        try {
            asio::io_context bContext;
            asio::ip::tcp::acceptor bAcceptor(bContext);
            //Configure the settings for the acceptor
            configureServerSettings(bContext, bAcceptor);
            std::cout << "Server listening on port: " << this->mPort << " on Address: " << this->mIpAddress << std::endl;
            //Start accepting connections asynchronously
            startAccepts(bAcceptor, bContext);

            bContext.run();
        } catch (const boost::system::system_error &bException) {
            std::cerr << "In Init server: Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
        }
    }

    /**
     * @brief Create a server endpoint and bind the acceptor
     * 
     * @param context 
     * @param acceptor 
     */
    void HTTPServer::configureServerSettings(asio::io_context &context, asio::ip::tcp::acceptor &acceptor) {
        try {
            //Some basic asio configurations, nothing too special, basic syntax and semantics for creating the endpoint and acceptor
            asio::ip::tcp::endpoint sEndpoint(asio::ip::address_v4(), this->mPort);
            acceptor.open(sEndpoint.protocol());
            acceptor.bind(sEndpoint);
            acceptor.listen();
        } catch (const boost::system::system_error &bException) {
            std::cerr << "In Configure Server Settings: Boost System Error: " << " (" << bException.code() << ") -> " << bException.what() << std::endl;
        }
    }
}