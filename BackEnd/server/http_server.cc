#include "http_server.h"
#include "function_map.h"

#include <iostream>


namespace Capstone {

    HTTPServer::HTTPServer(int port, std::string ipAddress) 
        : mPort(port), mIpAddress(ipAddress){ }


    HTTPServer::~HTTPServer() {
        
    }

    /**
     * @brief Handle the client once the connection is accepted and asynchronously read the data being inputted into a buffer
     * and then that buffer will be parsed and processed.
     * 
     * @param bSocket 
     */
    void HTTPServer::handleClients(std::shared_ptr<boost::asio::ip::tcp::socket> bSocket) {
        std::vector<char> received_data(1024);
        bSocket->async_read_some(asio::buffer(received_data),
            [this, bSocket, &received_data](const boost::system::error_code& error, std::size_t bytesRead) {
                if (!error) {
                    std::string clientData(received_data.begin(), received_data.begin() + bytesRead);
                    std::cout << clientData << "\n";
                    handleClients(bSocket);
                } else if (error == asio::error::eof){
                    std::cout << "Client disconnected\n";
                } else {
                    std::cerr << "Error reading from client: " << error.message() << "\n";
                }
        });
    }

    /**
     * @brief Start asynchronously accepting connections, once a connection is accepted, we want to handle it. If their is an error binding this
     * bSocket to the accept, then an error is thrown. Keep calling start accepts
     * 
     */
    void HTTPServer::startAccepts(asio::ip::tcp::acceptor &acceptor, asio::io_context &context) {
        try {
            std::shared_ptr<boost::asio::ip::tcp::socket> bSocket(new boost::asio::ip::tcp::socket(context));
            std::cout << "Waiting for client" << std::endl;
            acceptor.async_accept(*bSocket, [this, bSocket, &acceptor, &context](const boost::system::error_code& bErrorCode) {
                if (!bErrorCode) {
                    std::cout << "Client connected" << std::endl;
                    handleClients(bSocket);
                } else {
                    std::cout << "Boost Error: " << " (" << bErrorCode.value() << ") " << bErrorCode.message() << std::endl;
                }
                startAccepts(acceptor, context);
            });
        } catch (const boost::system::system_error& bException) {
            std::cout << "Boost Error: " << bException.what() << std::endl;
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

            configureServerSettings(bContext, bAcceptor);

            startAccepts(bAcceptor, bContext);

            bContext.run();
        } catch (const boost::system::system_error &bException) {
            std::cout << "Boost system error: " << bException.what() << std::endl;
        }
    }

    /**
     * @brief Create a server endpoint and bind the acceptor to that endpoint
     * 
     * @param context 
     * @param acceptor 
     */
    void HTTPServer::configureServerSettings(asio::io_context &context, asio::ip::tcp::acceptor &acceptor) {
        try {
            asio::ip::tcp::endpoint sEndpoint(asio::ip::address_v4(), this->mPort);
            acceptor.open(sEndpoint.protocol());
            acceptor.bind(sEndpoint);
            acceptor.listen();
        } catch (const boost::system::system_error &bException) {
            std::cout << "Boost system error: " << bException.what() << std::endl;
        }
    }



    








}