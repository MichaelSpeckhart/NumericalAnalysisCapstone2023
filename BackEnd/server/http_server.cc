#include "http_server.h"

#include <sstream>
#include <iostream>
#include <vector>

const int BUFFER_SIZE = 30720;
 

namespace http {
    /**
     * @brief Construct a new HTTPServer::HTTPServer object and initialize the port, ip, sock struct
     * 
     * @param ipAddress 
     * @param port 
     */
    HTTPServer::HTTPServer(std::string ipAddress, int port) : mIpAddress(ipAddress), mPort(port), 
        mSocketAddrLen(sizeof(mSocketAddress)) {
        mSocketAddress.sin_family = AF_INET;
        mSocketAddress.sin_port = htons(mPort);
        mSocketAddress.sin_addr.s_addr = inet_addr(mIpAddress.c_str());

        if (startServer() != 0) {
            std::cout << "Failed to start server with PORT: " << ntohs(mSocketAddress.sin_port) << "\n";
        }
    }

    HTTPServer::~HTTPServer() {
        closeServer();
    }

    /**
     * @brief Start the server by seeing if the socket can be created with the appropriate address, handle 
     * errors if neccesary
     * 
     * @return int 
     */
    int HTTPServer::startServer() {
        mSocket = socket(AF_INET, SOCK_STREAM, 0);
        if (mSocket < 0) {
            std::cout << "Cannot create socket\n";
            return 1;
        }

        if (bind(mSocket, (sockaddr *)&mSocketAddress, mSocketAddrLen) < 0) {
            std::cout << "Cannot connect socket to address\n";
            return 1;
        }

        return 0;
    }

    void HTTPServer::closeServer() {
        close(mSocket);
        close(mNewSocket);
        exit(0);
    }

    /**
     * @brief Start listening on the socket for clients, max 20 clients in the queue at a time
     * 
     */
    void HTTPServer::startListen() {
        if (listen(mSocket, 20) < 0) {
            std::cout << "Socket listen failed\n";
        }

        std::cout << " Listening on ADDRESS " << inet_ntoa(mSocketAddress.sin_addr) 
            << " PORT: " << ntohs(mSocketAddress.sin_port) << "\n";
        ssize_t bytesReceived;

        while (true) {
            acceptConnection(mNewSocket);

            std::vector<char> buffer(BUFFER_SIZE, 0);
            bytesReceived = read(mNewSocket, buffer.data(), BUFFER_SIZE);
            if (bytesReceived < 0) {
               std::cerr << "Failed to read bytes from client socket connection\n";
               exit(0);
            }

            std::cout <<  "Received Request from client\n";

            printClientMessage(buffer);
            mServerMsg = handleResponse();
            sendResponse();

            close(mNewSocket);
        }
    }

    /**
     * @brief Print the message from the client, will maybe log it or something later on
     * 
     * @param clientData 
     */
    void HTTPServer::printClientMessage(std::vector<char> clientData) {
        for (size_t i = 0; i < clientData.size(); ++i) {
            std::cout << clientData[i];
        }

        std::cout << "\n";
    }

    /**
     * @brief Accept the connection, handle the error and print it if neccessary
     * 
     * @param newSocket 
     */
    void HTTPServer::acceptConnection(int &newSocket) {
        newSocket = accept(mSocket, (sockaddr *)&mSocketAddress, &mSocketAddrLen);
        if (newSocket < 0) {
             std::cerr << "Server failed to accept incoming connection from ADDRESS: " << inet_ntoa(mSocketAddress.sin_addr) 
                << "; PORT: " << ntohs(mSocketAddress.sin_port) << "\n";
        }
    }

    /**
     * @brief Construct a response, for right now this is just a simple hello back
     * 
     * @return std::string 
     */
    std::string HTTPServer::handleResponse() {
        std::string response = "Hello from server";

        return response;
    }

    /**
     * @brief Write the response to the socket and send it back to client, handle errors
     * 
     */
    void HTTPServer::sendResponse() {
        ssize_t bytesSent;

        bytesSent = write(mNewSocket, mServerMsg.c_str(), mServerMsg.size());

        if (bytesSent == mServerMsg.size()) {
            std::cout << "Server Response sent to client\n";
        }
        else {
            std::cerr << "Error sending response to client\n";
        }
    }



}