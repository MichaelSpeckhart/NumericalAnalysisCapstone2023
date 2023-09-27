#pragma once

#ifndef INCLUDED_HTTP_SERVER
#define INCLUDED_HTTP_SERVER

#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include <algorithm>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <memory>
#include <string>
#include <vector>

#include <boost/beast/http.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/beast/version.hpp>

namespace beast = boost::beast;
namespace http = beast::http;
namespace asio = boost::asio;

namespace Capstone {
    class HTTPServer {
        public:
            HTTPServer(int port, std::string ipAddress);
            ~HTTPServer();
            void startListen();
            void init();

        private:
            int mPort;
            std::string mIpAddress;
            long mIncomingMessage;
            std::string mServerMsg;

            int startServer();
            void closeServer();
            void acceptConnection(int &new_socket);
            std::string handleResponse();
            void sendResponse();
            std::string constructResponse(std::string responseBody);
            void run(asio::io_context& context);
            void configureServerSettings(asio::io_context &context, asio::ip::tcp::acceptor &acceptor);
            void startAccepts(asio::ip::tcp::acceptor &acceptor, asio::io_context &context);
            void handleClients(std::shared_ptr<boost::asio::ip::tcp::socket> sSocket);
            void printClientMessage(std::vector<char> clientData);
    };
}

#endif