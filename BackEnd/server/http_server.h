#ifndef INCLUDED_HTTP_SERVER
#define INCLUDED_HTTP_SERVER

#include <algorithm>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>

#include <string>
#include <vector>

namespace http {
    class HTTPServer {

        public:
            HTTPServer(std::string ipAddress, int port);
            ~HTTPServer();
            void startListen();

        private:
            std::string mIpAddress;
            int mPort;
            int mSocket;
            int mNewSocket;
            long mIncomingMessage;
            struct sockaddr_in mSocketAddress;
            unsigned int mSocketAddrLen;
            std::string mServerMsg;


            int startServer();
            void closeServer();
            void acceptConnection(int &new_socket);
            std::string handleResponse();
            void sendResponse();
            void printClientMessage(std::vector<char> clientData);

    };
}

#endif