#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>

const int PORT = 8080;
const int MAX_CONNECTIONS = 5;

int main() {
    // Create socket
    int server_socket = socket(AF_INET, SOCK_STREAM, 0);
    if (server_socket == -1) {
        std::cerr << "Error: failed to create socket\n";
        return EXIT_FAILURE;
    }

    // Bind socket to port
    struct sockaddr_in server_address;
    server_address.sin_family = AF_INET;
    server_address.sin_addr.s_addr = INADDR_ANY;
    server_address.sin_port = htons(PORT);
    if (bind(server_socket, (struct sockaddr *)&server_address, sizeof(server_address)) < 0) {
        std::cerr << "Error: failed to bind socket to port " << PORT << "\n";
        return EXIT_FAILURE;
    }

    // Listen for connections
    if (listen(server_socket, MAX_CONNECTIONS) < 0) {
        std::cerr << "Error: failed to listen for connections\n";
        return EXIT_FAILURE;
    }

    std::cout << "Server started and listening on port " << PORT << "\n";

    // Accept connections and handle requests
    struct sockaddr_in client_address;
    socklen_t client_address_size = sizeof(client_address);
    while (true) {
        int client_socket = accept(server_socket, (struct sockaddr *)&client_address, &client_address_size);
        if (client_socket < 0) {
            std::cerr << "Error: failed to accept connection\n";
            continue;
        }

        std::cout << "Client connected\n";

        // Read request from client
        char buffer[1024] = {0};
        int bytes_received = read(client_socket, buffer, sizeof(buffer));
        if (bytes_received < 0) {
            std::cerr << "Error: failed to read request\n";
            close(client_socket);
            continue;
        }

        // Process request
        std::string request(buffer, bytes_received);
        std::cout << "Received request: " << request << "\n";
        std::string response = "Hello, client!";
        write(client_socket, response.c_str(), response.length());

        // Close connection
        close(client_socket);
        std::cout << "Client disconnected\n";
    }

    // Close server socket
    close(server_socket);

    return EXIT_SUCCESS;
}
