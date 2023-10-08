import json
import socket

def main():
    server_ip = "127.0.0.1" 
    server_port = 1025  

    client_data = '{ "operation": 0x10, "args": 1, "matrixData":"2,2 \r\n 1,2,3,4"}'

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client_socket:
        client_socket.connect((server_ip, server_port))
        print("Connected to the server")

        while True:
            if not client_data:
                break
            client_socket.send(client_data.encode('utf-8'))
            response = client_socket.recv(1024).decode()
            print("Server response:", response)

        print("Connection closed")

if __name__ == "__main__":
    main()

