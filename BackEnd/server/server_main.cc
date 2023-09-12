#include "http_server.h"

#include <iostream>


void printHelpArguments() {
    std::cout << "Usage";


}


int main(int argc, char **argv) {

    using namespace Capstone;

 

    HTTPServer server = HTTPServer(8080, "127.0.0.1");

    server.init();


    return 0;


}