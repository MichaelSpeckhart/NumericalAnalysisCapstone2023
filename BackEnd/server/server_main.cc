#include "http_server.h"

int main(int argc, char **argv) {

    using namespace http;

    HTTPServer server = HTTPServer("127.0.0.1", 8080);

    server.startListen();

    return 0;


}