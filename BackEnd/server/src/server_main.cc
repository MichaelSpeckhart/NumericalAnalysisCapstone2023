#include "http_server.h"

#include <iostream>
#include <string>


struct arg_t {
  int port;                    // The port on which to listen
  std::string ipAddress;       //IP Address for connection

  /// Construct an arg_t from the command-line arguments to the program
  ///
  /// @param argc The number of command-line arguments passed to the program
  /// @param argv The list of command-line arguments
  ///
  /// @throw An intmd5eger exception (1) if an invalid argument is given, or if
  ///        `-h` is passed in
  arg_t(int argc, char **argv) {
    long opt;
    while ((opt = getopt(argc, argv, "p:f:k:ht:b:i:u:d:r:o:a:")) != -1) {
      switch (opt) {
        case 'p':
          port = atoi(optarg);
          break;
        case 'a':
          ipAddress = std::string(optarg);
          break;
        default: 
          throw 1;
          return;
      }
    }
  }

  static void usage(char *progname) {
    std::cout << basename(progname) << ": company user directory server\n"
         << "  -p [int]    Port on which to listen for incoming connections\n"
         << "  -a [string] IP Address that the server is connected to\n"
         << "  -h          Print help (this message)\n";
  }
};


int main(int argc, char **argv) {
    using namespace Capstone;
    arg_t *args;
    try {
        args = new arg_t(argc, argv);
    } catch (int i) {
        arg_t::usage(argv[0]);
        return 1;
    }
 

    HTTPServer server = HTTPServer(args->port, args->ipAddress);

    server.init();


    return 0;


}