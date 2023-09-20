#include <unordered_map>
#include <cstdlib>
#include <functional>
#include <iostream>


namespace Capstone {
    /**
    * @brief Match the function ID given in the JSON to one of the matrix functions
    * 
    * @param funcID 
    */
   void map_id_to_function(size_t funcID) {
       switch (funcID) {
            case 0x10: std::cout << "Adding Matrices Together\n"; break;
            case 0x11: std::cout << "";
            case 0x12: std::cout << "";
            case 0x13: std::cout << "";
            case 0x14: std::cout << "";
            case 0x15: std::cout << "";
       }
    }   
}   

