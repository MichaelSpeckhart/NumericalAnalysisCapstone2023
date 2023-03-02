# build script --> needs modification for ubuntu environment
#
#
c++ -O3 -Wall -shared -undefined dynamic_lookup -std=c++11 $(python3 -m pybind11 --includes) Backend/functions.cc -o build/functions$(python3-config --extension-suffix)

#modified compilation line for ubuntu
# g++ -O3 -Wall -shared -std=c++17 -fPIC `python3 -m pybind11 --includes` functions.cc -o functions$(python3-config --extension-suffix) 
