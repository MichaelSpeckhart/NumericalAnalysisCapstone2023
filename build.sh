# build script --> needs modification for ubuntu environment
#
#
#
#c++ -O3 -Wall -shared -std=c++11 $(python3 -m pybind11 --includes) functions.cc -o functions$(python3-config --extension-suffix)

#modified compilation line for ubuntu
g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` functions.cc -o functions$(python3-config --extension-suffix)
