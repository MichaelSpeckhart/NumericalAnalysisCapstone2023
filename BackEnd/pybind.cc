#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void init_functions(py::module &);
void init_functionsCOO(py::module &);
void init_functionsCSR(py::module &);
void init_functionsCSC(py::module &);

PYBIND11_MODULE(functions, m) {
  init_functions(m);
  // init_functionsCOO(m);
  // init_functionsCSR(m);
  // init_functionsCSC(m);
}