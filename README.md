# A Programming Platform for Numerical Analysis - Backend

The goal of this project is to learn Numerical analysis through programming. We developed a library of matrix operations and solvers in C++ with no dependencies. We aimed for accessibility and high performance by making a front-end website to demonstrate the operations as well as developing parallel versions of the algorithms.

## Getting Started
Each set of operations can be found in ``BackEnd/src`` which has files for both dense and sparse data structures as well as parallel versions of the operations.

Examples of how to use the operations can be found in ``Backend/src/UnitTesting``

Our server code can be found in ``Backend/server``

Our Benchmarking code can be found in ``benchmarking``

## Running the unit tests

1. ``cd Backend/src/UnitTesting``
2. ``make`` to build the unit tests
3. ``./buildUnitTest/doctest`` to run the tests
