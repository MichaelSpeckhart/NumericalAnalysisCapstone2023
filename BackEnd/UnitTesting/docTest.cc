#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
//This file just defines are Doctest config

int factorial(int number) { return number > 1 ? factorial(number - 1) * number : 1; }

TEST_CASE("testing the factorial function") {
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}