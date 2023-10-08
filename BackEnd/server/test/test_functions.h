#pragma once

#include <vector>

using namespace std;

vector<vector<double>> transpose(vector<vector<double>> m1)
{
    const int d1 = static_cast<int>(m1.size()),
              d2 = static_cast<int>(m1[0].size());

    vector<vector<double>> m2(d2);

//#pragma omp parallel
    for (auto i = 0; i < d2; ++i)
        m2[i] = vector<double>(d1);

//#pragma omp parallel
    for (auto i = 0; i < d1; ++i)
        for (auto j = 0; j < d2; ++j)
            m2[j][i] = m1[i][j];

    return m2;
}

vector<vector<double>> sum_matrix(vector<vector<double>> m1,
                                  const vector<vector<double>> m2)
{
    //   NB: in the future it may be important to use the compiler defined sizes
    //   for platform portability
    size_t d1 = static_cast<size_t>(m1.size()),
           d2 = static_cast<size_t>(m1[0].size());

    // If the matrices do not have the same dimensions we return a dimension error
    // for the GUI to handle
    // if (d1 != static_cast<size_t>(m2.size()) || d2 != static_cast<size_t>(m2[0].size()))
    //     return d_err;

#//pragma omp parallel
    for (size_t i = 0; i < d1; ++i)
        for (size_t j = 0; j < d2; ++j)
            m1[i][j] -= m2[i][j];

    return m1;
}