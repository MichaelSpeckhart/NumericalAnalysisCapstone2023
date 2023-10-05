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