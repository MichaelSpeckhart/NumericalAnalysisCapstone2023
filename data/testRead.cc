#include <fstream>
#include <iostream>
#include <algorithm>

int main() {
    std::ifstream file("matrices/will199.mtx");
    int num_row = 0, num_col = 0, num_lines = 0;

    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows, columns, and non-zero values
    file >> num_row >> num_col >> num_lines;

    // Create 2D array and fill with zeros
    double* matrix = new double[num_row * num_col]();
    // Note: use () instead of {} to zero-initialize the array.

    // fill the matrix with data
    for (int l = 0; l < num_lines; l++)
    {
        double data;
        int row, col;
        file >> row >> col >> data;
        matrix[(row - 1) + (col - 1) * num_row] = data;
    }

    file.close();
    std::fstream key;
    key.open("test.txt",std::ios::out);

    // Print the matrix
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            key << matrix[i + j * num_row] << " ";
        }
        key << std::endl;
    }
    key.close();

    // Deallocate memory
    delete[] matrix;

    return 0;
}
