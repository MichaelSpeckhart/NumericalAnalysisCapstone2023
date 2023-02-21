#include <stdlib.h>
#include <vector>
#include <iostream>

//global variables
//rows
//columns
//number of non-zero integers

struct CompressedCoordinate {
    size_t row;
    size_t column;
    double value;
};


class COO {
    
    //save position and value of non-zero elements
    public:
        std::vector<int> rowCoord;
        std::vector<int> colCoord;
        std::vector<int> values;
        size_t numRows;
        size_t numCols;
        size_t nnz;
        COO(size_t nonz) { 
            rowCoord.reserve(nonz);
            colCoord.reserve(nonz);
            values.reserve(nonz);
        }            

            

};

/**
 * @brief Get the Value COO object
 * 
 * @param compressedCoord 
 * @param row 
 * @param col 
 * @return double 
 */
double getValueCOO(COO compressedCoord, int row, int col) {
    if (compressedCoord.numRows < row) {
        throw std::invalid_argument("Row is not within range\n");
    }
    if (compressedCoord.numCols < col) {
        throw std::invalid_argument("Column is not within range\n");
    }

    for (int i = 0; i < compressedCoord.values.size(); ++i) {
        if (row == compressedCoord.rowCoord[i] && col == compressedCoord.colCoord[i]) {
            return double(compressedCoord.values[i]);
        }
    }

    return -1;
}


/**
 * @brief 
 * 
 * @param denseMatrix 
 * @param rows 
 * @param cols 
 * @param nonz 
 * @return COO 
 */
COO convertDenseToCOO(std::vector<std::vector<int>> denseMatrix, size_t rows, size_t cols, size_t nonz) {
    size_t nnz_id = 0;
    COO coo(nonz);
    for (int i = 0; i < denseMatrix.size(); ++i) {
        for (int j = 0; j < denseMatrix[i].size(); ++j) {
            if (denseMatrix[i][j] != 0) {
                coo.rowCoord.push_back(i);
                coo.colCoord.push_back(j);
                coo.values.push_back(denseMatrix[i][j]);
                nnz_id++;
            }
        }
    }
    std::cout << "Values = [";
    for (int i = 0; i < coo.colCoord.size(); ++i) {
        std::cout << coo.values[i] << " ";
    }
    std::cout << "]\n";

    std::cout << "Row Coordinates = [";
    for (int i = 0; i < coo.colCoord.size(); ++i) {
        std::cout << coo.rowCoord[i] << " ";
    }
    std::cout << "]\n";

     std::cout << "Column Coordinates = [";
    for (int i = 0; i < coo.colCoord.size(); ++i) {
        std::cout << coo.colCoord[i] << " ";
    }
    std::cout << "]\n";

    return coo;
}

/**
 * @brief Convert a Compressed Coordinate Structure to a Dense Matrix
 * 
 * @param compressedCoord 
 * @return std::vector<std::vector<int>> 
 */
std::vector<std::vector<int>> convertCOOtoDense(COO compressedCoord) {
    std::vector<std::vector<int>> dense;
    int nnz_id = 0;

    for (int i = 0; i < compressedCoord.numRows; ++i) {
        std::vector<int> rowVector;
        for (int j = 0; j < compressedCoord.numCols; ++j) {
            if (i == compressedCoord.rowCoord[nnz_id] && j == compressedCoord.colCoord[nnz_id]) {
                rowVector.push_back(compressedCoord.values[nnz_id]);
                nnz_id++;
            } else {
                rowVector.push_back(0);
            }
        }
        dense.push_back(rowVector);
    }

    for (int i = 0; i < compressedCoord.numRows; ++i) {
        for (int j = 0; j < compressedCoord.numCols; ++j) {
            std::cout << dense[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    return dense;
}

COO multCOO(COO &compressedCoord1, COO compressedCoord2);

COO scalarMultCOO(COO compressedCoord, int scalar) {
    if (scalar == 0) {
        throw std::invalid_argument("Error: cannot zero out matrix\n");
    }

    for (int i = 0; i < compressedCoord.values.size(); ++i) {
        compressedCoord.values.at(i) = compressedCoord.values.at(i) * scalar;
    }

    return compressedCoord;
}

COO addCOO();

int main(int argc, char** argv) {
    std::vector<std::vector<int>> sample = {{2,0,0,2,0}, {3,4,2,5,0}, {5,0,0,8,17}, {0,0,10,16,0}, {0,0,0,0,14}};
    COO coordinate = convertDenseToCOO(sample, sample.size(), sample[0].size(), 3);
    coordinate.numRows = 5;
    coordinate.numCols = 5;
    convertCOOtoDense(coordinate);


}