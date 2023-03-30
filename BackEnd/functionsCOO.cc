#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


//global variables
//rows
//columns
//number of non-zero integers

//Rather than create three separate vectors, could store one vector with a struct that consists of the
//coordinates and the values all as one.

template <typename T>
    class COOMatrix {
        struct CompressedCoordinate {
            size_t row;
            size_t column;
            double value;
        };

    //save position and value of non-zero elements
    public:
        std::vector<int> rowCoord;
        std::vector<int> colCoord;
        std::vector<int> values;
        size_t numRows;
        size_t numCols;
        size_t nnz;
        COOMatrix(size_t nonz) { 
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
template<typename T>
    T get_matrixCOO(COOMatrix<T> compressedCoord, int row, int col) {
        if (compressedCoord.numRows < row) {
            throw std::invalid_argument("Row is not within range\n");
        }
        if (compressedCoord.numCols < col) {
            throw std::invalid_argument("Column is not within range\n");
        }

        for (size_t i = 0; i < compressedCoord.values.size(); ++i) {
            if (row == compressedCoord.rowCoord[i] && col == compressedCoord.colCoord[i]) {
                return compressedCoord.values.at(i);
            }
        }

        return 0;
    }


/**
 * @brief Convert a dense matrix to a sparse compressed column object. Any non-zero values are added 
 *        along with their coordinates
 * 
 * @param denseMatrix 
 * @param rows 
 * @param cols 
 * @param nonz 
 * @return COOMatrix 
 */
template<typename T>
    COOMatrix<T> from_vector(std::vector<std::vector<double>> denseMatrix, size_t rows, size_t cols, size_t nonz) {
        size_t nnz_id = 0;
        COOMatrix coo(nonz);
        for (size_t i = 0; i < denseMatrix.size(); ++i) {
            for (size_t j = 0; j < denseMatrix[i].size(); ++j) {
                if (denseMatrix[i][j] != 0) {
                    coo.rowCoord.push_back(i);
                    coo.colCoord.push_back(j);
                    coo.values.push_back(denseMatrix[i][j]);
                    nnz_id++;
                }
            }
        }

        return coo;
    }

/**
 * @brief Convert a Compressed Coordinate Structure to a Dense Matrix
 * 
 * @param compressedCoord 
 * @return std::vector<std::vector<int>> 
 */
template<typename T>
    std::vector<std::vector<T>> convertCOOtoDense(COOMatrix<T> compressedCoord) {
        std::vector<std::vector<int>> dense;
        int nnz_id = 0;

        for (size_t i = 0; i < compressedCoord.numRows; ++i) {
            std::vector<int> rowVector;
            for (size_t j = 0; j < compressedCoord.numCols; ++j) {
                if (i == compressedCoord.rowCoord[nnz_id] && j == compressedCoord.colCoord[nnz_id]) {
                    rowVector.push_back(compressedCoord.values[nnz_id]);
                    nnz_id++;
                } else {
                    rowVector.push_back(0);
                }
            }
            dense.push_back(rowVector);
        }

        return dense;
    }

//@Todo
template<typename T>
    COOMatrix<T> multCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2) {
        if (compressedCoord1.numRows != compressedCoord2.numCols) {
            throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.\n");
        }
        COOMatrix<T> returnMatrix;
        returnMatrix.numRows = compressedCoord1.numRows;
        returnMatrix.numCols = compressedCoord2.numCols;
        int k = 0;
        for (size_t i = 0; i < compressedCoord1.numRows; ++i) {
            for (size_t j = 0; j < compressedCoord2.numCols; ++j) {
                T sum = 0;

            }
        }
    }

/**
 * @brief Scale up a matrix by a certain scalar, throws error if scalar is zero
 * 
 * @param compressedCoord 
 * @param scalar 
 * @return COO 
 */
template<typename T>
    COOMatrix<T> scalarMultCOO(COOMatrix<T> &compressedCoord, int scalar) {
        if (scalar == 0) {
            throw std::invalid_argument("Error: cannot zero out matrix\n");
        }

        const int valueSize = static_cast<int>(compressedCoord.values.size());
        for (int i = 0; i < valueSize; ++i) {
            compressedCoord.values.at(i) = compressedCoord.values.at(i) * scalar;
        }

        return compressedCoord;
    }

/**
 * @brief Scale down a matrix by a certain scalar, throws error if scalar is zero
 * 
 * @param compressedCoord 
 * @param scalar 
 * @return COO 
 */
template<typename T>
    COOMatrix<T> scalarDivCOO(COOMatrix<T> &compressedCoord, int scalar) {
        if (scalar == 0) {
            throw std::invalid_argument("Error: cannot divide by zero\n");
        }

        for (int i = 0; i < compressedCoord.values.size(); ++i) {
            compressedCoord.values.at(i) = compressedCoord.values.at(i) / scalar;
        }

        return compressedCoord;
    }

/**
 * @brief 
 * 
 * @param compressedCoord 
 * @param scalar 
 */
template<typename T>
    void scalarAddCOO(COOMatrix<T> &compressedCoord, int scalar) {

    }

/**
 * @brief Add two Compressed Coordinate matrices together, both matrices must be of the same size
 * 
 * @param compressedCoord1 
 * @param compressedCoord2 
 * @return COO 
 */
template<typename T>
    COOMatrix<T> addCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2) {
        if (compressedCoord1.numRows != compressedCoord2.numRows) {
            throw std::invalid_argument("Error: Matrices do not have the same number of rows\n");
        }
        if (compressedCoord1.numCols != compressedCoord2.numCols) {
            throw std::invalid_argument("Error: Matrices do not have the same number of columns\n");
        }

        COOMatrix<T> returnMatrix;
        returnMatrix.numRows = compressedCoord1.numRows;
        returnMatrix.numCols = compressedCoord1.numCols;
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        while (i < compressedCoord1.numRows && j < compressedCoord2.numRows) {
            if (compressedCoord1.rowCoord.at(i) < compressedCoord2.rowCoord.at(j) 
                || compressedCoord1.rowCoord.at(i) == compressedCoord2.rowCoord.at(j) 
                && compressedCoord1.colCoord.at(i) < compressedCoord2.colCoord.at(j)) {
                    returnMatrix.values.at(k) = compressedCoord1.values.at(i);
                    i++;
                    k++;
            } else if (compressedCoord1.rowCoord.at(i) > compressedCoord2.rowCoord.at(j)
                        || compressedCoord1.rowCoord.at(i) == compressedCoord2.rowCoord.at(j)
                        && compressedCoord1.colCoord.at(i) > compressedCoord2.colCoord.at(j)) {
                    returnMatrix.values.at(k) = compressedCoord2.values.at(j);
                    j++;
                    k++;
            } else {
                returnMatrix.rowCoord.at(k) = compressedCoord1.rowCoord.at(i);
                returnMatrix.colCoord.at(k) = compressedCoord1.colCoord.at(i);
                returnMatrix.values.at(k) = compressedCoord1.values.at(i) + compressedCoord2.values.at(j);
                i++;
                j++;
                k++;
            }
        }

        while (i < compressedCoord1.numRows) {
            returnMatrix.values.at(k) = compressedCoord1.values.at(i);
            i++;
            k++;
        }
        while (j < compressedCoord2.numRows) {
            returnMatrix.values.at(k) = compressedCoord2.values.at(j);
            j++;
            k++;
        }

        return returnMatrix;
    }

// int main(int argc, char** argv) {
//     std::vector<std::vector<int>> sample = {{2,0,0,2,0}, {3,4,2,5,0}, {5,0,0,8,17}, {0,0,10,16,0}, {0,0,0,0,14}};
//     COO coordinate = convertDenseToCOO(sample, sample.size(), sample[0].size(), 3);
//     coordinate.numRows = 5;
//     coordinate.numCols = 5;
//     convertCOOtoDense(coordinate);


// }
