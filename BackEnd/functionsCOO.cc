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

namespace COO {

    /**
     * @brief Compressed Coordinate Structure implementation
     * 
     * @tparam T 
     */
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
        COOMatrix<T> from_vector(std::vector<std::vector<T>> denseMatrix) {
            size_t nnz_id = 0;
            COOMatrix<T> coo;
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
            for (size_t i = 0; i < valueSize; ++i) {
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

            for (size_t i = 0; i < compressedCoord.values.size(); ++i) {
                compressedCoord.values.at(i) = compressedCoord.values.at(i) / scalar;
            }

            return compressedCoord;
        }

    /**
     * @brief Takes a COO matrix structure and a scalar, modifies the matrix by that scalar using pass by reference
     * 
     * @tparam T (double,int,float)
     * @param compressedCoord 
     * @param scalar 
     */
    template<typename T>
        void scalar_add_matrixCOO(COOMatrix<T> &compressedCoord, int scalar) {
            if (scalar == 0) {
                return;
            }

            for (size_t i = 0; i < compressedCoord.values.size(); ++i) {
                compressedCoord.values.at(i) = compressedCoord.values.at(i) + scalar;
            }
        }

    /**
     * @brief Takes a COO matrix structure and a scalar, modifies the matrix by that scalar using pass by reference
     * 
     * @tparam T (double,int,float)
     * @param compressedCoord 
     * @param scalar 
     */
    template<typename T>
        void scalar_sub_matrixCOO(COOMatrix<T> &compressedCoord, int scalar) {
            //No modification, return nothing
            if (scalar == 0) {
                return;
            }

            for (size_t i = 0; i < compressedCoord.values.size(); ++i) {
                compressedCoord.values.at(i) = compressedCoord.values.at(i) - scalar;
            }
        }

    /**
     * @brief Add two Compressed Coordinate matrices together, both matrices must be of the same size
     * 
     * @param compressedCoord1 
     * @param compressedCoord2 
     * @return COO 
     */
    template<typename T>
        COOMatrix<T> add_matrixCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2) {
            //Arguments do not pass the required checks for matrix addition, must both be same size
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
                        returnMatrix.values.push_back(compressedCoord1.values.at(i));
                        i++;
                        k++;
                } else if (compressedCoord1.rowCoord.at(i) > compressedCoord2.rowCoord.at(j)
                            || compressedCoord1.rowCoord.at(i) == compressedCoord2.rowCoord.at(j)
                            && compressedCoord1.colCoord.at(i) > compressedCoord2.colCoord.at(j)) {
                        returnMatrix.values.push_back(compressedCoord2.values.at(j));
                        j++;
                        k++;
                } else {
                    returnMatrix.rowCoord.push_back(compressedCoord1.rowCoord.at(i));
                    returnMatrix.colCoord.push_back(compressedCoord1.colCoord.at(i));
                    returnMatrix.values.push_back(compressedCoord1.values.at(i) + compressedCoord2.values.at(j));
                    i++;
                    j++;
                    k++;
                }
            }

            while (i < compressedCoord1.numRows) {
                returnMatrix.values.push_back(compressedCoord1.values.at(i));
                i++;
                k++;
            }
            while (j < compressedCoord2.numRows) {
                returnMatrix.values.push_back(compressedCoord2.values.at(j));
                j++;
                k++;
            }

            return returnMatrix;
        }

        template<typename T>
            COOMatrix<T> multiply_matrixCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2) {
                //Arguments do not pass requirements for matrix multiplication, must be pxn * mxp
                if (compressedCoord1.numRows != compressedCoord2.numCols) {
                    throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.\n");
                }
                COOMatrix<T> returnMatrix;
                returnMatrix.numRows = compressedCoord1.numRows;
                returnMatrix.numCols = compressedCoord2.numCols;

                for (int i = 0; i < compressedCoord1.numRows; ++i) {
                    std::vector<T> row_data(returnMatrix.numCols, 0);

                    for (int j = compressedCoord1.rowCoord[i]; j < compressedCoord1.rowCoord[i+1]; ++j) {
                        int colA = compressedCoord1.colCoord[j];
                        T valA = compressedCoord1.values[j];

                        for (int k = compressedCoord2.rowCoord[colA]; k < compressedCoord2.rowCoord[colA + 1]; ++k) {
                            int colB = compressedCoord2.colCoord[k];
                            T valB = compressedCoord2.values[k];

                            row_data[colB] += valA * valB;
                        }
                    }

                    for (int j = 0; j < returnMatrix.numCols; ++j) {
                        if (row_data[j] != 0) {
                            returnMatrix.values.push_back(row_data[j]);
                            returnMatrix.colCoord.push_back(j);
                        }
                    }

                    returnMatrix.rowCoord.push_back(returnMatrix.values.size());
                }

                return returnMatrix;
            }

        /**
         * @brief 
         * 
         * @tparam T 
         * @param compressedCoord 
         * @return COOMatrix<T> 
         */
        template<typename T>
            COOMatrix<T> transpose_matrixCOO(COOMatrix<T> compressedCoord) {
                COOMatrix<T> returnMatrix;
                returnMatrix.numRows = compressedCoord.numRows;
                returnMatrix.numCols = compressedCoord.numCols;
                returnMatrix.nnz = compressedCoord.nnz;

                std::vector<std::vector<T>> returnValues(returnMatrix.numRows, std::vector<T>(returnMatrix.numCols, 0));

                for (int i = 0; i < compressedCoord.values.size(); ++i) {
                    int row = compressedCoord.rowCoord[i];
                    int col = compressedCoord.colCoord[i];
                    T value = compressedCoord.values[i];
                    returnValues[row][col] = value;
                }

                for (int i = 0; i < compressedCoord.numRows; ++i) {
                    for (int j = 0; j < compressedCoord.numCols; ++j) {
                        if (returnValues[i][j] != 0) {
                            returnMatrix.rowCoord.push_back(i);
                            returnMatrix.colCoord.push_back(j);
                            returnMatrix.values.push_back(returnValues[i][j]);
                        }
                    }
                }

                return returnMatrix;
            }
}

