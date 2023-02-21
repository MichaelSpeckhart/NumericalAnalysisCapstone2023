#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

template<typename T>
	class CSRMatrix{
        struct {
            int col;
            T value;
        }colValStruct;
        public:
        int numRows, numColumns;
        vector<colValStruct*> rowPointersVector;
        //fist value is col_ind and the second value is the value
        vector<colValStruct> columnValueVector;
    };
    //loadfile
    //savefile

    double get_matrixCSR(CSRMatrix<double> m1, int row, int col){
        if(m1.numRows< row){
            //throw error
        }
        if(m1.numColumns < col){
            //throw error
        }
        for(auto it = m1.rowPointersVector.at(row); it != m1.rowPointersVector.at(row+1);it++){
            if(it->col == col){
                return it->colValStruct.value;
            }
        }
        return -1;

    }

    CSRMatrix<double> add_matrixCSR(CSRMatrix<double> m1, CSRMatrix<double> m2){
        if(m1.numRows!= m2.numRows){
            //throw error
        }
        if(m1.numColumns!= m2.numColumns){
            //throw error
        }
        CSRMatrix<double> returnMatrix;
        for(size_t i =0; i < m2.rowPointersVector.size();i++){

        }

        return m1;

    }
    //TODO: Multiply and transpose