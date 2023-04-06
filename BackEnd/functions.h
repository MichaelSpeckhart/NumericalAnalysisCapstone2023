using namespace std;

vector<vector<vector<double> > > read_file(char *filename);
vector<vector<double> > sum_matrix(vector<vector<double> > m1, const vector<vector<double> > m2);
vector<vector<double> > sub_matrix(vector<vector<double> > m1, const vector<vector<double> > m2);
vector<vector<double> > generate_random_matrix(const int d1, const int d2, const double min, const double max);
vector<vector<double> > mult_matrix(const vector<vector<double> > m1, const vector<vector<double> > m2);
vector<vector<double> > scale_up(vector<vector<double> > m1, const double s);
vector<vector<double> > scale_down(vector<vector<double> > m1, const double s);
vector<vector<double> > transpose(const vector<vector<double> > m1);
bool save_file(const vector<vector<vector<double> > > matrices, const char *filename);
void swap_row(vector<vector<double> > &m, const int i, const int j);
int forward_elimination(vector<vector<double> > &m);
vector<double> backward_substitution(vector<vector<double> > &m);
vector<vector<double> > gaussian_elimination(vector<vector<double> > &m);

template <typename T>
    class COOMatrix;
template<typename T>
    T get_matrixCOO(COOMatrix<T> compressedCoord, int row, int col);
template<typename T>
    COOMatrix<T> from_vector(std::vector<std::vector<double>> denseMatrix, size_t rows, size_t cols, size_t nonz);
template<typename T>
    std::vector<std::vector<T>> convertCOOtoDense(COOMatrix<T> compressedCoord);
template<typename T>
    COOMatrix<T> multCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2);
template<typename T>
    COOMatrix<T> scalarMultCOO(COOMatrix<T> &compressedCoord, int scalar);
template<typename T>
    COOMatrix<T> scalarDivCOO(COOMatrix<T> &compressedCoord, int scalar);
template<typename T>
    void scalarAddCOO(COOMatrix<T> &compressedCoord, int scalar);
template<typename T>
    COOMatrix<T> addCOO(COOMatrix<T> compressedCoord1, COOMatrix<T> compressedCoord2);

template<typename T>
	class CSCMatrix;
template<typename T>
    T get_matrixCSC(CSCMatrix<T> m1, size_t row, size_t col);
template<typename T>
    CSCMatrix<T> from_vector(vector<vector<T> > &array);
template<typename T>
    void print_matrixCSC(CSCMatrix<T> m1);
template<typename T>
    CSCMatrix<T> add_matrixCSC(CSCMatrix<T> m1, CSCMatrix<T> m2);
template<typename T>
    CSCMatrix<T> transpose_matrixCSC(CSCMatrix<T> m1);
template<typename T>
    CSCMatrix<T> multiply_matrixCSC(CSCMatrix<T> m1, CSCMatrix<T> m2);

template<typename T>
	class CSRMatrix;
template<typename T>
    T get_matrixCSR(CSRMatrix<T> m1, size_t row, size_t col);
template<typename T>
    CSRMatrix<T> from_vector(vector<vector<T> > &array);
template<typename T>
    void print_matrixCSR(CSRMatrix<T> m1);
template<typename T>
    CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2);
template<typename T>
    CSRMatrix<T> transpose_matrixCSR(CSRMatrix<T> m1);
template<typename T>
    CSRMatrix<T> multiply_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2);
vector<vector<double> > load_fileCSR(string fileName);



