#ifndef MATRIX_CPP
#define MATRIX_CPP

#include "matrix.h"
#include <cmath>
#include <bits/stdc++.h>
#include <stdexcept>
#include <iomanip>

void Matrix::init(int rows, int columns)
{
    this -> rows = rows;
    this -> columns = columns;
    matrix.resize(rows, std::vector<double>(columns));
}

void Matrix::resizeMatrix(int rows, int columns)
{
    this -> rows = rows;
    this -> columns = columns;
    matrix.resize(rows);
    for (int i(0); i < rows; i++)
        matrix[i].resize(columns);
}

Matrix::Matrix(const Matrix &m)
{
    init(m.rows, m.columns);
    for(int i(0); i < rows; i++)
        for(int j(0); j < columns; j++)
            matrix[i][j] = m.matrix[i][j];
}

Matrix::Matrix(Matrix &&rhs) noexcept
{
    init(rhs.rows, rhs.columns);
    for(int i(0); i < rows; i++)
        matrix[i] = std::move(rhs.matrix[i]);
}

Matrix::Matrix(std::vector<std::vector<double>> &v)
{
    init(v.size(), v[0].size());
    for (int i(0); i < rows; i++)
    {
        if (v[0].size() != v[i].size())
            throw std::invalid_argument("Invalid matrix format");

        for (int j(0); j < columns; j++)
            matrix[i][j] = v[i][j];
    }
}

Matrix& Matrix::operator=(const Matrix &m)
{
    resizeMatrix(m.matrix.size(), m.matrix[0].size());
    for(int i(0); i < rows; i++)
        for(int j(0); j < columns; j++)
            matrix[i][j] = m.matrix[i][j];
    return *this;
}

Matrix& Matrix::operator=(Matrix &&rhs) noexcept
{
    rows = rhs.rows;
    columns = rhs.columns;
    matrix.swap(rhs.matrix);

    return *this;
}

Matrix operator+(const Matrix &m1, const Matrix &m2)
{
    if (!equalSize(m1, m2))
        throw std::invalid_argument("Only matrices of same dimensions can be added");

    Matrix result(m1);

    for (int i(0); i < m1.rows; i++)
        for (int j(0); j < m1.columns; j++)
            result[i][j] = m1[i][j] + m2[i][j];

    return result;
}

Matrix operator-(const Matrix &m1, const Matrix &m2)
{
    if (!equalSize(m1, m2))
        throw std::invalid_argument("Only matrices of same dimensions can be subtracted");

    Matrix result(m1);

    for (int i(0); i < m1.rows; i++)
        for (int j(0); j < m1.columns; j++)
            result[i][j] = m1[i][j] - m2[i][j];

    return result;
}

int findNearestPowerOfTwo(int n)
{
    int a(log2(n));

    if (std::pow(2, a) == n)
        return n;
    return std::pow(2, a + 1);
}

Matrix operator*(const Matrix &m1, const Matrix &m2)
{
    if (m1.columns != m2.rows)
        throw std::logic_error("Cannot multiply matrices with incompatible dimensions");

    Matrix m1_copy(m1);
    Matrix m2_copy(m2);

    // Size of the biggest side of both matrices
    int dimension(getBiggestSide(m1_copy, m2_copy));

    if (dimension <= 2)
        return multiplyBruteForce(m1, m2);

    // The size isn't a power of two
    if (ceil(log2(dimension)) != floor(log2(dimension)))
        dimension = findNearestPowerOfTwo(dimension);

    // Number of rows and columns that need to be added to make the matrices square
    int m1_rows(dimension - m1_copy.rows);
    int m1_columns(dimension - m1_copy.columns);
    int m2_rows(dimension - m2_copy.rows);
    int m2_columns(dimension - m2_copy.columns);

    m1_copy.upSize(m1_rows, m1_columns);
    m2_copy.upSize(m2_rows, m2_columns);

    Matrix result(multiplyStrassen(m1_copy, m2_copy));
    result.downSize(dimension - m1.rows, dimension - m2.columns);

    return result;
}

Matrix operator*(const Matrix &m, double n)
{
    Matrix copy(m);

    for( int i(0); i < copy.rows; i++)
        for (int j(0); j < copy.columns; j++)
            copy[i][j] = m[i][j] * n;

    return copy;
}

Matrix operator*(double n, const Matrix &m) { return m * n; }

Matrix operator/(const Matrix &m, double n)
{
    if (n == 0)
        throw std::invalid_argument("Cannot divide a matrix with zero");

    Matrix copy(m);

    for (int i(0); i < copy.rows; i++)
        for (int j(0); j < copy.columns; j++)
            copy[i][j] = m[i][j] / n;

    return copy;
}

Matrix operator/(double n, const Matrix &m)
{
    if (n == 0)
        throw std::invalid_argument("Cannot divide a matrix with zero");

    Matrix copy(m);

    for (int i(0); i < copy.rows; i++)
        for (int j(0); j < copy.columns; j++)
            copy[i][j] = n / m[i][j];

    return copy;
}

Matrix operator/(const Matrix &m1, const Matrix &m2)
{
    Matrix copy(m2);
    copy.inverse();

    return m1 * copy;
}

// Upsizes matrix for a given number of rows and columns
void Matrix::upSize(int rows, int columns)
{
    for (int i(0); i < this -> rows; i++)
        for (int k(0); k < columns; k++)
            matrix[i].push_back(0);

    for (int i(0); i < rows; i++)
        matrix.push_back(std::vector<double>(this -> columns + rows, 0));

    this -> rows += rows;
    this -> columns += columns;
}

// Downsizes matrix for a given number of rows and columns
void Matrix::downSize(int rows, int columns)
{
    this -> rows -= rows;
    this -> columns -= columns;
    matrix.resize(this -> rows);

    for (int i(0); i < this -> rows; i++)
        matrix[i].resize(this -> columns);
}

// Gets the biggest size of two matrices
int getBiggestSide(const Matrix &m1, const Matrix &m2)
{
    int max_m1(std::max(m1.rows, m1.columns));
    int max_m2(std::max(m2.rows, m2.columns));

    return std::max(max_m1, max_m2);
}

// Checks if the matrices are of equal size
bool equalSize(const Matrix &m1, const Matrix &m2)
{
    return m1.rows == m2.rows && m1.columns == m2.columns;
}

// Raises the matrix to a power of n
void Matrix::pow(int n)
{
    if (!isSquare())
        throw std::logic_error("Only a square matrix can be raised to a power or have an inverse");

    if (n == 0)
        identity();

    else if (n == -1)
        inverse();

    else if (n > 0)
    {
        Matrix copy(*this);
        for (int i(0); i < n-1; i++)
            *this *= copy;
    }
    else
    {
        pow(-n);
        inverse();
    }
}

// Transposes the matrix
void Matrix::transpose()
{
    std::vector<std::vector<double>> temp(columns, std::vector<double>(rows));

    for (int i(0); i < rows; i++)
        for (int j(0); j < columns; j++)
            temp[j][i] = matrix[i][j];
    matrix.swap(temp);

    int rows_size(rows);
    rows = columns;
    columns = rows_size;
}

// Transforms the matrix to an identity matrix, if it's square
void Matrix::identity()
{
    if (!isSquare())
        throw std::logic_error("Only a square matrix can have an identity matrix");

    for(int i(0); i < rows; i++)
        for(int j(0); j < columns; j++)
            if (i == j)
                matrix[i][j] = 1;
            else matrix[i][j] = 0;
}

// Turns the matrix to its inverse
void Matrix::inverse()
{
    if (!isSquare())
        throw std::logic_error("Only a square matrix can have it's inverse");

    double determinant = getDeterminant();

    if (determinant == 0)
        throw std::logic_error("A matrix cannot have an inverse if it's determinant is zero");

    adjugate();
    *this /= determinant;
}

// Calculates the adjugate matrix needed for finding the inverse
void Matrix::adjugate()
{
    Matrix adj(rows, columns);

    for (int i(0); i < rows; i++)
        for (int j(0); j < columns; j++)
        {
            Matrix sub(findSubMatrix(i, j));
            double cofactor = sub.getDeterminant();
            if ((i+j) % 2 == 1)
                cofactor = -cofactor;
            adj[i][j] = cofactor;
        }
    adj.transpose();
    *this = adj;
}

// For a matrix of size nxn, returns a matrix of size (n-1)x(n-1)
Matrix Matrix::findSubMatrix(int row_num, int col_num) const
{
    Matrix sub_matrix(rows-1, columns-1);
    bool row = false;
    for(int i(0); i < rows; i++)
    {
        if (i == row_num)
            row = true;
        bool column = false;
        for(int j(0); j < columns; j++)
        {
            if (j == col_num)
                column = true;
            if (i != row_num && j != col_num)
            {
                int index_i = i;
                int index_j = j;
                if (row) index_i--;
                if (column) index_j--;
                sub_matrix[index_i][index_j] = matrix[i][j];
            }
        }
    }
    return sub_matrix;
}

// The determinant is calculated using the Laplace expansion
double Matrix::determinant() const
{
    int matrix_size(matrix.size());

    if (matrix_size == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double determinant(0);
    int sign(1);

    for(int j(0); j < matrix_size; j++)
    {
        Matrix sub_matrix = this -> findSubMatrix(0, j);

        determinant += this -> matrix[0][j] * sub_matrix.determinant() * sign;
        sign = -sign;
    }
    return determinant;
}

double Matrix::getDeterminant() const
{
    if (matrix.size() == 1)
        return matrix[0][0];
    return determinant();
}

// Multiplies a matrix in n^3 time
Matrix multiplyBruteForce(const Matrix &m1, const Matrix &m2)
{
    Matrix result(m1.rows, m2.columns);

    for (int i(0); i < m1.rows; i++)
        for (int j(0); j < m2.columns; j++)
        {
            int product(0);
            for (int k(0); k < m1.columns; k++)
                product += m1[i][k] * m2[k][j];
            result[i][j] = product;
        }
    return result;
}

// Copies a specified amount of the matrices elements into the container matrix
// Used for splitting a matrix into 4 sub-matrices of equal size
void splitMatrix(const Matrix &m, Matrix &container, int row, int col, int dim)
{
    for (int i1 = 0, i2 = row; i1 < dim; i1++, i2++)
        for (int j1 = 0, j2 = col; j1 < dim; j1++, j2++)
            container[i1][j1] = m[i2][j2];
}

// Similar to the previous, but used to merge 4 sub-matrices of equal size into a single one
void joinMatrix(const Matrix &m, Matrix &container, int row, int col, int dim)
{
    for (int i1 = 0, i2 = row; i1 < dim; i1++, i2++)
        for (int j1 = 0, j2 = col; j1 < dim; j1++, j2++)
            container[i2][j2] = m[i1][j1];
}

// Strassen's algorithm for matrix multiplication. Time complexity O(n^log7)
Matrix multiplyStrassen(const Matrix &A, const Matrix &B)
{
    if (A.columns <= 2)
        return multiplyBruteForce(A, B);

    //Size of the sub-matrices
    int new_d(A.columns / 2);

    Matrix A11(new_d, new_d);
    Matrix A12(new_d, new_d);
    Matrix A21(new_d, new_d);
    Matrix A22(new_d, new_d);
    Matrix B11(new_d, new_d);
    Matrix B12(new_d, new_d);
    Matrix B21(new_d, new_d);
    Matrix B22(new_d, new_d);

    //Split matrices A and B into 8 sub-matrices
    splitMatrix(A, A11, 0 , 0, new_d);
    splitMatrix(A, A12, 0 , new_d, new_d);
    splitMatrix(A, A21, new_d, 0, new_d);
    splitMatrix(A, A22, new_d, new_d, new_d);
    splitMatrix(B, B11, 0 , 0, new_d);
    splitMatrix(B, B12, 0 , new_d, new_d);
    splitMatrix(B, B21, new_d, 0, new_d);
    splitMatrix(B, B22, new_d, new_d, new_d);

    Matrix P1 = multiplyStrassen(A11+A22, B11+B22);
    Matrix P2 = multiplyStrassen(A22, B21-B11);
    Matrix P3 = multiplyStrassen(A11+A12, B22);
    Matrix P4 = multiplyStrassen(A12-A22, B21+B22);
    Matrix P5 = multiplyStrassen(A11, B12-B22);
    Matrix P6 = multiplyStrassen(A21+A22, B11);
    Matrix P7 = multiplyStrassen(A11-A21, B11+B12);

    Matrix C11 = P1 + P2 - P3 + P4;
    Matrix C12 = P5 + P3;
    Matrix C21 = P6 + P2;
    Matrix C22 = P5 + P1 - P6 - P7;

    //Merge the sub-matrices
    Matrix C(A.columns, A.columns);

    joinMatrix(C11, C, 0 , 0, new_d);
    joinMatrix(C12, C, 0 , new_d, new_d);
    joinMatrix(C21, C, new_d, 0, new_d);
    joinMatrix(C22, C, new_d, new_d, new_d);

    return C;
}

// Outputs the matrix to the console
std::ostream& operator<<(std::ostream &out, const Matrix &m)
{
    for(int i(0); i < m.rows; i++)
    {
        for(int j(0); j < m.columns; j++)
        {
            out << std::fixed
                << std::setprecision(2)
                << std::setw(10)
                << m[i][j];
        }
        out << std::endl;
    }
    return out;
}

// Inputs the matrix from the console into a variable
std::istream& operator>>(std::istream &in, Matrix &m)
{
    if (in.get() != '[')
        throw std::invalid_argument("Invalid matrix format");

    std::vector<std::vector<double>> matrix;
    std::vector<double> column;

    while (in.peek() != ']')
    {
        bool is_negative(false);
        double n;

        if (in.peek() == ' ')
            in.get();

        if (in.peek() == '-')
        {
            in.get();
            is_negative = true;

            if (in.peek() < '0' && in.peek() > '9')
                throw std::invalid_argument("Invalid matrix format");
        }

        if (in.peek() >= '0' && in.peek() <= '9')
        {
            in >> n;

            if (is_negative)
                n = -n;

            column.push_back(n);
        }
        else if (in.peek() == ';')
        {
            in.get();
            matrix.push_back(column);
            column.clear();
        }
        else throw std::invalid_argument("Illegal character");
    }
    in.get();

    matrix.push_back(column);

    for (int i(1); i < matrix.size(); i++)
        if (matrix[0].size() != matrix[i].size())
            throw std::invalid_argument("Invalid matrix format");

    m.rows = matrix.size();
    m.columns = matrix[0].size();
    m.matrix = std::move(matrix);

    return in;
}

#endif // MATRIX_CPP
