#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

class Matrix
{
    std::vector<std::vector<double>> matrix;
    int rows, columns;

    void init(int, int);
    void resizeMatrix(int, int);
    void downSize(int, int);
    void upSize(int, int);
    double determinant() const;
    Matrix findSubMatrix(int, int) const;
    friend Matrix multiplyBruteForce(const Matrix&, const Matrix&);
    friend Matrix multiplyStrassen(const Matrix&, const Matrix&);
    friend void splitMatrix(const Matrix&, Matrix&, int, int, int);
    friend void joinMatrix(const Matrix&, Matrix&, int, int, int);
    friend int findNearestPowerOfTwo(int);
    friend int getBiggestSide(const Matrix&, const Matrix&);
    friend bool equalSize(const Matrix&, const Matrix&);

public:

    Matrix() { init(0, 0); }
    Matrix(int rows, int columns) { init(rows, columns); }
    Matrix(const Matrix&);
    Matrix(Matrix&&) noexcept;
    Matrix(std::vector<std::vector<double>>&);
    Matrix& operator=(const Matrix&);
    Matrix& operator=(Matrix&&) noexcept;
    void transpose();
    void inverse();
    void adjugate();
    void identity();
    void pow(int);
    double getDeterminant() const;
    bool isSquare() const { return rows == columns; }
    int getRows() const { return rows; }
    int getColumns() const { return columns; }
    std::vector<double> operator[](int index) const { return matrix[index]; }
    std::vector<double>& operator[](int index) { return matrix[index]; }
    friend std::ostream& operator<<(std::ostream&, const Matrix&);
    friend std::istream& operator>>(std::istream&, Matrix&);
    friend Matrix operator+(const Matrix&, const Matrix&);
    friend Matrix operator-(const Matrix&, const Matrix&);
    friend Matrix operator*(const Matrix&, const Matrix&);
    friend Matrix operator*(const Matrix&, double);
    friend Matrix operator*(double, const Matrix&);
    friend Matrix operator/(const Matrix&, const Matrix&);
    friend Matrix operator/(const Matrix&, double);
    friend Matrix operator/(double, const Matrix&);
    void operator+=(const Matrix &m) { *this = *this + m; }
    void operator-=(const Matrix &m) { *this = *this - m; }
    void operator*=(const Matrix &m) { *this = *this * m; }
    void operator/=(const Matrix &m) { *this = *this / m; }
    void operator*=(double n) { *this = *this * n; }
    void operator/=(double n) { *this = *this / n; }
};

#endif // MATRIX_H
