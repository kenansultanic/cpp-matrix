#ifndef TEST_H
#define TEST_H

#include <vector>
#include "matrix.h"

using matrix_elements = std::vector<std::vector<double>>;

void test1()
{
    matrix_elements me1 {{1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9}};
    matrix_elements me2 {{4, 7}, {5, -4}, {9, 0}, {-8, 3}};
    matrix_elements me3 {{13, 4}, {0, -1}, {0, 0}, {1, 1}};

    Matrix matrix1(me1);
    Matrix matrix2(me2);
    Matrix matrix3(me3);

    Matrix matrix4(matrix1*matrix2);

    matrix4 *= 1.5;
    matrix4 -= matrix3;

    std::cout << matrix4 << std::endl;
}

void test2()
{
    matrix_elements me {{1, 3, 5, 9}, {1, 3, 1, 7}, {4, 3, 9, 7}, {5, 2, 0, 9}};

    Matrix matrix(me);

    std::cout << matrix << std::endl;

    std::cout << "The determinant is equal to: \n\n";

    std::cout << matrix.getDeterminant() << std::endl << std::endl;

    std::cout << "If we transpose the matrix, we get: \n\n";

    matrix.transpose();

    std::cout << matrix << std::endl;

    std::cout << "If we power it by two we get: \n\n";

    matrix.pow(2);

    std::cout << matrix << std::endl;

    std::cout << "Its inverse matrix is: \n\n";

    matrix.inverse();

    std::cout << matrix << std::endl;

}

#endif // TEST_H
