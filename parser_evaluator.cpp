#ifndef PARSER_EVALUATOR_CPP
#define PARSER_EVALUATOR_CPP

#include "parser_evaluator.h"
#include <iostream>
#include <cmath>

int priority(char operation)
{
    if (operation == '+' || operation == '-')
        return 1;
    if (operation == '*' || operation == '/')
        return 2;
    if (operation == '^')
        return 3;
    return 0;
}

void performBinary(std::stack<Matrix> &matrices,
                   std::stack<double> &numbers,
                   std::stack<char> &chars,
                   std::stack<Symbol> &order)
{
    if (chars.empty() || priority(chars.top()) == 0 || order.size() < 2)
        throw std::invalid_argument("Invalid expression format");

    char operation(chars.top());
    chars.pop();

    // Operands
    Symbol operand2(order.top());
    order.pop();
    Symbol operand1(order.top());
    order.pop();

    if (operand1 == matrix && operand2 == matrix && matrices.size() < 2)
        throw std::invalid_argument("Invalid expression format");

    else if (operand1 == number && operand2 == number && numbers.size() < 2)
        throw std::invalid_argument("Invalid expression format");

    else if (operand1 != operand2 && (matrices.empty() || numbers.empty()))
        throw std::invalid_argument("Invalid expression format");

    // The operation is either + or -
    if (priority(operation) == 1)
    {
        if (operand1 != operand2)
            throw std::invalid_argument("Illegal operation");

        // Both operands are matrices
        if (operand1 == matrix)
        {
            Matrix m2(matrices.top());
            matrices.pop();
            Matrix m1(matrices.top());
            matrices.pop();

            if (operation == '+')
                matrices.push(m1 + m2);
            else matrices.push(m1 - m2);
            order.push(matrix);
        }
        // Both operands are numbers
        else
        {
            double n2(numbers.top());
            numbers.pop();
            double n1(numbers.top());
            numbers.pop();

            if (operation == '+')
                numbers.push(n1 + n2);
            else numbers.push(n1 - n2);
            order.push(number);
        }
    }
    // The operation is either * or /
    else if (priority(operation) == 2)
    {
        if (operand1 == matrix)
        {
            Matrix m2(matrices.top());
            matrices.pop();

            // Both operands are matrices
            if (operand2 == matrix)
            {
                Matrix m1(matrices.top());
                matrices.pop();

                if (operation == '*')
                    matrices.push(m1 * m2);
                else matrices.push(m1 / m2);
            }
            // First operand is a matrix, second is a number
            else
            {
                double n1(numbers.top());
                numbers.pop();

                if (operation == '*')
                    matrices.push(n1 * m2);
                else matrices.push(n1 / m2);
            }
            order.push(matrix);
        }
        else
        {
            double n2(numbers.top());
            numbers.pop();

            // Second operand is a matrix, first is a number
            if (operand2 == matrix)
            {
                Matrix m1(matrices.top());
                matrices.pop();

                if (operation == '*')
                    matrices.push(m1 * n2);
                else matrices.push(m1 / n2);
                order.push(matrix);
            }
            // Both operands are numbers
            else
            {
                double n1(numbers.top());
                numbers.pop();

                if (operation == '*')
                    numbers.push(n1 * n2);
                else numbers.push(n1 / n2);
                order.push(number);
            }
        }
    }
}

void performUnary(std::stack<Matrix> &matrices, int n)
{
    if (matrices.empty())
        throw std::invalid_argument("Invalid expression format");

    Matrix m(matrices.top());
    matrices.pop();
    m.pow(n);
    matrices.push(m);
}

void parse_evaluate()
{
    std::stack<char> chars;
    std::stack<Matrix> matrices;
    std::stack<double> numbers;
    std::stack<Symbol> order;
    Symbol previous(openBracket);

    if (std::cin.peek() == '\n')
    {
        std::cin.get();
        return;
    }

    try
    {
        while (std::cin.peek() != '\n')
        {
            while (std::cin.peek() == ' ')
                std::cin.get();

            if (std::cin.peek() == '(')
            {
                if (!(previous == binaryOperation || previous == openBracket))
                    throw std::invalid_argument("Only an operation can come before an open bracket");

                chars.push('(');
                previous = openBracket;
                std::cin.get();
            }
            else if (std::cin.peek() == ')')
            {
                if (previous == binaryOperation)
                    throw std::invalid_argument("A binary operation cannot come after a closed bracket");

                // Performs all operations inside the brackets
                while (!chars.empty() && chars.top() != '(')
                {
                    performBinary(matrices, numbers, chars, order);
                    if (chars.empty())
                        throw std::invalid_argument("Missing open bracket");
                }
                chars.pop();
                previous = closedBracket;
                std::cin.get();
            }
            else if (std::cin.peek() >= '0' && std::cin.peek() <= '9')
            {
                if (!(previous == binaryOperation || previous == unaryOperation || previous == openBracket))
                    throw std::invalid_argument("A number can only come after an operation or open bracket");

                double n;
                std::cin >> n;

                numbers.push(n);
                order.push(number);
                previous = number;
            }
            else if (std::cin.peek() == '[')
            {
                if (!(previous == openBracket || previous == binaryOperation))
                    throw std::invalid_argument("A matrix can only come after an open bracket or operation");

                Matrix m;
                std::cin >> m;

                matrices.push(m);
                order.push(matrix);
                previous = matrix;
            }
            else if (std::cin.peek() == 'I')
            {
                std::cin.get();

                if (!(std::cin.peek() >= '0' && std::cin.peek() <= '9'))
                    throw std::invalid_argument("Invalid expression format");

                int n;
                std::cin >> n;

                Matrix m(n, n);
                m.identity();

                matrices.push(m);
                order.push(matrix);
                previous = matrix;
            }
            else
            {
                char character(std::cin.peek());

                if (priority(character) == 0)
                    throw std::invalid_argument("Illegal character");

                // The operation is either +, -, * or /
                else if (priority(character) <= 2)
                {
                    std::cin.get();

                    // Negative number
                    if (character == '-' && previous == openBracket)

                        if (std::cin.peek() == '[' || (std::cin.peek() >= '0' && std::cin.peek() <= '9'))
                        {
                            numbers.push(-1);
                            order.push(number);
                            chars.push('*');
                            previous = binaryOperation;
                            continue;
                        }
                        else throw std::invalid_argument("Invalid expression format");

                    else if (previous == openBracket || previous == binaryOperation)
                        throw std::invalid_argument("A binary operation cannot come after an open bracket or another operation");

                    // Performs previous operations with a higher or equal priority
                    while (!chars.empty() && priority(chars.top()) >= priority(character))
                        performBinary(matrices, numbers, chars, order);

                    chars.push(character);
                    previous = binaryOperation;
                }
                // The operation is an exponent
                else if (priority(character) == 3)
                {
                    if (previous == openBracket || previous == binaryOperation || previous == unaryOperation)
                        throw std::invalid_argument("A unary operation cannot come after an open bracket or another operation");

                    std::cin.get();
                    char next(std::cin.peek());

                    // Transpose the matrix
                    if (next == 'T')
                    {
                        if (order.top() == number)
                            throw std::invalid_argument("Cannot transpose a number");

                        Matrix m(matrices.top());
                        matrices.pop();
                        m.transpose();
                        matrices.push(m);
                        std::cin.get();
                        continue;
                    }

                    char is_negative(false);

                    // The number is negative
                    if (next == '-')
                    {
                        is_negative = true;
                        std::cin.get();
                        next = std::cin.peek();
                    }

                    // Multiply the matrix to the power of n
                    if (next >= '0' && next <= '9')
                    {
                        int n;
                        std::cin >> n;

                        if (is_negative)
                            n = -n;

                        // Calculating the matrix exponential
                        if (order.top() == matrix)
                            performUnary(matrices, n);

                        // Calculating the number exponential
                        else
                        {
                            double number(numbers.top());
                            numbers.pop();
                            numbers.push(std::pow(number, n));
                        }
                        continue;
                    }
                    else throw std::invalid_argument("Illegal character");
                }
                else throw std::invalid_argument("Illegal character");
            }
        }

        while (!chars.empty())
        {
            if (chars.top() == '(')
                throw std::invalid_argument("Missing closing bracket");

            performBinary(matrices, numbers, chars, order);
        }

        // Outputs the final result

        if (matrices.size() == 1 && numbers.empty())
            std::cout << std::endl << matrices.top() << std::endl;

        else if (matrices.empty() && numbers.size() == 1)
            std::cout << std::endl << numbers.top() << std::endl << std::endl;

        else throw std::logic_error("Expression not balanced");
    }

    // Catch any errors
    catch (const std::invalid_argument &e)
    {
        std::cerr << "\nException: " << e.what() << std::endl << std::endl;
    }

    catch (const std::logic_error &e)
    {
        std::cerr << "\nException: " << e.what() << std::endl << std::endl;
    }

    // Unknown error occurred
    catch (...)
    {
        std::cerr << "\nException: An unknown error occurred" << std::endl << std::endl;
    }
    std::cin.ignore(1000, '\n');
}

#endif // PARSER_EVALUATOR_CPP
