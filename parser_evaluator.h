#ifndef PARSER_EVALUATOR_H
#define PARSER_EVALUATOR_H

#include "matrix.h"
#include <stack>

enum Symbol {openBracket, closedBracket, unaryOperation, binaryOperation, number, matrix};

int priority(char);
void performBinary(std::stack<Matrix>&, std::stack<double>&, std::stack<char>&, std::stack<Symbol>&);
void performUnary(std::stack<Matrix>&, int);
void parse_evaluate();

#endif // PARSER_EVALUATOR_H
