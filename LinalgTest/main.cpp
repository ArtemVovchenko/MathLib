//
//  main.cpp
//  TestCPP
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#include <iostream>
#include "linalg.h"

template<typename T>
void printAdditionNumericTask(linMatrix<T> left, linMatrix<T> right) {
    double time;
    std::cout << "----------------- Addition Numeric -----------------" << std::endl;
    
    std::cout << "Matrices shape:\t("<< left.length << ", " << left.heigth << ")" << std::endl;
    
    linMatrix<int> res = numericMatricesAddition(left, right, &time);
    std::cout << "Stdandart algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = optimizedNumericMatricesAddition(left, right, &time);
    std::cout << "Optimized (AVX) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelNumericMatricesAddition(left, right, &time);
    std::cout << "Asynchronous (OMP) algorithm time:\t" << time << " sec." << std::endl;
    res.dealloc();

    res = parallelOptimizedNumericMatricesAddition(left, right, &time);
    std::cout << "Asynchronous optimized (OMP + AVX) algorithm time:\t" << time << " sec." << std::endl;
    res.dealloc();
    
    std::cout << "----------------------------------------------------" << std::endl << std::endl;
}


template<typename T>
void printMultiplicationNumericTask(linMatrix<T> left, linMatrix<T> right) {
    double time;
    std::cout << "-------------- Multiplication Numeric --------------" << std::endl;
    
    std::cout << "Matrices shape:\t("<< left.length << ", " << left.heigth << ")" << std::endl;
    
    linMatrix<int> res = numericMatricesMultiplication(left, right, &time);
    std::cout << "Stdandart algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = optimizedNumericMatricesMultiplication(left, right, &time);
    std::cout << "Optimized algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelNumericMatricesMultiplication(left, right, &time);
    std::cout << "Asynchronous (OMP) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelOptimizedNumericMatricesMultiplication(left, right, &time);
    std::cout << "Asynchronous optimized algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();
    
    std::cout << "----------------------------------------------------" << std::endl << std::endl;
}


template<typename T>
void printAdditionBitTask(linMatrix<T> left, linMatrix<T> right) {
    
    toBitMatrix(left);
    toBitMatrix(right);
    
    double time;
    std::cout << "------------------- Addition Bit -------------------" << std::endl;
    
    std::cout << "Matrices shape:\t("<< left.length << ", " << left.heigth << ")" << std::endl;
    
    linMatrix<int> res = bitMatricesAddition(left, right, &time);
    std::cout << "Stdandart algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = optimizedBitMatricesAddition(left, right, &time);
    std::cout << "Optimized (AVX) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelBitMatricesAddition(left, right, &time);
    std::cout << "Asynchronous (OMP) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelOptimizedBitMatricesAddition(left, right, &time);
    std::cout << "Asynchronous optimized (OMP + AVX) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();
    
    std::cout << "----------------------------------------------------" << std::endl << std::endl;
}


template<typename T>
void printMultiplicationBitTask(linMatrix<T> left, linMatrix<T> right) {
    
    toBitMatrix(left);
    toBitMatrix(right);
    
    double time;
    std::cout << "---------------- Multiplication Bit ----------------" << std::endl;
    
    std::cout << "Matrices shape:\t("<< left.length << ", " << left.heigth << ")" << std::endl;
    
    linMatrix<int> res = bitMatricesMultiplication(left, right, &time);
    std::cout << "Stdandart algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = optimizedBitMatricesMultiplication(left, right, &time);
    std::cout << "Optimized algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelBitMatricesMultiplication(left, right, &time);
    std::cout << "Asynchronous (OMP) algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();

    res = parallelOptimizedBitMatricesMultiplication(left, right, &time);
    std::cout << "Asynchronous optimized algorithm time:\t"<< time << " sec." << std::endl;
    res.dealloc();
    
    std::cout << "----------------------------------------------------" << std::endl << std::endl;
}


int main(int argc, const char * argv[]) {
    linMatrix<int> a = linMatrix<int>(1000, 1000);
    fillMatrixWithRandomData(a);
    linMatrix<int> b = linMatrix<int>(1000, 1000);
    fillMatrixWithRandomData(b);

    printAdditionNumericTask(a, b);
    printMultiplicationNumericTask(a, b);
    printAdditionBitTask(a, b);
    printMultiplicationBitTask(a, b);

    a.dealloc();
    b.dealloc();
    
    return 0;
}
