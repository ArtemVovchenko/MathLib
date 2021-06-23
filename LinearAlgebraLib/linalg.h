//
//  linalg.h
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#ifndef linalg_h
#define linalg_h

#include <iostream>
#include <string>

#pragma mark - Structures

template <typename T>
struct linMatrix {
    unsigned int length;
    unsigned int heigth;
    T* startPtr;

    linMatrix(unsigned int matrixLength, unsigned int matrixHeigth) {
        length = matrixLength;
        heigth = matrixHeigth;
        unsigned long matrixLinearRepresentationLength = length * heigth;
        startPtr = (T*) malloc(sizeof(T) * matrixLinearRepresentationLength);//new T[matrixLinearRepresentationLength];
        memset(startPtr, 0, (size_t)length * (size_t)heigth);
    }
    
    void dealloc() {
        delete [] startPtr;
        startPtr = NULL;
    }
    
    void realloc() {
        if (startPtr != NULL) {
            //TODO: Update Exception Type (Declare Own Exception)
            throw std::invalid_argument("The data in struct is not deallocated");
        }
        startPtr = new T[length * heigth];
        memset(startPtr, 0, (size_t)length * (size_t)heigth);
    }
    
    void print() {
        for (int i = 0; i < heigth; ++i) {
            for (int j = 0; j < length; ++j) {
                std::cout << startPtr[i * length + j] << '\t';
            }
            std::cout << std::endl;
        }
    }
};

#pragma mark - Standart Algorithms

template<typename T>
linMatrix<T> numericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> numericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> bitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> bitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);


#pragma mark - Optimized Algorithms

template<typename T>
linMatrix<T> optimizedNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> optimizedNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> optimizedBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> optimizedBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);


#pragma mark - Parallel Algorithms

template<typename T>
linMatrix<T> parallelNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelOptimizedNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelOptimizedNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelOptimizedBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);

template<typename T>
linMatrix<T> parallelOptimizedBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time);




#pragma mark - Random Generators

template <typename T>
void fillMatrixWithRandomData(linMatrix<T> matrix);

void toBitMatrix(linMatrix<int> matrix);

#endif /* linalg_h */
