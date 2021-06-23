//
//  linalg_matrix_standart.cpp
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 22.12.2020.
//

#include "linalg.h"
#include "linalg_matrix_standart.h"


#pragma mark - Template Functions Defenitions
template<typename T>
INTERNAL void _numericMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    for (unsigned int i = 0; i < commonMatrixArrayRepresentationLength; ++i)
    {
        resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
    }
}


template<typename T>
INTERNAL void _numericMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    for (unsigned int i = 0; i < leftMatrixWidth; ++i)
    {
        for (unsigned int j = 0; j < rightMatrixLength; ++j)
        {
            for (unsigned int k = 0; k < commonShapeNumber; ++k) {
                resultMatrix[i * rightMatrixLength + j] += leftMatrix[i * commonShapeNumber + k] * rightMatrix[k * rightMatrixLength + j];
            }
        }
    }
}


template<typename T>
INTERNAL void _bitMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    for (unsigned int i = 0; i < commonMatrixArrayRepresentationLength; ++i)
    {
        resultMatrix[i] = leftMatrix[i] | rightMatrix[i];
    }
}


template<typename T>
INTERNAL void _bitMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    for (unsigned int i = 0; i < leftMatrixWidth; ++i)
    {
        for (unsigned int j = 0; j < rightMatrixLength; ++j)
        {
            for (unsigned int k = 0; k < commonShapeNumber; ++k) {
                resultMatrix[i * rightMatrixLength + j] ^= leftMatrix[i * commonShapeNumber + k] & rightMatrix[k * rightMatrixLength + j];
            }
        }
    }
}


#pragma mark - External Functions

template<typename T>
EXPORT linMatrix<T> numericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _numericMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> numericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _numericMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> bitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _bitMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}

template<typename T>
EXPORT linMatrix<T> bitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _bitMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


#pragma mark - Templates

template linMatrix<signed char> numericMatricesAddition(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> numericMatricesAddition(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> numericMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> numericMatricesAddition(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> numericMatricesAddition(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> numericMatricesAddition(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);

template linMatrix<signed char> numericMatricesMultiplication(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> numericMatricesMultiplication(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> numericMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> numericMatricesMultiplication(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> numericMatricesMultiplication(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> numericMatricesMultiplication(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);

template linMatrix<int> bitMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);

template linMatrix<int> bitMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
