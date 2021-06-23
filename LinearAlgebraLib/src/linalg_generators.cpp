//
//  linalg_generators.cpp
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#include "linalg.h"
#include "linalg_types.h"
#include "linalg_attributes.h"
#include "linalg_generators.h"

void _int8Generator(signed char *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_int_distribution<> distribution_func(RANDOM_INT_8_MIN, RANDOM_INT_8_MAX);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (signed char)distribution_func(generator);
    }
}

void _int16Generator(short *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_int_distribution<> distribution_func(RANDOM_INT_16_MIN, RANDOM_INT_16_MAX);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (short)distribution_func(generator);
    }
}

void _intGenerator(int *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_int_distribution<> distribution_func(RANDOM_INT_32_MIN, RANDOM_INT_32_MAX);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (int)distribution_func(generator);
    }
}

void _longGenerator(long *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_int_distribution<> distribution_func(RANDOM_INT_32_MIN, RANDOM_INT_32_MIN);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (long)distribution_func(generator);
    }
}

void _floatGenerator(float *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_real_distribution<> distribution_func(RATIONAL_MIN, RATIONAL_MAX);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (float)distribution_func(generator);
    }
}

void _doubleGenerator(double *adrArrayStartPointer, size_t arraySize) {
    std::random_device seed;
    std::mt19937 generator(seed());
    std::uniform_real_distribution<> distribution_func(RATIONAL_MIN, RATIONAL_MAX);
    
    for (int i = 0; i < arraySize; ++i) {
        adrArrayStartPointer[i] = (double)distribution_func(generator);
    }
}

template <typename T>
void _genericGenerator(T *adrArrayStartPointer, size_t arraySize) {
    
    int elementTypeCode = _getElementType(adrArrayStartPointer);
    switch (elementTypeCode) {
        case SUPPORTED_TYPES::INT_8:
            _int8Generator((signed char*) adrArrayStartPointer, arraySize);
            break;
            
        case SUPPORTED_TYPES::INT_16:
            _int16Generator((short*) adrArrayStartPointer, arraySize);
            break;
            
        case SUPPORTED_TYPES::INT:
            _intGenerator((int*) adrArrayStartPointer, arraySize);
            break;
            
        case SUPPORTED_TYPES::LONG:
            _longGenerator((long*) adrArrayStartPointer, arraySize);
            break;
            
        case SUPPORTED_TYPES::FLOAT:
            _floatGenerator((float*) adrArrayStartPointer, arraySize);
            break;
            
        case SUPPORTED_TYPES::DOUBLE:
            _doubleGenerator((double*) adrArrayStartPointer, arraySize);
            break;
            
        default:
            break;
    }
}

template <typename T>
void fillMatrixWithRandomData(linMatrix<T> matrix) {
    size_t matrixElementsNumber = (size_t)matrix.heigth * (size_t)matrix.length;
    T *linearMatrixArrayStartAddress = matrix.startPtr;
    
    _genericGenerator(linearMatrixArrayStartAddress, matrixElementsNumber);
}

void toBitMatrix(linMatrix<int> matrix) {
    for (unsigned long i = 0; i < matrix.length * matrix.heigth; ++i) {
        matrix.startPtr[i] %= 2 * (matrix.startPtr[i] < 0 ? -1 : 1);
    }
}

template void fillMatrixWithRandomData<signed char>(linMatrix<signed char>);
template void fillMatrixWithRandomData<short>(linMatrix<short>);
template void fillMatrixWithRandomData<int>(linMatrix<int>);
template void fillMatrixWithRandomData<long>(linMatrix<long>);
template void fillMatrixWithRandomData<float>(linMatrix<float>);
template void fillMatrixWithRandomData<double>(linMatrix<double>);

