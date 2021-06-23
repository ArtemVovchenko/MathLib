//
//  linalg_matrix_optimized.cpp
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 22.12.2020.
//

#include "linalg.h"
#include "linalg_matrix_optimized.h"


#pragma mark - Internal Functions
template<typename T>
INTERNAL void _optimizedNumericMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    
    unsigned int tail, avxVectorLength;
    __m256i *y_vector_i, *z_vector_i, *out_vector_i;
    __m256  *y_vector, *z_vector, *out_vector;
    __m256d *y_vector_d, *z_vector_d, *out_vector_d;
    
    
    switch (_getElementType(resultMatrix)) {
    case SUPPORTED_TYPES::INT_8:
        tail = commonMatrixArrayRepresentationLength % 32;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m256i*)leftMatrix;
        z_vector_i = (__m256i*)rightMatrix;
        out_vector_i = (__m256i*)resultMatrix;
            
        for (unsigned int i = 0; i < avxVectorLength / 32; ++i) {
            out_vector_i[i] = _mm256_add_epi8(y_vector_i[i], z_vector_i[i]);
        }
            
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    case SUPPORTED_TYPES::INT_16:
        tail = commonMatrixArrayRepresentationLength % 16;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m256i*)leftMatrix;
        z_vector_i = (__m256i*)rightMatrix;
        out_vector_i = (__m256i*)resultMatrix;
            
        for (unsigned int i = 0; i < avxVectorLength / 16; ++i) {
            out_vector_i[i] = _mm256_add_epi16(y_vector_i[i], z_vector_i[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    case SUPPORTED_TYPES::INT:
        tail = commonMatrixArrayRepresentationLength % 8;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m256i*)leftMatrix;
        z_vector_i = (__m256i*)rightMatrix;
        out_vector_i = (__m256i*)resultMatrix;
        
        for (unsigned int i = 0; i < avxVectorLength / 8; ++i) {
            out_vector_i[i] = _mm256_add_epi32(y_vector_i[i], z_vector_i[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    case SUPPORTED_TYPES::LONG:
        tail = commonMatrixArrayRepresentationLength % 4;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m256i*)leftMatrix;
        z_vector_i = (__m256i*)rightMatrix;
        out_vector_i = (__m256i*)resultMatrix;
            
        for (unsigned int i = 0; i < avxVectorLength / 4; ++i) {
            out_vector_i[i] = _mm256_add_epi64(y_vector_i[i], z_vector_i[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    case SUPPORTED_TYPES::FLOAT:
        tail = commonMatrixArrayRepresentationLength % 8;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector = (__m256*)leftMatrix;
        z_vector = (__m256*)rightMatrix;
        out_vector = (__m256*)resultMatrix;
        
        for (unsigned int i = 0; i < avxVectorLength / 8; ++i) {
            out_vector[i] = _mm256_add_ps(y_vector[i], z_vector[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    case SUPPORTED_TYPES::DOUBLE:
        tail = commonMatrixArrayRepresentationLength % 4;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_d = (__m256d*)leftMatrix;
        z_vector_d = (__m256d*)rightMatrix;
        out_vector_d = (__m256d*)resultMatrix;
        for (unsigned int i = 0; i < avxVectorLength / 4; ++i) {
            out_vector_d[i] = _mm256_add_pd(y_vector_d[i], z_vector_d[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
        }
        break;

    default:
        break;
    }
}


template<typename T>
INTERNAL void _optimizedNumericMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    for (int i = 0; i < leftMatrixWidth; ++i)
    {
        T* resultMatrixPointer = resultMatrix + i * rightMatrixLength;
        for (unsigned int j = 0; j < commonShapeNumber; ++j)
        {
            const T* rigthMatrixPointer = rightMatrix + j * rightMatrixLength;
            T leftMatrixPointer = leftMatrix[i * commonShapeNumber + j];
            for (unsigned int k = 0; k < rightMatrixLength; ++k) {
                resultMatrixPointer[k] += leftMatrixPointer * rigthMatrixPointer[k];
            }
        }
    }
}


template<typename T>
INTERNAL void _optimizedBitMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix, unsigned int commonMatrixArrayRepresentationLength) {
    unsigned int tail, avxVectorLength;
    __m128i *y_vector_i, *z_vectror_i, *out_vector_i;
    
    switch (_getElementType(resultMatrix)) {
    case SUPPORTED_TYPES::INT:
        tail = commonMatrixArrayRepresentationLength % 4;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m128i*)leftMatrix;
        z_vectror_i = (__m128i*)rightMatrix;
        out_vector_i = (__m128i*)resultMatrix;
        
        for (unsigned int i = 0; i < avxVectorLength / 4; ++i) {
            out_vector_i[i] = _mm_or_si128(y_vector_i[i], z_vectror_i[i]);
        }
        for (unsigned int i = avxVectorLength; i < commonMatrixArrayRepresentationLength; ++i) {
            resultMatrix[i] = leftMatrix[i] | rightMatrix[i];
        }
        break;

    default:
        break;
    }
}


template<typename T>
INTERNAL void _optimizedBitMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    for (unsigned int i = 0; i < leftMatrixWidth; ++i)
    {
        T* resultMatrixPointer = resultMatrix + i * rightMatrixLength;
        for (unsigned int j = 0; j < commonShapeNumber; ++j)
        {
            const T* rigthMatrixPointer = rightMatrix + j * rightMatrixLength;
            T leftMatrixPointer = leftMatrix[i * commonShapeNumber + j];
            if (leftMatrixPointer) {
                for (unsigned int k = 0; k < rightMatrixLength; ++k) {
                    resultMatrixPointer[k] ^= rigthMatrixPointer[k];
                }
            }
        }
    }
}



#pragma mark - External Functions

template<typename T>
EXPORT linMatrix<T> optimizedNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _optimizedNumericMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> optimizedNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _optimizedNumericMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> optimizedBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _optimizedBitMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}

template<typename T>
EXPORT linMatrix<T> optimizedBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _optimizedBitMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


#pragma mark - Templates

template linMatrix<signed char> optimizedNumericMatricesAddition(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> optimizedNumericMatricesAddition(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> optimizedNumericMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> optimizedNumericMatricesAddition(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> optimizedNumericMatricesAddition(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> optimizedNumericMatricesAddition(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);

template linMatrix<signed char> optimizedNumericMatricesMultiplication(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> optimizedNumericMatricesMultiplication(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> optimizedNumericMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> optimizedNumericMatricesMultiplication(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> optimizedNumericMatricesMultiplication(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> optimizedNumericMatricesMultiplication(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);

template linMatrix<int> optimizedBitMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);

template linMatrix<int> optimizedBitMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
