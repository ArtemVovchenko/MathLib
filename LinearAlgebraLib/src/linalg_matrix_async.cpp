//
//  linalg_matrix_async.cpp
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 22.12.2020.
//

#include "linalg.h"

#include "linalg.h"
#include "linalg_matrix_async.h"


#pragma mark - Internal Functions

template<typename T>
void _parallelNumericMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    #pragma omp parallel for
    for (unsigned int i = 0; i < commonMatrixArrayRepresentationLength; ++i)
    {
        resultMatrix[i] = leftMatrix[i] + rightMatrix[i];
    }
}

template<typename T>
void _parallelNumericMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    #pragma omp parallel for
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
void _parallelOptimizedNumericMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {

    unsigned int tail;
    unsigned int avxVectorLength;

    __m256i* y_vector_i;
    __m256i* z_vector_i;
    __m256i* out_vector_i;

    __m256* y_vector;
    __m256* z_vector;
    __m256* out_vector;

    __m256d* y_vector_d;
    __m256d* z_vector_d;
    __m256d* out_vector_d;

    switch (_getElementType(resultMatrix)) {
    case SUPPORTED_TYPES::INT_8:
        tail = commonMatrixArrayRepresentationLength % 32;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector_i = (__m256i*)leftMatrix;
        z_vector_i = (__m256i*)rightMatrix;
        out_vector_i = (__m256i*)resultMatrix;
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
void _parallelOptimizedNumericMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    #pragma omp parallel for
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
void _parallelBitMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    #pragma omp parallel for
    for (unsigned int i = 0; i < commonMatrixArrayRepresentationLength; ++i)
    {
        resultMatrix[i] = leftMatrix[i] | rightMatrix[i];
    }
}

template<typename T>
void _parallelBitMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    #pragma omp parallel for
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

template<typename T>
void _parallelOptimizedBitMatricesAddition(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int commonMatrixArrayRepresentationLength) {
    switch (_getElementType(resultMatrix)) {

        unsigned int tail;
        unsigned int avxVectorLength;
        __m128i* y_vector;
        __m128i* z_vector;
        __m128i* out_vector;

    case SUPPORTED_TYPES::INT:
        tail = commonMatrixArrayRepresentationLength % 4;
        avxVectorLength = commonMatrixArrayRepresentationLength - tail;
        y_vector = (__m128i*)leftMatrix;
        z_vector = (__m128i*)rightMatrix;
        out_vector = (__m128i*)resultMatrix;
        #pragma omp parallel for
        for (unsigned int i = 0; i < avxVectorLength / 4; ++i) {
            out_vector[i] = _mm_or_si128(y_vector[i], z_vector[i]);
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
void _parallelOptimizedBitMatricesMultiplication(const T* leftMatrix, const T* rightMatrix, T* resultMatrix,
    unsigned int leftMatrixWidth, unsigned int rightMatrixLength, unsigned int commonShapeNumber) {
    #pragma omp parallel for
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
EXPORT linMatrix<T> parallelNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelNumericMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> parallelOptimizedNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelOptimizedNumericMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> parallelBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelBitMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}

template<typename T>
EXPORT linMatrix<T> parallelBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelBitMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> parallelOptimizedNumericMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelOptimizedNumericMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> parallelNumericMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelNumericMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}


template<typename T>
EXPORT linMatrix<T> parallelOptimizedBitMatricesAddition(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if ((leftMatrix.length != rightMatrix.length) || (leftMatrix.heigth != rightMatrix.heigth)) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    unsigned int resultMatrixArrayLength = resultMatrix.heigth * resultMatrix.heigth;
    
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelOptimizedBitMatricesAddition(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, resultMatrixArrayLength);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}

template<typename T>
EXPORT linMatrix<T> parallelOptimizedBitMatricesMultiplication(const linMatrix<T> leftMatrix, const linMatrix<T> rightMatrix, double *time) {
    if (leftMatrix.length != rightMatrix.heigth) {
        throw std::invalid_argument("Invalid matricies shapes");
    }
    
    linMatrix<T> resultMatrix = linMatrix<T>(leftMatrix.length, leftMatrix.heigth);
    double executionStartTime, executionEndTime;
    
    executionStartTime = omp_get_wtime();
    _parallelOptimizedBitMatricesMultiplication(leftMatrix.startPtr, rightMatrix.startPtr, resultMatrix.startPtr, leftMatrix.heigth, rightMatrix.length, leftMatrix.length);
    executionEndTime = omp_get_wtime();
    
    if (time != NULL) {
        *time = executionEndTime - executionStartTime;
    }
    return resultMatrix;
}




#pragma mark - Templates

template linMatrix<signed char> parallelNumericMatricesAddition(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> parallelNumericMatricesAddition(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> parallelNumericMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> parallelNumericMatricesAddition(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> parallelNumericMatricesAddition(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> parallelNumericMatricesAddition(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);


template linMatrix<signed char> parallelOptimizedNumericMatricesAddition(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> parallelOptimizedNumericMatricesAddition(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> parallelOptimizedNumericMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> parallelOptimizedNumericMatricesAddition(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> parallelOptimizedNumericMatricesAddition(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> parallelOptimizedNumericMatricesAddition(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);


template linMatrix<signed char> parallelNumericMatricesMultiplication(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> parallelNumericMatricesMultiplication(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> parallelNumericMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> parallelNumericMatricesMultiplication(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> parallelNumericMatricesMultiplication(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> parallelNumericMatricesMultiplication(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);


template linMatrix<signed char> parallelOptimizedNumericMatricesMultiplication(const linMatrix<signed char> leftMatrix, const linMatrix<signed char> rightMatrix, double *time);
template linMatrix<short> parallelOptimizedNumericMatricesMultiplication(const linMatrix<short> leftMatrix, const linMatrix<short> rightMatrix, double *time);
template linMatrix<int> parallelOptimizedNumericMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
template linMatrix<long> parallelOptimizedNumericMatricesMultiplication(const linMatrix<long> leftMatrix, const linMatrix<long> rightMatrix, double *time);
template linMatrix<float> parallelOptimizedNumericMatricesMultiplication(const linMatrix<float> leftMatrix, const linMatrix<float> rightMatrix, double *time);
template linMatrix<double> parallelOptimizedNumericMatricesMultiplication(const linMatrix<double> leftMatrix, const linMatrix<double> rightMatrix, double *time);

template linMatrix<int> parallelBitMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);

template linMatrix<int> parallelBitMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);

template linMatrix<int> parallelOptimizedBitMatricesAddition(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);

template linMatrix<int> parallelOptimizedBitMatricesMultiplication(const linMatrix<int> leftMatrix, const linMatrix<int> rightMatrix, double *time);
