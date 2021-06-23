//
//  linalg_types.cpp
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#include "linalg_types.h"

template <typename T>
INTERNAL int _getElementType(T element) {
    
    if (typeid(T) == typeid(signed char*)) {
        return SUPPORTED_TYPES::INT_8;
    }
    if (typeid(T) == typeid(short*)) {
        return SUPPORTED_TYPES::INT_16;
    }
    if (typeid(T) == typeid(int*)) {
        return SUPPORTED_TYPES::INT;
    }
    if (typeid(T) == typeid(float*)) {
        return SUPPORTED_TYPES::FLOAT;
    }
    if (typeid(T) == typeid(double*)) {
        return SUPPORTED_TYPES::DOUBLE;
    }
    
    std::stringstream errorMessageFormater;
    errorMessageFormater << "Unsupported type of agrument: typeid(arg).name() == \""
    << typeid(T).name() << "\"";
    throw std::invalid_argument(errorMessageFormater.str());
}

template int _getElementType<signed char*>(signed char* element);
template int _getElementType<short*>(short* element);
template int _getElementType<int*>(int* element);
template int _getElementType<long*>(long* element);
template int _getElementType<float*>(float* element);
template int _getElementType<double*>(double* element);

