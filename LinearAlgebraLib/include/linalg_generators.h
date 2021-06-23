//
//  linalg_generators.h
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#ifndef linalg_generators_h
#define linalg_generators_h

#define RANDOM_INT_8_MIN    0x40
#define RANDOM_INT_8_MAX    0x3F

#define RANDOM_INT_16_MIN   0x40
#define RANDOM_INT_16_MAX   0x3FF

#define RANDOM_INT_32_MIN   0x40
#define RANDOM_INT_32_MAX   0x3FF

#define RANDOM_INT_64_MIN   0x400
#define RANDOM_INT_64_MAX   0x3FF

#define RATIONAL_MIN -100
#define RATIONAL_MAX 100

#include <random>
#include <string>
#include <unordered_map>

#include "linalg_attributes.h"
#include "linalg_types.h"


template <typename T>
void fillMatrixWithRandomData(linMatrix<T> matrix);


#endif /* linalg_generators_h */
