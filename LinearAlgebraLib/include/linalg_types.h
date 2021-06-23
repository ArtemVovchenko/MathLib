//
//  linalg_types.h
//  LinearAlgebra
//
//  Created by Artem Vovchenko on 21.12.2020.
//

#ifndef linalg_types_h
#define linalg_types_h

#include "linalg.h"
#include "linalg_attributes.h"

#include <sstream>
#include <stdexcept>

typedef enum {
    INT_8,
    INT_16,
    INT,
    LONG,
    FLOAT,
    DOUBLE
} SUPPORTED_TYPES;

template <typename T>
INTERNAL int _getElementType(T element);


#endif /* linalg_types_h */
