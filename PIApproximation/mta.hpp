//
//  mta.hpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#ifndef mta_hpp
#define mta_hpp

#include <stdio.h>
#include <iostream>
#include "StaticFuncs.hpp"
#include "omp.h"

#define LEIBNIZ_ASYNC_ITER 300000000
#define LEIBNIZ_COUNT_TRUE_NUMS 8
#define CHUDNOVSKY_ALG_ITERATIONS 10

void startLeibnizOMPSeries(unsigned digitsNumAfterPoint);

double startChudnovskiyOMPAlgorithm();

void startOMPPlaffAlgorithm(unsigned int digitsNumAfterPoint);

#endif /* mta_hpp */
