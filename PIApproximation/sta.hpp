//
//  sta.hpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#ifndef sta_hpp
#define sta_hpp

#include <stdio.h>
#include <iostream>
#include "omp.h"
#include "StaticFuncs.hpp"

#define LEIBNIZ_COUNT_PI_DEC_PLACES 1e-9
#define LEIBNIZ_COUNT_TRUE_NUMS 8
#define CHUDNOVSKY_ALG_ITERATIONS 10

void startLeibnizSeries(unsigned digitsNumAfterPoint);

void startChudnovskiyAlgorithm();

void startPlaffAlgorithm(unsigned int digitsNumAfterPoint);


#endif /* sta_hpp */
