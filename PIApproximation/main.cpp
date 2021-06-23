//
//  main.cpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#include <iostream>
#include <sstream>
#include <iomanip>

#include "omp.h"
#include "sta.hpp"
#include "mta.hpp"


int main(int argc, const char * argv[]) {
    double start, end;
    
    start = omp_get_wtime();
    startOMPPlaffAlgorithm(50000);
    end = omp_get_wtime();
    
    startChudnovskiyAlgorithm();
    
    std::cout << end - start <<std::endl;
    
    return 0;
}
