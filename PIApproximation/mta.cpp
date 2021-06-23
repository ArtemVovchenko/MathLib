//
//  mta.cpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#include "mta.hpp"


double leibnizOMPSeries() {
    double pi = 0.0;
    long serieElementsCount = LEIBNIZ_ASYNC_ITER;
    #pragma omp parallel default(none) shared(serieElementsCount, pi)
    {
        #pragma omp for reduction( + : pi )
        for (long i = 0; i < serieElementsCount; i += 2 ) {
            pi += 1.0 / (2 * i + 1);
        }
        
        #pragma omp for reduction( + : pi )
        for (long i = 1; i < serieElementsCount; i += 2) {
            pi += -1.0 / (2 * i + 1);
        }
    }
    
    return pi * 4;
}


void startLeibnizOMPSeries(unsigned digitsNumAfterPoint) {
    double piLeibnizValue = leibnizOMPSeries();
    
    //TODO: Исправить вывод
    std::stringstream converterToString;
    converterToString << std::dec << piLeibnizValue;
    std::string piString = converterToString.str();
    if (LEIBNIZ_COUNT_TRUE_NUMS <= digitsNumAfterPoint) {
        std::cout << piString << std::endl;
        return;
    }
    std::cout << piString.substr(0, digitsNumAfterPoint + 2);
}

double startChudnovskiyOMPAlgorithm() {
    double serieSum = 0.0;
    int loopIterations = CHUDNOVSKY_ALG_ITERATIONS;
    
    #pragma omp parallel default(none) shared(loopIterations, serieSum)
    {
        #pragma omp for reduction(+ : serieSum)
        for (int k = 0; k < loopIterations; k += 2) {
            serieSum += (1.0 * findFactorialOfNumber(6 * k) * (13591409 + (545140134 * k)))
            / (findFactorialOfNumber(3 * k) * pow(findFactorialOfNumber(k), 3.0) * pow(640320.0, 3.0 * k + 3.0/2.0));
        }
        
        #pragma omp for reduction(+ : serieSum)
        for (int k = 1; k < loopIterations; k += 2) {
            serieSum += (-1.0 * findFactorialOfNumber(6 * k) * (13591409 + (545140134 * k)))
            / (findFactorialOfNumber(3 * k) * pow(findFactorialOfNumber(k), 3.0) * pow(640320.0, 3.0 * k + 3.0/2.0));
        }
    }
    
    return 1.0 / (12.0 * serieSum);
}

void startOMPPlaffAlgorithm(unsigned int digitsNumAfterPoint) {
    unsigned int nessesaryPackagesCount = (unsigned int) (digitsNumAfterPoint / 16) + 1;
    unsigned int numbersOffseet = 16 - (nessesaryPackagesCount * 16 - digitsNumAfterPoint);
    std::vector<std::string> packedPiDigits(nessesaryPackagesCount);
    
    #pragma omp parallel for
    for (unsigned int i = 0; i < nessesaryPackagesCount; ++i) {
        packedPiDigits[i] = generateDigitsOfPiFromPosition(i * 16);
    }
    
    packedPiDigits[packedPiDigits.size() - 1] = packedPiDigits[packedPiDigits.size() - 1].substr(0, numbersOffseet);
    
    auto resPiString = createPiFromGeneratedStrings(packedPiDigits);
    std::cout << resPiString <<std::endl;
}
