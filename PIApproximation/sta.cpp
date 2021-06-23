//
//  sta.cpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#include "sta.hpp"

double leibnizSeries() {
    double pi = 0.;
    double sign = 1.;
    double i = 0.;
    do
    {
        pi += sign / (2. * i + 1.);
        sign = -sign;
        ++i;
    } while (1. / (2. * i + 1.) > LEIBNIZ_COUNT_PI_DEC_PLACES);
    pi *= 4;
    return pi;
}


void startLeibnizSeries(unsigned digitsNumAfterPoint) {
    double piLeibnizValue = leibnizSeries();
    
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


double calculatePiWithChudnovskiyEquastion() {
    double serieSum = 0.0;
    for (unsigned int k = 0; k < CHUDNOVSKY_ALG_ITERATIONS; ++k) {
        serieSum += ((k % 2 == 1 ? -1.0: 1.0) * findFactorialOfNumber(6 * k) * (13591409 + (545140134 * k)))
        / (findFactorialOfNumber(3 * k) * pow(findFactorialOfNumber(k), 3.0) * pow(640320.0, 3.0 * k + 3.0/2.0));
    }
    return 1.0 / (12.0 * serieSum);
}


void startChudnovskiyAlgorithm() {
    double time = 0.0;
    time = omp_get_wtime();
    double pi = calculatePiWithChudnovskiyEquastion();
    time = omp_get_wtime() - time;
    
    std::cout << std::string("PI with Chudnovskiy Algorithm: ")
    << std::setprecision(6)
    << pi
    << std::endl;
    
    std::cout << std::string("Full PI with Chudnovskiy Algorithm: ")
    << std::setprecision(std::numeric_limits<long double>::digits10)
    << pi
    << std::endl;
    
    std::cout << std::string("Time for getting the limit PI value: ")
    << time
    << std::endl;
}



void startPlaffAlgorithm(unsigned int digitsNumAfterPoint) {
    unsigned int nessesaryPackagesCount = (unsigned int) (digitsNumAfterPoint / 16) + 1;
    unsigned int numbersOffseet = 16 - (nessesaryPackagesCount * 16 - digitsNumAfterPoint);
    std::vector<std::string> packedPiDigits(nessesaryPackagesCount);
    
    for (unsigned int i = 0; i < nessesaryPackagesCount; ++i) {
        packedPiDigits[i] = generateDigitsOfPiFromPosition(i * 16);
    }
    
    packedPiDigits[packedPiDigits.size() - 1] = packedPiDigits[packedPiDigits.size() - 1].substr(0, numbersOffseet);
    
    auto resPiString = createPiFromGeneratedStrings(packedPiDigits);
    std::cout << resPiString <<std::endl;
}



