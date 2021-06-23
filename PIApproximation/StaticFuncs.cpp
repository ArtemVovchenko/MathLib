//
//  StaticFuncs.cpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#include "StaticFuncs.hpp"

unsigned long calculateProductFactorialTreeForNumbers(unsigned long leftNumber, unsigned long rightNumber) {
    if (leftNumber > rightNumber) {
        return 1;
    }
    if (leftNumber == rightNumber) {
        return leftNumber;
    }
    if (rightNumber - leftNumber == 1) {
        return leftNumber * rightNumber;
    }
    
    unsigned long numbersMiddle = (leftNumber + rightNumber) / 2;
    return calculateProductFactorialTreeForNumbers(leftNumber, numbersMiddle) * calculateProductFactorialTreeForNumbers(numbersMiddle + 1, rightNumber);
}

unsigned long findFactorialOfNumber(unsigned int number) {
    if (number < 0) {
        return 0;
    }
    if (number == 0) {
        return 1;
    }
    if (number == 1 || number == 2) {
        return number;
    }
    return calculateProductFactorialTreeForNumbers(2, number);
}

double modularExponent(double power, double base) {
    static double twoPowersTable[25];
    static int powersTableIsFilled = 0;
    double p1, pt, r;
    
    if (!powersTableIsFilled) {
        powersTableIsFilled = 1;
        twoPowersTable[0] = 1.;

        for (int i = 1; i < 25; ++i){
            twoPowersTable[i] = 2. * twoPowersTable[i-1];
        }
    }
    
    if (base == 1) {
        return 0.0;
    }
    
    int closestTwoPowerToInteredPower;
    for (closestTwoPowerToInteredPower = 0; closestTwoPowerToInteredPower < 25; ++closestTwoPowerToInteredPower) {
        if (twoPowersTable[closestTwoPowerToInteredPower] > power) {
            break;
        }
    }
    
    pt = twoPowersTable[closestTwoPowerToInteredPower - 1];
    p1 = power;
    r = 1.0;
    
    for (unsigned int j = 1; j <= closestTwoPowerToInteredPower; j++){
        if (p1 >= pt) {
            r = 16. * r;
            r = r - (int) (r / base) * base;
            p1 = p1 - pt;
        }
        pt = 0.5 * pt;
        if (pt >= 1.) {
            r = r * r;
            r = r - (int) (r / base) * base;
        }
    }
    
    return r;
}


double getSeries(unsigned int denominatorAddCoeficient, unsigned int digitsStartPosition) {
    double epsilon = 1e-17;
    double denominator, p, t;
    double serieSum = 0.0;
    
    for (long k = 0; k < digitsStartPosition; ++k) {
        denominator = 8 * k + denominatorAddCoeficient;
        p = digitsStartPosition - k;
        t = modularExponent(p, denominator);
        serieSum = serieSum + t / denominator;
        serieSum = serieSum - (int) serieSum;
    }
    
    for (long k = digitsStartPosition; k <= digitsStartPosition + 100; k++){
        denominator = 8 * k + denominatorAddCoeficient;
        t = pow (16.0, (double) ((long)digitsStartPosition - k)) / denominator;
        if (t < epsilon) break;
        serieSum = serieSum + t;
        serieSum = serieSum - (int) serieSum;
    }
    
    return serieSum;
}

std::string toHexPackage (double x) {
    int i;
    double y;
    char hexDigits[] = "0123456789ABCDEF";
    char hexAnswer[17];
    y = abs(x);
    for (i = 0; i < 16; i++){
        y = 16.0 * (y - floor(y));
        hexAnswer[i] = hexDigits[(int) y];
    }
    hexAnswer[16] = '\0';
    std::stringstream hexToDecimalConverter;
    std::string hexAnswerString = std::string(hexAnswer);
    return hexAnswerString;
}

std::string generateDigitsOfPiFromPosition(unsigned int position) {
    double firstSerie = getSeries(1, position);
    double secondSerie = getSeries(4, position);
    double thirdSerie = getSeries(5, position);
    double fourthSerie = getSeries(6, position);
    
    double pid = 4. * firstSerie - 2. * secondSerie - thirdSerie - fourthSerie;
    pid = pid - (int) pid + 1.0;
    return toHexPackage(pid);
}

std::string createPiFromGeneratedStrings(std::vector<std::string> piDigitsStrings) {
    std::stringstream piDigitsConcatenater;
    piDigitsConcatenater << std::hex <<"0x3,";
    for (std::string digitsPackage : piDigitsStrings) {
        piDigitsConcatenater << std::hex << digitsPackage;
    }
    return piDigitsConcatenater.str();
}
