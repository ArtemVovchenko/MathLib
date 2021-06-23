//
//  StaticFuncs.hpp
//  Lab4
//
//  Created by Artem Vovchenko on 20.12.2020.
//

#ifndef StaticFuncs_hpp
#define StaticFuncs_hpp

#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

unsigned long findFactorialOfNumber(unsigned int number);

std::string generateDigitsOfPiFromPosition(unsigned int position);

std::string createPiFromGeneratedStrings(std::vector<std::string> piDigitsStrings);

#endif /* StaticFuncs_hpp */
