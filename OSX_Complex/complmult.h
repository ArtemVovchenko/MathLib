//
//  complmult.h
//  OSX_Complex
//
//  Created by Artem on 19.10.2020.
//

#ifndef complmult_h
#define complmult_h

#include <stdio.h>
#include <immintrin.h>
#include <complex.h>
#include <time.h>

#define ARRAY_SIZE 16777216


void _generate_complex_array(complex float *__dst, unsigned __dst_size);


void _mult_complex_floats(complex float *__x,
                          complex float *__y,
                          complex float *__dst,
                          unsigned __dst_size,
                          double *time);


void _mult_complex_float_SSE(complex float *__x,
                             complex float *__y,
                             complex float *__out,
                             unsigned __out_size,
                             double *time);


void _mult_complex_float_AVX(complex float *__x,
                             complex float *__y,
                             complex float *__out,
                             unsigned __out_size,
                             double *time);


#endif /* complmult_h */
