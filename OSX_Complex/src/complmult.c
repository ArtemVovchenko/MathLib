//
//  complmult.c
//  OSX_Complex
//
//  Created by Artem on 19.10.2020.
//

#include "complmult.h"

void _generate_complex_array(complex float *__dst,
                             unsigned __dst_size) {
    complex float generated;
    for (unsigned i = 0 ; i < __dst_size; ++i) {
        double real = (double)rand()/RAND_MAX;
        double imag = (double)rand()/RAND_MAX;
        if (!(i % 2)) {
            real *= -1;
        }
        if (!(i % 3) && !(i % 2)) {
            imag *= -1;
        }
        generated = (float)real + (double)imag * I;
        __dst[i] = generated;
    }
}


void _mult_complex_floats(complex float *__x,
                          complex float *__y,
                          complex float *__dst,
                          unsigned __dst_size,
                          double *time) {
    clock_t start, finish;
    start = clock();
    for (unsigned i = 0; i < __dst_size; ++i) {
        __dst[i] = __x[i] * __y[i];
    }
    finish = clock();
    
    *time = (double)(finish - start) / CLOCKS_PER_SEC;
}


    
void _mult_complex_float_SSE(complex float *__x,
                             complex float *__y,
                             complex float *__out,
                             unsigned __out_size,
                             double *time) {
    
    __m128 *x_vec = (__m128*)__x;
    __m128 *y_vec = (__m128*)__y;
    __m128 *out_vec = (__m128*)__out;
    __m128 temp1_vec, temp2_vec;
    
    clock_t start, finish;
    start = clock();
    
    for (unsigned int i = 0 ; i < __out_size / 2; ++i) {
        temp1_vec = _mm_shuffle_ps(x_vec[i], x_vec[i], 0xA0);
        temp2_vec = _mm_shuffle_ps(x_vec[i], x_vec[i], 0xF5);
        temp1_vec = _mm_mul_ps(temp1_vec, y_vec[i]);
        temp2_vec = _mm_mul_ps(temp2_vec, y_vec[i]);
        temp2_vec = _mm_shuffle_ps(temp2_vec, temp2_vec, 0xB1);
        out_vec[i] = _mm_addsub_ps(temp1_vec, temp2_vec);
    }
    
    finish = clock();
    *time = (double)(finish - start) / CLOCKS_PER_SEC;
}


#ifdef __AVX__
void _mult_complex_float_AVX(complex float *__x,
                             complex float *__y,
                             complex float *__out,
                             unsigned __out_size,
                             double *time) {
    
    __m256 *x_vec = (__m256*)__x;
    __m256 *y_vec = (__m256*)__y;
    __m256 *out_vec = (__m256*)__out;
    __m256 temp1_vec, temp2_vec;
    
    clock_t start, finish;
    start = clock();
    
    for (unsigned i = 0 ; i < __out_size / 4; ++i) {
        temp1_vec = _mm256_shuffle_ps(x_vec[i], x_vec[i], 0xA0);
        temp2_vec = _mm256_shuffle_ps(x_vec[i], x_vec[i], 0xF5);
        temp1_vec = _mm256_mul_ps(temp1_vec, y_vec[i]);
        temp2_vec = _mm256_mul_ps(temp2_vec, y_vec[i]);
        temp2_vec = _mm256_shuffle_ps(temp2_vec, temp2_vec, 0XB1);
        out_vec[i] = _mm256_addsub_ps(temp1_vec, temp2_vec);
    }
    
    finish = clock();
    *time = (double)(finish - start) / CLOCKS_PER_SEC;
}
#endif //__AVX__
