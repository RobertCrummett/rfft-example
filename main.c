#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define FOURIER_PREFIX(name) name
#include "fourier.h"

#define DATA_SIZE 8
static double my_data[DATA_SIZE] = {1, 2, 3, 4, 5, 6, 7, 8};

void auxspace(size_t size, complex_t** aux) {
    size_t larger_size = 1;
    while (larger_size >> 1 < size) larger_size <<= 1;
    *aux = calloc(size + 2 * larger_size + larger_size / 2, sizeof(**aux));
}

void rfft(double *data, size_t size) {
    //
    // Sources:
    // [1] https://www.robinscheibler.org/2013/02/13/real-fft.html
    // [2] https://kovleventer.com/blog/fft_real/
    //
    size_t half_size = size/2;
    complex_t *cdata = (complex_t *)data;

    complex_t* aux = NULL;
    auxspace(half_size, &aux); // for fft and precomputing twiddles
    if (aux == NULL) exit(1);
    fft(cdata, aux, half_size, 1);

    // Unmixing twiddles
    for (size_t i = 1; i < half_size/2+1; i++) {
        size_t j = half_size - i;
        aux[i] = complex_exp((complex_t){0.0, -2.0*FOURIER_PI*j/size});
    }
    // Note:
    // aux[0] = (complex_t){1.0, 0.0};
    // aux[half_size] = (complex_t){0.0, -1.0};
 
    complex_t Zx = complex_add(cdata[0], complex_conj(cdata[0]));
    Zx = complex_mul(Zx, (complex_t){0.5, 0.0});

    complex_t Zy = complex_sub(cdata[0], complex_conj(cdata[0]));
    Zy = complex_mul(Zy, (complex_t){0.0, -0.5});

    //
    // Packing without reallocation requires DC and Nyquist
    // be pushed together. This only works for series of even
    // length (half_size is even)
    //
    // TODO:
    // Organizing packing for series of odd length (half_size
    // is odd)
    //
    cdata[0] = (complex_t){
        complex_real(complex_add(Zx, Zy)), // DC
        complex_real(complex_sub(Zx, Zy))  // Nyquist
    };

    for (size_t i = 1; i < half_size/2+1; i++) {
        size_t j = half_size - i;

        Zx = complex_add(cdata[i], complex_conj(cdata[j]));
        Zx = complex_mul(Zx, (complex_t){0.5, 0.0});

        Zy = complex_sub(cdata[i], complex_conj(cdata[j]));
        Zy = complex_mul(Zy, (complex_t){0.0, -0.5});

        if (i == half_size-i)
            cdata[i] = complex_add(Zx, complex_mul(aux[i], complex_conj(Zy)));
        else {
            complex_t WZy = complex_mul(aux[i], complex_conj(Zy));
            cdata[i] = complex_add(Zx, WZy);
            cdata[half_size-i] = complex_conj(complex_sub(Zx, WZy));
        }
    }

    free(aux);
}

int main(void) {
    printf("Input:\n");
    for (int i = 0; i < DATA_SIZE; i++)
        printf("%lf\n", my_data[i]);

    rfft(my_data, DATA_SIZE);

    printf("\n");
    printf("Output:\n");
    for (int i = 0; i < DATA_SIZE; i+=2)
        printf("%lf %lf\n", my_data[i], my_data[i+1]);

    return 0;
}
