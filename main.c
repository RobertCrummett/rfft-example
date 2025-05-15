#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define FOURIER_PREFIX(name) name
#include "fourier.h"

#define DATA_SIZE 8
static double my_data[DATA_SIZE] = {1, 2, 3, 4, 5, 6, 7, 8};
static complex_t my_cdata[DATA_SIZE] = {
    {1, 0}, {2, 0}, {3, 0}, {4, 0},
    {5, 0}, {6, 0}, {7, 0}, {8, 0}
};

void auxspace(size_t size, complex_t** aux) {
    size_t larger_size = 1;
    while (larger_size >> 1 < size) larger_size <<= 1;
    *aux = calloc(size + 2 * larger_size + larger_size / 2, sizeof(**aux));
}

void rfft(double *data, complex_t *aux, size_t size, size_t stride) {
    //
    // Sources:
    // [1] https://www.robinscheibler.org/2013/02/13/real-fft.html
    // [2] https://kovleventer.com/blog/fft_real/
    //
    if (size == 0 || size == 1)
        return;

    if (size & 1) {
        fprintf(stderr, "error, the size of the data (%zu) must be even to use rfft; otherwise use fft\n", size);
        exit(1);
    }

    size_t half_size = size / 2;
    complex_t *cdata = (complex_t *)data;

    fft(cdata, aux, half_size, stride);

    // Unmixing twiddles
    for (size_t i = 1; i < half_size; i++)
        aux[i] = complex_exp((complex_t){0.0, -2.0*FOURIER_PI*i/size});

    // Note:
    // aux[0]         = (complex_t){1.0,  0.0};
    // aux[half_size] = (complex_t){0.0, -1.0};

    //
    // Packing without reallocation requires DC and Nyquist
    // be pushed together. This only works for series of even
    // length (size is even)
    //
    cdata[0] = (complex_t){
        complex_real(cdata[0])+complex_imag(cdata[0]),
        complex_real(cdata[0])-complex_imag(cdata[0]),
    };

    // complex_t Zx = complex_add(cdata[0], complex_conj(cdata[0]));
    // Zx = complex_mul(Zx, (complex_t){0.5, 0.0});

    // complex_t Zy = complex_sub(cdata[0], complex_conj(cdata[0]));
    // Zy = complex_mul(Zy, (complex_t){0.0, -0.5});
    //
    // cdata[0] = (complex_t){
    //     complex_real(complex_add(Zx, Zy)), // DC
    //     complex_real(complex_sub(Zx, Zy))  // Nyquist
    // };

    for (size_t i = 1; i < half_size/2+1; i++) {
        size_t j = half_size - i;

        complex_t Ci = cdata[i*stride];
        complex_t Cj = cdata[j*stride];
        complex_t Wi = aux[i];

        complex_t Zx = complex_add(Ci, complex_conj(Cj));
        Zx = complex_mul(Zx, (complex_t){0.5, 0.0});

        complex_t Zy = complex_sub(Ci, complex_conj(Cj));
        Zy = complex_mul(Zy, (complex_t){0.0, -0.5});

        if (i == j)
            cdata[i*stride] = complex_add(Zx, complex_mul(Wi, Zy));
        else {
            complex_t WiZy = complex_mul(Wi, Zy);
            cdata[i*stride] = complex_add(Zx, WiZy);
            cdata[j*stride] = complex_conj(complex_sub(Zx, WiZy));
        }
    }
}

complex_t *packed_ref(double *data, size_t size, size_t stride, size_t index) {
    size_t half_size = size / 2;
    if (index == 0) // DC FREQUENCY
        return (complex_t *)&data[0];
    else if (index == half_size) // NYQUIST FREQUENCY
        return (complex_t *)&data[1];
    else if (index < half_size) // POSITIVE FREQUENCY
        return (complex_t *)&data[2*index];
    else // NEGATIVE FREQUENCY
        return (complex_t *)&data[2*(size-index)];
}

complex_t rfft_index(double *data, size_t size, size_t stride, size_t index) {
    size_t half_size = size / 2;
    if (index == 0) // DC FREQUENCY
        return (complex_t){data[0], 0.0};
    else if (index == half_size) // NYQUIST FREQUENCY
        return (complex_t){data[1], 0.0};
    else if (index < half_size) // POSITIVE FREQUENCY
        return (complex_t){data[2*index], data[2*index+1]};
    else // NEGATIVE FREQUENCY
        return (complex_t){data[2*(size-index)], -data[2*(size-index)+1]};
}

int main(void) {
    complex_t* aux = NULL;
    auxspace(DATA_SIZE/2, &aux);
    if (aux == NULL)
        exit(1);

#ifdef PRINTOUT
    printf("Input:\n");
    for (int i = 0; i < DATA_SIZE; i++)
        printf("%lf\n", my_data[i]);
#endif // PRINTOUT

    rfft(my_data, aux, DATA_SIZE, 1);

#ifdef PRINTOUT
    printf("\n");
    printf("RFFT Output:\n");
    for (int i = 0; i < DATA_SIZE; i++)
        printf("%lf %lf\n", 
               complex_real(rfft_index(my_data, DATA_SIZE, 1, i)),
               complex_imag(rfft_index(my_data, DATA_SIZE, 1, i)));
#endif // PRINTOUT

    free(aux);
    aux = NULL;
    auxspace(DATA_SIZE, &aux);
    if (aux == NULL)
        exit(1);

    fft(my_cdata, aux, DATA_SIZE, 1);

#ifdef PRINTOUT
    printf("\n");
    printf("FFT Output:\n");
    for (int i = 0; i < DATA_SIZE; i++)
        printf("%lf %lf\n", complex_real(my_cdata[i]), complex_imag(my_cdata[i]));
#endif // PRINTOUT

    free(aux);
    return 0;
}
