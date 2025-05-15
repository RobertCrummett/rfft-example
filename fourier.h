#ifndef FOURIER_H
#define FOURIER_H

#include <stddef.h>

#ifdef FOURIER_STATIC
#define FOURIER_PUBLICDEC static
#define FOURIER_PUBLICDEF static
#else
#define FOURIER_PUBLICDEC extern
#define FOURIER_PUBLICDEF
#endif // FOURIER_STATIC

#ifdef COMPLEX_STATIC
#define COMPLEX_PUBLICDEC static
#define COMPLEX_PUBLICDEF static
#else
#define COMPLEX_PUBLICDEC extern
#define COMPLEX_PUBLICDEF
#endif // COMPLEX_STATIC

#ifndef FOURIER_PREFIX
#define FOURIER_PREFIX(name) fourier_##name
#endif // FOURIER_PREFIX

typedef struct complex_t {
    double real;
    double imag;
} complex_t;

/* Expected size of the auxillary space for the chirp Z transform and internal radix FOURIER calls */

// size_t larger_size = 1;
// while (larger_size >> 1 < size) larger_size <<= 1;
//
// complex_t* aux = calloc(size + 2 * larger_size + larger_size / 2, sizeof(*aux));
//
FOURIER_PUBLICDEC void FOURIER_PREFIX(fft)(complex_t *data, complex_t *aux, size_t size, size_t stride);
FOURIER_PUBLICDEC void FOURIER_PREFIX(ifft)(complex_t *data, complex_t *aux, size_t size, size_t stride);

COMPLEX_PUBLICDEC complex_t complex_add(complex_t a, complex_t b);
COMPLEX_PUBLICDEC complex_t complex_sub(complex_t a, complex_t b);
COMPLEX_PUBLICDEC complex_t complex_mul(complex_t a, complex_t b);
COMPLEX_PUBLICDEC complex_t complex_div(complex_t a, complex_t b);

COMPLEX_PUBLICDEC double complex_real(complex_t z);
COMPLEX_PUBLICDEC double complex_imag(complex_t z);
COMPLEX_PUBLICDEC double complex_abs(complex_t z);
COMPLEX_PUBLICDEC double complex_arg(complex_t z);

COMPLEX_PUBLICDEC complex_t complex_conj(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_proj(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_exp(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_log(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_pow(complex_t base, complex_t exponent);
COMPLEX_PUBLICDEC complex_t complex_sqrt(complex_t z);

COMPLEX_PUBLICDEC complex_t complex_sin(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_cos(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_tan(complex_t z);

COMPLEX_PUBLICDEC complex_t complex_sinh(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_cosh(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_tanh(complex_t z);

COMPLEX_PUBLICDEC complex_t complex_asin(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_acos(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_atan(complex_t z);

COMPLEX_PUBLICDEC complex_t complex_asinh(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_acosh(complex_t z);
COMPLEX_PUBLICDEC complex_t complex_atanh(complex_t z);

#endif // FOURIER_H

#ifdef FOURIER_IMPLEMENTATION

#include <math.h>
#include <stddef.h>

#ifndef FOURIER_PI
#define FOURIER_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#endif // FOURIER_PI

#ifndef COMPLEX_PI
#define COMPLEX_PI 3.1415926535897932384626433832795028841971693993751058209749445923078164062
#endif // COMPLEX_PI

#ifndef FOURIER_I
#define FOURIER_I ((complex_t){0.0, 1.0})
#endif // FOURIER_I

static void *fourier__memset(void *s, int c, size_t n) {
    unsigned char *dst = s;
    while (n > 0) {
        *dst = (unsigned char)c;
        dst++;
        n--;
    }
    return s;
}

static void fourier__swap(complex_t* a, complex_t* b) {
	complex_t temp = *a;
	*a = *b;
	*b = temp;
}

static void fourier__shuffle(complex_t* data, size_t size, size_t stride) {
	size_t position, mask, target = 0;
	for (position = 0; position < size; position++) {
		if (target > position)
			fourier__swap(&data[position * stride], &data[target * stride]);

		mask = size >> 1;

		while (target & mask) {
			target &= ~mask;
			mask >>= 1;
		}
		target |= mask;
	}
}

static void fourier__radixfft(complex_t *data, complex_t *aux, size_t size, size_t stride) {
	size_t i, step, jump, group, pair, match, halfsize = size / 2;
	complex_t twiddle, product;

	fourier__shuffle(data, size, stride);

	for (i = 0; i < halfsize; i++)
		aux[i] = complex_exp(complex_mul(FOURIER_I, (complex_t){-FOURIER_PI * i / halfsize, 0.0}));

	for (step = 1; step < size; step <<= 1) {
		jump = step << 1;
		twiddle = (complex_t){1.0, 0.0};
		for (group = 0; group < step; group++) {
			for (pair = group; pair < size; pair += jump) {
				match = pair + step;
				product = complex_mul(twiddle, data[match * stride]);
				data[match * stride] = complex_sub(data[pair * stride], product);
				data[pair * stride] = complex_add(data[pair * stride], product);
			}

			if (group + 1 == step)
				continue;

			twiddle = aux[(group + 1) * halfsize / step];
		}
	}
}

static void fourier__radixifft(complex_t *data, complex_t *aux, size_t size, size_t stride) {
	size_t i, step, jump, group, pair, match, halfsize = size / 2;
	complex_t twiddle, product;

	fourier__shuffle(data, size, stride);

	for (i = 0; i < halfsize; i++)
		aux[i] = complex_exp(complex_mul(FOURIER_I, (complex_t){FOURIER_PI * i / halfsize, 0.0}));


	for (step = 1; step < size; step <<= 1) {
		jump = step << 1;
		twiddle = (complex_t){1.0, 0.0};
		for (group = 0; group < step; group++) {
			for (pair = group; pair < size; pair += jump) {
				match = pair + step;
				product = complex_mul(twiddle, data[match * stride]);
				data[match * stride] = complex_sub(data[pair * stride], product);
				data[pair * stride] = complex_add(data[pair * stride], product);
			}

			if (group + 1 == step)
				continue;

			twiddle = aux[(group + 1) * halfsize / step];
		}
	}
}

FOURIER_PUBLICDEF void FOURIER_PREFIX(fft)(complex_t *data, complex_t *aux, size_t size, size_t stride) {
	size_t larger_size = 1, i;
	while (larger_size >> 1 < size)
		larger_size <<= 1;

	for (i = 0; i < size; i++)
		aux[i] = complex_exp(complex_mul(FOURIER_I, (complex_t){-FOURIER_PI * i * i / size, 0.0}));

	for (i = 0; i < size; i++)
		aux[size + i] = complex_mul(aux[i], data[i * stride]);

	fourier__memset(aux + size + larger_size, 0, larger_size * sizeof(*aux));
	aux[size + larger_size] = aux[0];
	for (i = 1; i < size; i++) {
		aux[size + larger_size + i]     = complex_conj(aux[i]);
		aux[size + 2 * larger_size - i] = complex_conj(aux[i]);
	}

	fourier__radixfft(aux + size, aux + size + 2 * larger_size, larger_size, 1);
	fourier__radixfft(aux + size + larger_size, aux + size + 2 * larger_size, larger_size, 1);

	for (i = 0; i < larger_size; i++)
		aux[size + i] = complex_mul(aux[size + i], aux[size + larger_size + i]);

	fourier__radixifft(aux + size, aux + size + 2 * larger_size, larger_size, 1);

	for (i = 0; i < size; i++)
		data[i * stride] = complex_div(complex_mul(aux[size + i], aux[i]), (complex_t){larger_size, 0.0});
}

FOURIER_PUBLICDEF void FOURIER_PREFIX(ifft)(complex_t *data, complex_t *aux, size_t size, size_t stride) {
	size_t larger_size = 1, i;
	while (larger_size >> 1 < size)
		larger_size <<= 1;

	for (i = 0; i < size; i++)
		aux[i] = complex_exp(complex_mul(FOURIER_I, (complex_t){FOURIER_PI * i * i / size, 0.0}));

	for (i = 0; i < size; i++)
		aux[size + i] = complex_mul(aux[i], data[i * stride]);

	fourier__memset(aux + size + larger_size, 0, larger_size * sizeof(*aux));
	aux[size + larger_size] = aux[0];
	for (i = 1; i < size; i++) {
		aux[size + larger_size + i]     = complex_conj(aux[i]);
		aux[size + 2 * larger_size - i] = complex_conj(aux[i]);
	}

	fourier__radixfft(aux + size, aux + size + 2 * larger_size,larger_size, 1);
	fourier__radixfft(aux + size + larger_size, aux + size + 2 * larger_size, larger_size, 1);

	for (i = 0; i < larger_size; i++)
		aux[size + i] = complex_mul(aux[size + i], aux[size + larger_size + i]);

	fourier__radixifft(aux + size, aux + size + 2 * larger_size, larger_size, 1);

	for (i = 0; i < size; i++)
		data[i * stride] = complex_div(complex_mul(aux[size + i], aux[i]), (complex_t){larger_size, 0.0});
}

COMPLEX_PUBLICDEF complex_t complex_add(complex_t a, complex_t b) {
    return (complex_t){complex_real(a)+complex_real(b), complex_imag(a)+complex_imag(b)};
}

COMPLEX_PUBLICDEF complex_t complex_sub(complex_t a, complex_t b) {
    return (complex_t){complex_real(a)-complex_real(b), complex_imag(a)-complex_imag(b)};
}

COMPLEX_PUBLICDEF complex_t complex_mul(complex_t a, complex_t b) {
    return (complex_t){complex_real(a)*complex_real(b)-complex_imag(a)*complex_imag(b),
                       complex_real(a)*complex_imag(b)+complex_imag(a)*complex_real(b)};
}

COMPLEX_PUBLICDEF complex_t complex_div(complex_t a, complex_t b) {
    double b2 = complex_real(b)*complex_real(b) + complex_imag(b)*complex_imag(b);
    return (complex_t){(complex_real(a)*complex_real(b)+complex_imag(a)*complex_imag(b))/b2,
                       (complex_imag(a)*complex_real(b)-complex_real(a)*complex_imag(b))/b2};
}

COMPLEX_PUBLICDEF double complex_real(complex_t z) {
    return z.real;
}

COMPLEX_PUBLICDEF double complex_imag(complex_t z) {
    return z.imag;
}

COMPLEX_PUBLICDEF double complex_abs(complex_t z) {
    // TODO: remove sqrt
    return sqrt(complex_real(z)*complex_real(z) + complex_imag(z)*complex_imag(z));
}

COMPLEX_PUBLICDEF double complex_arg(complex_t z) {
    // TODO: remove atan2
    return atan2(complex_imag(z), complex_real(z));
}

COMPLEX_PUBLICDEF complex_t complex_conj(complex_t z) {
    return (complex_t){complex_real(z), -complex_imag(z)};
}

COMPLEX_PUBLICDEF complex_t complex_proj(complex_t z) {
    // TODO: Remove INFINITY, copysign, isinf
    if (isinf(complex_real(z)) || isinf(complex_imag(z)))
        return (complex_t){INFINITY, copysign(0.0, complex_imag(z))};
    else
        return z;
}

COMPLEX_PUBLICDEF complex_t complex_exp(complex_t z) {
    // TODO: Remove exp, sin, cos
    double a = exp(complex_real(z));
    double x = cos(complex_imag(z));
    double y = sin(complex_imag(z));
    return (complex_t){a*x, a*y};
}

COMPLEX_PUBLICDEF complex_t complex_log(complex_t z) {
    // TODO: Remove log
    return (complex_t){log(complex_abs(z)), complex_arg(z)};
}

COMPLEX_PUBLICDEF complex_t complex_pow(complex_t base, complex_t exponent) {
    return complex_exp(complex_mul(complex_log(base), exponent));
}

COMPLEX_PUBLICDEF complex_t complex_sqrt(complex_t z) {
    const complex_t half = {0.5, 0.0};
    return complex_exp(complex_mul(complex_log(z), half));
}

COMPLEX_PUBLICDEF complex_t complex_sin(complex_t z) {
    double x = complex_real(z);
    double y = complex_imag(z);
    // TODO: Remove cos, sin, cosh, sinh
    return (complex_t){sin(x)*cosh(y), cos(x)*sinh(y)};
}

COMPLEX_PUBLICDEF complex_t complex_cos(complex_t z) {
    double x = complex_real(z);
    double y = complex_imag(z);
    // TODO: Remove cos, sin, cosh, sinh
    return (complex_t){cos(x)*cosh(y), -sin(x)*sinh(y)};
}

COMPLEX_PUBLICDEF complex_t complex_tan(complex_t z) {
    return complex_div(complex_sin(z), complex_cos(z));
}

COMPLEX_PUBLICDEF complex_t complex_sinh(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    const complex_t ineg = {0.0,-1.0};
    return complex_mul(ineg, complex_sin(complex_mul(ipos, z)));
}

COMPLEX_PUBLICDEF complex_t complex_cosh(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    return complex_cos(complex_mul(ipos, z));
}

COMPLEX_PUBLICDEF complex_t complex_tanh(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    const complex_t ineg = {0.0,-1.0};
    return complex_mul(ineg, complex_tan(complex_mul(ipos, z)));
}

COMPLEX_PUBLICDEF complex_t complex_asin(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    const complex_t ineg = {0.0,-1.0};
    return complex_mul(ineg, complex_asinh(complex_mul(ipos, z)));
}

COMPLEX_PUBLICDEF complex_t complex_acos(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    const complex_t one = {1.0, 0.0};
    const complex_t half_pi = {COMPLEX_PI / 2.0, 0.0};
    
    complex_t a = complex_sqrt(complex_sub(one, complex_mul(z, z)));
    complex_t b = complex_log(complex_add(complex_mul(ipos, z), a));

    return complex_add(half_pi, complex_mul(ipos, b));
}

COMPLEX_PUBLICDEF complex_t complex_atan(complex_t z) {
    const complex_t ipos = {0.0, 1.0};
    const complex_t ineg = {0.0,-1.0};
    return complex_mul(ineg, complex_atanh(complex_mul(ipos, z)));
}

COMPLEX_PUBLICDEF complex_t complex_asinh(complex_t z) {
    const complex_t one = {1.0, 0.0};

    complex_t a = complex_add(one, complex_mul(z, z));
    complex_t b = complex_add(z, complex_sqrt(a));

    return complex_log(b);
}

COMPLEX_PUBLICDEF complex_t complex_acosh(complex_t z) {
    const complex_t one = {1.0, 0.0};
    
    complex_t a = complex_sqrt(complex_add(z, one));
    complex_t b = complex_sqrt(complex_sub(z, one));
    complex_t c = complex_add(z, complex_mul(a, b));

    return complex_log(c);
}

COMPLEX_PUBLICDEF complex_t complex_atanh(complex_t z) {
    const complex_t one = {1.0, 0.0};
    const complex_t two = {2.0, 0.0};

    complex_t a = complex_log(complex_add(one, z));
    complex_t b = complex_log(complex_sub(one, z));
    complex_t c = complex_sub(a, b);
    
    return complex_div(c, two);
}

#endif // FOURIER_IMPLEMENTATION
