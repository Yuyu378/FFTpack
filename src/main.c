// 
// main.c
// 
//      Main file - Usage Example
// 

#include <math.h>
#include <time.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft.h"
#include "random.h"

static void output_complex_value(dcomplex value);
static void output_complex_array(dcomplex* array, unsigned len);

int main(void) {

    pcg32_init((uint64_t)time(NULL));                           // initialize pcg
    unsigned N = randint(1, 20);                                // generate a random length
    dcomplex* data = (dcomplex*)calloc(N, sizeof(dcomplex));    // allocate memory
    if (data == NULL) {                                         // check if allocation was succeed
        printf("Allocation failed.\n");
        printf("\nPress any key to continue . . . ");
        if (_getch()) printf("\n");
        exit(1);                                                // output error msg and exit if unsuccessful
    }
    for (unsigned i = 0; i < N; i++) {
        data[i] = uniform(-20, 20) + I * uniform(-20, 20);      // fill data array with random numbers
        output_complex_value(data[i]);
    }
    printf("\n");                                               // output complex value

    dcomplex* fftres = fft(data, N, backward);                  // Compute FFT
    if (fftres == NULL) {                                       // raise error if calculation fails
        printf("Calculate fft failed.\n");
        printf("\nPress any key to continue . . . ");
        if (_getch()) printf("\n");
        exit(1);
    }
    output_complex_array(fftres, N);                            // output FFT result

    dcomplex* ifftres = ifft(fftres, N, backward);              // perform backward calculation
    if (ifftres == NULL) {                                      // raise error if calculation fails
        printf("Calculate ifft failed.\n");
        printf("\nPress any key to continue . . . ");
        if (_getch()) printf("\n");
        exit(1);
    }
    output_complex_array(ifftres, N);                           // output IFFT result, should be same as data array

    free(data);                                                 // free memory
    free(fftres);                                               // Important!! fft() & ifft() will allocate memory,
    free(ifftres);                                              // should be freed by the caller

    printf("\nPress any key to continue . . . ");
    if (_getch()) printf("\n");
    return 0;
}

static void output_complex_value(dcomplex value) {
    #ifdef __GNUC__
    double re = creal(value);
    double im = cimag(value);
    #endif
    #ifdef _MSC_VER
    double re = value._Val[0];
    double im = value._Val[1];
    #endif
    #ifdef ARDUINO
    double re = value.re;
    double im = value.im;
    #endif
    printf("    %.6lf", re);
    if (im < 0) {
        printf(" - ");
    }
    else {
        printf(" + ");
    }
    printf("%.6lfi\n", fabs(im));
    return;
}

static void output_complex_array(dcomplex* array, unsigned len) {
    for (unsigned i = 0; i < len; i++) {
        output_complex_value(array[i]);
    }
    printf("\n");
    return;
}
