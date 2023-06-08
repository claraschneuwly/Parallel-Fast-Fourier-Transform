// Implementation of a parallel radix-2 algorithm folloing Paper (3) : "A parallel FFT on an MIMD machine" by Amir AVERBUCH
#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <thread>

typedef std::complex<double> Complex;

void fft(Complex  *FFT_transformed, Complex  *x){
    const size_t N = x.size();
    // Trivial size-1 DFT base case
    if (N == 1){
        for (int i = 0; i < N; i++){
            FFT_transformed[i] = x[i];
        }
        return;
    }

    // Divide x into even and odd indices
    Complex x_even[N/2];
    Complex x_odd[N/2];
    Complex evenTransformed[N/2];
    Complex oddTransformed[N/2];

    for (size_t i = 0; i < N; i += 2) {
        x_even.push_back(x[i]);
        x_odd.push_back(x[i+1]);
    }

    // Recursion
    fft(evenTransformed, x_even);
    fft(oddTransformed, x_odd);

    // Combine the even and odd arrays
    for (int k = 0; k < N / 2; k++) {
        Complex factor = std::polar(1.0, -2.0 * M_PI * k / N) * oddTransformed[k];
        FFT_transformed[k] = evenTransformed[k] + factor;
        FFT_transformed[k + N / 2] = evenTransformed[k] - factor;
    }
}


Complex bit_reversal(double x, int n) {
    Complex res[n]; 
    for (int i=0; i<n/2; i++) {
            res[i] = x[n-i-1]; 
    } 
    return res; 
}
