// Implement polynomial multiplication using the FFT algorithm 

#include "fft.cpp"
#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <thread>
#include <fstream>


void multi_poly(complex<double>* res, complex<double>* p1, complex<double>* p2, int n){
    complex<double> poly1[2 * n]
    complex<double> poly2[2 * n]
    complex<double> result_poly[2 * n];
    #Apply fft to both polynomials 
    parallel_fft(poly1, p1, 2 * n);                    
    parallel_fft(poly2, p2, 2 * n); 
    #multiply the coefficients of the transformed polynomials
    for (i = 0; i<n; i++){
        result_poly[i] = poly1[i]*poly2[i]; 
    }
    #Apply the inverse fft to the multiplied polynomial to obtain the final polynomial 
    inverse_fft(res, result_poly, 2 * n); 
}; 

int main(){
    int n = 4;
    complex<double> poly1[2 * n] = {1,2,3,4}
    complex<double> poly2[2 * n] = {1,2,3,4}
    complex <double> res[2 * n];
    poly_mult(res, poly1, poly2, 4);
    for (int i = 0; i < 2 * n; i ++) {
        std::cout << res[i] << std::endl;
    }
}
