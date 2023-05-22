#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

typedef std::complex<double> Complex;
typedef std::vector<Complex> CArray;

// Compute the DFT of x
void cooley_tukey(CArray& x) {
    const size_t N = x.size();

    // Trivial size-1 DFT base case
    if (N == 1) {
        return;
    } 
    
    else {
        // Divide x into even and odd indices
        CArray x_even;
        CArray x_odd;
        for (size_t i = 0; i < N; i += 2) {
            x_even.push_back(x[i]);
            x_odd.push_back(x[i+1]);
        }
        
        // Recursion
        cooley_tukey(x_even);
        cooley_tukey(x_odd);

        // Combine the even and odd arrays
        for (size_t j = 0; j < N/2; ++j) {
            Complex q = std::exp((-2 * M_PI * j) / N) * x_odd[j];
            x[j] = x_even[j] + q;
            x[j + (N/2)] = x_even[j] - q;
        }
    }
}

// Test the algorithm
int main() {
    // Test with values 1,2,3,4
    CArray x; 
    const Complex temp[] = {1.0,2.0,3.0,4.0};
    for (size_t i = 0; i < 4; i += 1) {
            x.push_back(temp[i]);
    }

    // Print test values
    std::cout << "Input sequence: " << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << x[i] << std::endl;
    }

    // Call algorithm
    cooley_tukey(x);

    std::cout << "FFT result: " << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << x[i] << std::endl;
    }

    return 0;
}