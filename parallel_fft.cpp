// Parallel Cooley-Tukey radix-2 algorithm

#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

typedef Complex Complex;

// Compute the DFT of x#include <complex>
#include <cmath>
#include <vector>
#include <thread>

typedef std::complex<double> Complex;

void cooley_tukey_parallel(std::vector<Complex>& x) {
    const size_t N = x.size();

    // Trivial size-1 DFT base case
    if (N == 1) {
        return;
    }

    // Divide x into even and odd indices
    std::vector<Complex> x_even;
    std::vector<Complex> x_odd;
    for (size_t i = 0; i < N; i += 2) {
        x_even.push_back(x[i]);
        x_odd.push_back(x[i + 1]);
    }

    // Create threads for the recursive calls on even and odd subarrays
    std::thread even_thread(cooley_tukey_parallel, std::ref(x_even));
    std::thread odd_thread(cooley_tukey_parallel, std::ref(x_odd));

    // Wait for the threads to finish
    even_thread.join();
    odd_thread.join();

    // Combine the even and odd arrays
    for (size_t j = 0; j < N / 2; ++j) {
        Complex W(std::cos(2 * M_PI * j / N), std::sin(2 * M_PI * j / N));
        x[j] = x_even[j] + W * x_odd[j];
        x[j + (N / 2)] = x_even[j] - W * x_odd[j];
    }
}


// Test the algorithm
int main() {
    // Test with values 1,2,3,4
    std::vector<Complex> x; 
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
    
    // Print results
    std::cout << "FFT result: " << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << x[i] << std::endl;
    }

    return 0;
}
