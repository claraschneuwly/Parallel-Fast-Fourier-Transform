#include <iostream>
#include <complex>
#include <thread>
#include <cmath>
#include <vector>

//using namespace std;

typedef std::complex<double> Complex;

std::vector<Complex> fft_parallel(const std::vector<Complex>& x) {
    int N = x.size();

    // Trivial size-1 DFT base case
    if (N <= 1) {
        return x;
    }

    // Divide x into even and odd indices
    std::vector<Complex> even;
    std::vector<Complex> odd;

    for (size_t i = 0; i < N; i += 2) {
        even.push_back(x[i]);
        odd.push_back(x[i + 1]);
    }

    // Recursion
    std::vector<Complex> transformed(N);
    std::vector<Complex> evenTransformed = fft_parallel(even);
    std::vector<Complex> oddTransformed = fft_parallel(odd);

    // Combine the even and odd arrays
    for (int k = 0; k < N / 2; k++) {
        Complex factor = std::polar(1.0, -2.0 * M_PI * k / N) * oddTransformed[k];
        transformed[k] = evenTransformed[k] + factor;
        transformed[k + N / 2] = evenTransformed[k] - factor;
    }

    return transformed;
}

std::vector<Complex> inverse_fft_parallel(const std::vector<Complex>& x) {
    int N = x.size();

    // Trivial size-1 DFT base case
    if (N <= 1) {
        return x;
    }

    // Divide x into even and odd indices
    std::vector<Complex> even;
    std::vector<Complex> odd;

    for (size_t i = 0; i < N; i += 2) {
        even.push_back(x[i]);
        odd.push_back(x[i + 1]);
    }

    // Recursion
    std::vector<Complex> transformed(N);
    std::vector<Complex> evenTransformed = inverse_fft_parallel(even);
    std::vector<Complex> oddTransformed = inverse_fft_parallel(odd);

    // Combine the even and odd arrays
    for (int k = 0; k < N / 2; k++) {
        Complex factor = std::polar(1.0, 2.0 * M_PI * k / N) * oddTransformed[k];
        transformed[k] = evenTransformed[k] + factor;
        transformed[k + N / 2] = evenTransformed[k] - factor;
    }

    return transformed;
}

void parallel_multiply(const std::vector<Complex>& a, const std::vector<Complex>& b, std::vector<Complex>& result, int start, int end) {
    for (int i = start; i < end; i++) {
        result[i] = a[i] * b[i];
    }
}

std::vector<double> parallel_fft_multiply(const std::vector<double>& res, const std::vector<double>& coefficients) {
    // Final function: perform the forward FFT computation on the vector of input values and the coefficients, multiply  
    // the transformed arrays in parallel, and apply the inverse FFT to obtain the final result.

    int N = res.size();

    std::vector<Complex> transformed_res(N);
    std::vector<Complex> transformed_coefficients(N);
    std::vector<Complex> multiplied_transform(N);
    std::vector<Complex> multiplied_res(N);

    // Convert input signals to complex numbers
    for (int i = 0; i < N; i++) {
        transformed_res[i] = Complex(res[i], 0.0);
        transformed_coefficients[i] = Complex(coefficients[i], 0.0);
    }

    // Perform FFT on both signal and coefficients
    transformed_res = fft_parallel(transformed_res);
    transformed_coefficients = fft_parallel(transformed_coefficients);

    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers(num_threads);
    int block_size = N / num_threads;

    // Multiply transformed arrays in parallel
    for (int t = 0; t < num_threads; t++) {
        int start = t * block_size;
        int end;
        if (t == num_threads - 1) {
            end = N;
        } else {
            end = start + block_size;
        }
        workers[t] = std::thread(&parallel_multiply, std::ref(transformed_res), std::ref(transformed_coefficients), std::ref(multiplied_transform), start, end);
    }

    // Wait for all threads to finish
    for (int t = 0; t < num_threads; t++) {
        workers[t].join();
    }

    // Perform inverse FFT on the multiplied transform
    multiplied_res = inverse_fft_parallel(multiplied_transform);

    // Extract the real part of the result
    std::vector<double> result(N);
    for (int i = 0; i < N; i++) {
        //result[i] = multiplied_res[i].real();
        result[i] = std::real(multiplied_res[i]);
    }

    return result;
}

int main() {
    // Test with values 1,2,3,4
    std::vector<double> signal; 
    const double temp1[] = {1.0,2.0,3.0,4.0};
    for (size_t i = 0; i < 4; i += 1) {
            signal.push_back(temp1[i]);
    } 
    std::vector<double> coefficients; // Coefficient array
    const double temp2[] = {1, 2, 1, 0};
    for (size_t i = 0; i < 4; i += 1) {
            coefficients.push_back(temp2[i]);
    }

    std::vector<double> result = parallel_fft_multiply(signal, coefficients);

    for (double value : result) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
