#include <iostream>
#include <complex>
#include <thread>
#include <cmath>
#include <vector>

//using namespace std;

typedef std::complex<double> Complex;

std::vector<Complex> fft_forward(const std::vector<Complex>& x) {
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
    std::vector<Complex> evenTransformed = fft_forward(even);
    std::vector<Complex> oddTransformed = fft_forward(odd);

    // Combine the even and odd arrays
    for (int k = 0; k < N / 2; k++) {
        Complex coeff (cos(2*M_PI*k/N), sin(2*M_PI*k/N));
        Complex q = coeff * oddTransformed[k];
        transformed[k] = evenTransformed[k] + q;
        transformed[k + (N/2)] = evenTransformed[k] - q;
    }

    return transformed;
}

std::vector<Complex> inverse_fft(const std::vector<Complex>& x) {
    std::vector<Complex> conj_x(x.size());

    for (int i = 0; i < x.size(); i++) {
        conj_x[i] = std::conj(x[i]);
    }

    std::vector<Complex> transformed = fft_forward(conj_x);

    int N = transformed.size();
    for (int i = 0; i < transformed.size(); i++) {
        transformed[i] = std::conj(transformed[i])/(Complex)N;
    }

    return transformed;   
}

void parallel_multiply(const std::vector<Complex>& a, const std::vector<Complex>& b, std::vector<Complex>& result, int start, int end) {
    for (int i = start; i < end; i++) {
        result[i] = a[i] * b[i];
    }
}

std::vector<Complex> parallel_fft_multiply(std::vector<Complex>& x, const std::vector<Complex>& coefficients) {
    // Final function: perform the forward FFT computation on the vector of input values and the coefficients, multiply  
    // the transformed arrays in parallel, and apply the inverse FFT to obtain the final result.

    int N = x.size();

    std::vector<Complex> res(N);
    // std::vector<Complex> transformed_coefficients(N);
    std::vector<Complex> multiplied_transform(N);
    std::vector<Complex> multiplied_res(N);

    // Perform FFT on both signal and coefficients
    res = fft_forward(x);

    // PARALLEL VERSION
    // int num_threads = 4;
    // std::vector<std::thread> workers(num_threads);
    // int block_size = N / num_threads;

    // NOTE: even without parallelization it does not work
    for (int i = 0; i < res.size(); i++) {
        multiplied_transform[i] = res[i] * coefficients[i];
    }

    // PARALLEL VERSION 
    // for (int t = 0; t < num_threads; t++) {
    //     int start = t * block_size;
    //     int end;
    //     if (t == num_threads - 1) {
    //         end = N;
    //     } else {
    //         end = start + block_size;
    //     }
    //     std::cout << "Thread " << t <<  "start: " <<  start << "end: " << end << std::endl;
    //     workers[t] = std::thread(&parallel_multiply, std::cref(res), std::cref(coefficients), std::ref(multiplied_transform), start, end);
    // }

    // // Wait for all threads to finish
    // for (int t = 0; t < num_threads; t++) {
    //     workers[t].join();
    // }

    // Perform inverse FFT on the multiplied transform
    multiplied_res = inverse_fft(multiplied_transform);

    return multiplied_res;
}

int main() {
    // Test with values 1,2,3,4
    std::vector<Complex> x; 
    const double temp1[] = {1.0, 2.0, 3.0, 4.0};
    for (size_t i = 0; i < 4; i += 1) {
            x.push_back(temp1[i]);
    } 
    std::vector<Complex> coefficients; // Coefficient array
    
    for (size_t i = 0; i < 4; i += 1) {
        double theta = -2.0 * M_PI * i / 4;
        coefficients.push_back(Complex(std::cos(theta), std::sin(theta)));
    }

    std::vector<Complex> result = parallel_fft_multiply(x, coefficients);

    for (Complex value : result) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
