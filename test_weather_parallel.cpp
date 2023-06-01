#include <iostream>
#include <complex>
#include <thread>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES
#include <fstream>
#include <string>


//using namespace std;

typedef std::complex<double> Complex;
typedef std::vector<Complex> CArray;

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
        // Open historical weather data
    std::ifstream inputFile("weather_data_clean.csv");
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return 1;
    }
    // Put data in a CArray
    CArray numbers;
    std::string line;
    while (std::getline(inputFile, line)) {
        int number = std::stoi(line);
        numbers.push_back(number);
    }

    inputFile.close();


    std::cout << "Input sequence first 6 numbers: " << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << numbers[i] << std::endl;
    }
    // Call algorithm
    fft_parallel(numbers);

    // Open output file
    std::ofstream outputFile("output_weather_par_test.csv");
    if (!outputFile.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return 1;
    }

    // Write the output data to the CSV output file
    for (const auto& n : numbers) {
        outputFile << n << std::endl;
    }

    outputFile.close();

    // Print output
    std::cout << "FFT first 6 result: " << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << numbers[i] << std::endl;
    }

    return 0;
}
