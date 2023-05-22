// Test the Cooley-Tukey algorithm with historical CO2 data (weather_data_celan.csv) and write results to output_weather_test.csv

#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <string>

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
    cooley_tukey(numbers);

    // Open output file
    std::ofstream outputFile("output_weather_test.csv");
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
