// Implementation of a parallel radix-2 algorithm folloing Paper (3) : "A parallel FFT on an MIMD machine" by Amir AVERBUCH
#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <thread>
#include <chrono>

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

void fft_aux(int begin, int end, Complex*  FFT_transformed, Complex*  x){
    // Aux function for parallel FFT computation
    int N = end - begin;
    Complex x_res[N], fft_transformed_res[N];
    for (int i = 0; i < N; i++){
        x_res[i] = x[i + begin];
    }

    fft(fft_transformed_res, x_res, N);

    for (int i = 0; i < N; i++){
        FFT_transformed[i + begin] = fft_transformed_res[i];
    }
}

void order_fft(Complex  FFT_ordered[], Complex  x[], int N, int num_threads){
    // Recursive function to order the input vector x in an ascending order for parallel FFT

    // Base case: if size is 1, set FFT_ordered to x
    if (num_threads == 1){
        for (int i = 0; i < N; i++){
            FFT_ordered[i] = x[i];
        }
        return;
    }

    // Divide x into even and odd indices
    Complex x_even[N/2];
    Complex x_odd[N/2];
    Complex evenOrdered[N/2];
    Complex oddOrdered[N/2];

    for (int i=0; i<N/2; i++){
        x_even[i] = x[2*i];
        x_odd[i] = x[2*i+1];
    }

    // Recursively order x_even and x_odd
    order_fft(evenOrdered, x_even, N/2, num_threads/2);
    order_fft(oddOrdered, x_odd, N/2, num_threads/2);

    // Combine the even and odd arrays
    for (int k = 0; k < N/2; k++){
        FFT_ordered[k] = evenOrdered[k];
        FFT_ordered[k + N / 2] = oddOrdered[k];
    }
}

void parallel_fft(Complex*  FFT_transformed, Complex* x, int N, int num_threads){
    // Main function for parallel FFT computation

    // 1) Order the input array x using order_fft function.
    Complex  FFT_ordered[N];
    order_fft(FFT_ordered, x, N, num_threads);

    // Creates threads to compute the FFT in parallel
    int block_size = N / num_threads;
    std::vector<std::thread> workers(num_threads - 1);

    int start_block = 0;
    
    for (int i = 0; i < num_threads-1; i++){
        int end_block = start_block + block_size;
        workers[i] = std::thread(&fft_aux, start_block, end_block, &FFT_transformed[0], &FFT_ordered[0]);
        start_block = end_block;
    }

    // Join threads
    for (int i = 0; i < num_threads - 1; i++){
        workers[i].join();
    }

    fft_aux(start_block, start_block + block_size, &FFT_transformed[0], &FFT_ordered[0]);

    // Combine the results from all threads to get the final FFT result
    int m = N / num_threads * 2;
    for (int t = 0; t < log2(num_threads); t++){
        int begin = 0;
        int end = m;
        while (end <= N){
            Complex  evenTransformed[m/2], oddTransformed[m/2];
            for (int i = 0; i < m/2; i++){
                evenTransformed[i] = FFT_transformed[begin + i];
                oddTransformed[i] = FFT_transformed[begin + i+ m/2];
            }

            for (int i = 0; i < m/2; i++){
                Complex factor = std::polar(1.0, -2 * M_PI * i / N) * oddTransformed[i];
                FFT_transformed[begin+i] = evenTransformed[i] + factor;
                FFT_transformed[begin+i + m/2] = evenTransformed[i] - factor;
            }

            begin += m;
            end += m;
        }
        m *= 2;
    }

}

int main(){
    int N = 8;

    int num_thread = 2;

    Complex p0[8]{Complex(1,0),Complex(2,0),Complex(3,0),Complex(4,0), Complex(5,0),Complex(6,0),Complex(7,0),Complex(8,0)};

    Complex FFT_transformed[N],x[N],pFFT_transformed[N];

    for (int i=0; i<N; i++){
        x[i] = p0[i];
    }

    Complex FFT_transformed[N],pFFT_transformed[N];
    
    auto start = std::chrono::steady_clock::now();
    fft(FFT_transformed, x, N);
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for fft is " << elapsed << " microseconds" << std::endl;

    // Print results
    std::cout << "FFT result: " << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << FFT_transformed[i] << std::endl;
    }

    start = std::chrono::steady_clock::now();
    parallel_fft(pFFT_transformed, x, N, num_thread);
    finish = std::chrono::steady_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for pfft is " << elapsed << " microseconds" << std::endl;

    // Print results
    std::cout << "pFFT result: " << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << pFFT_transformed[i] << std::endl;
    }
    
    return 0;
}
