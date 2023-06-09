#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <thread>
#include <fstream>

#define TEST 0

typedef std::complex<double> Complex;

void fft_rec(std::vector<Complex> &FFT_transformed, std::vector<Complex> &x, int N){
    // Recursive function for computing the FFT

    // Trivial size-1 DFT base case: set FFT_transformed to x
    if (N == 1){
        for (int i = 0; i < N; i++){
            FFT_transformed[i] = x[i];
        }
        return;
    }

    // Divide x into even and odd indices
    std::vector<Complex> x_even(N/2), x_odd(N/2), evenTransformed(N/2), oddTransformed(N/2);
    // Complex x_even[N/2];
    // Complex x_odd[N/2];
    // Complex evenTransformed[N/2];
    // Complex oddTransformed[N/2];

    for (int i=0; i<N/2; i++){
        x_even[i] = x[2*i];
        x_odd[i] = x[2*i+1];
    }

    // Recursion
    fft_rec(evenTransformed, x_even, N/2);
    fft_rec(oddTransformed, x_odd, N/2);

    // Combine the even and odd arrays
    for (int k = 0; k < N / 2; k++) {
        Complex factor = std::polar(1.0, -2.0 * M_PI * k / N) * oddTransformed[k];
        FFT_transformed[k] = evenTransformed[k] + factor;
        FFT_transformed[k + N / 2] = evenTransformed[k] - factor;
    }
}

unsigned int bit_reversal(unsigned int n, int s){
    unsigned int res = 0;
    for (int i = 0; i < s; i++){
        res <<= 1;
        res |= (n & 1);
        n >>= 1;
    }
    return res;
}


void butterfly(std::vector<Complex> &x, int start, int num_blocks, int len) {
    int N = x.size();
    double ang = 2 * M_PI / len;
    Complex wlen(cos(ang), sin(ang));
    for (int i = start; i < num_blocks; i += len) {
        Complex W(1);
        for (int j = 0; j < len / 2 && i + j + len/2 < N; j++) {
            Complex u = x[i + j];
            Complex v = x[i + j + len / 2] * W;
            x[i + j] = u + v;
            x[i + j + len / 2] = u - v;
            W *= wlen;
        }
    }
}

void p_fft_iter(std::vector<Complex>& x, int original_num_threads) {
    int N = x.size();
    int log2N = log2(N);

    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            std::swap(x[i], x[j]);
    }

    for (int len = 2; len <= N; len <<= 1) {
        int num_threads = original_num_threads;
        // Execute each level of butterfly operations in parallel

        int num_blocks = N / len;
        if (N % len != 0)
            ++num_blocks;

        if (num_threads > num_blocks) //if more threads than blocks
            num_threads = num_blocks;

        int blocks_per_thread = num_blocks / N ;
        int threads_with_extra_blocks = num_blocks - blocks_per_thread * num_threads;
        std::vector<std::thread> workers(num_threads - 1);

        int start = 0;
        //launch threads at each iteration level
        for (int i = 0; i < num_threads-1; i++){
            std::cout << "start: " << start << std::endl;
            if (threads_with_extra_blocks > 0){
                workers[i] = std::thread(&butterfly, std::ref(x), start, blocks_per_thread + 1, len);
                start = start + (blocks_per_thread + 1) * len; 
                threads_with_extra_blocks--;
            }
            else{
                std::cout << "BLOCKS " << blocks_per_thread << " " << len << std::endl;
                workers[i] = std::thread(&butterfly, std::ref(x), start, blocks_per_thread, len);
                start = start + blocks_per_thread*len;
            }
        }
        std::cout << "start : " << start << std::endl;
        if (threads_with_extra_blocks > 0)
            butterfly(x, start, blocks_per_thread+1, len);
        else   
            butterfly(x, start, blocks_per_thread, len);


        //join threads
        for (int i = 0; i < num_threads - 1; i++){
            workers[i].join();
        }
        std::cout << x[0] << std::endl;
    }
}


void fft_iter(std::vector<Complex> &x){ 
    int N = x.size();
    int log2N = log2(N);
    
    for (int i = 1, j = 0; i < N; i++) { //NOTE: change this code if possible 
        int bit = N >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            std::swap(x[i], x[j]);
    }

    for (int len = 2; len <= N; len <<= 1) {
        std::cout << "len: " << len << std::endl;
        double ang = 2 * M_PI / len;
        Complex wlen(cos(ang), sin(ang));
        for (int i = 0; i < N; i += len) {
            Complex W(1);
            for (int j = 0; j < len / 2; j++) {
                Complex u = x[i+j];
                Complex v = x[i+j+len/2] * W;
                x[i+j] = u + v;
                x[i+j+len/2] = u - v;
                W *= wlen;
            }
        }
        std::cout << x[0] << std::endl;
    }
}

void fft_iter_aux(int begin, int end, std::vector<Complex> &FFT_transformed, std::vector<Complex> &x, int len){
    // Aux function for parallel FFT computation
    size_t chunk_size = end - begin;
    std::vector<Complex> x_res(chunk_size), fft_transformed_res(chunk_size);
    // Complex x_res[chunk_size], fft_transformed_res[chunk_size];

    for (size_t i = 0; i < chunk_size; i++){
        x_res[i] = x[i + begin];
    }

    fft_rec(fft_transformed_res, x_res, chunk_size);

    for (int i = 0; i < chunk_size; i++){
        FFT_transformed[i + begin] = fft_transformed_res[i];
    }
}

void p_ifft(std::vector<Complex> &FFT_inverse, std::vector<Complex>& x, int num_threads){
    int N = x.size();

    for (auto &val : x) {val = std::conj(val);}

    p_fft_iter(x, num_threads);

    for (size_t i = 0; i < N; i++) {
        FFT_inverse[i] = std::conj(x[i])/(Complex)N;
    }
}

void inverse_fft(std::vector<Complex> &FFT_inverse, std::vector<Complex>& x, bool recursive) {
    int N = x.size();

    for (auto &val : x) {val = std::conj(val);}

    if (recursive){
        fft_rec(FFT_inverse, x, N);

        for (size_t i = 0; i < N; i++) {
            FFT_inverse[i] = std::conj(FFT_inverse[i])/(Complex)N;
        }
    }
    else {
        fft_iter(x);

        for (size_t i = 0; i < N; i++) {
            FFT_inverse[i] = std::conj(x[i])/(Complex)N;
        }
    }
}

void fft_aux(int begin, int end, std::vector<Complex> &FFT_transformed, std::vector<Complex> &x){
    // Aux function for parallel FFT computation
    size_t chunk_size = end - begin;
    std::vector<Complex> x_res(chunk_size), fft_transformed_res(chunk_size);
    // Complex x_res[chunk_size], fft_transformed_res[chunk_size];

    for (size_t i = 0; i < chunk_size; i++){
        x_res[i] = x[i + begin];
    }

    fft_rec(fft_transformed_res, x_res, chunk_size);

    for (int i = 0; i < chunk_size; i++){
        FFT_transformed[i + begin] = fft_transformed_res[i];
    }
}


void order_fft(std::vector<Complex> & FFT_ordered, std::vector<Complex> &x, int N, int num_threads){
    // Recursive function to order the input vector x in an ascending order for parallel FFT

    // Base case: if size is 1, set FFT_ordered to x
    if (num_threads == 1){
        for (int i = 0; i < N; i++){
            FFT_ordered[i] = x[i];
        }
        return;
    }

    // Divide x into even and odd indices
    std::vector<Complex> x_even(N/2), x_odd(N/2), evenOrdered(N/2), oddOrdered(N/2);
    // Complex x_even[N/2];
    // Complex x_odd[N/2];
    // Complex evenOrdered[N/2];
    // Complex oddOrdered[N/2];

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

void parallel_fft(std::vector<Complex> &FFT_transformed, std::vector<Complex> &x, int N, int num_threads){
    int log2N = log2(N);
    

    // Order the input array x using order_fft function.
    std::vector<Complex> FFT_ordered(N);
    order_fft(FFT_ordered, x, N, num_threads);

    // Create threads to compute FFT in parallel
    int block_size = N / num_threads;
    std::vector<std::thread> workers(num_threads - 1);

    int start_block = 0;

    for (int i = 0; i < num_threads-1; i++){
        int end_block = start_block + block_size;
        workers[i] = std::thread(&fft_aux, start_block, end_block, std::ref(FFT_transformed), std::ref(FFT_ordered));
        start_block = end_block;
    }
    fft_aux(start_block, start_block + block_size, FFT_transformed, FFT_ordered);

    for (int i = 0; i < num_threads-1; i++){
        workers[i].join();
    }

    // join threads


    int m = N / num_threads * 2; //NOTE: could be improved
    for (int t = 0; t < log2(num_threads); t++){
        int begin = 0;
        int end = m;
        while (end <= N){
            std::vector<Complex> evenTransformed(m/2), oddTransformed(m/2);
            // Complex  evenTransformed[m/2], oddTransformed[m/2];
            for (int i = 0; i < m/2; i++){
                evenTransformed[i] = FFT_transformed[begin + i];
                oddTransformed[i] = FFT_transformed[begin + i+ m/2];
            }

            for (int i = 0; i < m/2; i++){
                Complex factor = std::polar(1.0, -2 * M_PI * i / N) * oddTransformed[i];
                FFT_transformed[begin+i] = evenTransformed[i] + factor ;
                FFT_transformed[begin+i + m/2] = evenTransformed[i] - factor;
            }

            begin += m;
            end += m;
        }
        m *= 2;
    }
}



// Test the algorithm
int main() {
    #if TEST == 0 //implement custom code 

    int N = 8;

    int num_thread = 4;

    // Complex input[2] = {Complex(1,0), Complex(2,0)};
    Complex input[8] = {Complex(1,0),Complex(2,0),Complex(3,0),Complex(4,0), Complex(5,0),Complex(6,0),Complex(7,0),Complex(8,0)};

    std::vector<Complex> FFT_transformed(N), pFFT_transformed(N);
    
    std::vector<Complex> x(N);
    // Complex FFT_transformed[N], x[N], pFFT_transformed[N];

    for (int i=0; i<N; i++){
        x[i] = input[i];
    }

    auto start = std::chrono::steady_clock::now();
    p_fft_iter(x, num_thread);
    // fft_iter(*x);
    std::cout << "p_fft_iter done" << std::endl;
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for fft is " << elapsed << " microseconds" << std::endl;

    // Print results
    std::cout << "FFT result: " << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << x[i] << std::endl;
    }

    // start = std::chrono::steady_clock::now();
    // parallel_fft(pFFT_transformed, x, N, num_thread);
    // finish = std::chrono::steady_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    // std::cout << "Time for parallel fft is " << elapsed << " microseconds" << std::endl;

    // // Print results
    // std::cout << "Parallel FFT result: " << std::endl;
    // for (int i = 0; i < 8; ++i) {
    //     std::cout << pFFT_transformed[i] << std::endl;
    // }
    
    // start = std::chrono::steady_clock::now();
    // p_fft_iter(x, 4);
    // // parallel_fft(pFFT_transformed, x, N, num_thread);
    // finish = std::chrono::steady_clock::now();
    // elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    // std::cout << "Time for sequential fft is " << elapsed << " microseconds" << std::endl;
    // std::cout << "Sequential FFT result: " << std::endl;
    // for (int i = 0; i < 8; ++i) {
    //     std::cout << x[i] << std::endl;
    // }


    #elif TEST == 1 // Test on weather data
    // std::cout << 8/3 << std::endl;
    // Open historical weather data
    std::ifstream inputFile("test_data/weather_data_clean.csv");
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return 1;
    }
    // Store data in a complex array
    std::vector<Complex> weather_data, FFT_transformed;

    int num_frequencies = 100;

    std::string line;

    while (std::getline(inputFile, line)) {
        double number = std::stoi(line);
        weather_data.push_back(Complex(number, 0));
        // FFT_transformed.push_back(Complex(number, 0));
    }
    inputFile.close();
    std::vector<Complex> FFT_inverse(weather_data.size());


    // Print datapoints
    std::cout << "FFT first 6 data-points: " << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << weather_data[i] << std::endl;
    }

    // Compute FFT of weather - sequential (iterative)
    auto start = std::chrono::steady_clock::now();
    fft_iter(weather_data);
    if (weather_data.size() > num_frequencies) { // iterative
        std::fill(weather_data.begin() + num_frequencies, weather_data.end() + 1, 0); 
    }
    std::cout << weather_data.size() << std::endl;

    // inverse_fft(FFT_inverse, weather_data, false); //sequential
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time for fft is " << elapsed << " microseconds" << std::endl;


    // fft_rec(FFT_transformed, weather_data, weather_data.size());
    // parallel_fft(FFT_transformed, weather_data, weather_data.size(), 4);
    
    // Truncate the DFT 
    // if (FFT_transformed.size() > num_frequencies) { // recursive 
    //     std::fill(FFT_transformed.begin() + num_frequencies, FFT_transformed.end(), Complex(0, 0)); 
    // }


    // Compute inverse FFT
    
    // inverse_fft(FFT_transformed, weather_data, true);
    
    // Create output file
    std::ofstream outputFile("test_data/output_weather_test_simple.csv");
    if (!outputFile.is_open()) {
        std::cout << "Failed to open the file." << std::endl;
        return 1;
    }

    // Write the output data to the CSV output file
    for (const auto& n : FFT_inverse) {
        outputFile << n << std::endl;
    }

    outputFile.close();

    // Print output
    std::cout << "FFT first 6 result: " << std::endl;
    for (int i = 0; i < 6; ++i) {
        std::cout << FFT_inverse[i] << std::endl;
    }
   
    #endif
    return 0;
}
