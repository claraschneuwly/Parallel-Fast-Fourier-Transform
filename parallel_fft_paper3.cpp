// Implementation of a parallel radix-2 algorithm folloing Paper (3) : "A parallel FFT on an MIMD machine" by Amir AVERBUCH

Complex bit_reversal(double x, int n) {
    Complex res[n]; 
    for (int i=0; i<n/2; i++) {
            res[i] = x[n-i-1]; 
    } 
    return res; 
}
