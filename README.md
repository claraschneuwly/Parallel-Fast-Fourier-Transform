# Concurrent-Project: Fast Fourier Transform (FFT)

In this project, the goal will be to implement and then parallelize computation of the discrete Fourier transform (DFT). Although DFT looks “mathish”, it and its variants are at the heart of many important algorithms (ranging from jpeg and mp3 compressions to polynomial and number multiplication)1. For basics on DFT, we refer to the book by Jeff Erickson2.

Here is the approach:
1. compute the DFT [x ̃1, . . . , x ̃N ] of the original data;
2. keep only few largest terms, and make others to be zeroes.

Tentative plan for the project is:

1. Parallelize DFT:
(a) implement the standard Cooley-Tukey radix-2 algorithm (see sequential_fft.cpp);
(b) take historical weather data (or any other data, but periodicity would be great) and test the quality of approximation by the algorithm;
(c) implement a parallel radix-2 algorithm (following this paper: https://doi.org/10.1016/0167-8191(90)90031- 4 and this paper: https://doi.org/10.1109/SUPERC.1994.344263);
(d) Perform a detailed comparison of the resulting algorithms and their versions on the
weather data and generated benchmarks.

2. Used the implemented algorithms for some other application of FFT of your choice, and try to parallelize the whole code as much as possible. Potential applications are:
• implement an analogue of Wolfram picture-curves, see5.
• implement a polynomial multiplication algorithm for polynomials with integer coeffi- cients (so you adapt your DFT from complex number to numbers modulo prime!) and compare with straightforward algorithms, see exercises to Chapter 32 of CLRS.
