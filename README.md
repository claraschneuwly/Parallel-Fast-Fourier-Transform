# Concurrent-Project: Fast Fourier Transform (FFT)
Course Project for CSE305 - Concurrent and Distributed Computing at Ecole Polytechnique

In this project, the goal will be to implement and then parallelize computation of the discrete Fourier transform (DFT). Although DFT looks “mathish”, it and its variants are at the heart of many important algorithms (ranging from jpeg and mp3 compressions to polynomial and number multiplication)1. For basics on DFT, we refer to the book by Jeff Erickson(1).

Here is the approach:
1. compute the DFT [x ̃1, . . . , x ̃N ] of the original data;
2. keep only few largest terms, and make others to be zeroes.

Tentative plan for the project is:

1. Parallelize DFT:<br>
(a) implement the standard Cooley-Tukey radix-2 algorithm (see sequential_fft.cpp); <br>
(b) take historical weather data (or any other data, but periodicity would be great) and test the quality of approximation by the algorithm;<br>
(c) implement a parallel radix-2 algorithm (following paper [3] (see parallel_fft_paper3.cpp) and paper [4] (see parallel_fft_paper4.cpp);<br>
(d) Perform a detailed comparison of the resulting algorithms and their versions on the weather data and generated benchmarks (test_weather.cpp and test_weather_parallel.cpp).

3. Used the implemented algorithms for some other application of FFT of your choice, and try to parallelize the whole code as much as possible. Potential applications are: <br>
• implement an analogue of Wolfram picture-curves, see5. <br>
• implement a polynomial multiplication algorithm for polynomials with integer coefficients (so you adapt your DFT from complex number to numbers modulo prime!) and compare with straightforward algorithms, see exercises to Chapter 32 of CLRS (see mult_poly.cpp).

[1] https://doi.org/10.1119/1.3254017  <br>
[2] http://jeffe.cs.illinois.edu/teaching/algorithms/notes/A-fft.pdf <br>
[3] https://doi.org/10.1016/0167-8191(90)90031-4 <br>
[4] https://doi.org/10.1109/SUPERC.1994.344263 <br>
