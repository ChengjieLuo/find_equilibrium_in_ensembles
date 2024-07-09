#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>
#include <iomanip>
#include <fstream>
#include <tuple>
#include <string>
#include <fftw3.h>

// Typedefs for clarity
typedef std::vector<double> Vec;
typedef std::vector<std::vector<double>> Vec2D;
typedef std::vector<std::vector<std::vector<double>>> Vec3D;

typedef std::complex<double> Complex;
typedef std::vector<Complex> ComplexArray;

void print(int a)
{
    std::cout << a << std::endl;
    return;
}

void print(double a)
{
    std::cout << std::setprecision(14);

    std::cout << a << std::endl;
    return;
}

void print(const char *a)
{
    std::cout << a << std::endl;
    return;
}

void print(char a[])
{
    std::cout << a << std::endl;
    return;
}

void print(const std::string a)
{
    std::cout << a << std::endl;
    return;
}

void print(Vec a)
{
    std::cout << std::setprecision(14);
    std::cout << "[ ";
    for (const auto &val : a)
    {
        std::cout << val << ", ";
    }
    std::cout << "]" << std::endl;
    return;
}

void print(Vec2D a)
{
    std::cout << std::setprecision(14);
    std::cout << "[\n";
    for (const auto &val : a)
    {
        std::cout << "[ ";
        for (const auto &val1 : val)
        {
            std::cout << val1 << ", ";
        }
        std::cout << "],";
        std::cout << std::endl;
    }
    std::cout << "]" << std::endl;
    return;
}

void print(Vec3D a)
{
    std::cout << std::setprecision(14);
    std::cout << "[\n";
    for (const auto &val : a)
    {
        std::cout << "[\n";
        for (const auto &val1 : val)
        {
            std::cout << "[ ";
            for (const auto &val2 : val1)
            {
                std::cout << val2 << ", ";
            }
            std::cout << "],";
            std::cout << std::endl;
        }
        std::cout << "]," << std::endl;
    }
    std::cout << "]" << std::endl;
    return;
}

// Function to perform 1D Fourier Transform using FFTW
ComplexArray FourierTransform1D(const ComplexArray &input)
{
    int N = input.size();
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    // Initialize input array
    for (int i = 0; i < N; ++i)
    {
        in[i][0] = input[i].real();
        in[i][1] = input[i].imag();
    }

    // Create plan for forward DFT
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the plan
    fftw_execute(plan);

    // Store the result in a vector
    ComplexArray result(N);
    for (int i = 0; i < N; ++i)
    {
        result[i] = Complex(out[i][0], out[i][1]);
    }

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}

// Function to perform 1D Inverse Fourier Transform using FFTW
ComplexArray InverseFourierTransform1D(const ComplexArray &input)
{
    int N = input.size();
    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    // Initialize input array
    for (int i = 0; i < N; ++i)
    {
        in[i][0] = input[i].real();
        in[i][1] = input[i].imag();
    }

    // Create plan for backward DFT
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the plan
    fftw_execute(plan);

    // Store the result in a vector
    ComplexArray result(N);
    for (int i = 0; i < N; ++i)
    {
        result[i] = Complex(out[i][0] / N, out[i][1] / N); // Normalize the result
    }

    // Clean up
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return result;
}

// Function to perform convolution using Fourier transforms
ComplexArray Convolution(const ComplexArray &X, const ComplexArray &K)
{
    int N = X.size();

    // Perform Fourier Transform of X
    ComplexArray X_fft = FourierTransform1D(X);

    // Element-wise multiplication of X_fft and K
    ComplexArray Y_fft(N);
    for (int i = 0; i < N; ++i)
    {
        Y_fft[i] = X_fft[i] * K[i];
    }

    // Perform Inverse Fourier Transform of the result
    ComplexArray Y = InverseFourierTransform1D(Y_fft);

    return Y;
}

Vec Convolution(const Vec &Xr, const Vec &Kr)
{
    int N = Xr.size();
    // print("N=");
    // print(N);
    ComplexArray X(N), K(N), Y(N);
    for (int i = 0; i < N; ++i)
    {
        X[i].real(Xr[i]);
        X[i].imag(0.0);
        K[i].real(Kr[i]);
        K[i].imag(0.0);
    }
    Y = Convolution(X, K);

    Vec Yr(N, 0);
    for (int i = 0; i < N; ++i)
    {
        Yr[i] = Y[i].real();
    }
    return Yr;
}

// Function to calculate frequency values
std::vector<double> calculateFrequencies(double samplingRate, int numDataPoints)
{
    double frequencyResolution = samplingRate / numDataPoints;

    std::vector<double> frequencies(numDataPoints);
    if (numDataPoints % 2 == 0)
    {
        for (int i = 0; i < numDataPoints / 2; ++i)
        {
            frequencies[i] = i * frequencyResolution;
        }
        for (int i = numDataPoints / 2; i < numDataPoints; ++i)
        {
            frequencies[i] = (i - numDataPoints) * frequencyResolution;
        }
    }
    else
    {
        for (int i = 0; i < (numDataPoints + 1) / 2; ++i)
        {
            frequencies[i] = i * frequencyResolution;
        }
        for (int i = (numDataPoints + 1) / 2; i < numDataPoints; ++i)
        {
            frequencies[i] = (i - numDataPoints) * frequencyResolution;
        }
    }

    return frequencies;
}

int find_period(const Vec &phi, int num_coord)
{
    int N = phi.size();
    if (num_coord != N)
    {
        print("num_coord=");
        print(num_coord);
        print("N=");
        print(N);
        print("error in the shape of phi");
    }
    ComplexArray X(N), Y(N);
    for (int i = 0; i < N; ++i)
    {
        X[i].real(phi[i]);
        X[i].imag(0.0);
        // K[i].real(Kr[i]);
        // K[i].imag(0.0);
    }

    Y = FourierTransform1D(X);
    Vec absY = Vec(N, 0.0);
    for (int i = 0; i < N; ++i)
    {

        absY[i] = std::sqrt(Y[i].real() * Y[i].real() + Y[i].imag() * Y[i].imag());
        // print(absY[i]);
        // print(Y[i].real());
        // print(Y[i].imag());
    }

    int index_max1 = 0;
    int index_max2 = 0;
    double maxabsY1 = absY[index_max1];
    double maxabsY2 = absY[index_max2];

    for (int i = 0; i < N; ++i)
    {
        if (absY[i] >= maxabsY1)
        {
            maxabsY1 = absY[i];
            index_max2 = index_max1;
            index_max1 = i;
        }
    }

    Vec fre = calculateFrequencies(1.0, num_coord);
    // print("fre=");
    // print(fre);
    // print("index2=");
    // print(index_max2);

    double periods = N * std::abs(fre[index_max2]);
    // print("periods=");
    // print(periods);
    return int(periods + 0.0000001);
}

void rescale(Vec &phi, int periods, int num_coord)
{
    Vec oldphi = phi;

    int Nold = num_coord / periods;

    for (int icoord = 0; icoord < Nold; icoord++)
    {
        for (int index = 0; index < periods; index++)
        {
            phi[periods * icoord + index] = oldphi[icoord] + (oldphi[icoord + 1] - oldphi[icoord]) / periods * index;
        }
    }
}

// // Example input arrays
// ComplexArray X = {Complex(1, 0), Complex(1, 0), Complex(1, 0), Complex(1, 0),
//                   Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};
// ComplexArray K = {Complex(1, 0), Complex(0.5, 0), Complex(0.25, 0), Complex(0, 0),
//                   Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};

// // Assuming K is already in the frequency domain, otherwise you would need to Fourier Transform K
// // K = FourierTransform1D(K);

// // Perform convolution
// ComplexArray Y = Convolution(X, K);

// std::cout << "X :" << std::endl;
// for (const auto& val : X) {
//     std::cout << val << std::endl;
// }
// std::cout << "K :" << std::endl;
// for (const auto& val : K) {
//     std::cout << val << std::endl;
// }
// // Print the result
// std::cout << "Convolution result:" << std::endl;
// for (const auto& val : Y) {
//     std::cout << val << std::endl;
// }

// double samplingRate = 1;
// int numDataPoints = 8;
// std::vector<double> freq = calculateFrequencies(samplingRate, numDataPoints);
// for (const auto &val : freq)
// {
//     std::cout << val << std::endl;
// }

// #include <iostream>
// #include <fftw3.h>

// int main() {
//     const int N = 8; // Size of the array

//     // Allocate input and output arrays
//     fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//     fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

//     // Initialize input array with some data
//     for (int i = 0; i < N; ++i) {
//         in[i][0] = i;   // Real part
//         in[i][1] = 0.0; // Imaginary part
//     }

//     // Create a plan for forward DFT
//     fftw_plan plan_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

//     // Execute the plan
//     fftw_execute(plan_forward);

//     // Print the output of forward transform
//     std::cout << "Forward Transform:" << std::endl;
//     for (int i = 0; i < N; ++i) {
//         std::cout << "(" << out[i][0] << ", " << out[i][1] << ")" << std::endl;
//     }

//     // Create a plan for backward (inverse) DFT
//     fftw_plan plan_backward = fftw_plan_dft_1d(N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

//     // Execute the plan
//     fftw_execute(plan_backward);

//     // Normalize the output of the inverse transform
//     for (int i = 0; i < N; ++i) {
//         in[i][0] /= N;
//         in[i][1] /= N;
//     }

//     // Print the output of backward transform
//     std::cout << "Backward Transform (after normalization):" << std::endl;
//     for (int i = 0; i < N; ++i) {
//         std::cout << "(" << in[i][0] << ", " << in[i][1] << ")" << std::endl;
//     }

//     // Clean up
//     fftw_destroy_plan(plan_forward);
//     fftw_destroy_plan(plan_backward);
//     fftw_free(in);
//     fftw_free(out);

//     return 0;
// }