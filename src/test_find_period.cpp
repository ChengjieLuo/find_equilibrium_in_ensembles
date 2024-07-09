#include "gibbs_dynamics.h"

// Example usage of the Convolution function
int main(int argc, char *argv[])
{
    int N = 128;
    Vec phi(N, 0.0);
    for (int i = 0; i < N; i++)
    {
        phi[i] = std::sin(2.0 * M_PI / 32. * double(i));
        // print(phi[i]);
        // print("end phi");
    }
    int period = find_period(phi, N);

    Vec phiall=phi;
    print(phiall);
    rescale(phiall,period,N);
    print(phiall);




    return 0;
}