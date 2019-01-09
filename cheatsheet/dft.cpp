// is in O(n^2)
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <iostream>

using namespace Eigen;
using namespace std;

int main()
{
}

VectorXcd pconvfft(const VectorXcd &u, const VectorXcd &x)
{                           // complex v e c t o r s
    Eigen::FFT<double> fft; // Eigen b uil t −i n d a t a type f o r FFT.
    VectorXcd tmp = (fft.fwd(u)).cwiseProduct(fft.fwd(x));
    return fft.inv(tmp);
}