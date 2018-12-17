#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
    VectorXd test(5);
    test << 1, 2, 3, 4, 5;
    cout << test << endl;

    return 0;
}