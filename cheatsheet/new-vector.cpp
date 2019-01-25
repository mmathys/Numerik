#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

int main()
{
    VectorXd test(5);
    test << 1, 2, 3, 4, 5;
    cout << test << endl;

    /**
     * eigen vector from std vector
     */
    const std::vector<double> c = {-0.90617984593, -0.53846931010, 0.0, 0.53846931010, 0.90617984593};
    VectorXd C = VectorXd::Map(c.data(), c.size());

    /**
     * std vector from eigen vector
     */
    C.data();

    return 0;
}