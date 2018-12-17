#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

int main();
void debug(string);

int main()
{
    VectorXd test(5);
    test << 1, 2, 3, 4, 5;

    int localRef = 0;
    test = test.unaryExpr([&localRef](double item) {
        localRef++;
        return item + 1;
    });

    cout << "local ref: " << localRef << endl;
    cout << "vector:" << endl;
    cout << test << endl;

    debug("something");

    return 0;
}

void debug(string msg)
{
    cout << "DEBUG" << endl;
    cout << msg << endl;
    cout << "END DEBUG" << endl;
}