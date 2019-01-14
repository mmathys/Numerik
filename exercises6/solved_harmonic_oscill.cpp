#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <iostream>
#include "writer.hpp"
#include <assert.h>

/// Uses the explicit Euler method to compute y from time 0 to time T
/// where y is a 2x1 vector solving the linear system of ODEs as in the exercise
///
/// @param[out] yT at the end of the call, this will have vNext
/// @param[in] y0 the initial conditions
/// @param[in] zeta the damping factor (see exercise)
/// @param[in] h the step size
/// @param[in] T the final time at which to compute the solution.
///
/// The first component of y (the position) will be stored to y1, the second component (the velocity) to y2. The i-th entry of y1 (resp. y2) will contain the first (resp. second) component of y at time i*h.
///

using namespace std;
using namespace Eigen;

//----------------explicitEulerBegin----------------
void explicitEuler(std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &time,
				   const Eigen::Vector2d &y0,
				   double zeta, double h, double T)
{
	int steps = T / h;
	int size = steps + 1;
	VectorXd Time = VectorXd::LinSpaced(size, 0, T);
	assert(Time.size() == size);
	MatrixXd Y(size, 2);
	Y(0, 0) = y0(0);
	Y(0, 1) = y0(1);
	MatrixXd A(2, 2);
	A << 0, 1, -1, -2 * zeta;

	for (int i = 1; i < size; i++)
	{
		VectorXd res = Y.row(i - 1).transpose() + h * A * Y.row(i - 1).transpose();
		Y.row(i) = res.transpose();
	}

	// convert eigen vectors into std vectors
	for (int i = 0; i < size; i++)
	{
		time.push_back(Time(i));
		y1.push_back(Y(i, 0));
		y2.push_back(Y(i, 1));
	}
}
//----------------explicitEulerEnd----------------

// Implements the implicit Euler. Analogous to explicit Euler, same input and output parameters
//----------------implicitEulerBegin----------------
void implicitEuler(std::vector<double> &y1, std::vector<double> &y2, std::vector<double> &time,
				   const Eigen::Vector2d &y0,
				   double zeta, double h, double T)
{
	int steps = T / h;
	int size = steps + 1;
	VectorXd Time = VectorXd::LinSpaced(size, 0, T);
	assert(Time.size() == size);
	MatrixXd Y(size, 2);
	Y(0, 0) = y0(0);
	Y(0, 1) = y0(1);
	MatrixXd A(2, 2);
	A << 0, 1, -1, -2 * zeta;

	for (int i = 1; i < size; i++)
	{
		MatrixXd B = MatrixXd::Identity(2, 2) - h * A;
		// B * y_k+1 = y_k
		VectorXd prev = Y.row(i - 1).transpose();
		VectorXd res = B.fullPivLu().solve(prev);
		Y.row(i) = res.transpose();
	}

	time.resize(0);
	y1.resize(0);
	y2.resize(0);
	// convert eigen vectors into std vectors
	for (int i = 0; i < size; i++)
	{
		time.push_back(Time(i));
		y1.push_back(Y(i, 0));
		y2.push_back(Y(i, 1));
	}
}
//----------------implicitEulerEnd----------------

//----------------energyBegin----------------
// Energy computation given the velocity. Assume the energy vector to be already initialized with the correct size.
void Energy(const std::vector<double> &v, std::vector<double> &energy)
{
	assert(v.size() == energy.size());
	// TODO: Task (e)
	// ...
	// ...
}
//----------------energyEnd----------------

int main()
{

	double T = 20.0;
	double h = 0.5; // Change this for explicit / implicit time stepping comparison
	const Eigen::Vector2d y0(1, 0);
	double zeta = 0.2;
	std::vector<double> y1;
	std::vector<double> y2;
	std::vector<double> time;
	explicitEuler(y1, y2, time, y0, zeta, h, T);
	writeToFile("position_expl.txt", y1);
	writeToFile("velocity_expl.txt", y2);
	writeToFile("time_expl.txt", time);
	std::vector<double> energy(y2.size());
	//Energy(y2, energy);
	writeToFile("energy_expl.txt", energy);

	y1.assign(y1.size(), 0);
	y2.assign(y2.size(), 0);
	time.assign(time.size(), 0);
	energy.assign(energy.size(), 0);
	implicitEuler(y1, y2, time, y0, zeta, h, T);
	writeToFile("position_impl.txt", y1);
	writeToFile("velocity_impl.txt", y2);
	writeToFile("time_impl.txt", time);
	//Energy(y2, energy);
	writeToFile("energy_impl.txt", energy);

	return 0;
}
