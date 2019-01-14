#include <iostream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <chrono>

using namespace Eigen;

MatrixXd quadraticSpline(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	VectorXd dx = x.tail(N) - x.head(N);
	VectorXd dx2 = dx.cwiseProduct(dx);
	MatrixXd A = MatrixXd::Zero(3*N, 3*N);

	// require f(x_i) = y_i
	for(int i=0; i<N; i++){
		A(i, 2*N+i) = 1;
	}

	// require f(x_{i+1}) = y_{i+1}
	for(int i=0; i<N; i++){
		A(N+i, i) = dx2(i);
		A(N+i, N+i) = dx(i);
		A(N+i, 2*N+i) = 1;
	}

	// require f' to be continuous at x_i
	for(int i=0; i<N; i++){
		A(2*N+i, i) = 2*dx(i);
		A(2*N+i, N+i) = 1;
		A(2*N+i, N + ((i+1) % N)) = -1;
	}

	VectorXd rhs(3*N);
	rhs << y.head(N), y.tail(N), VectorXd::Zero(N);

	VectorXd out = A.fullPivLu().solve(rhs);
	return Map<MatrixXd>(out.data(), N, 3);
}

MatrixXd quadraticSplineFast(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	VectorXd dx = x.tail(N) - x.head(N);
	VectorXd dx2 = dx.cwiseProduct(dx);
	VectorXd dy = y.tail(N) - y.head(N);
	MatrixXd out(N,3);

	VectorXd b = dy.cwiseQuotient(dx);
	int sign = 1;

	for(int i=0; i<N-1; i++){
		sign *= -1;
		b(N-1) += sign*b(i);
	}

	for(int i=N-2; i>=0; i--){
		b(i) = -b(i+1) + 2*dy(i)/dx(i);
	}

	VectorXd a = - b.cwiseQuotient(dx) + dy.cwiseQuotient(dx2);

	out << a, b, y.head(N);
	return out;
}

MatrixXd quadraticSplineFast_SparseLU(const VectorXd &x, const VectorXd &y) {
	int N = x.size()-1;
	VectorXd dx = x.tail(N) - x.head(N);
	VectorXd dx2 = dx.cwiseProduct(dx);
	SparseMatrix<double> A(3*N, 3*N);
	A.reserve(VectorXd::Constant(3*N,3));

	// require f(x_i) = y_i
	for(int i=0; i<N; i++){
		A.insert(i, 2*N+i) = 1;
	}

	// require f(x_{i+1}) = y_{i+1}
	for(int i=0; i<N; i++){
		A.insert(N+i, i) = dx2(i);
		A.insert(N+i, N+i) = dx(i);
		A.insert(N+i, 2*N+i) = 1;
	}

	// require f' to be continuous at x_i
	for(int i=0; i<N; i++){
		A.insert(2*N+i, i) = 2*dx(i);
		A.insert(2*N+i, N+i) = 1;
		A.insert(2*N+i, N + ((i+1) % N)) = -1;
	}
	A.makeCompressed();

	VectorXd b(3*N);
	b << y.head(N), y.tail(N), VectorXd::Zero(N);

	SparseLU<SparseMatrix<double>> solver;
	solver.compute(A);
	VectorXd out = solver.solve(b);

	return Map<MatrixXd>(out.data(), N, 3);
}



/***** TESTING ******/

int main() {
	
	/*
	std::cout<<std::fixed;
	std::cout.precision(6);
	int max=14;
	std::cout<<max<<std::endl;
	int N=3;
	for(int j=0; j<max; j++){
		std::cout<<N<<std::endl;
		VectorXd x = VectorXd::LinSpaced(N+1,0,100);
		VectorXd y = VectorXd::Random(N+1);
		std::cout<<x<<std::endl;
		std::cout<<y<<std::endl;
		std::cout<<quadraticSplineFast(x,y)<<std::endl;
		int flag=0;
		std::cout<<flag<<std::endl;
		N<<=1;
		N++;
	}
	*/


	const double TOL = 10e-5;
	int testnum;
	std::cin >> testnum;
	for(int i=1; i<testnum; i++){
		int N;
		std::cin >> N;
		VectorXd x(N+1), y(N+1);
		MatrixXd ans(N,3);

		for(int j=0; j<N+1; j++)
			std::cin >> x(j);
		for(int j=0; j<N+1; j++)
			std::cin >> y(j);


		for(int j=0; j<N; j++){
			std::cin >> ans(j,0);
			std::cin >> ans(j,1);
			std::cin >> ans(j,2);
		}
		
		int flag;
		std::cin >> flag;

		auto start = std::chrono::steady_clock::now();	
		MatrixXd test;
		if(flag)
			test = quadraticSpline(x,y);
		else 
			test = quadraticSplineFast(x,y);
		auto end = std::chrono::steady_clock::now();
		double difftime = std::chrono::duration<double,std::milli>(end-start).count();

		if((test-ans).norm() > TOL * ans.norm()){
			std::cout << "Test " << i << " INCORRECT.\n";
			std::cerr << "Your answer:\n" << test << "\n\n" << "Correct answer:\n" << ans << "\n\n";
			break;
		}
		else
			std::cout << "Test " << i << " with N=" << N << " correct in " << difftime << "ms.\n";
	}
}
