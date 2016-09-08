// ConsoleApplication6.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>

class Solver {
	// Functor
private:
	std::vector<double> A, B, C, D;
	int n;

	

public:
	Solver(Eigen::VectorXd x0, Eigen::VectorXd data) {

		int count=0;

		n = x0.size();
		A.resize(n);
		B.resize(n*n);
		C.resize(n*n);
		D.resize(n*n);

		
		for (count; count < n; count++)
			A[count] = data[count];
		for (count; count < n*n+n; count++)
			B[count-n] = data[count];
		for (count; count < 2*n*n+n; count++)
			C[count-n*n-n] = data[count];
		for (count; count < 3*n*n+n; count++)
			D[count-2*n*n-n] = data[count];

	};
	

	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;
	int inputs() const;// inputs is the dimension of x.
	int values() const; // "values" is the dimension of F 

};

int Solver::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
	
	double temp, temp1, temp2;
	assert(fvec.size() == n);
	//Implementation of Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional sistems

	
	for (int k = 0; k < n; k++){

		temp = 0;
		temp = A[k];
		for (int j = 0; j<n; j++)
		
			temp = temp + B[j + k * n] * x[j] / abs(x[j]);
		
		temp1 = 0;
		for (int j=0; j<n;j++)
			temp1 = temp1+ C[j + k * n] * x[j];


		temp2 = 0;
		for (int j = 0; j<n; j++)
			temp2 = temp2 + D[j + k * n] * abs(x[j])*x[j];

		fvec[k] = temp + temp1 + temp2;
		}
		return 0;
	}



int Solver::df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const {

	assert(fjac.rows() == n);
	assert(fjac.cols() == n);
	//Implementation of Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional sistems Jacobian
	for (int k = 0; k < n; k++)
		for (int j = 0; j < n; j++) {
	
		fjac(k, j) = 2 * D[n*k+j]*x[j] + C[n*k+j];
	}
	return (0);

}
int Solver::inputs() const { return n; } // inputs is the dimension of x.
int Solver::values() const { return n; } // "values" is the number of f_i 



std::vector <double> tests(std::vector<double> A, std::vector<double>  B, std::vector<double> C, std::vector<double> D)
{
	/*
	A is an independent constant vector
	B is a vector dependent on the solution sing
	C ia a vector containing the lineal coefficients
	D ia a vector containing the quadratic coefficients	
	Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional
	*/

	int n;
	n = A.size();
	std::vector <double> _results;
	int info;
	Eigen::VectorXd x(n);
	Eigen::VectorXd data(10 * n);
	int count = 0;


	//Starting values 	
	x.setConstant(n, 5);

	x[0] = 0.001;
	x[1] = -0.001;
	


	//The matrix an vector's coeficients are stored in a single Eigen vector
	for (count; count < n; count++)
		data[count] = A[count];
	for (count; count < n*4; count++)
		data[count] = B[count-n];
	for (count; count < n * 7; count++)
		data[count] = C[count - n*4];
	for (count;  count < n * 10; count++)
		data[count] = D[count - n*7];

	for (int j = 0; j < x.size(); j++)
		if (abs(x[j]) < 1E-9)
			x[j] = 1E-9;



	//Computation
	Solver functor ( x, data);								
	/*Eigen::LevenbergMarquardt<Solver,double> lm(functor);
	lm.parameters.xtol = 1e-16;
	lm.parameters.ftol = 1e-16;
	info = lm.minimize(x);*/
	Eigen::HybridNonLinearSolver<Solver> solver(functor);
	solver.diag.setConstant(n, 1.);
	solver.useExternalScaling = true;
	info = solver.hybrj1(x);


	//The results are stored in the vector _results
	_results.resize(x.size());
	for (int count = 0; count < x.size(); count++)
		_results.at(count) = x(count);


	return _results;

}

int main()
{
	std::vector <double> results;
	std::vector<double> A = { 0.001,2.120926,2.120926 }, B = { 0,0,0,0,0,0,0,0,0 }, C = { 0,1,1,-10000,0,0,-10000,0,0 }, D = { 0,0,0,0,-50000,0,0,0,-100000 };
	results = tests(A,B,C,D);

	return 0;
}

