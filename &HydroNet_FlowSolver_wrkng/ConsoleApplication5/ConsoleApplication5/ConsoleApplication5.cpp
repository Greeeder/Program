// ConsoleApplication5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>

#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>

struct MyFunctor
{
	int operator()(const Eigen::VectorXf &x, Eigen::VectorXf &fvec) const
	{
		// Implement y = x^2
		fvec(0) = x(0)*x(0) - 5.0 * x(0);
		return 0;
	}

	int df(const Eigen::VectorXf &x, Eigen::MatrixXf &fjac) const
	{
		// Implement dy/dx = 2*x
		fjac(0) = 2.0f * x(0)+5.0;
		return 0;
	}

	int inputs() const { return 1; }
	int values() const { return 1; } // number of constraints
};


int main(int argc, char *argv[])
{
	Eigen::VectorXf x(1);
	x(0) = 4.0;
	std::cout << "x: " << x << std::endl;

	MyFunctor functor;
	Eigen::LevenbergMarquardt<MyFunctor, float> lm(functor);
	lm.minimize(x);

	std::cout << "x that minimizes the function: " << x << std::endl;

	return 0;
}