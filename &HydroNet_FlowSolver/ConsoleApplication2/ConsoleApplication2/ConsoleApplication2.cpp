// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <Eigen/NonLinearOptimization>
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include<Eigen/SparseQR>



template < typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic >
	struct Functor {
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		OutputsAtCompileTime = NY
		 
	};

	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<Scalar, OutputsAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, OutputsAtCompileTime> JacobianType;

	int m_inputs, m_values;

	Functor() : m_inputs(In), m_values(ValuesAtCompileTime) {}

	Functor(int inputs, int values) : m_inputs(inputs),	m_values(values)
	{}
	int inputs() const {
		return m_inputs;
	}

	int values() const {
		return m_values;
	}
	
	
};

struct MyFunctor: Functor <double> {

	MyFunctor(void):	Functor <double> (4, 4) {}

	int operator () (const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {

		//std::cout << "x:\n" << x << "\n" <<std::endl;
	//	std::cout << "fvec:\n" << fvec << " | norm: " << fvec.norm() << "\n" << std::endl;

		// Random equations
		fvec[0] = x[0] - x[1];
		fvec[1] = 2 * x[1];
		fvec[2] = 0.5 * x[2] - x[0] + 5;
		fvec[3] = x[4] - 2 * x[1];
		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac)
	{
		const int n = x.size();
		assert(fjac.rows() == n);
		assert(fjac.cols() == n);
		for (int k = 0; k < n; k++)
		{
			for (int j = 0; j < n; j++)
				fjac(k, j) = 0.;
		}


		return 0;
	}
};




int main()
{
	int info;
	Eigen::VectorXd x = { -1.0,6.5,-3.2,4.1 };
	
	
	//wpr << -2.0, 5.6;
	//wpr << -3.5, 2.0, -25.0;
	//std::cout << "Initial guess of the solution:\n" << wpr << std::endl;  // The answer to the problem is: -3, 5, -16

	MyFunctor functor;  
	Eigen::HybridNonLinearSolver <MyFunctor>  solver(functor);  
	info = solver.solve(x);
	
	
	//int ret = lm.minimize(wpr);  // Minimize using LM algorithm; initial guess is "x"
	//std::cout << "\nIterations: " << lm.iterations << std::endl;

	//std::cout << "\nx that minimizes the function:\n" << wpr << std::endl;  // x has been minimized

    return 0;
}