// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>


template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Eigen::Matrix <Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix <Scalar, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix <Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

	const int m_inputs, m_values;

	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	int inputs() const { return m_inputs; }
	int values() const { return m_values; }

	// you should define that in the subclass :
	//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

struct hybrj_functor : Functor <double>
{
	hybrj_functor(void) : Functor <double>(3,3) {}

	


	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
	{
		double temp, temp1, temp2;
		const int n = x.size();
		const double A[] = { 21.0,-10.0, 9.0 }, B[] = { -5.0,-1.0,0.0 }, C[] = { -3.5,1.2,- 3.0 };

		assert(fvec.size() == n);
		for (int k = 0; k < n; k++)
		{
			
			temp = A[k];
			temp1 = x[k] * B[k];
			temp2 = x[k] * x[k] * C[k];
			fvec[k] = temp + temp1 + temp2;

		}
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
			fjac(k, k) = 3. - 4.*x[k];
			if (k) fjac(k, k - 1) = -1.;
			if (k != n - 1) fjac(k, k + 1) = -2.;
		}
		return (0);
	}
};


void tests()
{
	const int n =3;
	int info;
	double c;
	Eigen::VectorXd x(n);

	/* the following starting values provide a rough fit. */
	x.setConstant(n,-1.);


	// do the computation
	hybrj_functor functor;
	Eigen::HybridNonLinearSolver<hybrj_functor> solver(functor);
	solver.diag.setConstant(n, 1.);
	solver.useExternalScaling = true;
	info = solver.hybrd1(x);

	c = x(0);
	c = x(1);
	c = x[2];
	
	
}

int main()
{
	tests();

    return 0;
}

