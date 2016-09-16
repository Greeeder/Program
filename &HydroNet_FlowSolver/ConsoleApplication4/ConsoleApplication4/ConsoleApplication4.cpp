// ConsoleApplication4.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include<unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>



template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
	typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

	const int m_inputs, m_values;

	Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
	Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

	int inputs() const { return m_inputs; }
	int values() const { return m_values; }

	// you should define that in the subclass :
	//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

struct lmder_functor : Functor<double>
{
	lmder_functor(void) : Functor<double>(3, 3) {}
	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
	{
		
		for (int i=0; i < 15; i++) {
			fvec(0) = 30 - x(0) + x(1)*x(1) - sqrt(x(2));
			fvec(1) = 10 - x(0) / x(2);
			fvec(2) = 0;
		}
		return 0;
	}

	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
	{
		
		
			
			
		
		return 0;
	}
};

void test()
{
	int n = 3, info;
	double c;

	Eigen::VectorXd x;

	/* the following starting values provide a rough fit. */
	x.setConstant(n, 8.);

	// do the computation
	lmder_functor functor;
	Eigen::LevenbergMarquardt<lmder_functor> lm(functor);
	info = lm.minimize(x);

	c = x(0);
	c = x(1);
	c = x(2);



}


int main()
{
	test(); 

    return 0;
}

