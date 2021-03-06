// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>
#include <vector>


template<typename _Scalar,int DATA = Eigen:Dynamic, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
	typedef _Scalar Scalar;
	enum {
		
		InputsAtCompileTime = NX,
		ValuesAtCompileTime = NY
	};
	typedef Eigen::Matrix <Scalar, data, 1> InputType;
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
		/*const*/ double A[] = { 0,0,0 } /*= { 3.0, 10.0, 7.0 }*/, B[] = { 0,0,0,0,0,0,0,0,0 } /*= { -1.0,-2.0,-1.0,2.2,3.2,-5.2,-8.0,-1.2,1.0 }*/, C[] = { 0,0,0,0,0,0,0,0,0 }  /*= { -.05,-1.3,2.0,-2.0,-1.0,2.2,3.2,-5.2,0.2 }*/;
		int count = 0;
		
		for (count; count < n; count++)
			 A[count]= data [count];
		for (count; n < count < n * 4; count++)
			B[count - n] = data [count];
		for (count; n * 4 < count < n * 7; count++)
			C[count - 4 * n]= data [count] ;

		assert(fvec.size() == n);
		for (int k = 0,K=0; k < n; k++,K=K+3)
		{
			if (k < 4) {
				temp = A[k];
				if (k <( 8 - 2)) {
					temp1 = x[0] * B[K] + x[1]*B[K+1]+x[2]*B[K+2];
					temp2 = x[0] * x[0] * C[K] + x[1] * x[1] * C[K+1] + x[2] * x[2] * C[K+2];
				}
				fvec[k] = temp + temp1 + temp2;
			}
			//if (k >= 4)
				/*if ((k - 4 < 1) || (k - 4 > 2))
					fvec[k] = x[k - 6] - x[k - 6 + 1];
				else if (k - 6 == 1)
					fvec[k] = x[k - 6 - 1] - (x[k - 6] + x[k - 6 + 1]);
				else */
				///	fvec[k] = x[k - 4 +2] - (x[k - 4+1] + x[k - 4 ]);
					
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
	std::vector <double> c;
	Eigen::VectorXd x(n);
	Eigen::VectorXd data(7*n);
	double A[] = { 3.0,10.0,7.0 }, B[] = { -1.0,-2.0,-1.0,2.2,3.2,-5.2,-8.0,-1.2,1.0 }, C[] = { -.05,-1.3,2.0,-2.0,-1.0,2.2,3.2,-5.2,0.2 };
	int count = 0;


	/* the following starting values provide a rough fit. */
	x.setConstant(n,0.0);


		for (count; count < n; count++)
			data[count] = A[count];
		for (count; n < count < n * 4; count++)
			data[count] = B[count - n];
		for (count; n * 4 < count < n * 7; count++)
			data[count] = C[count-4*n];





	// do the computation
	hybrj_functor functor ;
	Eigen::HybridNonLinearSolver<hybrj_functor> solver(functor);
	solver.diag.setConstant(n, 1.);
	solver.useExternalScaling = true;
	info = solver.hybrj1(x);

	c.resize(x.size());

	for (int count = 0; count < x.size(); count++)
		c.at(count) = x(count);

	
	
}

int main()
{
	tests();

    return 0;
}

