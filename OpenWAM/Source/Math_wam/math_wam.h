/* --------------------------------------------------------------------------------*\
==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 \\ |  X  | //  W ave     |
 \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*-------------------------------------------------------------------------------- */

/*!
 * \file Math_wam.h
 * \author Francisco Jose Arnau <farnau@mot.upv.es>
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 *
 * \section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenWAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \section DESCRIPTION
 * This file declares several auxiliary math functions, as well as some
 * typedefs.
 */

// ---------------------------------------------------------------------------
#ifndef Math_wamH
#define Math_wamH
// ---------------------------------------------------------------------------

#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <iostream>
#include <Eigen/Dense>
#include "Wrappers/Exception.hpp"

using namespace std;
using namespace Eigen;

typedef unsigned int Uint; ///< Unsigned integer
typedef std::vector<double> dVector; ///< Double vector
typedef std::vector<std::vector<double> > dMatrix; ///< 2-dimensional double matrix
typedef std::vector<int> iVector; ///< Integer vector
typedef std::vector<std::vector<int> > iMatrix; ///< 2-dimensional integer matrix
typedef std::vector<bool> bVector; ///< Boolean vector
typedef std::vector<std::vector<bool> > bMatrix; ///< 2-dimensional boolean matrix
typedef Array<double, Dynamic, 1, ColMajor, Dynamic, 1> ColVector; ///< Column vector.
typedef Array<double, 1, Dynamic, RowMajor, 1, Dynamic> RowVector; ///< Row vector.
typedef Array<double, Dynamic, Dynamic, ColMajor, Dynamic, Dynamic> ColArray; ///< Column major array.
typedef Array<double, Dynamic, Dynamic, RowMajor, Dynamic, Dynamic> RowArray; ///< Row major array.
typedef RowArray dArray;

/*!
 * \brief Computes a square root of a double.
 *
 * Computes @f$ \sqrt{x} @f$, x being a double.
 *
 * \param x The value.
 * \return @f$ x ^ 2 @f$
 */
inline double Sqrt(double x) {
#ifdef usetry
	if(x < 0.) {
		throw(Exception("Sqrt of negative number."));
	}
#endif
	return std::sqrt(x);
}

/*!
 * \fn	inline ArrayXd Sqrt(ArrayXd x)
 *
 * \brief	Computes a square root of a double.
 *
 * Computes @f$ \sqrt{x} @f$, x being an array.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	18/02/2016
 *
 * \param	x	The array.
 *
 * \return	@f$ x ^ 2 @f$.
 */

inline ArrayXd Sqrt(ArrayXd x) {
	return x.sqrt();
}

/*!
 * \brief Computes a square root of a float.
 *
 * Computes @f$ \sqrt{x} @f$, x being a float.
 *
 * \param x The value.
 * \return @f$ x ^ 2 @f$
 */
inline float Sqrt(float x) {
#ifdef usetry
	if(x < 0.) {
		throw(Exception("Sqrt of negative number."));
	}
#endif
	return std::sqrt(x);
}

/*!
 * \brief Computes a square root of a long double.
 *
 * Computes @f$ \sqrt{x} @f$, x being a long double.
 *
 * \param x The value.
 * \return @f$ x ^ 2 @f$
 */
inline long double Sqrt(long double x) {
#ifdef usetry
	if(x < 0.) {
		throw(Exception("Sqrt of negative number."));
	}
#endif
	return std::sqrt(x);
}

/*!
 * \brief Computes a square root of an integer type.
 *
 * Computes @f$ \sqrt{x} @f$, x being of integer type.
 *
 * \param x The value.
 * \return @f$ x ^ 2 @f$
 */
template<class T>
inline double Sqrt(T x) {
#ifdef usetry
	if(x < 0.) {
		throw(Exception("Sqrt of negative number."));
	}
#endif
	return std::sqrt(x);
}

template<class T>

/*!
 * \fn	void DimdMatrix(T& x, int n, int m)
 *
 * \brief	Add dimensions to a matrix std::vector<std::vector<T>>.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param [in,out]	x	The matrix.
 * \param	n		 	The number of rows.
 * \param	m		 	The number of columns.
 */

void DimdMatrix(T& x, int n, int m) {
	x.resize(n);
	for(int i = 0; i < n; i++)
		x[i].resize(m);
}

double Interpola(double vizq, double vder, double axid, double xif);

/*!
 * \brief Returns x to the power of 2.
 *
 * Computes @f$ x ^ 2 @f$
 *
 * \param x The value.
 * \return @f$ x ^ 2 @f$
 */
template<class T>
inline T pow2(T x) {
	return x * x;
}

/*!
 * \brief Returns x to the power of 3.
 *
 * Computes @f$ x ^ 3 @f$
 *
 * \param x The value.
 * \return @f$ x ^ 3 @f$
 */
template<class T>
inline T pow3(T x) {
	return x * x * x;
}

/*!
 * \brief Returns x to the power of 4.
 *
 * Computes @f$ x ^ 4 @f$
 *
 * \param x The value.
 * \return @f$ x ^ 4 @f$
 */
template<class T>
inline T pow4(T x) {
	return x * x * x * x;
}

/*!
 * \brief Returns x to the power of 0.25.
 *
 * Computes @f$ x ^ {0.25} @f$
 *
 * \param x The value.
 * \return @f$ x ^ {0.25} @f$
 */
template<class T>
inline T pow025(T x) {
	return Sqrt(Sqrt(x));
}

/*!
 * \brief Returns x to the power of 1.5.
 *
 * Computes @f$ x ^ {1.5} @f$
 *
 * \param x The value.
 * \return @f$ x ^ {1.5} @f$
 */
template<class T>
inline T pow150(T x) {
	return Sqrt(pow3(x));
}

/*!
 * \brief Returns x to the power of 0.75.
 *
 * Computes @f$ x ^ {0.75} @f$
 *
 * \param x The value.
 * \return @f$ x ^ {0.75} @f$
 */
template<class T>
inline T pow075(T x) {
	return Sqrt(pow150(x));
}

template<class T, class U>

/*!
 * \fn	inline T poww(T x, U y)
 *
 * \brief	Computes @f$ x ^ {y} @f$.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	x	The value.
 * \param	y	The exponent.
 *
 * \return	@f$ x ^ {y} @f$ if x is higher of the precision limit of the computer else return 0.
 */

inline T poww(T x, U y) {
	return abs(x) > std::numeric_limits<T>::epsilon() ? pow(x, y) : 0;
}

template<class T>

/*!
 * \fn	inline T Sqrtw(T x)
 *
 * \brief	Computes @f$ x ^ {0.5} @f$.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	x	The value.
 *
 * \return	@f$ x ^ {0.5} @f$ if x is higher of the precision limit of the computer else return 0.
 */

inline T Sqrtw(T x) {
	return abs(x) > std::numeric_limits<T>::epsilon() ? Sqrt(x) : 0;
}

#ifdef __BORLANDC__
template<class T>
inline T cbrt(T x) {
	return pow(x, 1 / 3) > std::numeric_limits<T>::epsilon() ? pow(x, 1 / 3) : 0;
}
#endif

template<class T>

/*!
 * \fn	inline bool DoubEqZero(T x)
 *
 * \brief	Compare if a double is equal to zero.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	x	The double.
 *
 * \return	true if it succeeds, false if it fails.
 */

inline bool DoubEqZero(T x) {
	return abs(x) > std::numeric_limits<T>::epsilon() ? false : true;
}

/*!
 * \fn	double QuadraticEqP(double A, double B, double C);
 *
 * \brief	Computes the positive root of a quadratic equation
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	A	The coeficient of \f$x^2\f$.
 * \param	B	The coeficient of \f$x\f$.
 * \param	C	The independent coeficient.
 *
 * \return	The positive root.
 */

double QuadraticEqP(double A, double B, double C);

/*!
 * \fn	inline double QuadraticEqN(double A, double B, double C);
 *
 * \brief	Computes the negative root of a quadratic equation
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	A	The coeficient of \f$x^2\f$.
 * \param	B	The coeficient of \f$x\f$.
 * \param	C	The independent coeficient.
 *
 * \return	The negative root.
 */

inline double QuadraticEqN(double A, double B, double C);

/*!
 * \fn	void Poly2(dVector &X, dVector &Y, dVector &A);
 *
 * \brief	Compute the coeficients of a second order polinomy that cross \f$(x, y)\f$.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param [in,out]	X	The x values.
 * \param [in,out]	Y	The y values.
 * \param [in,out]	A	The coefficients of the polinomy.
 */

void Poly2(dVector &X, dVector &Y, dVector &A);

template<class T>

/*!
 * \fn	inline const T &Max(const T &a, const T &b)
 *
 * \brief	Determines the maximum between a and b.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	a	The a parameter.
 * \param	b	The b parameter.
 *
 * \return	The maximum value.
 */

inline const T &Max(const T &a, const T &b) {
	return b > a ? (b) : (a);
}

template<class T>

/*!
 * \fn	inline const T &Min(const T &a, const T &b)
 *
 * \brief	Determines the minimum between a and b.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	a	The a parameter.
 * \param	b	The b parameter.
 *
 * \return	The minimum value.
 */

inline const T &Min(const T &a, const T &b) {
	return b < a ? (b) : (a);
}

/*!
 * \fn	inline double Mean(dVector a)
 *
 * \brief	Determines the mean value of the elements of a dVector.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	a	The dVector.
 *
 * \return	The mean value.
 */

inline double Mean(dVector a) {
	double val = 0.;
	for(int i = 0; i < a.size(); i++) {
		val += a[i];
	}
	return val / (double) a.size();
}

template<class T>

/*!
 * \fn	inline T Sign(const T &a, const T &b)
 *
 * \brief	Return a multiplied by the sign of b.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/03/2016
 *
 * \param	a	The a.
 * \param	b	The b.
 *
 * \return	a multiplied by the sign of b.
 */

inline T Sign(const T &a, const T &b) {
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline float Sign(const float &a, const double &b);

inline float Sign(const double &a, const float &b);

template<class T>
inline void Swap(T &a, T &b) {
	T dum = a;
	a = b;
	b = dum;
}

template<class T>
inline T MaxComponent(std::vector<T> &x) {
	T max = x[0];
	for(Uint i = 1; i < x.size(); i++) {
		if(x[i] > max)
			max = x[i];
	}
	return max;
}

template<class T>
inline T MinComponent(std::vector<T> &x) {
	T min = x[0];
	for(Uint i = 1; i < x.size(); i++) {
		if(x[i] < min)
			min = x[i];
	}
	return min;
}


struct stPolar {
	double Ang;
	double Mod;

	stPolar();

	stPolar(double X, double Y);

	void operator()(double X, double Y);
};

struct stRectan {
	double X;
	double Y;

	stRectan();

	stRectan(double Ang, double Mod);

	void operator()(double Ang, double Mod);
};

struct Base_interp {
	int n, mm, jsav, cor, dj;
	const double *xx, *yy;

	Base_interp();

	Base_interp(dVector &x, const double *y, int m);

	double interp(double x);

	int locate(const double x);
	int hunt(const double x);

	double virtual rawinterp(int jlo, double x) = 0;
};

struct Linear_interp: Base_interp {
	Linear_interp();

	Linear_interp(dVector &xv, dVector &yv);

	void operator()(dVector & xv, dVector & yv);

	double rawinterp(int j, double x);
};

struct Hermite_interp: Base_interp {
	dVector y2;

	Hermite_interp();

	Hermite_interp(dVector &xv, dVector &yv);

	void operator()(dVector & xv, dVector & yv);

	void sety2(const double *xv, const double *yv);
	double rawinterp(int j, double x);

};

struct Step_interp: Base_interp {
	Step_interp();

	Step_interp(dVector &xv, dVector &yv);

	void operator()(dVector & xv, dVector & yv);

	double rawinterp(int j, double x);
};

/*!
 * \brief Solves a function using Brent's method.
 *
 * Finds the root of a function using Brent's method.  It uses a combination of
 * the bisection method, the secant method and inverse quadratic interpolation,
 * and finds the root in an interval.  The function must have different signs
 * at each bound of the interval.  If the signs are the same, returns one of
 * the bounds, whichever is minimum.
 *
 * \param func The function.
 * \param x1 The lower bound of the interval.
 * \param x2 The upper bound of the interval.
 * \param tol The tolerance of the method.
 *
 * \return The interpolated value.
 */
template<class T>
inline double zbrent(T& func, const double& x1, const double& x2, const double& tol) {
	const int ITMAX = 100;
	const double EPS = std::numeric_limits<double>::epsilon();
	double a = x1;
	double b = x2;
	double c = x2;
	double d;
	double e;
	double fa = func(a);
	double fb = func(b);
	double fc;
	double p;
	double q;
	double r;
	double s;
	double tol1;
	double xm;

	if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
		if(fabs(fa) < fabs(fb)) {
			fa = func(a);
			return a;
		} /* Condicion original if((fabs(fa) < fabs(fb)) && fabs(fa) < tol) */
		else if(fabs(fa) > fabs(fb)) {
			fb = func(b);
			return b;
		} /* Original if((fabs(fa) > fabs(fb)) && fabs(fb) < tol) */
		// return 0;
		/* throw("Root must be bracketed in zbrent"); */
	}

	fc = fb;
	for(int iter = 0; iter < ITMAX; iter++) {
		if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c = a;
			fc = fa;
			e = d = b - a;
		}
		if(fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
		xm = 0.5 * (c - b);
		if(fabs(xm) <= tol1 || fb == 0.0) {
			/* std::cout << iter << std::endl; */
			return b;
		}
		if(fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s = fb / fa;
			if(a == c) {
				p = 2.0 * xm * s;
				q = 1.0 - s;
			} else {
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			if(p > 0.0)
				q = -q;
			p = fabs(p);
			double min1 = 3.0 * xm * q - fabs(tol1 * q);
			double min2 = fabs(e * q);
			if(2.0 * p < (min1 < min2 ? min1 : min2)) {
				e = d;
				d = p / q;
			} else {
				d = xm;
				e = d;
			}

		} else {
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if(fabs(d) > tol1)
			b += d;
		else
			b += Sign(tol1, xm);
		fb = func(b);
	}
	return b;
}

/*!
 * \brief Finds the root of a function.
 *
 * Finds the root of a function using Brent's method.  It uses a combination of
 * the bisection method, the secant method and inverse quadratic interpolation,
 * and finds the root in an interval.  The function must have different signs
 * at each bound of the interval.  If the signs are the same, returns one of
 * the bounds, whichever is minimum.  It uses a default tolerance of 1e-10.
 *
 * \param func The function.
 * \param x1 The lower bound of the interval.
 * \param x2 The upper bound of the interval.
 *
 * \return The interpolated value.
 */
template<class T>
inline double FindRoot(T & func, const double x1, const double x2) {
	return zbrent(func, x1, x2, 1e-10);
}

template<class T>

/*!
 * \brief	Given a function or functor func and an initial guessed range x1 to x2, the routine expands
 *			the range geometrically until a root is bracketed by the returned values x1 and x2 (in which
 *			case zbrac returns true) or until the range becomes unacceptably large (in which case zbrac
 *			returns false).
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param [in,out]	x1  	The initial range limit.
 * \param [in,out]	x2  	The initial range limit.
 *
 * \return	true if it succeeds, false if it fails.
 */

inline bool zbrac(T & func, double & x1, double & x2) {
	const int NTRY = 200;
	const double FACTOR = 0.1;
	if(x1 == x2)
		throw("Bad initial range in zbrac");
	double f1 = func(x1);
	double f2 = func(x2);
	for(int j = 0; j < NTRY; j++) {
		if(f1 * f2 < 0.0)
			return true;
		if(fabs(f1) < fabs(f2))
			f1 = func(x1 += FACTOR * (x1 - x2));
		else
			f2 = func(x2 += FACTOR * (x2 - x1));
	}
}

template<class T>

/*!
 * \brief	Given a function or functor func and an initial guessed range x1 to x2, the routine expands
 *			the range geometrically until a root is bracketed by the returned values x1 and x2 (in which
 *			case zbrac returns true) or until the range becomes unacceptably large or the values reach the
 *			maximum or minimum values (in which case zbrac returns false).
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param [in,out]	x1  	The initial range limit.
 * \param [in,out]	x2  	The initial range limit.
 * \param	min				The minimum value of the range.
 * \param	max				The maximum value of the range.
 *
 * \return	true if it succeeds, false if it fails.
 */

inline bool zbrac2(T & func, double & x1, double & x2, const double & min, const double & max) {
	const int NTRY = 200;
	// const double FACTOR=0.1;
	if(x1 == x2)
		throw("Bad initial range in zbrac");
	double x1lim = min;
	double x2lim = max;
	if(x1 > x2) {
		x1lim = max;
		x2lim = min;
	}
	double f1 = func(x1);
	double f2 = func(x2);
	for(int j = 0; j < NTRY; j++) {
		if(f1 * f2 < 0.0)
			return true;
		if(fabs(f1) < fabs(f2))
			f1 = func(x1 = (0.9 * x1 + 0.1 * x1lim));
		else
			f2 = func(x2 = (0.9 * x2 + 0.1 * x2lim));
	}
	std::cout << "Do not exist solution" << std::endl;
	std::cout << "x1 :" << x1 << " f1: " << f1 << " x1max: " << x1lim << std::endl;
	std::cout << "x2 :" << x2 << " f2: " << f2 << " x2max: " << x2lim << std::endl;

	return false;

}

template<class T>

/*!
 * \brief	Given a function or functor fx defined on the interval [x1,x2], subdivide the interval into
 *			n equally spaced segments, and search for zero crossings of the function. nroot will be set
 *			to the number of bracketing pairs found. If it is positive, the arrays xb1[0..nroot-1] and
 *			xb2[0..nroot-1] will be filled sequentially with any bracketing pairs that are found. On input,
 *			these vectors may have any size, including zero; they will be resized to 
 *			nroot.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	fx   	The effects.
 * \param	x1			 	The low interval limit.
 * \param	x2			 	The high interval limit.
 * \param	n			 	The number of segments.
 * \param [in,out]	xb1  	The low limit of every root.
 * \param [in,out]	xb2  	The high limit of every root.
 * \param [in,out]	nroot	The number of roots.
 */

inline void zbrak(T & fx, const double x1, const double x2, const int n, dVector & xb1, dVector & xb2, int & nroot) {
	int nb = 20;
	xb1.resize(nb);
	xb2.resize(nb);
	nroot = 0;
	double dx = (x2 - x1) / n;
	double x = x1;
	double fp = fx(x += dx);
	for(int i = 0; i < n; i++) {
		double fc = fx(x += dx);
		if(fc * fp <= 0.0) {
			xb1[nroot] = x - dx;
			xb2[nroot++] = x;
			if(nroot == nb) {
				dVector tempvec1(xb1), tempvec2(xb2);
				xb1.resize(2 * nb);
				xb2.resize(2 * nb);
				for(int j = 0; j < nb; j++) {
					xb1[j] = tempvec1[j];
					xb2[j] = tempvec2[j];
				}
				nb *= 2;
			}
		}
		fp = fc;
	}
}

template<class T>

/*!
 * \brief	Using bisection, return the root of a function or functor func known to lie between x1 and x2.
 *			The root will be refined until its accuracy is xacc.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param	x1				The low limit of the range.
 * \param	x2				The high limit of the range.
 * \param	xacc			The accuracy.
 *
 * \return	.
 */

inline double rtbis(T & func, const double x1, const double x2, const double xacc) {
	const int JMAX = 50;

	double dx, xmid, rtb;

	double f = func(x1);
	double fmid = func(x2);
	if(f * fmid >= 0.0) {
		if((fabs(f) < fabs(fmid)) && fabs(f) < xacc) {
			return x1;
		} else if((fabs(f) > fabs(fmid)) && fabs(fmid) < xacc)
			return x2;
		throw("Root must be bracketed for bisection in rtbis");
	}
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for(int j = 0; j < JMAX; j++) {
		fmid = func(xmid = rtb + (dx *= 0.5));
		if(fmid <= 0.0)
			rtb = xmid;
		if(fabs(dx) < xacc || fmid == 0.0) {
			/* std::cout << j << std::endl; */return rtb;
		}
	}
	throw("Too many bisections in rtbis");
}

template<class T>

/*!
 * \brief	Using the false-position method, return the root of a function or functor func known to lie
 *			between x1 and x2. The root is refined until its accuracy is xacc.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param	x1				The low limit of the range.
 * \param	x2				The high limit of the range.
 * \param	xacc			The accuracy.
 *
 * \return	.
 */

inline double rtflsp(T & func, const double x1, const double x2, const double xacc) {
	const int MAXIT = 1000;

	double xl, xh, del;

	double fl = func(x1);
	double fh = func(x2);
	if(fl * fh > 0.0)
		throw("Root must be bracketed in rtflsp");
	if(fl < 0.0) {
		xl = x1;
		xh = x2;
	} else {
		xl = x2;
		xh = x1;
		Swap(fl, fh);
	}
	double dx = xh - xl;
	for(int j = 0; j < MAXIT; j++) {
		double rtf = xl + dx * fl / (fl - fh);
		double f = func(rtf);
		if(f < 0.0) {
			del = xl - rtf;
			xl = rtf;
			fl = f;
		} else {
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if(fabs(del) < xacc || f == 0.0) {
			/* std::cout << j << std::endl; */return rtf;
		}
	}
	throw("Maximum number of iterations exceeded in rtflsp");
}

template<class T>

/*!
 * \brief	Using the secant method, return the root of a function or functor func thought to lie between
 *			x1 and x2. The root is refined until its accuracy is xacc.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param	x1				The low limit of the range.
 * \param	x2				The high limit of the range.
 * \param	xacc			The accuracy.
 *
 * \return	.
 */

inline double rtsec(T & func, const double x1, const double x2, const double xacc) {
	const int MAXIT = 100;

	double xl, rts;

	double fl = func(x1);
	double f = func(x2);
	if(fabs(fl) < fabs(f)) {
		rts = x1;
		xl = x2;
		Swap(fl, f);
	} else {
		xl = x1;
		rts = x2;
	}
	for(int j = 0; j < MAXIT; j++) {
		double dx = (xl - rts) * f / (f - fl);
		xl = rts;
		fl = f;
		rts += dx;
		f = func(rts);
		if(fabs(dx) < xacc || f == 0.0) {
			/* std::cout << j << std::endl; */return rts;
		}
	}
	throw("Maximum number of iterations exceede in rtsec");
}

template<class T>

/*!
 * \brief	Using Ridders’ method, return the root of a function or functor func known to lie between x1
 *			and x2. The root will be refined to an approximate accuracy xacc.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	func	The function.
 * \param	x1				The low limit of the range.
 * \param	x2				The high limit or the range.
 * \param	xacc			The accuracy.
 *
 * \return	.
 */

inline double zriddr(T & func, const double x1, const double x2, const double xacc) {
	const int MAXIT = 100;
	double fl = func(x1);
	double fh = func(x2);
	if((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		double xl = x1;
		double xh = x2;
		double ans = -9.99e99;
		for(int j = 0; j < MAXIT; j++) {
			double xm = 0.5 * (xl + xh);
			double fm = func(xm);
			double s = Sqrt(fm * fm - fl * fh);
			if(s == 0.0) {
				/* std::cout << j << std::endl; */return ans;
			}
			double xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
			if(fabs((double)(xnew - ans <= xacc))) {
				/* std::cout << j << std::endl; */return ans;
			}
			ans = xnew;
			double fnew = func(ans);
			if(fnew == 0) {
				/* std::cout << j << std::endl; */return ans;
			}
			if(Sign(fm, fnew) != fm) {
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
			} else if(Sign(fl, fnew) != fl) {
				xh = ans;
				fh = fnew;
			} else if(Sign(fh, fnew) != fh) {
				xl = ans;
				fl = fnew;
			} else
				throw("never get here.");
			if(fabs(xh - xl) <= xacc) {
				/* std::cout << j << std::endl; */return ans;
			}
		}
		throw("zriddr exceed maximum iterations");
	} else {
		if(fl == 0.0)
			return x1;
		if(fh == 0.0)
			return x2;
		throw("root must be bracketed in zriddr.");
	}
}

template<class T>

/*!
 * \brief	Using the Newton-Raphson method, return the root of a function known to lie in the interval
 *			[x1; x2]. The root will be refined until its accuracy is known within xacc. funcd is a usersupplied
 *			struct that returns the function value as a functor and the first derivative of the function
 *			at the point x as the function df.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	funcd	The function.
 * \param	x1			 	The low limit of the range.
 * \param	x2			 	The high limit of the range.
 * \param	xacc		 	The accuracy.
 *
 * \return	.
 */

inline double rtnewt(T & funcd, const double x1, const double x2, const double xacc) {
	const int JMAX = 20.;
	double rtn = 0.5 * (x1 + x2);
	for(int j = 0; j < JMAX; j++) {
		double f = funcd(rtn);
		double df = funcd.df(rtn);
		double dx = f / df;
		rtn -= dx;
		if((x1 - rtn) * (rtn - x2) < 0.0)
			throw("Jumped out of brackets in rtnewt");
		if(fabs(dx) < xacc)
			return rtn;
	}
	throw("Maximum number of iterations exceede in rtnewt");
}

template<class T>

/*!
 * \brief	Using a combination of Newton-Raphson and bisection, return the root of a function bracketed
 *			between x1 and x2. The root will be refined until its accuracy is known within xacc. funcd
 *			is a user-supplied struct that returns the function value as a functor and the first derivative of
 *			the function at the point x as the function df (see text).
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	16/03/2016
 *
 * \param [in,out]	funcd	The function.
 * \param	x1			 	The low limit of the range.
 * \param	x2			 	The high limit of the range.
 * \param	xacc		 	The accuracy.
 *
 * \return	.
 */

inline double rtsafe(T & funcd, const double x1, const double x2, const double xacc) {
	const int MAXIT = 100;

	double xh, xl;

	double fl = funcd(x1);
	double fh = funcd(x2);
	if((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		throw("Root mus be bracketed in rt safe");
	}
	if(fl == 0.0)
		return x1;
	if(fh == 0.0)
		return x2;
	if(fl < 0.0) {
		xl = x1;
		xh = x2;
	} else {
		xh = x1;
		xl = x2;
	}
	double rts = 0.5 * (x1 + x2);
	double dxold = fabs(x2 - x1);
	double dx = dxold;
	double f = funcd(rts);
	double df = funcd.df(rts);
	for(int j = 0; j < MAXIT; j++) {
		if((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (fabs(2.0 * f) > fabs(dxold * df))) {
			dxold = dx;
			dx = 0.5 * (xh - xl);
			rts = xl + dx;
			if(xl == rts)
				return rts;
		} else {
			dxold = dx;
			dx = f / df;
			double temp = rts;
			rts -= dx;
			if(temp == rts)
				return rts;
		}
		if(fabs(dx) < xacc)
			return rts;
		double f = funcd(rts);
		double df = funcd.df(rts);
		if(f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	throw("Maximum number of iterations exceede in rtsafe");
}

template<typename T>
vector<size_t> sort_indexes_up(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	for(size_t i = 0; i != idx.size(); ++i)
		idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
		return v[i1] < v[i2];
	});

	return idx;

}

template<typename T>
vector<size_t> sort_indexes_down(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	for(size_t i = 0; i != idx.size(); ++i)
		idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
		return v[i1] > v[i2];
	});

	return idx;

}

void polcoe(dVector &x, dVector &y, dVector &Coef);

struct LUdcmp {
	int n;
	dMatrix lu;
	iVector indx;
	double d;

	LUdcmp(dMatrix &a);
	void solve(dVector &b, dVector &x);
	void inverse(dMatrix &ainv);

	dMatrix &aref;
};

///*!
// * \brief Computes pow element-wise.
// *
// * Computes pow element-wise for ArrayXd inputs.
// *
// * \param X The base.
// * \param Y The exponent.
// * \return The value of computing @f$ X ^ Y @f$ element-wise.
// */
//ArrayXd pow(const ArrayXd &X, const ArrayXd &Y);

/*!
 * \brief One-dimensional linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The point where a value is required.
 * \return The interpolated value.
 */
double linear_interp(const dVector& X, const dVector& Y, double x);

/*!
 * \brief One-dimensional linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points. RowVector version.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The point where a value is required.
 * \return The interpolated value.
 */
double linear_interp(const RowVector & X, const RowVector & Y, double x);

/*!
 * \brief One-dimensional linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points. RowVector version, with several
 * points to get the interpolated value.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The points where a value is required.
 * \return The interpolated values.
 */
RowVector linear_interp(const RowVector& X, const RowVector& Y, const RowVector& x);

/*!
 * \brief One-dimensional periodic linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points.  The function is supposed to
 * behave periodic.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The point where a value is required.
 * \return The interpolated value.
 */
double periodic_linear_interp(const dVector& X, const dVector& Y, double x);

/*!
 * \brief One-dimensional periodic linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points. RowVector version.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The point where a value is required.
 * \return The interpolated value.
 */
double periodic_linear_interp(const RowVector& X, const RowVector& Y, double x);

/*!
 * \brief One-dimensional periodic linear interpolator.
 *
 * Returns the one-dimensional piecewise linear interpolation to a function
 * with given values at discrete data-points. RowVector version, with several
 * points to get the interpolated value.
 *
 * \param X The known data points.
 * \param Y The values for the known data points.
 * \param x The points where a value is required.
 * \return The interpolated values.
 */
RowVector periodic_linear_interp(const RowVector& X, const RowVector& Y, const RowVector& x);

/*!
 * \brief Copy the data from an Eigen ArrayXd to a dVector.
 *
 * Copies the data from an Eigen ArrayXd to a dVector.  Useful for using
 * an Eigen array with the dVector interpolators.
 *
 * \param X The ArrayXd to copy.
 * \return dVector with the copied data.
 */
dVector CopyArrayXdTodVector(const ArrayXd & X);

/*!
 * \brief Copy the data from a dVector to an Eigen ArrayXd.
 *
 * Copies the data from a dVector to an Eigen ArrayXd.  Useful for using
 * an Eigen array with the dVector interpolators.
 *
 * \param X The dVector to copy.
 * \return Eigen ArrayXd with the copied data.
 */
ArrayXd CopydVectorToArrayXd(const dVector & X);

/*!
 * \brief Copy the data from a raw array to an Eigen ArrayXd.
 *
 * Copies the data from a raw array to an Eigen ArrayXd.
 *
 * \param X The raw array to copy.
 * \param size Array size.
 * \return Eigen ArrayXd with the copied data.
 */
ArrayXd CopyRawArrayToArrayXd(double *X, int size);

/*!
 * \brief Copy the data from a raw array to a ColVector.
 *
 * Copies the data from a raw array to a ColVector.
 *
 * \param X The raw array to copy.
 * \param size Array size.
 * \return ColVector with the copied data.
 */
ColVector CopyRawArrayToColVector(double *X, int size);

/*!
 * \brief Copy the data from a raw array to a RowVector.
 *
 * Copies the data from a raw array to RowVector.
 *
 * \param X The raw array to copy.
 * \param size Array size.
 * \return RowVector with the copied data.
 */
RowVector CopyRawArrayToRowVector(double *X, int size);

/*!
 * \brief Copy the data from a raw array to an Eigen ArrayXXd.
 *
 * Copies the data from a raw array to an Eigen ArrayXd.
 *
 * \param X The raw array to copy.
 * \param x_size First axis size.
 * \param y_size Second axis size.
 * \return Eigen ArrayXXd with the copied data.
 */
ArrayXXd CopyRawArrayToArrayXXd(double **X, int x_size, int y_size);

/*!
 * \brief Copy the data from a raw array to a dArray.
 *
 * Copies the data from a raw array to a dArray.
 *
 * \param X The raw array to copy.
 * \param x_size First axis size.
 * \param y_size Second axis size.
 * \return dArray with the copied data.
 */
dArray CopyRawArrayTodArray(double **X, int x_size, int y_size);

/*!
 * \brief A linear interpolator function typedef.
 */
typedef double (*linear_interp_fun)(const RowVector&, const RowVector&, double);

#endif
