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

/**
 * @file Math_wam.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 *
 * @section LICENSE
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
 * @section DESCRIPTION
 * This file defines several auxiliary math functions.
 */

// ---------------------------------------------------------------------------
#include "Math_wam.h"

double Interpola(double vizq, double vder, double axid, double xif) {
	return vizq + ((vder - vizq) / axid) * xif;
}

double QuadraticEqP(double A, double B, double C) {
	return (-B + Sqrt(B * B - 4 * A * C)) / 2 / A;
}

double QuadraticEqN(double A, double B, double C) {
	return (-B - Sqrt(B * B - 4 * A * C)) / 2 / A;
}

void Poly2(dVector &X, dVector &Y, dVector &A) {
	double x01 = X[0] - X[1];
	double x02 = X[0] - X[2];
	double x12 = X[1] - X[2];

	double a0 = Y[0] / (x01 * x02);
	double a1 = Y[1] / (-x01 * x12);
	double a2 = Y[2] / (x02 * x12);

	A[0] = a0 + a1 + a2;
	A[1] = -(a0 * (X[1] + X[2]) + a1 * (X[0] + X[2]) + a2 * (X[0] + X[1]));
	A[2] = a0 * X[1] * X[2] + a1 * X[0] * X[2] + a2 * X[0] * X[1];
}

//double Max(const double &a, const float &b) {
//	return b > a ? (b) : double(a);
//}
//
//float Max(const float &a, const double &b) {
//	return b > a ? float(b) : (a);
//}
//
//double Min(const double &a, const float &b) {
//	return b < a ? (b) : double(a);
//}
//
//float Min(const float &a, const double &b) {
//	return b < a ? float(b) : (a);
//}

float Sign(const float &a, const double &b) {
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

float Sign(const double &a, const float &b) {
	return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

void polcoe(dVector &x, dVector &y, dVector &Coef) {

	dMatrix A;
	DimdMatrix(A, x.size(), x.size());
	for(int i = 0; i < x.size(); i++) {
		for(int j = 0; j < x.size(); j++) {
			A[i][j] = pow(x[i], j);
		}
	}

	LUdcmp CoefSolver(A);

	CoefSolver.solve(y, Coef);

}

stPolar::stPolar() {
}

stPolar::stPolar(double X, double Y) {
	Mod = Sqrt(X * X + Y * Y);
	Ang = atan(Y / X);
}

void stPolar::operator()(double X, double Y) {
	Mod = Sqrt(X * X + Y * Y);
	Ang = atan(Y / X);
}

stRectan::stRectan() {
}

stRectan::stRectan(double Ang, double Mod) {
	X = Mod * cos(Ang);
	Y = Mod * sin(Ang);
}

void stRectan::operator()(double Ang, double Mod) {
	X = Mod * cos(Ang);
	Y = Mod * sin(Ang);
}

Base_interp::Base_interp() :
	n(0), mm(0), jsav(0), cor(0), xx(), yy() {
}

Base_interp::Base_interp(dVector &x, const double *y, int m) :
	n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
	dj = Max(1, (int) pow025((double) n));
}

double Base_interp::interp(double x) {
	int jlo = cor ? hunt(x) : locate(x);
	return rawinterp(jlo, x);
}

int Base_interp::locate(const double x) {
	int ju = 0, jm = 0, jl = 0;
	if(n < 2 || mm < 2 || mm > n)
		throw("locate size error");
	bool ascnd = (xx[n - 1] >= xx[0]);
	jl = 0;
	ju = n - 1;
	while(ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if((x >= xx[jm]) == ascnd) {
			jl = jm;
		} else {
			ju = jm;
		}
	}
	cor = abs(jl - jsav) > dj ? 0 : 1;
	jsav = jl;
	return Max(0, Min(n - mm, jl - ((mm - 2) >> 1)));

}

int Base_interp::hunt(const double x) {
	int jl = jsav, jm, ju, inc = 1;
	if(n < 2 || mm < 2 || mm > n)
		throw("locate size error");
	bool ascnd = (xx[n - 1] >= xx[0]);
	if(jl < 0 || jl > n - 1) {
		jl = 0;
		ju = n - 1;
	} else {
		if((x >= xx[jl]) == ascnd) {
			for(;;) {
				ju = jl + inc;
				if(ju >= n) {
					ju = n - 1;
					break;
				} else if((x < xx[ju]) == ascnd)
					break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for(;;) {
				jl = jl - inc;
				if(jl <= 0) {
					jl = 0;
					break;
				} else if((x >= xx[jl]) == ascnd)
					break;
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while(ju - jl > 1) {
		jm = (ju + jl) >> 1;
		if((x >= xx[jm]) == ascnd) {
			jl = jm;
		} else {
			ju = jm;
		}
	}
	cor = abs(jl - jsav) > dj ? 0 : 1;
	jsav = jl;
	return Max(0, Min(n - mm, jl - ((mm - 2) >> 1)));
}

Linear_interp::Linear_interp() :
	Base_interp() {
}

Linear_interp::Linear_interp(dVector &xv, dVector &yv) :
	Base_interp(xv, &yv[0], 2) {
}

void Linear_interp::operator()(dVector & xv, dVector & yv) {
	xx = &xv[0];
	yy = &yv[0];
	n = xv.size();
	mm = 2;
	jsav = 0;
	cor = 0;
	dj = Max(1, (int) pow025((double) n));
}

double Linear_interp::rawinterp(int j, double x) {
	if(xx[j] == xx[j + 1]) {
		return yy[j];
	} else {
		return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j])) * (yy[j + 1] - yy[j]);
	}
}

Hermite_interp::Hermite_interp() :
	Base_interp(), y2() {
}
;

Hermite_interp::Hermite_interp(dVector &xv, dVector &yv) :
	Base_interp(xv, &yv[0], 2), y2(xv.size()) {
	sety2(&xv[0], &yv[0]);
}

void Hermite_interp::operator()(dVector & xv, dVector & yv) {
	xx = &xv[0];
	yy = &yv[0];
	n = xv.size();
	mm = 2;
	jsav = 0;
	cor = 0;
	y2.resize(n);
	dj = Max(1, (int) pow025((double) n));
	sety2(&xv[0], &yv[0]);
}

void Hermite_interp::sety2(const double *xv, const double *yv) {
	double DeltaK = 0., AlphaK = 0., BetaK = 0., TauK = 0.;

	for(int i = 1; i < n - 1; ++i) {
		y2[i] = (yv[i] - yv[i - 1]) / 2. / (xv[i] - xv[i - 1]) + (yv[i + 1] - yv[i]) / 2. / (xv[i + 1] - xv[i]);
	}
	y2[0] = (yv[1] - yv[0]) / (xv[1] - xv[0]);
	y2[n - 1] = (yv[n - 1] - yv[n - 2]) / (xv[n - 1] - xv[n - 2]);

	for(int i = 0; i < n - 1; i++) {
		DeltaK = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]);
		if(DeltaK == 0) {
			y2[i] = 0;
			y2[i + 1] = 0;
		} else {
			AlphaK = y2[i] / DeltaK;
			BetaK = y2[i + 1] / DeltaK;
			if(BetaK * BetaK + AlphaK * AlphaK > 9) {
				TauK = 3 / Sqrt(BetaK * BetaK + AlphaK * AlphaK);
				y2[i] = TauK * AlphaK * DeltaK;
				y2[i + 1] = TauK * BetaK * DeltaK;
			}
		}
	}
}

double Hermite_interp::rawinterp(int j, double x) {
	double ret_val, h00, h10, h01, h11, t2, t3, t, h;
	// int k=0;

	if(x <= xx[j]) {
		ret_val = yy[j];
	} else if(x >= xx[j + 1]) {
		ret_val = yy[j + 1];
	} else {
		h = (xx[j + 1] - xx[j]);
		t = (x - xx[j]) / h;
		t2 = t * t;
		t3 = t2 * t;
		h00 = 2 * t3 - 3 * t2 + 1;
		h10 = t3 - 2 * t2 + t;
		h01 = -2 * t3 + 3 * t2;
		h11 = t3 - t2;
		ret_val = h00 * yy[j] + h * h10 * y2[j] + h01 * yy[j + 1] + h * h11 * y2[j + 1];
	}
	return ret_val;
}

Step_interp::Step_interp() :
	Base_interp() {
}

Step_interp::Step_interp(dVector &xv, dVector &yv) :
	Base_interp(xv, &yv[0], 2) {
}

void Step_interp::operator()(dVector & xv, dVector & yv) {
	xx = &xv[0];
	yy = &yv[0];
	n = xv.size();
	mm = 2;
	jsav = 0;
	cor = 0;
	dj = Max(1, (int) pow025((double) n));
}

double Step_interp::rawinterp(int j, double x) {
	if(xx[j] == xx[j + 1]) {
		return yy[j];
	} else {
		return yy[j];
	}
}

LUdcmp::LUdcmp(dMatrix &a) :
	n(a.size()), lu(a), aref(a), indx(n) {

	const double TINY = 1.0e-40;
	int i, imax, j, k;
	double big, temp;
	dVector vv(n);
	d = 1.0;
	for(i = 0; i < n; i++) {
		big = 0.0;
		for(j = 0; j < n; j++)
			if((temp = fabs(lu[i][j])) > big)
				big = temp;
		if(big == 0.0)
			throw("Singular matrix in LUdcmp");
		vv[i] = 1.0 / big;
	}
	for(k = 0; k < n; k++) {
		big = 0.0;
		for(i = k; i < n; i++) {
			temp = vv[i] * fabs(lu[i][k]);
			if(temp > big) {
				big = temp;
				imax = i;
			}
		}
		if(k != imax) {
			for(j = 0; j < n; j++) {
				temp = lu[imax][j];
				lu[imax][j] = lu[k][j];
				lu[k][j] = temp;
			}
			d = -d;
			vv[imax] = vv[k];
		}
		indx[k] = imax;
		if(lu[k][k] == 0.0)
			lu[k][k] = TINY;
		for(i = k + 1; i < n; i++) {
			temp = lu[i][k] /= lu[k][k];
			for(j = k + 1; j < n; j++)
				lu[i][j] -= temp * lu[k][j];

		}
	}
}

void LUdcmp::solve(dVector &b, dVector &x) {
	int i, ii = 0, ip, j;
	double sum = 0.;
	if(b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for(i = 0; i < n; i++)
		x[i] = b[i];
	for(i = 0; i < n; i++) {
		ip = indx[i];
		sum = x[ip];
		x[ip] = x[i];
		if(ii != 0)
			for(j = ii - 1; j < i; j++)
				sum -= lu[i][j] * x[j];
		else if(sum != 0.0)
			ii = i + 1;
		x[i] = sum;

	}
	for(i = n - 1; i >= 0; i--) {
		sum = x[i];
		for(j = i + 1; j < n; j++)
			sum -= lu[i][j] * x[j];
		x[i] = sum / lu[i][i];
	}

}

//ArrayXd pow(const ArrayXd &X, const ArrayXd &Y) {
//	ArrayXd Z = ArrayXd(X.size());
//	for(ArrayXd::Index i = 0; i < X.size(); i++) {
//		Z(i) = pow(X(i), Y(i));
//	}
//	return Z;
//}

double linear_interp(const dVector& X, const dVector& Y, double x) {
	Uint i = 0, i_0 = 0, i_1 = 0;
	double x_0 = 0, x_1 = 0, y_0 = 0, y_1 = 0;

	i_1 = std::distance(X.begin(), std::upper_bound(X.begin(), X.end(), x));
	i_0 = i_1 - 1;
	y_0 = Y[i_0];
	y_1 = Y[i_1];
	x_0 = X[i_0];
	x_1 = X[i_1];
	return y_0 + (x - x_0) / (x_1 - x_0) * (y_1 - y_0);
}

double linear_interp(const RowVector& X, const RowVector& Y, double x) {
	if (x <= X(0)) {
		return Y(0);
	} else if (x >= X.tail(1)(0)) {
		return Y.tail(1)(0);
	}
	Uint n, i_0 = 0, i_1 = 0;
	double x_0 = 0, x_1 = 0, y_0 = 0, y_1 = 0;
	n = X.size();
	i_1 = std::distance(&X(0), std::upper_bound(&X(0), &X(n - 1), x));
	i_0 = i_1 - 1;
	y_0 = Y[i_0];
	y_1 = Y[i_1];
	x_0 = X[i_0];
	x_1 = X[i_1];
	return y_0 + (x - x_0) / (x_1 - x_0) * (y_1 - y_0);
}

RowVector linear_interp(const RowVector& X, const RowVector& Y, const RowVector& x) {
	RowVector y(x.size());
	for(auto i = 0; i < x.size(); i++) {
		y(i) = linear_interp(X, Y, x(i));
	}
	return y;
}

double periodic_linear_interp(const dVector& X, const dVector& Y, double x) {
	double x_0 = X.front();
	double x_end = X.back();
	double dx_1 = x - x_0;
	double dx_2 = x_end - x_0;
	double x_ = x - std::floor(dx_1 / dx_2) * dx_2;
	return linear_interp(X, Y, x_);
}

double periodic_linear_interp(const RowVector& X, const RowVector& Y, double x) {
	double x_0 = X(0);
	double x_end = X.tail(1)(0);
	double dx_1 = x - x_0;
	double dx_2 = x_end - x_0;
	double x_ = dx_1;
	if (dx_2 > 0) {
		x_ = x - std::floor(dx_1 / dx_2) * dx_2;
	}
	return linear_interp(X, Y, x_);
}

RowVector periodic_linear_interp(const RowVector& X, const RowVector& Y, const RowVector& x) {
	RowVector y(x.size());
	for(auto i = 0; i < x.size(); i++) {
		y(i) = periodic_linear_interp(X, Y, x(i));
	}
	return y;
}

dVector CopyArrayXdTodVector(const ArrayXd & X) {
	dVector output(X.size());
	Map<ArrayXd> Y(output.data(), X.size());
	Y = X;
	return output;
}

ArrayXd CopydVectorToArrayXd(const dVector & X) {
	ArrayXd output = Map<const ArrayXd>(X.data(), X.size());
	return output;
}

ArrayXd CopyRawArrayToArrayXd(double *X, int size) {
	ArrayXd output = Map<ArrayXd>(X, size);
	return output;
}

ColVector CopyRawArrayToColVector(double *X, int size) {
	ColVector output = Map<ColVector>(X, size);
	return output;
}

RowVector CopyRawArrayToRowVector(double *X, int size) {
	RowVector output = Map<RowVector>(X, size);
	return output;
}

ArrayXXd CopyRawArrayToArrayXXd(double **X, int x_size, int y_size) {
	ArrayXXd output(x_size, y_size);
	for(int i = 0; i < x_size; i++) {
		output.row(i) = Map<ArrayXd>(X[i], y_size);
	}
	return output;
}

dArray CopyRawArrayTodArray(double **X, int x_size, int y_size) {
	dArray output(x_size, y_size);
	for(int i = 0; i < x_size; i++) {
		output.row(i) = Map<RowVector>(X[i], y_size);
	}
	return output;
}
