/*--------------------------------------------------------------------------------*\
 *==========================|
 *\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 * \\ |  X  | //  W ave     |
 *  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 *   \\/   \//    M odel    |
 *----------------------------------------------------------------------------------
 *License
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
 * ExtrapMapFunctions.h
 *
 * Created on: 29/7/2015
 *     Author: farnau
 *
 \*--------------------------------------------------------------------------------*/

#ifndef SOURCE_TURBOCOMPRESSOR_EXTRAPMAPFUNCTIONS_H_
#define SOURCE_TURBOCOMPRESSOR_EXTRAPMAPFUNCTIONS_H_

#include "Globales.h"
#include "Correlations.h"

struct stElipse: stCorrelation {

	stElipse() :
		stCorrelation() {

	}

	void Fun(double & x, dVector & a, double&y, dVector & dyda) {

		if(x >= a[1]) {
			y = 0.;

			dyda[0] = 0.;
			dyda[1] = 0.;
			dyda[2] = 0.;
		} else {
			y = val(x, a);

			dyda[0] = -(a[2] * ((2 * a[0] - 2 * x) / pow2(a[0] - a[1]) - (2 * pow2(a[0] - x)) / pow3(a[0] - a[1]))) / (2 * sqrt(
						  1 - pow2(a[0] - x) / pow2(a[0] - a[1])));
			dyda[1] = -(a[2] * pow2(a[0] - x)) / (pow3(a[0] - a[1]) * sqrt(1 - pow2(a[0] - x) / pow2(a[0] - a[1])));
			dyda[2] = sqrt(1 - pow2(fabs(x - a[0]) / (a[1] - a[0])));
		}
	}

	double val(double & x, dVector & a) {
		if(x > a[1])
			return 0.;
		else
			return a[2] * sqrt(1 - pow2(fabs(x - a[0]) / (a[1] - a[0])));
	}

};

struct stWmax: stCorrelation {

	stWmax() {

		p.resize(2);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = val(x, a);

		dyda[0] = 1.;
		dyda[1] = x * 1e-5;

	}

	double val(double & x, dVector &a) {

		return a[0] + a[1] * x * 1e-5;
	}

	double val(double & x) {

		return val(x, p);
	}

};

struct stWzs1: stCorrelation {

	stWzs1() {

		p.resize(2);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = val(x, a);

		dyda[0] = pow(x * 1e-5, a[1]);
		dyda[1] = y * log(x * 1e-5);

	}

	double val(double & x, dVector &a) {

		return a[0] * pow(x * 1e-5, a[1]);
	}

	double val(double & x) {

		return val(x, p);
	}

};

struct stFparm: stCorrelation {

	stFparm() {

		p.resize(3);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = val(x, a);

		if(x <= 0) {
			dyda[0] = 0.;
			dyda[1] = 0.;
			dyda[2] = 0.;
		} else {
			double kk5 = pow(x * 1e-5 / a[1], a[2]);
			dyda[0] = (1 - y) * kk5;
			dyda[1] = -a[0] * a[2] * x * 1e-5 * (1 - y) * kk5 / (x * 1e-5 / a[1]) / pow2(a[1]);
			dyda[2] = a[0] * (1 - y) * log(x * 1e-5 / a[1]) * kk5;
		}

	}

	double val(double & x, dVector &a) {

		return 1 - exp(-a[0] * pow(x * 1e-5 / a[1], a[2]));
	}

	double val(double & x) {

		return val(x, p);
	}

};

struct stRCzs1: stCorrelation {

	stRCzs1() {

		p.resize(2);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = val(x, a);

		dyda[0] = pow(x * 1e-5, a[1]);
		dyda[1] = y * log(x * 1e-5);

	}

	double val(double & x, dVector &a) {

		return 1 + a[0] * pow(x * 1e-5, a[1]);
	}

	double val(double & x) {

		return val(x, p);
	}

};

struct stLeufven: stCorrelation {

	stLeufven() {

		p.resize(11);
	}

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {

		double n = x[0] * 1e-5;

		double c1 = a[0] + a[1] * n;
		double c2 = a[2] + a[3] * pow(n, a[4]);
		double wmax = a[5] + a[6] * n;
		if(x[1] >= wmax) {
			y = 0;
			for(int i = 0; i < 10; i++) {
				dyda[i] = 0.;
			}
		} else {
			double wzsl = a[7] * pow(n, a[8]);
			double rczsl = 1 + a[9] * pow(n, a[10]);

			double dm = (x[1] - wzsl);
			if(dm < 0.)
				dm = 0.;
			if(dm > 1.)
				dm = 1.;

			double wrel = dm / (wmax - wzsl);
			double kk1 = pow(wrel, c1);
			double kk2 = pow(1. - kk1, 1 / c2 - 1);
			double kk3 = pow(1. - kk1, 1 / c2);
			double kk4 = pow(wrel, c1 - 1);

			y = rczsl * kk3;

			if(wrel <= 0) {
				dyda[0] = 0.;
				dyda[1] = 0.;
			} else {
				dyda[0] = -(log(wrel) * kk2 * rczsl * kk1) / c2;
				dyda[1] = -(n * log(wrel) * kk2 * rczsl * kk1) / c2;
			}
			if(kk1 >= 1.) {
				dyda[2] = 0.;
				dyda[3] = 0.;
				dyda[4] = 0.;
			} else {
				dyda[2] = -(log(1. - kk1) * kk3 * rczsl) / pow2(c2);
				dyda[3] = -(pow(n, a[4]) * log(1. - kk1) * kk3 * rczsl) / pow2(c2);
				dyda[4] = -(a[3] * pow(n, a[4]) * log(n) * log(1. - kk1) * kk3 * rczsl) / pow2(c2);
			}
			dyda[5] = (kk2 * rczsl * kk4 * dm * c1) / (c2 * pow2(wmax - wzsl));
			dyda[6] = (n * kk2 * rczsl * kk4 * dm * c1) / (c2 * pow2(wmax - wzsl));
			dyda[7] = (kk2 * rczsl * kk4 * (pow(n, a[8]) / (wmax - wzsl) - (pow(n, a[8]) * dm) / pow2(wmax - wzsl)) * c1) / c2;
			dyda[8] = (kk2 * ((wzsl * log(n)) / (wmax - wzsl) - (wzsl * log(n) * dm) / pow2(wmax - wzsl)) * rczsl * kk4 * c1) / c2;
			dyda[9] = pow(n, a[10]) * kk3;
			dyda[10] = a[9] * pow(n, a[10]) * log(n) * kk3;
		}

	}

	double val_rc(dVector & x, dVector &a) {

		double n = x[0] * 1e-5;

		double c1 = a[0] + a[1] * n;
		double c2 = a[2] + a[3] * pow(n, a[4]);

		double wmax = a[5] + a[6] * n;
		double wzsl = a[7] * pow(n, a[8]);
		double rczsl = 1 + a[9] * pow(n, a[10]);

		if(x[1] < wzsl)
			return rczsl;
		else if(x[1] > wmax)
			return 0.;
		else
			return rczsl * pow(1 - pow((x[1] - wzsl) / (wmax - wzsl), c1), 1 / c2);
	}

	double val_rc(dVector & x) {

		return val_rc(x, p);
	}

	double val_m(dVector & x, dVector &a) {

		double n = x[0] * 1e-5;

		double c1 = a[0] + a[1] * n;
		double c2 = a[2] + a[3] * pow(n, a[4]);

		double wmax = a[5] + a[6] * n;
		double wzsl = a[7] * pow(n, a[8]);
		double rczsl = 1 + a[9] * pow(n, a[10]);

		if(x[1] > rczsl)
			return 0.;
		else
			return wzsl + (wmax - wzsl) * pow(1 - pow(x[1] / rczsl, c2), 1 / c1);
	}

	double val_m(dVector & x) {

		return val_m(x, p);
	}

	void check_dif(double m, double n) {
		dVector c;
		dVector dyda;
		dVector dyda2;
		dVector in;
		double y;
		double y1, y2, x1, x2;

		in.resize(2);
		dyda.resize(p.size());

		in[0] = n;
		in[1] = m;

		c = p;

		Fun(in, c, y, dyda);
		for(int i = 0; i < p.size(); i++) {
			c = p;
			x1 = p[i] * 1.00001;
			x2 = p[i] * 0.99999;
			if(p[i] == 0) {
				x1 = p[i] + 0.00001;
				x2 = p[i] - 0.00001;
			}
			c[i] = x1;
			y1 = val_rc(in, c);
			c[i] = x2;
			y2 = val_rc(in, c);
			dyda2.push_back((y1 - y2) / (x1 - x2));

		}

	}

};

struct stLeufvenMod: stCorrelation {

	stLeufvenMod() {

		p.resize(12);
	}

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {

		double n = x[0] * 1e-5;

		double wmax = a[5] + a[6] * n;

		if(x[1] >= wmax) {
			y = 0;
			for(int i = 0; i < 10; i++) {
				dyda[i] = 0.;
			}
		} else {
			double c1 = a[0] + a[1] * n;
			double c2 = a[2] + a[3] * pow(n, a[4]);
			double F = 1 - exp(-a[7] * pow(n / a[8], a[9]));
			double rczsl = 1 + a[10] * pow(n, a[11]);

			double dm = (x[1] - wmax * F);
			if(dm < 0.)
				dm = 0.;
			if(dm > 1.)
				dm = 1.;

			double wrel = dm / (wmax * (1 - F));
			double kk1 = pow(wrel, c1);
			double kk2 = pow(1. - kk1, 1 / c2 - 1);
			double kk3 = pow(1. - kk1, 1 / c2);
			double kk4 = pow(wrel, c1 - 1);

			y = rczsl * kk3;

			if(wrel <= 0) {
				dyda[0] = 0.;
				dyda[1] = 0.;
			} else {
				dyda[0] = -(log(wrel) * kk2 * rczsl * kk1) / c2;
				dyda[1] = -(n * log(wrel) * kk2 * rczsl * kk1) / c2;
			}
			if(kk1 >= 1.) {
				dyda[2] = 0.;
				dyda[3] = 0.;
				dyda[4] = 0.;
			} else {
				dyda[2] = -(log(1. - kk1) * kk3 * rczsl) / pow2(c2);
				dyda[3] = -(pow(n, a[4]) * log(1. - kk1) * kk3 * rczsl) / pow2(c2);
				dyda[4] = -(a[3] * pow(n, a[4]) * log(n) * log(1. - kk1) * kk3 * rczsl) / pow2(c2);
			}
			dyda[5] = (kk2 * rczsl * kk4 * dm * c1) / (c2 * pow2(wmax * (1 - F)));
			dyda[6] = (n * kk2 * rczsl * kk4 * dm * c1) / (c2 * pow2(wmax * (1 - F)));
			if(n <= 0) {
				dyda[7] = 0.;
				dyda[8] = 0.;
				dyda[9] = 0.;
			} else {
				double dydF = rczsl * c1 * kk2 * kk4 * (wmax - x[1]) / (pow2(F - 1) * wmax);
				double kk5 = pow(n / a[8], a[9]);
				dyda[7] = dydF * (1 - F) * kk5;
				dyda[8] = -dydF * a[7] * a[9] * n * (1 - F) * kk5 / (n / a[8]) / pow2(a[8]);
				dyda[9] = dydF * a[7] * (1 - F) * log(n / a[8]) * kk5;
			}
			dyda[10] = pow(n, a[10]) * kk3;
			dyda[11] = a[9] * pow(n, a[10]) * log(n) * kk3;
		}

	}

	double val_rc(dVector & x, dVector &a) {

		double n = x[0] * 1e-5;

		double c1 = a[0] + a[1] * n;
		double c2 = a[2] + a[3] * pow(n, a[4]);

		double wmax = a[5] + a[6] * n;
		double wzsl = wmax * (1 - exp(-a[7] * pow(n / a[8], a[9])));
		double rczsl = 1 + a[10] * pow(n, a[11]);

		if(x[1] < wzsl)
			return rczsl;
		else if(x[1] > wmax)
			return 0.;
		else
			return rczsl * pow(1 - pow((x[1] - wzsl) / (wmax - wzsl), c1), 1 / c2);
	}

	double val_rc(dVector & x) {

		return val_rc(x, p);
	}

	double val_m(dVector & x, dVector &a) {

		double n = x[0] * 1e-5;

		double c1 = a[0] + a[1] * n;
		double c2 = a[2] + a[3] * pow(n, a[4]);

		double wmax = a[5] + a[6] * n;
		double wzsl = wmax * (1 - exp(-a[7] * pow(n / a[8], a[9])));
		double rczsl = 1 + a[10] * pow(n, a[11]);

		if(x[1] > rczsl)
			return 0.;
		else
			return wzsl + (wmax - wzsl) * pow(1 - pow(x[1] / rczsl, c2), 1 / c1);
	}

	double val_m(dVector & x) {

		return val_m(x, p);
	}

	void check_dif(double m, double n) {
		dVector c;
		dVector dyda;
		dVector dyda2;
		dVector in;
		double y;
		double y1, y2, x1, x2;

		in.resize(2);
		dyda.resize(p.size());

		in[0] = n;
		in[1] = m;

		c = p;

		Fun(in, c, y, dyda);
		for(int i = 0; i < p.size(); i++) {
			c = p;
			x1 = p[i] * 1.00001;
			x2 = p[i] * 0.99999;
			if(p[i] == 0) {
				x1 = p[i] + 0.00001;
				x2 = p[i] - 0.00001;
			}
			c[i] = x1;
			y1 = val_rc(in, c);
			c[i] = x2;
			y2 = val_rc(in, c);
			dyda2.push_back((y1 - y2) / (x1 - x2));

		}

	}

};

struct stMartin: stCorrelation {

	stMartin() {

		p.resize(3);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = (a[0] + a[1] * x) / (a[2] - x);

		dyda[0] = 1 / (a[2] - x);
		dyda[1] = x / (a[2] - x);
		dyda[2] = -(a[0] + a[1] * x) / pow2(a[2] - x);
	}

	double val_phi(double Psi, dVector &a) {
		return (a[2] * Psi - a[0]) / (a[1] + Psi);
	}

	double val_psi(double Phi, dVector &a) {
		return (a[0] + a[1] * Phi) / (a[2] - Phi);
	}
};

struct stParabola: stCorrelation {
	stParabola() {
		p.resize(3);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		y = a[0] * pow2(x) + a[1] * x + a[2];

		dyda[0] = pow2(x);
		dyda[1] = x;
		dyda[2] = 1.;
	}

	double val(double x) {
		return p[0] * pow2(x) + p[1] * x + p[2];
	}
};

struct stEffCurve: stCorrelation {
	stEffCurve() {
		p.resize(2);
	}

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {

		if(x > 1 - pow(0.5, a[0]) + 1e-5) {

			y = 1 - pow(1 - 2 * pow(1 - x, 1 / a[0]), a[1]);

			dyda[0] = -(2 * a[1] * log(1 - x) * pow(1 - 2 * pow(1 - x, 1 / a[0]), a[1] - 1) * pow(1 - x, 1 / a[0])) / pow2(a[0]);
			dyda[1] = -log(1 - 2 * pow(1 - x, 1 / a[0])) * pow(1 - 2 * pow(1 - x, 1 / a[0]), a[1]);
		} else if(x < 1 - pow(0.5, a[0]) - 1e-5) {
			y = 1 - pow(-1 + 2 * pow(1 - x, 1 / a[0]), a[1]);

			dyda[0] = (2 * a[1] * log(1 - x) * pow(-1 + 2 * pow(1 - x, 1 / a[0]), a[1] - 1) * pow(1 - x, 1 / a[0])) / pow2(a[0]);
			dyda[1] = -log(-1 + 2 * pow(1 - x, 1 / a[0])) * pow(-1 + 2 * pow(1 - x, 1 / a[0]), a[1]);
		} else {
			y = 1;

			dyda[0] = 0.;
			dyda[1] = 0.;
		}
	}

	double val(double &x, dVector &a) {

		if(x > 1 - pow(0.5, a[0])) {
			return 1 - pow(1 - 2 * pow(1 - x, 1 / a[0]), a[1]);
		} else {
			return 1 - pow(-1 + 2 * pow(1 - x, 1 / a[0]), a[1]);
		}
	}

};

#endif /* SOURCE_TURBOCOMPRESSOR_EXTRAPMAPFUNCTIONS_H_ */
