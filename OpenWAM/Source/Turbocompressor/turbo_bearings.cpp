#include "turbo_bearings.hpp"
/**
 * @file turbo_bearings.cpp
 * @author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * @version 0.3.5
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
 * The TurboBearings class represent the bearing system in a turbocharger.
 *
 * It has methods to compute the power losses in the bearing system. This file
 * has the implementation of TurboBearings.
 *
 */

using namespace std;

TurboBearings::TurboBearings() {

	_L_jb = 29.80E-3;
	_R_jb = 7.73E-3 / 2.;
	_h_jb = 7.5E-6;
	_k_jb = 0.4;
	_k_tb = 0.;
	_k_A_c = 1.;
	_k_A_t = 1.;
	_A_c = 1.;
	_A_t = 1.;
	_k_m = 0.33;
	_R_tb_min = 0.1;
	_R_tb_max = 1.;
}

TurboBearings::TurboBearings(TFluid *Oil, double L_jb, double R_jb, double h_jb, double k_jb, double k_A_c,
							 double k_A_t, double A_c, double A_t, double k_m, double R_tb_min, double R_tb_max,
							 double k_tb) {

	_Oil = Oil;
	_L_jb = L_jb;
	_R_jb = R_jb;
	_h_jb = h_jb;
	_k_jb = k_jb;
	_k_tb = k_tb;
	_k_A_c = k_A_c;
	_k_A_t = k_A_t;
	_A_c = A_c;
	_A_t = A_t;
	_k_m = k_m;
	_R_tb_min = R_tb_min;
	_R_tb_max = R_tb_max;
}

double TurboBearings::h_tb(double T) {

	double g = 0.;

	g = pow(_R_tb_max, 2) / 2. * (log(_R_tb_max) - 0.5) + pow(_R_tb_min, 2) / 2. * (log(_R_tb_min) - 0.5);

	double htb = pow(_k_m * _m * abs(g) * 12. * _Oil->FunVisc(T) / abs(_k_A_c * _A_c * (_p2 - _p1) / 4. - _k_A_t * _A_t *
					 (_p3 - _p4) / 2.) / _Oil->FunRHO(T), 1. / 3.);

	if(htb < 1e-6)
		htb = 1e-6;

	return htb;
}

double TurboBearings::P_jb(double T) {

	return 2. * 4. * atan(1.0) * pow(_R_jb, 3) * _L_jb / _h_jb * _k_jb * _Oil->FunVisc(T) * pow(_n, 2);
}

double TurboBearings::P_tb(double T) {

	return 4. * atan(1.0) * (pow(_R_tb_max, 2) - pow(_R_tb_min, 2)) * pow((_R_tb_max + _R_tb_min) / 2.,
			2) * _k_tb * _Oil->FunVisc(T) * pow(_n, 2) / h_tb(T);
}

double TurboBearings::get_T_oil_m() {

	auto loop = [&](double T) {
		return (P_jb(T) + P_tb(T)) / _m / _Oil->FunCp(T) * 0.75
			+ _T1 - T;
	};

	if (_T1 <= 523.15) {
		return zbrent(loop, _T1, 523.15, 0.1);
	}
	else {
		return 523.15;
	}
}

double TurboBearings::P_oil(double T1, double n, double p1, double p2, double p3, double p4, double m) {

	double Tm = 0.;

	_T1 = T1;
	_n = n;
	_p1 = p1;
	_p2 = p2;
	_p3 = p3;
	_p4 = p4;
	_m = m;
	if(_m < 0.00055556) {
		_m = 0.00055556;
	}
	Tm = get_T_oil_m();
	return P_jb(Tm) + P_tb(Tm);

}

