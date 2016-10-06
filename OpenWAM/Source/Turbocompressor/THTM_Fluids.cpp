/*--------------------------------------------------------------------------------*\
|==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 | License
 |
 |	This file is part of OpenWAM.
 |
 |	OpenWAM is free software: you can redistribute it and/or modify
 |	it under the terms of the GNU General Public License as published by
 |	the Free Software Foundation, either version 3 of the License, or
 |	(at your option) any later version.
 |
 |	OpenWAM is distributed in the hope that it will be useful,
 |	but WITHOUT ANY WARRANTY; without even the implied warranty of
 |	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 |	GNU General Public License for more details.
 |
 |	You should have received a copy of the GNU General Public License
 |	along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 |
 \*--------------------------------------------------------------------------------*/

#include "THTM_Fluids.h"

stHTMair::stHTMair() :
	stHTM_Fluid() {
	AF = 0.;
	Hum = 0.;
}

double stHTMair::fun_mu(double T) {
	dVector A;
	A.resize(7);
	A[0] = 0.0000997217;
	A[1] = 0.000000214683;
	A[2] = -0.000000000226615;
	A[3] = 1.1354e-13;
	A[4] = -0.0000218018;
	A[5] = -0.001168773;
	A[6] = 0;

	if(T > 1500.) {
		T = 1500.;
	}

	return Property(T, A);
}

double stHTMair::fun_k(double T) {
	dVector A;
	A.resize(7);
	A[0] = -0.028882591;
	A[1] = 0.0000746349;
	A[2] = -0.000000020381;
	A[3] = 6.66102E-13;
	A[4] = 0.005818895;
	A[5] = 0.407241757;
	A[6] = 0;

	if(T > 1500.) {
		T = 1500.;
	}

	return Property(T, A);
}

double stHTMair::fun_Cp(double T) {

	double Cp = 0., x = 0.;

	if(AF == 0.0) {
		x = Hum / (1 + Hum);
		Cp = (1 - x) * fun_Cp_AirS(T) + x * fun_Cp_W(T);
	} else {
		x = 1 / AF;
		Cp = (1 - x) * fun_Cp_AirS(T) + x * fun_Cp_Gas(T);
	}
	return Cp;
}

double stHTMair::fun_Cp_AirS(double T) {
	if(T > 1500.) {
		T = 1500.;
	}
	return -10.4199 * sqrt(T) + 2809.87 - 67227.1 / sqrt(T) + 917124.4 / T - 4174853.6 / pow150(T);
}

double stHTMair::fun_Cp_W(double T) {
	if(T > 1500.) {
		T = 1500.;
	}
	return -41.9055 * sqrt(T) + 10447.493 - 382002.49 / sqrt(T) + 6456647.7 / T - 37951136.5 / pow150(T);
}

double stHTMair::fun_Cp_Gas(double T) {
	if(T > 1500.) {
		T = 1500.;
	}
	return 926.554 + 0.43045 * T - 0.0001125 * pow2(T) + 0.000000008979 * pow3(T);
}

double stHTMair::fun_g(double T) {
	return fun_Cp(T) / (fun_Cp(T) - fun_R());
}

double stHTMair::fun_R() {
	double x = 0., R_Air = 0.;

	if(AF == 0.0) {
		x = Hum / (1 + Hum);
		R_Air = (1 - x) * 287 + x * 462.1762;
	} else {
		x = 1 / AF;
		R_Air = (1 - x) * 287 + x * 285.4;
	}
	return R_Air;
}

double stHTMair::fun_rho(double p, double T) {
	return (100000 * p) / (fun_R() * T);
}

double stHTMair::fun_Pr(double T) {
	return fun_mu(T) * fun_Cp(T) / fun_k(T);
}

void stHTMair::CalcProperties(double p, double T) {
	mu = fun_mu(T);
	R = fun_R();
	Cp = fun_Cp(T);
	g = fun_g(T);
	rho = fun_rho(p, T);
	k = fun_k(T);
	Pr = mu * Cp / k;
}

stHTMoil::stHTMoil() :
	stHTM_Fluid() {
	mu_c1 = 0.0115374825;
	mu_c2 = 62.1420218;
	mu_c3 = 275.888115;
}

double stHTMoil::fun_rho(double T) {
	dVector A;
	A.resize(7);
	A[0] = -16572.74497;
	A[1] = -38.71871842;
	A[2] = 0.054641551;
	A[3] = -3.45327E-05;
	A[4] = 4398.316892;
	A[5] = 0;
	A[6] = 0;

	if(T > 550.) {
		T = 550.;
	}

	return Property(T, A);
}

double stHTMoil::fun_Cp(double T) {
	dVector A;
	A.resize(7);
	A[0] = -34040.12013;
	A[1] = -71.24978467;
	A[2] = 0.107361063;
	A[3] = -6.63966E-05;
	A[4] = 8670.542844;
	A[5] = 0;
	A[6] = 0;

	if(T > 550.) {
		T = 550.;
	}

	return Property(T, A);
}

double stHTMoil::fun_mu(double T) {
	if(T > 550.) {
		T = 550.;
	}

	return mu_c1 * exp(mu_c2 / (T - mu_c3));
}

double stHTMoil::fun_k(double T) {
	dVector A;
	A.resize(7);
	A[0] = -13.78897843;
	A[1] = -0.028315003;
	A[2] = 3.80794E-05;
	A[3] = -2.25321E-08;
	A[4] = 3.437962885;
	A[5] = 0;
	A[6] = 0;

	if(T > 550.) {
		T = 550.;
	}

	return Property(T, A);
}

double stHTMoil::fun_Pr(double T) {
	return fun_mu(T) * fun_Cp(T) / fun_k(T);
}

void stHTMoil::CalcProperties(double p, double T) {
	mu = fun_mu(T);
	Cp = fun_Cp(T);
	rho = fun_rho(T);
	k = fun_k(T);
	Pr = mu * Cp / k;
}

