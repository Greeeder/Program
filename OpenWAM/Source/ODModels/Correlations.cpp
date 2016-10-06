/* --------------------------------------------------------------------------------*\
|==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
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

#include "Correlations.h"

stBMVaned::stBMVaned(double VGTtoAlpha1, double VGTtoAlpha2, double R2geom, double LTE, double Z0, double HBlade) :
	stBladeMechanism() {
	sVGTtoAlpha1 = VGTtoAlpha1;
	sVGTtoAlpha2 = VGTtoAlpha2;
	sR2geom = R2geom;
	sLTE = LTE;
	sZ0 = Z0;
	sGamma = __cons::Pi_x_2 / sZ0;
	sHBlade = HBlade;
}

double stBMVaned::Toz3(double Position) {
	double alpha = ToAlpha(Position);
	double r2 = sqrt(pow2(sR2geom - sLTE * cos(alpha)) + pow2(sLTE * sin(alpha)));
	return __cons::Pi_x_2 * r2 / ToL2(Position) / sZ0;

}

double stBMVaned::ToL2(double Position) {
	double alpha = ToAlpha(Position);
	double r2 = sqrt(pow2(sR2geom - sLTE * cos(alpha)) + pow2(sLTE * sin(alpha)));
	return 2 * r2 * sin(sGamma / 2) * cos(alpha);
}

double stBMVaned::ToAN(double Position) {
	return sZ0 * ToL2(Position) * sHBlade;
}

double stBMVaned::ToAlpha(double Position) {
	return sVGTtoAlpha1 * Position + sVGTtoAlpha2;
}

stBMVaneless::stBMVaneless(double R2geom, double HBlade) :
	stBladeMechanism() {
	sVGTtoAlpha1 = 0;
	sVGTtoAlpha2 = __cons::Pi_2;
	sR2geom = R2geom;
	sLTE = 0;
	sZ0 = 1;
	sGamma = __cons::Pi_x_2 / sZ0;
	sAN = __cons::Pi_x_2 * sR2geom * HBlade;
}

double stBMVaneless::Toz3(double Position) {
	return 1.0;
}

double stBMVaneless::ToL2(double Position) {
	return __cons::Pi_x_2 * sR2geom;
}

double stBMVaneless::ToAN(double Position) {
	return sAN;
}

double stBMVaneless::ToAlpha(double Position) {
	return sVGTtoAlpha2;
}

