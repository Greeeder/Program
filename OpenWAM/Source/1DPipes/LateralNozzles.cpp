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


 \*--------------------------------------------------------------------------------*/

/*!
 * \file LateralNozzles.cpp
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
 * This file defines a lateral nozzles source term for a TPipe object.
 * It also includes helper functions.
 */

#include "LateralNozzles.hpp"
#include "TPipe.hpp"


TLateralNozzles::
TLateralNozzles(const Pipe_ptr & pipe):TGenericSource(pipe) {
	FArea.setZero(1, FSource.cols());
	FOutletPressure = __units::BarToPa(__cons::PRef);
	FOutletTemperature = __cons::TRef;
}

RowArray TLateralNozzles::ComputeSource(double t, double dt) {
	for (auto i = 0; i < FSource.cols(); i++) {
		FSource(0, i) = -nozzle_mass_flow(FPipe->getTotalPressure()((unsigned int)i),
			FOutletPressure, FPipe->getTotalTemperature()((unsigned int)i),
			FOutletTemperature, FPipe->getR()((unsigned int)i),
			FPipe->getcp()((unsigned int)i), FPipe->getGamma()((unsigned int)i),
			FArea(i)
		);
	}
	FSource.row(1) = FSource.row(0) * FPipe->getSpeed();
	FSource.row(2) = FSource.row(0) * FPipe->getcp() * FPipe->getTotalTemperature();
	return FSource;
}

RowVector TLateralNozzles::getArea() const {
	return FArea;
}

double TLateralNozzles::getArea(unsigned int i) const {
	return FArea(i);
}

void TLateralNozzles::setArea(const RowVector & area) {
	FArea = area;
}

void TLateralNozzles::setArea(double area, unsigned int i) {
	FArea(i) = area;
}


void TLateralNozzles::setArea(double area) {
	FArea.setConstant(area);
}

void TLateralNozzles::setOutletConditions(double p, double T)
{
	FOutletPressure = p;
	FOutletTemperature = T;
}

double nozzle_mass_flow(double p_1, double p_2, double T_1, double T_2,
	double R, double cp, double g, double A) {
	double result = 0.;
	if (p_1 >= p_2) {
		double p_throat = p_1 / pow(2 / (g + 1.), g / (1. - g));
		if (p_2 <= p_throat) {p_2 = p_throat;};
		result = A * p_1 / R / T_1 * pow(p_2 / p_1, 1. / g) *
			Sqrt(2. * cp * T_1 * (1. - pow(p_2 / p_1, (g - 1.) / g)));
	} else {
// 		result = -A * p_2 / R / T_2 * pow(p_1 / p_2, 1. / g) *
// 			Sqrt(2. * cp * T_2 * (1. - pow(p_1 / p_2, (g - 1.) / g)));
		result = 0.;
	}
	return result;
}


void set_lateral_nozzles(const Pipe_ptr & pipe) {
	Source_ptr source(new TLateralNozzles(pipe));
	pipe->setExtraSource(move(source));
}
