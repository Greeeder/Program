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

/**
 * \file BoundaryCondition.cpp
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
 * This file defines a basic boundary condition.
 */

#include "BoundaryCondition.hpp"

TBoundaryCondition::TBoundaryCondition() {
	FCurrentTime = 0.;
	FPressure = 1E5;
	FSpeed = 0.;
	FSpeedOfSound = __cons::ARef;
	FTemperature = __cons::TRef;
	FFlux_0.setZero(3, 1);
}

TBoundaryCondition::TBoundaryCondition(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end):
	TBoundaryCondition() {
	FPipe_0 = pipe_0.get();
	FCurrentTime = pipe_0->getCurrentTime();
	FPipeEnd_0 = pipe_end;
	if(FPipeEnd_0 == nmLeft) {
		FCell_0 = 0;
	} else {
		FCell_0 = pipe_0->getPressure().cols() - 1;
	}
}

TBoundaryCondition::TBoundaryCondition(const TFlowObject_ptr & pipe_0, nmPipeEnd pipe_end) :
TBoundaryCondition() {
	FlowObj = pipe_0.get();
	FCurrentTime = pipe_0->getCurrentTime();
	FPipeEnd_0 = pipe_end;
	if (FPipeEnd_0 == nmLeft) {
		FCell_0 = 0;
	}
	else {
		FCell_0 = pipe_0->getNin() - 1;
	}
}

TBoundaryCondition::TBoundaryCondition(const TFlowObject_ptr & zerod, int con_id) :
TBoundaryCondition() {
	FlowObj = zerod.get();
	FCurrentTime = zerod->getCurrentTime();
	FCon_Index = con_id;
}

TBoundaryCondition::~TBoundaryCondition() {}

double TBoundaryCondition::getArea() const {
	return FArea;
}

double TBoundaryCondition::getBeta() const {
	return FBeta;
}

double TBoundaryCondition::getEnthalpyFlow() const {
	return FFlux_0(2) * getArea();
}

double TBoundaryCondition::getEntropy() const {
	return FEntropy;
}

double TBoundaryCondition::getLambda() const {
	return FLambda;
}

double TBoundaryCondition::getMassFlow() const {
	return FFlux_0(0) * getArea();
}

double TBoundaryCondition::getPressure() const {
	return FPressure;
}

double TBoundaryCondition::getSpeed() const {
	return FSpeed;
}

double TBoundaryCondition::getTemperature() const {
	return FTemperature;
}

void TBoundaryCondition::setBeta(double beta) {
	FBeta = beta;
}

void TBoundaryCondition::setEntropy(double entropy) {
	FEntropy = entropy;
}

void TBoundaryCondition::setLambda(double lambda) {
	FLambda = lambda;
}

void TBoundaryCondition::computePTUFromCharacteristics(double t, double dt) {
	FCurrentTime = t;
	updateBCCharacteristics(t, dt);
	/* TODO: Gamma(T)! */
	double g = __Gamma::G;
	FSpeedOfSound = (getLambda() + getBeta()) / 2. * __cons::ARef;
	FSpeed = (getLambda() - getBeta()) / __Gamma::G1(g) * __cons::ARef;
	FPressure = __units::BarToPa(pow(FSpeedOfSound / __cons::ARef / FEntropy, __Gamma::G4(g)));
	FTemperature = pow2(FSpeedOfSound) / g / __R::Air;
}

