/**
 * @file T1DAcousticTurbine.cpp
 * @author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * @author Francisco Jose Arnau Martinez <farnau@mot.upv.es>
 * @version 0.1
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
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 * The T1DAcousticTurbine class represents the acoustics of a classic
 * one-dimensional turbine.
 *
 * This file defines a T1DAcousticTurbine class.
 */

#include "T1DAcousticTurbine.h"

T1DAcousticTurbine::T1DAcousticTurbine(Pipe_ptr inlet, Volute_ptr volute,
	Pipe_ptr outlet) {
	FInlet = inlet;
	FOutlet = outlet;
	FVolute = volute;
	FMembers.push_back(inlet);
	FMembers.push_back(volute);
	FMembers.push_back(outlet);
}

T1DAcousticTurbine::T1DAcousticTurbine() {}

void T1DAcousticTurbine::setVoluteNozzleEffectiveArea(double A, unsigned int i) {
	dynamic_cast<TPressureBC*>(FVolute->getRightBC())->setCD(A
		/ dynamic_cast<TPressureBC*>(FVolute->getRightBC())->getArea());
}

void T1DAcousticTurbine::setVoluteOutletConditions(double p, double T) {
	dynamic_cast<TPressureBC*>(FVolute->getRightBC())->setPT(p, T);
}

RowVector T1DAcousticTurbine::VoluteOutletMassFlow() const {
	return RowVector::Constant(1, FVolute->getMassFlow(100.));
}

RowVector T1DAcousticTurbine::VoluteOutletp() const {
	return FVolute->getPressure().tail(1);
}

RowVector T1DAcousticTurbine::VoluteOutletp0() const {
	return FVolute->getTotalPressure().tail(1);
}

RowVector T1DAcousticTurbine::VoluteOutletT() const {
	return FVolute->getTemperature().tail(1);
}

RowVector T1DAcousticTurbine::VoluteOutletT0() const {
	return FVolute->getTotalTemperature().tail(1);
}
