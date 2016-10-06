/**
 * @file TQ2DAcousticTurbine.cpp
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
 * The TQ2DAcousticTurbine class represents the acoustics of a
 * quasi-two-dimensional turbine.
 *
 * This file defines a TQ2DAcousticTurbine class.
 */

#include "TQ2DAcousticTurbine.h"

TQ2DAcousticTurbine::TQ2DAcousticTurbine(Pipe_ptr inlet, Volute_ptr volute,
	Pipe_ptr outlet) {
	FInlet = inlet;
	FOutlet = outlet;
	FVolute = volute;
	FMembers.push_back(inlet);
	FMembers.push_back(volute);
	FMembers.push_back(outlet);
}

TQ2DAcousticTurbine::TQ2DAcousticTurbine() {}

double TQ2DAcousticTurbine::DIn(int i) const {
	return FInlet->getD()(0);
}

double TQ2DAcousticTurbine::DInTot() const {
	return FInlet->getD()(0);
}

double TQ2DAcousticTurbine::DOut() const {
	return FOutlet->getD().tail(1)(0);
}

double TQ2DAcousticTurbine::ExpRatio() const {
	return P30() / P4();
}

double TQ2DAcousticTurbine::MassIn(int i) const {
	return MassIn();
}

double TQ2DAcousticTurbine::MassIn() const {
	return FInlet->getLeftBC()->getMassFlow();
}

double TQ2DAcousticTurbine::MassOut() const {
	return FOutlet->getRightBC()->getMassFlow();
}

double TQ2DAcousticTurbine::OutletNozzlep() const {
	return FOutlet->getPressure((unsigned int)0);
}


double TQ2DAcousticTurbine::P3() const {
	return __units::PaToBar(FInlet->getPressure((unsigned int)0));
}

double TQ2DAcousticTurbine::P30() const {
	return __units::PaToBar(FInlet->getTotalPressure((unsigned int)(0)));
}

double TQ2DAcousticTurbine::P30(int i) const {
	return __units::PaToBar(FInlet->getTotalPressure((unsigned int)(0)));
}

double TQ2DAcousticTurbine::P4() const {
	return __units::PaToBar(FOutlet->getPressure().tail(1)(0));
}

double TQ2DAcousticTurbine::P40() const {
	return __units::PaToBar(FOutlet->getTotalPressure().tail(1)(0));
}

double TQ2DAcousticTurbine::R() const {
	return FInlet->getR((unsigned int)(0));
}

double TQ2DAcousticTurbine::RotorInletMassFlow() const {
	return FOutlet->getLeftBC()->getMassFlow();
}

void TQ2DAcousticTurbine::setRotorOutletEffectiveArea(double A) {
	dynamic_cast<TPressureBC*>(FOutlet->getLeftBC())->setCD(A
		/ dynamic_cast<TPressureBC*>(FOutlet->getLeftBC())->getArea());
}

void TQ2DAcousticTurbine::setRotorInletConditions(double p, double T) {
	dynamic_cast<TPressureBC*>(FOutlet->getLeftBC())->setPT(p, T);
}

void TQ2DAcousticTurbine::setVoluteNozzleEffectiveArea(double A, unsigned int i) {
	FVolute->setLateralNozzleArea(A, i);
}

void TQ2DAcousticTurbine::setVoluteOutletConditions(double p, double T) {
	dynamic_cast<TLateralNozzles*>(
		FVolute->getExtraSources())->setOutletConditions(p, T);
}

double TQ2DAcousticTurbine::SIn() const {
	return FInlet->getArea((unsigned int)(0));
}

double TQ2DAcousticTurbine::SOut() const {
	return FOutlet->getArea().tail(1)(0);
}

double TQ2DAcousticTurbine::T3() const {
	return FInlet->getTemperature((unsigned int)0);
}

double TQ2DAcousticTurbine::T30() const {
	return FInlet->getTotalTemperature((unsigned int)0);
}

double TQ2DAcousticTurbine::T30(int i) const {
	return FInlet->getTotalTemperature((unsigned int)0);
}

double TQ2DAcousticTurbine::T4() const {
	return FOutlet->getTemperature().tail(1)(0);
}

double TQ2DAcousticTurbine::T40() const {
	return FOutlet->getTotalTemperature().tail(1)(0);
}

RowVector TQ2DAcousticTurbine::VoluteOutletMassFlow() const {
	return FVolute->getLateralMassFlow();
}

RowVector TQ2DAcousticTurbine::VoluteOutletp() const {
	return FVolute->getPressure();
}

RowVector TQ2DAcousticTurbine::VoluteOutletp0() const {
	return FVolute->getTotalPressure();
}

RowVector TQ2DAcousticTurbine::VoluteOutletT() const {
	return FVolute->getTemperature();
}

RowVector TQ2DAcousticTurbine::VoluteOutletT0() const {
	return FVolute->getTotalTemperature();
}
