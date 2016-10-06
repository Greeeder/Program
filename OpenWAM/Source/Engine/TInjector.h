/**
 * @file TInjector.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 11 de mar. de 2016
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
 * TODO Insert here the description.
 */
#ifndef SOURCE_ENGINE_TINJECTOR_H_
#define SOURCE_ENGINE_TINJECTOR_H_

#include "Globales.h"
#include <memory>

class TInjector {

	friend shared_ptr<TInjector> createInjector(xml_node node, TComponentArray_ptr fluid);

private:

	double Diameter;
	bool EndInjection;
	bool IsInjecting;
	int OrificeNumber;
	double DischargeCoefficient;
	double FuelInjected;
	double InjectedMass;
	double FPeriod;
	double FPreviousAngle;
	double FPreviousTime;
	double FRate;
	double FCurrentRate;
	double FTimeBegin;
	double FTimeEnd;
	TFluid_ptr Fuel;
	double FSOI;


public:
	TInjector();

	virtual ~TInjector();

	TFluid_ptr getFluid(){
		return Fuel;
	}

	double getFuelInjected();

	double getFuelRate(double angle, double t);

	double getFuelRate() {
		return FCurrentRate;
	}

	double getSOI(){
		return FSOI;
	}

	void setFuel(TComponentArray_ptr workfluid);



};

shared_ptr<TInjector> createInjector(xml_node node, TComponentArray_ptr fluid);

typedef shared_ptr<TInjector> TInjector_ptr;

#endif /* SOURCE_ENGINE_TINJECTOR_H_ */
