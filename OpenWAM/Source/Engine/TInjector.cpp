/**
 * @file TInjector.cpp
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
#include "TInjector.h"

TInjector::TInjector() {
	// TODO Auto-generated constructor stub
	FuelInjected = 0.;
	IsInjecting = false;
	FPreviousTime = 0.;
	OrificeNumber = 1.;
	Diameter = 0.;
	DischargeCoefficient = 1.;
	FTimeBegin = 0.;
	FTimeEnd = 0.;
	InjectedMass = 0.;

}

TInjector::~TInjector() {
	// TODO Auto-generated destructor stub
}

double TInjector::getFuelInjected(){
	return FuelInjected;
}

double TInjector::getFuelRate(double angle, double t){

	double rate = 0.;

	if (angle < FPreviousAngle){
		IsInjecting = false;
		EndInjection = false;
	}
	else{
		if (!IsInjecting){
			if (angle > FSOI){
				FTimeBegin = t;
				FTimeEnd = t + FPeriod;
				IsInjecting = true;
			}	
		}else{
			double tstep = t - FPreviousTime;
			if (t < FTimeEnd){
				InjectedMass += FRate * tstep;
				rate = FRate;
			}
			else{
				if (!EndInjection){
					rate = (FuelInjected - InjectedMass) / tstep;
					EndInjection = true;
				}
			}

		}
	}
	FPreviousTime = t;
	FPreviousAngle = angle;
	FCurrentRate = rate;
	return rate;
}

void TInjector::setFuel(TComponentArray_ptr workfluid) {
	RowVector Y;
	Y.setZero(workfluid->size());
	Fuel = make_shared<TFluid>(workfluid);
	for (int i = 0; i < workfluid->size(); i++){
		if (workfluid->at(i)->getName() == "FUEL"){
			Y(i) = 1.;
		}
	}
	Fuel->SetComposition(Y);
}

shared_ptr<TInjector> createInjector(xml_node node, TComponentArray_ptr fluid){

	auto Injector = make_shared<TInjector>();

	Injector->FSOI = GetAttributeAsDouble(node, "SOI");
	Injector->FPeriod = GetAttributeAsDouble(node, "Period");
	Injector->FuelInjected = GetAttributeAsDouble(node, "FuelMass");
	Injector->FRate = Injector->FuelInjected / Injector->FPeriod;

	RowVector Y;

	Y.setZero(fluid->size());
	Injector->Fuel = make_shared<TFluid>(fluid);
	for (int i = 0; i < fluid->size(); i++){
		if (fluid->at(i)->getName() == "FUEL"){
			Y(i) = 1.;
		}
	}
	Injector->Fuel->SetComposition(Y);

	return Injector;
}

