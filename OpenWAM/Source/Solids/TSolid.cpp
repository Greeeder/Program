/**
 * @file TSolid.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 9 de mar. de 2016
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
 * Class for the calculation of any generic solid material.
 */
#include "TSolid.h"

TSolid::TSolid() {

}

TSolid::TSolid(const TSolid* obj) {

	Conductivity = obj->Conductivity;
	HeatCapacity = obj->HeatCapacity;
	Density = obj->Density;
	Elasticity = obj->Elasticity;
	CoefCond = obj->CoefCond;
	ExpCond = obj->ExpCond;
	CoefHeCap = obj->CoefHeCap;
	ExpHeCap = obj->ExpHeCap;
	CoefDens = obj->CoefDens;
	ExpDens = obj->ExpDens;
	CoefElas = obj->CoefElas;
	ExpElas = obj->ExpElas;
	Name = obj->Name;
}

TSolid::TSolid(string name) {
	
	Name = name;

}

TSolid::~TSolid() {
	// TODO Auto-generated destructor stub
}

double TSolid::getConductivity(double T) {
	return CoefCond.matrix().transpose() * pow(T, ExpCond).matrix();
}

double TSolid::getDensity(double T) {
	return CoefDens.matrix().transpose() * pow(T, ExpDens).matrix();
}

double TSolid::getElasticity(double T) {
	return CoefElas.matrix().transpose() * pow(T, ExpElas).matrix();
}

double TSolid::getHeatCapacity(double T) {
	return CoefHeCap.matrix().transpose() * pow(T, ExpHeCap).matrix();
}

void TSolid::setTemperature(double T){
	Conductivity = getConductivity(T);
	Density = getDensity(T);
	HeatCapacity = getHeatCapacity(T);
}
