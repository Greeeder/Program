/**
 * @file TWiebe.cpp
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
 * Calculates the heat released by means of a Wiebe function.
 */
#include "TWiebe.h"

TWiebe::TWiebe() : TCombustion() {
	
	Mode = "Wiebe";

}

TWiebe::TWiebe(TWiebe* obj) : TCombustion(obj) {

	Duration = obj->Duration;
	SOC = obj->SOC;
	m = obj->m;
	C = obj->C;
	HRL = obj->HRL;
	Mode = "Wiebe";

}

TWiebe::~TWiebe() {
	// TODO Auto-generated destructor stub
}

double TWiebe::getHRL(double angle){
	if(angle <= SOC) {
		HRL = 0.;
	} else {
		double xxx = C * pow((angle - SOC) / Duration, m + 1);
		if((xxx) < 20.) {
			HRL = 1. - 1. / exp(xxx);
		} else {
			HRL = 1.;
		}
	}
	return HRL;
}

void TWiebe::ReadCombustionData(xml_node node_comb){

	Duration = GetAttributeAsDouble(node_comb, "Duration");
	SOC = GetAttributeAsDouble(node_comb, "SOC");
	m = GetAttributeAsDouble(node_comb, "m");
	C = GetAttributeAsDouble(node_comb, "C");

}
