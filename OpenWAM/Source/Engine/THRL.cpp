/**
 * @file THRL.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 14 de mar. de 2016
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
 * Impose heat realeas by means of a 1D look-up table (Angle-HRL).
 */
#include "THRL.h"

THRL::THRL() {
	
	Mode = "HRL";

}

THRL::THRL(THRL* obj) {
	
	Angle = obj->Angle;
	HeatRL = obj->HeatRL;
	SOC = obj->SOC;
	interpHRL = new Hermite_interp(Angle, HeatRL);
	Mode = "HRL";

}

THRL::~THRL() {
	// TODO Auto-generated destructor stub
}

double THRL::getHRL(double angle){

	HRL = interpHRL->interp(angle - SOC);
	return HRL;
}

void THRL::ReadCombustionData(xml_node node_comb){

	SOC = GetAttributeAsDouble(node_comb, "SOC");

	double Delta;
	bool getDelta = true;

	for (xml_node node_p = GetNodeChild(node_comb, "HRLrow"); node_p;
		node_p = node_p.next_sibling("HRLrow")){

		if (getDelta){
			Delta = GetAttributeAsDouble(node_p, "Angle");
			getDelta = false;
		}

		Angle.push_back(GetAttributeAsDouble(node_p, "Angle") - Delta);
		HeatRL.push_back(GetAttributeAsDouble(node_p, "HRLvalue"));
	}

	interpHRL = new Hermite_interp(Angle, HeatRL);
}

