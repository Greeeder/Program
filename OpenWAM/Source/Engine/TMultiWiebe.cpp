/**
 * @file TMultiWiebe.cpp
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
 * Calculates heat released by means of a multi-Wiebe function.
 */
#include "TMultiWiebe.h"

TMultiWiebe::TMultiWiebe() {
	
	Mode = "MultiWiebe";

}

TMultiWiebe::TMultiWiebe(TMultiWiebe * obj){
	Wiebe.resize(obj->Wiebe.size());
	for (int i = 0; i < obj->Wiebe.size(); i++){
		Wiebe[i] = new TWiebe(obj->Wiebe[i]);
	}
	Beta = obj->Beta;
	HRL = obj->HRL;
	Mode = "MultiWiebe";
}

TMultiWiebe::~TMultiWiebe() {
	// TODO Auto-generated destructor stub
}

double TMultiWiebe::getHRL(double angle){

	HRL = Beta[0] * Wiebe[0]->getHRL(angle);
	for(int i = 1; i < Wiebe.size(); i++){
		HRL += Beta[i] * Wiebe[i]->getHRL(angle);
	}
	return HRL;
}

void TMultiWiebe::ReadCombustionData(xml_node node_comb){

	int nlaw = CountNodes(node_comb, "SingleLaw");
	Wiebe.resize(nlaw);
	Beta.resize(nlaw);

	int i = 0;
	for(xml_node node_sl = GetNodeChild(node_comb,"SingleLaw"); node_sl;
		node_sl = node_sl.next_sibling("SingleLaw")){
		Beta[i] = GetAttributeAsDouble(node_sl, "Beta");
		xml_node node_wiebe = GetNodeChild(node_sl, "Wiebe");
		Wiebe[i] = new TWiebe();
		Wiebe[i]->ReadCombustionData(node_wiebe);
		i++;
	}
}

