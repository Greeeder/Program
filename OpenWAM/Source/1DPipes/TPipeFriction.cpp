/**
 * @file TPipeFriction.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 feb. 2016
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
 * This class calculates the friction between the gas and the wall..
 */

#include "TPipeFriction.h"
#include "TTubo.h"

TPipeFriction::TPipeFriction() {
	// TODO Auto-generated constructor stub

}

TPipeFriction::~TPipeFriction() {
	// TODO Auto-generated destructor stub
}

void TPipeFriction::ReadFrictionData(xml_node node_pipe) {

	xml_node node_fr = GetNodeChild(node_pipe, "Pip_Friction");
	FrictionMultiplier = GetAttributeAsDouble(node_fr, "Multiplier");
	Rugosity = GetAttributeAsDouble(node_fr, "Rugosity");

}

ArrayXd TPipeFriction::FrictionCoefficient(const ArrayXd & Re) {

	return ArrayXd::Zero(Re.size());

}

