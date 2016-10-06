/**
 * @file TPipeIntHeat.cpp
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
 * This class calculates the heat transfer between the gas and the wall..
 */

#include "TPipeIntHeat.h"

TPipeIntHeat::TPipeIntHeat(double HeatMult, VectorXd D, double cellsize) {

	HeatMultiplier = HeatMult;
	Dint = D.array();
	Aint = D.array().pow(2) * __cons::Pi_4;

}

TPipeIntHeat::~TPipeIntHeat() {
	// TODO Auto-generated destructor stub
}

