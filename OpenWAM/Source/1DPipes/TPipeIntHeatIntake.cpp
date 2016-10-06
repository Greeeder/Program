/**
 * @file TPipeIntHeatIntake.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 de feb. de 2016
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
 * Object to calculate internal heat transfer in intake pipes.
 */
#include "TPipeIntHeatIntake.h"

TPipeIntHeatIntake::TPipeIntHeatIntake(int ncells, double HeatMult, VectorXd D, double cellsize) : 
	TPipeIntHeat(HeatMult, D, cellsize) {
	Cond.setZero(ncells);
	Qint.setZero(ncells);
	h.setZero(ncells);
}

TPipeIntHeatIntake::~TPipeIntHeatIntake() {
	// TODO Auto-generated destructor stub
}

ArrayXd TPipeIntHeatIntake::Heat(const ArrayXd &Tgas, const ArrayXd &Twall, const ArrayXd &Re, const std::vector<TFluid_ptr>& Gas) {

	if(HeatMultiplier > 0) {
		for(int i = 0; i < Cond.size(); i++) {
			Cond(i) = Gas[i]->getConductivity();
		}

		h = 0.0694 * Re.pow(0.75) * Cond / Dint;

		Qint = Aint * h * (Tgas - Twall) * HeatMultiplier;
	}

	return Qint;
}

