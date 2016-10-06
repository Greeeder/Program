/**
 * @file TPortIntHeatExhaust.cpp
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
 * Object to calculate internal heat transfer in intake ports.
 */
#include "TPortIntHeatExhaust.h"

TPortIntHeatExhaust::TPortIntHeatExhaust(int ncells, double HeatMult, VectorXd D, double cellsize) : 
	TPipeIntHeat(HeatMult, D, cellsize) {
	Cond.setZero(ncells);
	Visc.setZero(ncells);
	ViscWall.setZero(ncells);
	Qint.setZero(ncells);
	h.setZero(ncells);

}

TPortIntHeatExhaust::~TPortIntHeatExhaust() {
	// TODO Auto-generated destructor stub
}

ArrayXd TPortIntHeatExhaust::Heat(const ArrayXd &Tgas, const ArrayXd &Twall, const ArrayXd &Re, const std::vector<TFluid_ptr>& Gas) {

	if (HeatMultiplier > 0) {
		for (int i = 0; i < Cond.size(); i++) {
			Cond(i) = Gas[i]->getConductivity();
			Visc(i) = Gas[i]->getViscosity();
			ViscWall(i) = Gas[i]->FunVisc(Twall(i));
		}


		h = 0.1 * Re.pow(0.8) * 0.709 * (Visc / ViscWall).pow(0.14) * Cond / Dint;

		Qint = Aint * h * (Tgas - Twall) * HeatMultiplier;
	}

	return Qint;
}

