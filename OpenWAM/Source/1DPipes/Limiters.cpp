/* --------------------------------------------------------------------------------*\
==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 \\ |  X  | //  W ave     |
 \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*--------------------------------------------------------------------------------*/

/**
 * @file Limiters.hpp
 * @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
 * The functions in this file represent limiters that can be used
 * in MUSCL-like schemes.
 *
 * This file contains the implementation of such functions.
 */


#include "Limiters.hpp"

void MC(const RowArray &ratio, RowArray * phi) {
	for(int j = 0; j < ratio.cols(); j++) {
		for(int i = 0; i < ratio.rows(); i++) {
			(*phi)(i, j) = max(0.,
				min(min(2. * ratio(i, j), 0.5 * (1. + ratio(i, j))),
					2.));
		}
	}
}

void Minmod(const RowArray &ratio, RowArray * phi) {
	for(int j = 0; j < ratio.cols(); j++) {
		for(int i = 0; i < ratio.rows(); i++) {
			(*phi)(i, j) = max(0., min(1., ratio(i, j)));
		}
	}
}

void Superbee(const RowArray &ratio, RowArray * phi) {
	double x;
	for(int j = 0; j < ratio.cols(); j++) {
		for(int i = 0; i < ratio.rows(); i++) {
			x = max(min(2. * ratio(i, j), 1.), min(ratio(i, j), 2.));
			(*phi)(i, j) = max(0., x);
		}
	}
}

void VanLeer(const RowArray &ratio, RowArray * phi) {
	*phi = (ratio + ratio.abs()) / (1. + ratio.abs());
	for(int i = 0; i < ratio.rows() - 1; i++) {
		for(int j = 0; j < ratio.cols() - 1; j++) {
			if((*phi)(i, j) < 0.) {
				(*phi)(i, j) = 0.;
			}
		}
	}
}