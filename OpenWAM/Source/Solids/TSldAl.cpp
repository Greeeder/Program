/**
 * @file TSldAl.cpp
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
 * Class for the calculation of the Aluminium properties.
 */
#include "TSldAl.h"

TSldAl::TSldAl() : TSolid() {
	Name = "Aluminium";
	CoefCond.setZero(1);
	CoefCond(0) = 237;
	ExpCond.setLinSpaced(CoefCond.size(), 0, CoefCond.size() - 1);

	CoefHeCap.setZero(1);
	CoefHeCap(0) = 900;
	ExpHeCap.setLinSpaced(CoefHeCap.size(), 0, CoefHeCap.size() - 1);

	CoefDens.setZero(1);
	CoefDens(0) = 2700;
	ExpDens.setLinSpaced(CoefDens.size(), 0, CoefDens.size() - 1);
}

TSldAl::~TSldAl() {
	// TODO Auto-generated destructor stub
}

