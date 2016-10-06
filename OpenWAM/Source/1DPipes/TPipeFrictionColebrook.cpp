/**
 * @file TPipeFrictionColebrook.cpp
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

#include "TPipeFrictionColebrook.h"

TPipeFrictionColebrook::TPipeFrictionColebrook() {
	// TODO Auto-generated constructor stub

}

TPipeFrictionColebrook::~TPipeFrictionColebrook() {
	// TODO Auto-generated destructor stub
}

ArrayXd TPipeFrictionColebrook::FrictionCoefficient(const ArrayXd& Re) {

	if (FrictionMultiplier > 0){
		f = ArrayXd::Constant(Re.size(), 32.);
		for (int i = 0; i < Re.size(); i++) {
			if (Re(i) > 2000) {
				f(i) = 0.0625 / pow2(log10(RelativeRugosity(i) + 5.74 / pow(Re(i), 0.9)));
			}
			else if (Re(i) > 1) {
				f(i) = 32. / Re(i);
			}
		}
		f *= FrictionMultiplier;
	}
	else{
		f = ArrayXd::Zero(Re.size());
	}

	return f;
}

void TPipeFrictionColebrook::setRelativeRugosity(const RowVector& D){
	RelativeRugosity = Rugosity / 3700 / D;
}

