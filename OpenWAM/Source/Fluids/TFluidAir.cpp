/**
 * @file TFluidAir.cpp
 * @author F.J. Arnau <farnau@mot.upv.es>
 * @version 0.1.0
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
 * The TFluidAir class is used for the different fluid to be used in OpenWAM.
 *
 * It represent a fluid composed by air.
 *
 */

//#include "stdafx.h"
#include "TFluidAir.h"
#include "Math_wam.h"

TFluidAir::TFluidAir() {

	Y.resize(1, 1);
	comp.resize(1);
	comp[0] = make_shared<TChAir>("AIR", __R::Universal, __PM::Air);
	_Component = make_shared<TComponentArray_obj>(comp);
}

TFluidAir::~TFluidAir() {
	// TODO Auto-generated destructor stub
}

double TFluidAir::FunCp(double T) {

	Cp = _Component->at(0)->FunCp(T);;
	return Cp;

}

double TFluidAir::FunCv(double T) {

	return _Component->at(0)->FunCv(T);
}

double TFluidAir::FunR() {

	return _Component->at(0)->FunR();
}

double TFluidAir::FunRHO(double T) {

	return _Component->at(0)->FunRHO(T);
}

double TFluidAir::FunU(double T) {

	return _Component->at(0)->FunU(T);
}

double TFluidAir::FunH(double T) {

	return _Component->at(0)->FunH(T);
}

double TFluidAir::Funk(double T) {

	Conductivity = _Component->at(0)->Funk(T);
	return Conductivity;

}

double TFluidAir::FunVisc(double T) {

	Viscosity = _Component->at(0)->FunVisc(T);
	return Viscosity;

}

double TFluidAir::FunPr(double T) {

	return _Component->at(0)->FunPr(T);

}

