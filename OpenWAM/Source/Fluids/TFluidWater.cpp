/**
 * @file TFluidWater.cpp
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
 * The TFluidWater class is used for the different fluid to be used in OpenWAM.
 *
 * It represent a fluid composed by water.
 *
 */

//#include "stdafx.h"
#include "TFluidWater.h"
#include "Math_wam.h"

TFluidWater::TFluidWater() {

	Y.resize(1, 1);
	comp.resize(1);
	comp[0] = make_shared<TChH2Ol>("WATER", __R::Universal, __PM::H2O);
	_Component = make_shared<TComponentArray_obj>(comp);

}

TFluidWater::~TFluidWater() {
	// TODO Auto-generated destructor stub
}

double TFluidWater::FunCp(double T) {

	Cp = _Component->at(0)->FunCp(T);
	return Cp;

}

double TFluidWater::FunCv(double T) {

	return _Component->at(0)->FunCv(T);
}

double TFluidWater::FunR() {

	return _Component->at(0)->FunR();
}

double TFluidWater::FunRHO(double T) {

	return _Component->at(0)->FunRHO(T);
}

double TFluidWater::FunU(double T) {

	return _Component->at(0)->FunU(T);
}

double TFluidWater::FunH(double T) {

	return _Component->at(0)->FunH(T);
}

double TFluidWater::Funk(double T) {

	Conductivity = _Component->at(0)->Funk(T);
	return Conductivity;

}

double TFluidWater::FunVisc(double T) {

	Viscosity = _Component->at(0)->FunVisc(T);
	return Viscosity;

}

double TFluidWater::FunPr(double T) {

	return _Component->at(0)->FunPr(T);

}

