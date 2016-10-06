/**
 * @file TFluidOil.cpp
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
 * The TFluidOil class is used for the different fluid to be used in OpenWAM.
 *
 * It represent a fluid composed by oil.
 *
 */

//#include "stdafx.h"
#include "TFluidOil.h"
#include "Math_wam.h"

TFluidOil::TFluidOil() {

	Y.resize(1, 1);
	comp.resize(1);
	comp[0] = make_shared<TChOil>("OIL", __R::Universal, 300.);
	_Component = make_shared<TComponentArray_obj>(comp);

}

TFluidOil::~TFluidOil() {
	// TODO Auto-generated destructor stub
}

double TFluidOil::FunCp(double T) {

	return _Component->at(0)->FunCp(T);;

}

double TFluidOil::FunCv(double T) {

	return _Component->at(0)->FunCv(T);
}

double TFluidOil::FunR() {

	return _Component->at(0)->FunR();
}

double TFluidOil::FunRHO(double T) {

	return _Component->at(0)->FunRHO(T);
}

double TFluidOil::FunU(double T) {

	return _Component->at(0)->FunU(T);
}

double TFluidOil::FunH(double T) {

	return _Component->at(0)->FunH(T);
}

double TFluidOil::Funk(double T) {

	return _Component->at(0)->Funk(T);

}

double TFluidOil::FunVisc(double T) {

	return _Component->at(0)->FunVisc(T);

}

double TFluidOil::FunPr(double T) {

	return _Component->at(0)->FunPr(T);

}

