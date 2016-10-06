/**
 * @file TFluidIdealAir.cpp
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
 * The TFluidIdealAir class is used for a fluid composed only by ideal air.
 *
 * It represents a fluid composed by ideal air. This file includes the
 * definition of such class.
 *
 */

//#include "stdafx.h"
#include "TFluidIdealAir.h"
#include "Math_wam.h"

TFluidIdealAir::TFluidIdealAir() {

	Y.resize(1, 1);
	comp.resize(1);
	comp[0] = make_shared<TChAir>("AIR", __R::Universal, __PM::IdealAir);
	_Component = make_shared<TComponentArray_obj>(comp);
}

TFluidIdealAir::~TFluidIdealAir() {
	// TODO Auto-generated destructor stub
}

double TFluidIdealAir::FunCp(double T) {
	Cp = _Component->at(0)->FunCp(T);;
	return Cp;

}

double TFluidIdealAir::FunCv(double T) {
	return _Component->at(0)->FunCv(T);
}

double TFluidIdealAir::FunR() {
	return _Component->at(0)->FunR();
}

double TFluidIdealAir::FunRHO(double T) {
	return _Component->at(0)->FunRHO(T);
}

double TFluidIdealAir::FunU(double T) {
	return _Component->at(0)->FunU(T);
}

double TFluidIdealAir::FunH(double T) {
	return _Component->at(0)->FunH(T);
}

double TFluidIdealAir::Funk(double T) {
	Conductivity = _Component->at(0)->Funk(T);
	return Conductivity;
}

double TFluidIdealAir::FunVisc(double T) {
	Viscosity = _Component->at(0)->FunVisc(T);
	return Viscosity;

}

double TFluidIdealAir::FunPr(double T) {
	return _Component->at(0)->FunPr(T);
}
