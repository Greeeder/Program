/**
 * @file TFluid.cpp
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
 * The TFluid class is used for the different fluid to be used in OpenWAM.
 *
 * It represent a fluid composed by several chemical components and calculates
 * the thermodynamic properties depending on the composition.
 *
 */

//#include "stdafx.h"
#include "TFluid.h"
#include "Math_wam.h"

TFluid::TFluid() {

}

TFluid::TFluid(TFluid* fluid){
	Y = fluid->Y;

	_Component = (fluid->_Component);

	R = fluid->R;

	Viscosity = fluid->Viscosity;
	Conductivity = fluid->Conductivity;
	Cp = fluid->Cp;
}

TFluid::TFluid(TComponentArray_obj *componentes) {

	_Component = make_shared<TComponentArray_obj>(*componentes);

	Y.setZero(componentes->size());

}

TFluid::TFluid(TComponentArray_ptr componentes) {

	_Component = componentes;

	Y.setZero(componentes->size());

}

TFluid::TFluid(string name) {

}

TFluid::~TFluid() {
	// TODO Auto-generated destructor stub
}

void TFluid::AppendFluid(TFluid_ptr fluid1, double m1, double m0) {

	Y = (Y * m0 + fluid1->GetComposition() * m1) / (m0 + m1);
}

double TFluid::FunCp(double T) {

	Cp = Y(0) * _Component->at(0)->FunCp(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			Cp += Y(i) * _Component->at(i)->FunCp(T);
			_Y += Y(i);
		}
	}
	Cp /= _Y;
	return Cp;
}

double TFluid::FunCv(double T) {

	double Cv = Y(0) * _Component->at(0)->FunCv(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			Cv += Y(i) * _Component->at(i)->FunCv(T);
			_Y += Y(i);
		}
	}
	return Cv / _Y;
}

double TFluid::FunR() {

	R = Y(0) * _Component->at(0)->FunR();
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			R += Y(i) * _Component->at(i)->FunR();
			_Y += Y(i);
		}
	}
	return R / _Y;;
}

double TFluid::FunRHO(double T) {

	double RHO = Y(0) * _Component->at(0)->FunRHO(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			RHO += Y(i) * _Component->at(i)->FunRHO(T);
			_Y += Y(i);
		}
	}
	return RHO / _Y;
}

double TFluid::FunU(double T) {

	double U = Y(0) * _Component->at(0)->FunU(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			U += Y(i) * _Component->at(i)->FunU(T);
			_Y += Y(i);
		}
	}
	return U / _Y;
}

double TFluid::FunUsens(double T){

	double U = Y(0) * _Component->at(0)->FunCv(T) * T;
	double _Y = Y(0);

	for (unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			U += Y(i) * _Component->at(i)->FunCv(T) * T;
			_Y += Y(i);
		}
	}
	return U / _Y;
}

double TFluid::FunH(double T) {

	double H = Y(0) * _Component->at(0)->FunH(T);

	for(unsigned int i = 1; i < Y.size(); i++) {
		H += Y(i) * _Component->at(i)->FunH(T);
	}
	return H;
}

double TFluid::FunGamma(double T) {
	//gamma = 1 + r / Cv;
	return 1 + FunR() / FunCv(T);
}

TChemicalComponent_ptr TFluid::getChemicalComponent(std::string nombre) {

	for(int i = 0; i < _Component->size(); i++) {
		if(_Component->at(i)->getName() == nombre) {
			return _Component->at(i);
		}
	}

	cout << "ERROR: The chemical component " << nombre << " does not exist in the fluid" << endl;
	return NULL;
}

TComponentArray_ptr TFluid::GetComponents() const {
	return _Component;
}

double TFluid::getMolarFraction(std::string nombre){

	double top = 0;
	double bottom = 0;

	for (int i = 0; i < _Component->size(); i++) {
		bottom += Y(i) / _Component->at(i)->FunPM();
		if (_Component->at(i)->getName() == nombre) {
			top = Y(i) / _Component->at(i)->FunPM();
		}
	}
	
	return top / bottom;
}

double TFluid::getMolecularWeight(){

	double bottom = 0;

	for (int i = 0; i < _Component->size(); i++) {
		bottom += Y(i) / _Component->at(i)->FunPM();
	}

	return 1 / bottom;
}

double TFluid::GetY(int i) {
	return Y(i);
}

double TFluid::GetY(std::string nombre) {

	for(int i = 0; i < _Component->size(); i++) {
		if(_Component->at(i)->getName() == nombre) {
			return Y(i);
		}
	}

	cout << "ERROR: The chemical component " << nombre << " does not exist in the fluid" << endl;
	return 0.;
}

TChemicalComponent_ptr TFluid::GetComponent(int i) {

	return  _Component->at(i);
}

void TFluid::SetY(std::string nombre, double valor) {

	for(int i = 0; i < _Component->size(); i++) {
		if(_Component->at(i)->getName() == nombre) {
			Y(i) = valor;
			return;
		}
	}
	cout << "ERROR: The chemical component " << nombre << " cannot be set the mass fraction" << endl;
}

double TFluid::FunRComposicion(std::vector<double>* comp) {

	double result = comp->at(0) * _Component->at(0)->FunR();

	for(int i = 1; i < _Component->size(); i++) {
		result += comp->at(i) * _Component->at(i)->FunR();
	}

	return result;
}

double TFluid::FunUComposicion(std::vector<double>* comp, double T) {

	double result = comp->at(0) * _Component->at(0)->FunU(T);

	for(int i = 1; i < _Component->size(); i++) {
		result += comp->at(i) * _Component->at(i)->FunU(T);
	}

	return result;
}

double TFluid::FunCvComposicion(std::vector<double>* comp, double T) {
	double result = comp->at(0) * _Component->at(0)->FunCv(T);

	for(int i = 1; i < _Component->size(); i++) {
		result += comp->at(i) * _Component->at(i)->FunCv(T);
	}

	return result;
}

double TFluid::FunGammaComposicion(std::vector<double>* comp, double T) {
	//gamma = 1 + r / Cv;
	return 1 + FunRComposicion(comp) / FunCvComposicion(comp, T);
}

void TFluid::SetComponentes(TComponentArray_obj *componentes) {

	_Component = make_shared<TComponentArray_obj>(*componentes);

	Y.resize(componentes->size());
}

void TFluid::SetComponentes(TComponentArray_ptr componentes) {

	_Component = componentes;

	Y.resize(componentes->size());
}

double TFluid::getTemperature(double u) {
	double T = u / 717;

	stTforRealGas TForRealGas(this, u);

	return T = zbrent(TForRealGas, 200., 2000., 1e-5);
}

double TFluid::getTfromEqOfState(double rho, double p) {
	return p / rho / R;
}

double TFluid::getPfromEqOfState(double rho, double T) {
	return rho * R * T;
}

double TFluid::getRhofromEqOfState(double p, double T) {
	return p / R / T;
}

void TFluid::SetComposition(RowVector Ynew) {
	if (Ynew.size() == Y.size()){
		Y = Ynew;
		R = FunR();
	}
	else{
		cout << "ERROR: The size of the composition array is not coherent with the fluid" << endl;
	}

}

double TFluid::Funk(double T) {

	Conductivity = Y(0) * _Component->at(0)->Funk(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			Conductivity += Y(i) * _Component->at(i)->Funk(T);
			_Y += Y(0);
		}
	}
	return Conductivity / _Y;

}

double TFluid::FunVisc(double T) {

	Viscosity = Y(0) * _Component->at(0)->FunVisc(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			Viscosity += Y(i) * _Component->at(i)->FunVisc(T);
			_Y += Y(i);
		}
	}
	return Viscosity / _Y;

}

double TFluid::FunPr(double T) {

	double Pr = Y(0) * _Component->at(0)->FunPr(T);
	double _Y = Y(0);

	for(unsigned int i = 1; i < Y.size(); i++) {
		if (Y(i) > 0.01){
			Pr += Y(i) * _Component->at(i)->FunPr(T);
			_Y += Y(i);
		}
	}
	return Pr / _Y;

}

