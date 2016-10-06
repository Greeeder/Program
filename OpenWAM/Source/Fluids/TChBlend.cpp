/**
 * @file TChBlend.cpp
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
 * The TChGasolina class represents the chemical component fuel blend.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of a fuel blend.
 *
 */

//#include "stdafx.h"
#include "TChBlend.h"

TChBlend::TChBlend(std::string nombre, TComponentArray_obj Fuel, std::vector<double> Fraction,
				   double RU) {

	std::vector<double> n;
	std::vector<double> m;
	std::vector<double> p;

	_Fuel = Fuel;
	_YFuel = Fraction;
	setName(nombre);

	if(_Fuel.size() != _YFuel.size()) {
		std::cerr << "ERROR: The size of the components for the blend fuel is different to the size of the fractions vector" <<
				  std::endl;
	}

	double totalcomp = 0.;
	for(int i = 0; i < _YFuel.size(); i++) {
		totalcomp += _YFuel[i];
	}
	if(totalcomp > 1 + 1e-6 || totalcomp < 1 - 1e-6) {
		std::cerr << "ERROR: The total fraction of blend fluel is different to 1" << std::endl;
	}

	for(int i = 0; i < _Fuel.size(); i++) {
		n.push_back(_Fuel[i]->FunPM() / (12 + _Fuel[i]->FunRelacionHC() + 16 * _Fuel[i]->FunRelacionOC()));
		m.push_back(n[i] * _Fuel[i]->FunRelacionHC());
		p.push_back(n[i] * _Fuel[i]->FunRelacionOC());
	}
	//TO-DO Falta completar el cálculo de los diferentes parámetros.

	_PM = 0.;
	_HV = 0.;
	_DENSF = 0.;
	_CAHFL = 0.;
	_FE = 0.;
	double n_blend = 0;
	double m_blend = 0;
	double p_blend = 0;
	for(int i = 0; i < _Fuel.size(); i++) {
		_PM += _YFuel[i] / _Fuel[i]->FunPM();
		_DENSF += _YFuel[i] / _Fuel[i]->FunRHO();
		n_blend += _YFuel[i] / _Fuel[i]->FunPM() * n[i];
		m_blend += _YFuel[i] / _Fuel[i]->FunPM() * m[i];
		p_blend += _YFuel[i] / _Fuel[i]->FunPM() * p[i];
		_HV += _YFuel[i] * _Fuel[i]->FunHV();
		_CAHFL += _YFuel[i] * _Fuel[i]->Fun_CAHFL();
		_CBHFL += _YFuel[i] * _Fuel[i]->Fun_CBHFL();
		_FE += _YFuel[i] / _Fuel[i]->getFE();

	}
	_PM = 1 / _PM;
	_DENSF = 1 / _DENSF;
	_FE = 1 / _FE;
	n_blend = n_blend * _PM;
	m_blend = m_blend * _PM;
	p_blend = p_blend * _PM;
	_HC = m_blend / n_blend;
	_OC = p_blend / n_blend;
	_R = _RU / _PM;

}

// TChBlend::~TChBlend() {
// 	// TODO Auto-generated destructor stub
// }

double TChBlend::FunCv(double T) {
	double Cv = 0;
	for(int i = 0; i < _Fuel.size(); i++) {
		Cv += _Fuel[i]->FunCv(T) * _YFuel[i];
	}
	return Cv;
}

double TChBlend::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return  641.154 + T * (0.43045 + T * (-0.0001125 + T * 8.979e-9));

	return 0.;
}

double TChBlend::FunU(double T) {
	double U = 0;
	for(int i = 0; i < _Fuel.size(); i++) {
		U += _Fuel[i]->FunU(T) * _YFuel[i];
	}
	return U;
}

double TChBlend::FunH(double T) {
	//return _CAHFL + _CBHFL * (_TF + _C_TIY);

	return _CAHFL + _CBHFL * T;
}

