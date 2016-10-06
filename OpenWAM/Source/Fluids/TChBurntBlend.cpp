/**
 * @file TChBurntBlend.cpp
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
 * The TChBurntGasoil class represents the chemical component burnt gas produced by a
 * combustion with a fuel blend.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of a burnt gas.
 *
 */

//#include "stdafx.h"
#include "TChBurntBlend.h"

TChBurntBlend::TChBurntBlend(std::string nombre, TComponentArray_obj Fuel, std::vector<double> Fraction,
							 std::vector<double> PMQ, double RU, double XO2) {

	std::vector<double> n;
	std::vector<double> m;
	std::vector<double> XF;
	std::vector<double> aux;
	std::vector<double> XQ;

	setName(nombre);

	_BurntGas.resize(Fuel.size());

	double PMF = 0.;
	double PMSUM = 0.;
	for(int i = 0; i < Fuel.size(); i++) {
		if(Fuel[i]->getName() == "GASOLINA") {
			_BurntGas[i] = make_shared<TChBurntGasolina>("BURNT_GASOLINA", RU, PMQ[i]);
		} else if(Fuel[i]->getName() == "GASOIL") {
			_BurntGas[i] = make_shared<TChBurntGasoil>("BURNT_GASOIL", RU, PMQ[i]);
		}
		n.push_back(Fuel[i]->FunPM() / (12 + Fuel[i]->FunRelacionHC() + 16 * Fuel[i]->FunRelacionOC()));
		m.push_back(n[i] * Fuel[i]->FunRelacionHC());
		PMF += Fraction[i] / Fuel[i]->FunPM();
	}
	PMF = 1 / PMF;
	double moles;

	for(int i = 0; i < Fuel.size(); i++) {
		XF.push_back(Fraction[i] * PMF / Fuel[i]->FunPM());
		aux.push_back(n[1] + m[1] / 2 + n[1] * (1 - XO2) / XO2 * (1 + Fuel[i]->FunRelacionHC() / 4 - Fuel[i]->FunRelacionOC() /
					  2));
		moles += XF[i] * aux[i];
	}

	_PM = 0.;
	for(int i = 0; i < Fuel.size(); i++) {
		XQ.push_back(XF[i] * aux[i] / moles);
		_PM += XQ[i] * PMQ[i];
	}

	_R = _RU / _PM;

	for(int i = 0; i < Fuel.size(); i++) {
		_Yq_BLEND.push_back(XQ[i] * PMQ[i] / _PM);
	}

}

TChBurntBlend::~TChBurntBlend() {
	// TODO Auto-generated destructor stub
}

double TChBurntBlend::FunCv(double T) {
	double Cv = 0;
	for(int i = 0; i < _BurntGas.size(); i++) {
		Cv += _BurntGas[i]->FunCv(T) * _Yq_BLEND[i];
	}
	return Cv;
}

double TChBurntBlend::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return  641.154 + T * (0.43045 + T * (-0.0001125 + T * 8.979e-9));

	return 0.;
}

double TChBurntBlend::FunU(double T) {
	double U = 0;
	for(int i = 0; i < _BurntGas.size(); i++) {
		U += _BurntGas[i]->FunU(T) * _Yq_BLEND[i];
	}
	return U;
}

double TChBurntBlend::FunH(double T) {
	return FunU(T) + FunR() * T;
}

