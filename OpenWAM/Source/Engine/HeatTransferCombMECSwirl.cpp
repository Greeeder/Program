/**
* @file HeatTransferCombMECSwirl.cpp
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @date 30 de mar. de 2016
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
* This file solve heat transfer in CIE with swirl.
*/

#include <iostream>
#include <fstream>
#include "HeatTransferCombMECSwirl.h"
#include <sstream>
#include "Constantes.h"
#include "Math.h"

using namespace std;
//namespace CalculoCALMEC
//{


void HeatTransferCombMECSwirl::Calcula_Calor_Woschni(double T, double p)
{
	double x;
	double A;
	double b;
	double c;
	double factor;
	double d_piston;
	double radio;
	double parr;
	double COMB;
	double c1;
	double pot1;
	double pot2;
	double pot3;
	double pot4;
	double vg;

	try
	{
		if (_inicializarPorAngulo.ciclo_cerrado) {

			parr = _inicializarAlRCA.PCA * pow((_inicializarAlRCA.VCA / _inicializarPorAngulo.VCIL), 1.36);
			COMB = ((_inicializar.ConfMotor_BASE_VD * _inicializarAlRCA.TCA) / (_inicializarAlRCA.PCA * _inicializarAlRCA.VCA)) * (p - parr);
			if (COMB < 0)
				COMB = 0;

			x = Xp(_inicializarPorAngulo.ang, _inicializar.ratio_ctms);
			_CU = _inicializar.k_nodis * _inicializar.cu_prima * x;

			A = _inicializar.Constante_TRANS_CALOR_C_A_WOSCHNI;
			b = _inicializar.Constante_TRANS_CALOR_C_B_WOSCHNI;
			c = _inicializar.Constante_TRANS_CALOR_C_C_WOSCHNI;

			factor = A + b * pow(_inicializar.DatosCalcGeneral_CM, c);

			c1 = _inicializarAlRCA.CW1 + _inicializarAlRCA.CW2 * _CU / _inicializar.DatosCalcGeneral_CM;
		}
		else {

			COMB = 0;
			d_piston = crankMechanism->getVolumeDisplaced();
			//d_piston = Dist_Despl_inst(_inicializarPorAngulo.ang, _inicializar.ConfMotor_BASE_LM, _inicializar.ConfMotor_BASE_LB, _inicializar.ConfMotor_BASE_E);
			radio = sqrt(pow(_inicializar.ConfMotor_PISTON_DB, 2) * (1 - d_piston / _inicializar.ConfMotor_BASE_S) + ((d_piston / _inicializar.ConfMotor_BASE_S)) * pow(_inicializar.ConfMotor_BASE_D, 2)) / 2;
			_CU = _inicializar.ConfMotor_CULATA_CTM * (2 * __cons::Pi * _inicializar.VarMed_N) * radio * 0.75 * pow((_inicializar.ConfMotor_BASE_D / (2 * radio)), 2);

			factor = 1;

			c1 = _inicializar.ConfMotor_BASE_WR1A + _inicializar.ConfMotor_BASE_WR1B * _CU / _inicializar.DatosCalcGeneral_CM;
		}

		pot1 = pow(_inicializar.ConfMotor_BASE_D, (_inicializar.Constante_TRANS_CALOR_C_N - 1));
		pot3 = pow(p, _inicializar.Constante_TRANS_CALOR_C_N);
		vg = (c1 * _inicializar.DatosCalcGeneral_CM + _inicializar.ConfMotor_BASE_W2 * COMB);
		pot4 = pow(vg, _inicializar.Constante_TRANS_CALOR_C_N);
		pot2 = pow(T, (0.75 - 1.62 * _inicializar.Constante_TRANS_CALOR_C_N));

		_h = _inicializar.Constante_TRANS_CALOR_C_WH1 * pot1 * pot2 * pot3 * pot4;

		_QCUL = _inicializarPorAngulo.deltat * factor * _h * _inicializar.ConfMotor_CULATA_ACUL * (T - _inicializarAlRCA.TCUL);
		_QPIS = _inicializarPorAngulo.deltat * factor * _h * _inicializar.ConfMotor_PISTON_AP * (T - _inicializarAlRCA.TPIS);
		_QCIL = _inicializarPorAngulo.deltat * factor * _h * _inicializarPorAngulo.areacil * (T - _inicializarAlRCA.TCIL);
		_QW = _QCUL + _QPIS + _QCIL;
		_QRAD = 0;

	}
	catch (const std::exception& e)
	{
		throw e;
	}
}
