/**
* @file HeatTransfer.h
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
* This file solve heat transfer in the cylinder.
*/

#include <string>
#include <vector>
//#include "Comun.h"
#include "TCrankMechanism.h"
#include "BasicHeatTransfer.h"

//namespace CalculoCALMEC
//{
#ifndef __HEATTRANSFER_H
#define __HEATTRANSFER_H



class HeatTransfer : public BasicHeatTransfer
{
public:

	struct stInicializar
	{
		double ConfMotor_BASE_D;
		double ConfMotor_BASE_S;

		double ConfMotor_BASE_VD;
		double ConfMotor_BASE_WR1A;
		double ConfMotor_BASE_WR1B;
		double ConfMotor_BASE_W2; //Constante_TRANS_CALOR_C_W2 o si hay radiación ConfMotor_BASE_W2_CON_RAD;

		double ConfMotor_PISTON_DB;
		double ConfMotor_PISTON_AP;

		double ConfMotor_CULATA_ACUL;
		double ConfMotor_CULATA_CTM;

		double Constante_TRANS_CALOR_C_A_WOSCHNI;
		double Constante_TRANS_CALOR_C_B_WOSCHNI;
		double Constante_TRANS_CALOR_C_C_WOSCHNI;
		double Constante_TRANS_CALOR_C_N;
		double Constante_TRANS_CALOR_C_WH1;

		double VarMed_N;

		double DatosCalcGeneral_CM;

		double Meanflow_adm; // Gasto medio efectivo de aire [kg/cc] - nret = rendimiento de retencion utilizado para woschni Tumble
		double cu_prima;
		double k_nodis;
		double ratio_ctms;
	};

	struct stInicializarAlRCA
	{
		double TPIS;
		double TCUL;
		double TCIL;
		double PCA;
		double VCA;
		double TCA;
		double CW1;
		double CW2;

	};

	struct stInicializarPorAngulo
	{
		double deltat;
		double ang;
		double areacil;
		bool ciclo_cerrado;
		double VCIL;
	};

	void Inicializar(TCrankMechanism *CrankMechanism, stInicializar Inicializar);

	void InicializarAlRCA(stInicializarAlRCA InicializarAlRCA);

	void InicializarPorAngulo(stInicializarPorAngulo InicializarPorAngulo);

	virtual void Calcula_Calor_Woschni(double T, double p);

	virtual double Heat(double T, double p);

	virtual double ConvectiveHeat() { return _QW; };

	virtual double RadiativeHeat() {
		return _QRAD;
	};

	virtual void ReadData(xml_node node);

	double getQW() const { return _QW; }
	double getQRAD() const { return _QRAD; }
	void setQRAD(double value) { _QRAD = value; }
	double geth() const { return _h; }
	double getCU() const { return _CU; }
	double getQCUL() const { return _QCUL; }
	double getQPIS() const { return _QPIS; }
	double getQCIL() const { return _QCIL; }
	double getQPIPA_ADM() const { return _QPIPA_ADM; }
	double getQPIPA_ESC() const { return _QPIPA_ESC; }
	double Xp(double x, double ratio_ctms);

protected:

	stInicializar _inicializar;
	stInicializarAlRCA _inicializarAlRCA;
	stInicializarPorAngulo _inicializarPorAngulo;
	TCrankMechanism *crankMechanism;

	//resultados
	double _QW;
	double _QRAD;
	double _h; //coeficiente de película
	double _CU;
	double _QCUL;
	double _QPIS;
	double _QCIL;
	double _QPIPA_ADM;
	double _QPIPA_ESC;


};
#endif