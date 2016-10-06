/* --------------------------------------------------------------------------------*\
==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 \\ |  X  | //  W ave     |
 \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*-------------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------
#ifndef TIsoSpeedLineH
#define TIsoSpeedLineH

#include <vector>
#include <iostream>

// #include "interp_1d.h"
#include "Globales.h"
//#include "TTC_HTM.h"
#include "TTCHTM2.h"
#include "TEffTFunction.h"

struct stFunEffectSec {
	double Gamma;
	double Area;
	double R;
	bool CalculaGR;
	double Angle;
	double Diam1;
	double Diam2;
	double Diam3;
	double Diam4;
	double n_limit;
	double BladeH;
	double CritAngle;

	void fun(double& SES, double& RES, double MASS, double ER, double EFF, double Speed);

};

// ---------------------------------------------------------------------------

class TIsoSpeedLine {
  private:

	double FSpeed;
	double FERMax;
	double FERMin;
	dVector FExpansionRatio;
	dVector FExpansionRatioAdim;
	dVector FReducedAirMassFlow;
	dVector FReducedAirMassFlowAdim;
	dVector FEfficiency;
	dVector StatorEffectiveSection;
	dVector RotorEffectiveSection;
	dVector FPower;
	double FPowerMin;
	double FPowerMax;
	int FNumDatos;

	Base_interp *iStator;
	Base_interp *iRotor;
	Base_interp *iEfficiency;

	TEffTFunction* FEffTFun;

  public:

	TIsoSpeedLine();
	~TIsoSpeedLine();

	void ReadIsoSpeed(int point, FILE *Input);

	void AsignaValores(double ER, double MF, double EF);

	void EffectiveSection(double Area, bool CalculaGR, double Angle, double Diam1, double Diam2, double Diam3,
						  double n_limit);

	void CalculatePower(double Tin);

	void PutSpeed(double sp);

	double Speed() /* {return FSpeed;} */;

	double ERMax() /* {return FERMax;} */;

	double ERMin() /* {return FERMin;} */;

	double PowerMin() {
		return FPowerMin;
	}

	double PowerMax() {
		return FPowerMax;
	}

	double Stator(double ERAdim);

	double Rotor(double ERAdim);

	double Efficiency(double ERAdim);

	void Adimensionaliza();

//	void GetAdiabaticEfficiency(TTC_HTM *HTM, double TinT, double TinC, double MaxEff);

	void GetAdiabaticEfficiency(TTC_HTM2 *HTM, double TinT, double TinC, double MaxEff);

	void PrintEffectiveSection(FILE *fich);

	double MaximumEfficiency();

	int size() {
		return FNumDatos;
	}
	;

	double MapEff(int i) {
		return FEfficiency[i];
	}
	;

	double MapER(int i) {
		return FExpansionRatio[i];
	}
	;

	double MapMass(int i) {
		return FReducedAirMassFlow[i];
	}
	;

	double MapSpeed() {
		return FSpeed;
	}
	;

	void AsignEfficiencyData(double K1, double A0, double beta, double K3, double g, double Dwheel);

	void Fit_Efficiency(int step, TAeffFunction *AeffFun);

	double GetEffCoef(int i) {
		return FEffTFun->get_z(i);
	}

	void Fix_z3(double val);

	void Put_z(double z1, double z2, double z3);

};

#endif
