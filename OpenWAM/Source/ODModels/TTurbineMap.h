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
#ifndef TTurbineMapH
#define TTurbineMapH

#include "TTurbPosition.h"

class TTurbina;

struct stMapData {
	dVector Position;
	dVector Speed;
	dVector ERatio;
	dVector Mass;
	dVector Efficiency;
	dVector BSR;
};

class TTurbineMap {
  private:

	double FDiamIN;
	double FDiamOUT;
	int FNumPositions;

	stMapData FMapData;

	std::vector<TTurbPosition> FTurbPosition;

	bool FFixedTurbine;

	bool FIsAdiabatic;
	double FTempMeasure;

	bool FExtrapolate;

	double FStatorES;
	double FRotorES;
	double FEffTurb;

	dVector FPowerMin;
	dVector FPowerMax;

	dVector FOpening;
	stFunEffectSec FFunEffectSec;
	dVector FMaxEff;

	Base_interp *iMaxEff;

	stBladeMechanism *BladeMechanism;
	TAeffFunction_v2 *AeffFun;
	TEffTFunction_v2 *EffFun;

	double FMaxOpening;
	double FMinOpening;

  public:
	TTurbineMap();

	~TTurbineMap();

	void LoadTurbineMap(FILE *Input, TTurbina *Turbine);

	void LoadTurbineMapXML(xml_node node_turb);

	void CurrentEffectiveSection(double n, double er, double rack, double T10T00);

	double StatorEF() {
		return FStatorES;
	}
	;

	double RotorEF() {
		return FRotorES;
	}
	;

	double EffTurb() {
		return FEffTurb;
	}
	;

	double getTempMeasure() {
		return FTempMeasure;
	}
	;

	bool IsFixed() {
		return FFixedTurbine;
	}
	;

	void PrintFinalMap(FILE *fich);

//	void CalculateAdiabaticEfficiency(TTC_HTM *HTM, double TinC);

	void CalculateAdiabaticEfficiency(TTC_HTM2 *HTM, double TinC);

	void ProcessTurbineMap(TTurbina *Turbine);

	double TempMeasure() {
		return FTempMeasure;
	}
	;

	void SearchMaximumEfficiency();

	double MaxEfficiency(double rack);

	void InterpExtendedMap(double n, double er, double bsr, double rack);

	void CorrectRotorEffArea(double T10T00);

	void InterpolateTurbineMap(double n, double er, double bsr, double rack, double &mred, double &eff);

	void PrintExtrapolatedMap(double Dw);

	void ExtrapoltePoint(double n, double er, double bsr, double rack, double& aeff, double& eff);

	void ExtrapolateTurbineMap(TTurbina *Turbina);

	bool IsExtrapolated() {
		return FExtrapolate;
	}
	;

};
// ---------------------------------------------------------------------------
#endif
