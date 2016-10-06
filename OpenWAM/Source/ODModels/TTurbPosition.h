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


 \*--------------------------------------------------------------------------------*/

//---------------------------------------------------------------------------
#ifndef TTurbPositionH
#define TTurbPositionH

#include "TIsoSpeedLine.h"
#include "TAeffFunction.h"
// ---------------------------------------------------------------------------

class TTurbPosition {
  private:

	double FPosition;
	double FAngle;

	double FStatorSec;
	double FRotorSec;

	double FEfficiency;

	double FSpeedMin;
	double FSpeedMax;

	double FPowerMin;
	double FPowerMax;

	int FLines;

	std::vector<TIsoSpeedLine> FSpeedLine;

	dVector FERMax;
	dVector h_FERMax;
	dVector FERMin;
	dVector h_FERMin;
	dVector FSpeed;

	dVector MapSpeed;
	dVector MapER;
	dVector MapEff;
	dVector MapMass;

	dVector MapAeff;
	dVector MapBSR;
	dVector MapK3;

	TAeffFunction *AeffFun;
	TEffTFunction *EffTFun;

	Hermite_interp *z1_interp;
	Hermite_interp *z2_interp;
	double z3_mean;

	dMatrix z_Eff;

  public:

	TTurbPosition();
	~TTurbPosition();

	void ReadTurbinPosition(FILE *Input, int rows, double pos, double ang);

	void LoadTurbinPosition(dVector SP, dVector ER, dVector MF, dVector EF, double pos);

	void ReadTurbinPositionXML(xml_node node_map);

	void EffectiveArea(double Area, bool CalculaGR, double Diam1, double Diam2, double Diam3, double n_limit);

	void CalculatePower(double Tin);

	void PutPosition(double Pos);

	void InterpolaPosicion(double n, double er);

	void SearchMapLimits();

	double StatorSec() /* {return FStatorSec;} */;

	double RotorSec() /* {return FRotorSec;} */;

	double Rack() /* {return FPosition;} */;

	double Efficiency() /* {return FEfficiency;} */;

	double getAngle() {
		return FAngle;
	}
	;

	void putAngle(double val) {
		FAngle = val;
	}
	;

	double getSpeedMax() {
		return FSpeedMax;
	}
	;

	double getSpeedMin() {
		return FSpeedMin;
	}
	;

	void PrintTurbinePosition(FILE *Fich);

	double MinPowerLimit(double rtc);

	double MaxPowerLimit(double rtc);

//	void AdiabaticEfficiency(TTC_HTM *HTM, double TinT, double TinC);

	void AdiabaticEfficiency(TTC_HTM2 *HTM, double TinT, double TinC);

	double MaximumEfficiency();

	void ExtrapolateTurbinePosition(FILE *fich, double Dwheel, double BladeH, double Dout, double Dnut, double Din,
		double beta);

	void MapInterpolation(double Speed, double BSR, double ER, double& Aeff_interp, double &Eff_interp);

	TIsoSpeedLine getSpeedLine(int i) const {
		return FSpeedLine[i];
	}

	unsigned int size() const {
		return FSpeedLine.size();
	}
};

#endif
