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
#pragma hdrstop

#include "TTurbineMap.h"
#include "TTurbina.h"
//#include <cmath>

TTurbineMap::TTurbineMap() {
	FIsAdiabatic = true;
	FFixedTurbine = false;
	FExtrapolate = false;

	iMaxEff = NULL;
	FEffTurb = 0.6;

}

TTurbineMap::~TTurbineMap() {

	FTurbPosition.clear();

	delete iMaxEff;
}

//void TTurbineMap::LoadTurbineMap(FILE *Input, TTurbina *Turbine) {
//	int rows = 0, Adiab = 0;
//	double pos = 0., ang = 0.;
//	double Area = pow2(Turbine->getDiametroTurbinaIn()) * Pi / 4;
//	double n_limit = 1.165;
//	bool CalculaGR = false;
//
//	fscanf(Input, "%d ", &Adiab);
//	if (Adiab == 0) {
//		FIsAdiabatic = false;
//		fscanf(Input, "%lf ", &FTempMeasure);
//	}
//	int Extrapolate = 0;
//	fscanf(Input, "%d ", &Extrapolate);
//	if (Extrapolate == 1) {
//		FExtrapolate = true;
//	}
//
//	fscanf(Input, "%d ", &FNumPositions);
//	FTurbPosition.resize(FNumPositions);
//	for (int i = 0; i < FNumPositions; i++) {
//		fscanf(Input, "%d %lf %lf", &rows, &pos, &ang);
//		FTurbPosition[i].ReadTurbinPosition(Input, rows, pos, ang);
//		if (ang > Turbine->getCriticalAngle())
//			CalculaGR = false;
//		else
//			CalculaGR = true;
//
//		if (!FExtrapolate) {
//			FTurbPosition[i].EffectiveArea(Area, CalculaGR,
//					Turbine->getDiametroRodete(),
//					Turbine->getDiametroRodeteOut(),
//					Turbine->getDiametroTuerca(), n_limit);
//
//			FTurbPosition[i].SearchMapLimits();
//		} else {
//			FFunEffectSec.Diam1 = Turbine->getDiametroRodete();
//			FFunEffectSec.Diam2 = Turbine->getDiametroRodeteOut();
//			FFunEffectSec.Diam3 = Turbine->getDiametroTuerca();
//			FFunEffectSec.Diam4 = Turbine->getDiametroTurbinaIn();
//			FFunEffectSec.CalculaGR = CalculaGR;
//			FFunEffectSec.Area = Area;
//			FFunEffectSec.n_limit = n_limit;
//			FFunEffectSec.Gamma = 1.4;
//			FFunEffectSec.R = 287.;
//			FFunEffectSec.BladeH = Turbine->getBladeHeight();
//		}
//
//	}
//	if (FNumPositions == 1)
//		FFixedTurbine = true;
//}

void TTurbineMap::LoadTurbineMap(FILE *Input, TTurbina *Turbine) {
	int rows = 0, Adiab = 0;
	double pos = 0., ang = 0.;
	double sp = 0.;
	double er = 0.;
	double ms = 0.;
	double ef = 0.;
	dVector v_sp;
	dVector v_er;
	dVector v_ms;
	dVector v_ef;
	dVector v_po;

	fscanf(Input, "%d ", &Adiab);
	if(Adiab == 0) {
		FIsAdiabatic = false;
		fscanf(Input, "%lf ", &FTempMeasure);
	}
	int Extrapolate = 0;
	fscanf(Input, "%d ", &Extrapolate);
	if(Extrapolate == 1) {
		FExtrapolate = true;
	}

	fscanf(Input, "%d ", &FNumPositions);
	FTurbPosition.resize(FNumPositions);
	for(int i = 0; i < FNumPositions; i++) {
		fscanf(Input, "%d %lf", &rows, &pos);
		v_sp.resize(rows);
		v_er.resize(rows);
		v_ms.resize(rows);
		v_ef.resize(rows);
		v_po.resize(rows);

		if(i == 0) {
			FMaxOpening = pos;
			FMinOpening = pos;
		} else {
			if(pos > FMaxOpening)
				FMaxOpening = pos;
			if(pos < FMinOpening)
				FMinOpening = pos;
		}

		for(int j = 0; j < rows; j++) {
			fscanf(Input, "%lf %lf %lf %lf", &sp, &er, &ms, &ef);
			v_sp[j] = sp;
			v_er[j] = er;
			v_ms[j] = ms;
			v_ef[j] = ef;
			v_po[j] = pos;

		}
		FTurbPosition[i].LoadTurbinPosition(v_sp, v_er, v_ms, v_ef, pos);

		FMapData.Position.insert(FMapData.Position.end(), v_po.begin(), v_po.end());
		FMapData.Speed.insert(FMapData.Speed.end(), v_sp.begin(), v_sp.end());
		FMapData.ERatio.insert(FMapData.ERatio.end(), v_er.begin(), v_er.end());
		FMapData.Mass.insert(FMapData.Mass.end(), v_ms.begin(), v_ms.end());
		FMapData.Efficiency.insert(FMapData.Efficiency.end(), v_ef.begin(), v_ef.end());
	}
	if(FNumPositions == 1)
		FFixedTurbine = true;
}

void TTurbineMap::LoadTurbineMapXML(xml_node node_turb) {

	dVector v_sp;
	dVector v_er;
	dVector v_ms;
	dVector v_ef;
	dVector v_po;

	int rows;
	double pos;

	FIsAdiabatic = GetAttributeAsBool(node_turb, "AdiabaticMap");
	if(!FIsAdiabatic) {
		FTempMeasure = GetXMLTemperature(node_turb, "MapTemperature");
	}
	FExtrapolate = GetAttributeAsBool(node_turb, "Extrapolate");

	FNumPositions = CountNodes(node_turb, "Trb_TurbineMap");

	FTurbPosition.resize(FNumPositions);
	int id = 0;
	int i = 0;
	int j = 0;
	for(xml_node node_map = GetNodeChild(node_turb, "Trb_TurbineMap"); node_map;
		node_map = node_map.next_sibling("Trb_TurbineMap")) {

		pos = GetAttributeAsInt(node_map, "RackPos");
		id = GetAttributeAsInt(node_map, "MapID");

		if(i == 0) {
			FMaxOpening = pos;
			FMinOpening = pos;
		} else {
			if(pos > FMaxOpening)
				FMaxOpening = pos;
			if(pos < FMinOpening)
				FMinOpening = pos;
		}
		i++;

		FTurbPosition[id - 1].ReadTurbinPositionXML(node_map);

		rows = CountNodes(node_map, "Tmp_TurbineMapPoint");
		v_sp.resize(rows);
		v_er.resize(rows);
		v_ms.resize(rows);
		v_ef.resize(rows);
		v_po.clear();
		v_po.resize(rows, pos);

		std::string RedSpeedUnit = "";
		std::string RedMassUnit = "";

		if(node_map.child("Units")) {
			xml_node node_unit = GetNodeChild(node_map, "Units");
			RedSpeedUnit = node_unit.attribute("ReducedSpeed").value();
			RedMassUnit = node_unit.attribute("ReducedMassFlow").value();
		}

		j = 0.;
		for(xml_node node_mappoint = GetNodeChild(node_map, "Tmp_TurbineMapPoint"); node_mappoint;
			node_mappoint = node_mappoint.next_sibling("Tmp_TurbineMapPoint")) {

			v_sp[j] = GetXMLReducedSpeed(node_mappoint, "Speed", RedSpeedUnit);
			v_ms[j] = GetXMLReducedMassFlow(node_mappoint, "Mass", RedMassUnit);
			v_er[j] = GetAttributeAsDouble(node_mappoint, "ExpRatio");
			v_ef[j] = GetAttributeAsDouble(node_mappoint, "Efficiency");
			j++;
		}

		FMapData.Position.insert(FMapData.Position.end(), v_po.begin(), v_po.end());
		FMapData.Speed.insert(FMapData.Speed.end(), v_sp.begin(), v_sp.end());
		FMapData.ERatio.insert(FMapData.ERatio.end(), v_er.begin(), v_er.end());
		FMapData.Mass.insert(FMapData.Mass.end(), v_ms.begin(), v_ms.end());
		FMapData.Efficiency.insert(FMapData.Efficiency.end(), v_ef.begin(), v_ef.end());

	}
	if(FNumPositions == 1)
		FFixedTurbine = true;
}

void TTurbineMap::CurrentEffectiveSection(double n, double er, double rack, double T10T00) {
	if(FFixedTurbine) {
		FTurbPosition[0].InterpolaPosicion(n, er);
		FStatorES = FTurbPosition[0].StatorSec();
		FRotorES = FTurbPosition[0].RotorSec() * Sqrt(T10T00);
		FEffTurb = FTurbPosition[0].Efficiency();
	} else {
		int i = 0;
		while(rack > FTurbPosition[i].Rack() && i < FNumPositions) {
			++i;
		}
		if(i == 0) {
			FTurbPosition[i].InterpolaPosicion(n, er);
			FStatorES = FTurbPosition[i].StatorSec();
			FRotorES = FTurbPosition[i].RotorSec() * sqrt(T10T00);
			FEffTurb = FTurbPosition[i].Efficiency();
		} else {
			FTurbPosition[i].InterpolaPosicion(n, er);
			FTurbPosition[i - 1].InterpolaPosicion(n, er);
			double DeltaRack = (rack - FTurbPosition[i - 1].Rack()) / (FTurbPosition[i].Rack() - FTurbPosition[i - 1].Rack());
			FStatorES = FTurbPosition[i - 1].StatorSec() * (1 - DeltaRack) + FTurbPosition[i].StatorSec() * DeltaRack;
			FRotorES = (FTurbPosition[i - 1].RotorSec() * (1 - DeltaRack) + FTurbPosition[i].RotorSec() * DeltaRack) * Sqrt(T10T00);
			FEffTurb = FTurbPosition[i - 1].Efficiency() * (1 - DeltaRack) + FTurbPosition[i].Efficiency() * DeltaRack;
		}
	}
}
// ---------------------------------------------------------------------------

void TTurbineMap::PrintFinalMap(FILE *fich) {
	for(int i = 0; i < FNumPositions; ++i) {
		fprintf(fich, "%lf\n", FTurbPosition[i].Rack());
		FTurbPosition[i].PrintTurbinePosition(fich);
	}
}

//void TTurbineMap::CalculateAdiabaticEfficiency(TTC_HTM *HTM, double TinC) {
//
//	if(!FIsAdiabatic) {
//
//		FILE *fres = fopen("TurbRes.dat", "w");
//
//		fprintf(fres, "Turbocharger_Speed(rmp)\t");
//		fprintf(fres, "Turbine_Mass(kg/s)\t");
//		fprintf(fres, "Expansion_ratio(-)\t");
//		fprintf(fres, "Adiabatic_efficiency(-)\t");
//		fprintf(fres, "Efficiency(-)\t");
//		fprintf(fres, "Turbine_Heat_Flow(W)\t");
//		fprintf(fres, "Compressor_Heat_Flow(W)\t");
//		fprintf(fres, "Oil_Heat_Flow(W)\t");
//		fprintf(fres, "Water_Heat_Flow(W)\t");
//		fprintf(fres, "Ambient_Heat_Flow(W)\t");
//		fprintf(fres, "Temperature_G(K)\t");
//		fprintf(fres, "Temperature_A(K)\t");
//		fprintf(fres, "Temperature_W(K)\t");
//		fprintf(fres, "Temperature_T(K)\t");
//		fprintf(fres, "Temperature_H1(K)\t");
//		fprintf(fres, "Temperature_H2(K)\t");
//		fprintf(fres, "Temperature_H3(K)\t");
//		fprintf(fres, "Temperature_C(K)\t");
//		fprintf(fres, "Temperature_Ai(K)\t");
//		fprintf(fres, "Temperature_Ao(K)\t");
//		fprintf(fres, "Temperature_Amb(K)\t");
//		fprintf(fres, "Temperature_O1(K)\t");
//		fprintf(fres, "Temperature_O2(K)\t");
//		fprintf(fres, "Temperature_O3(K)\t");
//		fprintf(fres, "Compressor_Power(W)\t");
//		fprintf(fres, "Turbine_Power(W)\t");
//		fprintf(fres, "Mechanical_Efficiency(W)\t");
//		fprintf(fres, "Turbine_In_Enthalpy(W)\n");
//
//		fclose(fres);
//
//		for(int i = 0; i < FNumPositions; ++i) {
//			FTurbPosition[i].AdiabaticEfficiency(HTM, FTempMeasure, TinC);
//		}
//	}
//	SearchMaximumEfficiency();
//}

void TTurbineMap::CalculateAdiabaticEfficiency(TTC_HTM2 *HTM, double TinC) {

	if(!FIsAdiabatic) {

		FILE *fres = fopen("TurbRes.dat", "w");

		fprintf(fres, "Turbocharger_Speed(rmp)\t");
		fprintf(fres, "Turbine_Mass(kg/s)\t");
		fprintf(fres, "Expansion_ratio(-)\t");
		fprintf(fres, "Adiabatic_efficiency(-)\t");
		fprintf(fres, "Efficiency(-)\t");
		fprintf(fres, "Turbine_Heat_Flow(W)\t");
		fprintf(fres, "Compressor_Heat_Flow(W)\t");
		fprintf(fres, "Oil_Heat_Flow(W)\t");
		fprintf(fres, "Water_Heat_Flow(W)\t");
		fprintf(fres, "Ambient_Heat_Flow(W)\t");
		fprintf(fres, "Temp_node_T(K)\t");
		fprintf(fres, "Temp_node_H1(K)\t");
		fprintf(fres, "Temp_node_H2(K)\t");
		fprintf(fres, "Temp_node_H3(K)\t");
		fprintf(fres, "Temp_node_C(K)\t");
		fprintf(fres, "Temp_turbine_in(K)\t");
		fprintf(fres, "Temp_compressor_in(K)\t");
		fprintf(fres, "Temp_oil_in(K)\t");
		fprintf(fres, "Temp_water_in(K)\t");
		fprintf(fres, "Temp_ambient(K)\t");
		fprintf(fres, "Temp_oil_out(K)\t");
		fprintf(fres, "Temp_compressor_out(K)\t");
		fprintf(fres, "Temp_water_out(K)\t");
		fprintf(fres, "Compressor_Power(W)\t");
		fprintf(fres, "Turbine_Power(W)\t");
		fprintf(fres, "Mechanical_Efficiency(W)\t");
		fprintf(fres, "Turbine_In_Enthalpy(W)\n");

		fclose(fres);

		int k = 0;
		for(int i = 0; i < FNumPositions; ++i) {
			FTurbPosition[i].AdiabaticEfficiency(HTM, FTempMeasure, TinC);
			for (int n = 0; n < FTurbPosition[i].size(); ++n){
				for (int m = 0; m < FTurbPosition[i].getSpeedLine(n).size(); ++m){
					FMapData.Efficiency[k] = FTurbPosition[i].getSpeedLine(n).MapEff(m);
					k++;
				}
			}
		}
		
	}
	SearchMaximumEfficiency();
}

void TTurbineMap::ProcessTurbineMap(TTurbina *Turbine) {

	double Area = __geom::Circle_area(Turbine->getDiametroTurbinaIn());
	double n_limit = 1.165;
	bool CalculaGR = false;

	FFunEffectSec.Diam1 = Turbine->getDiametroRodete();
	FFunEffectSec.Diam2 = Turbine->getDiametroRodeteOut();
	FFunEffectSec.Diam3 = Turbine->getDiametroTuerca();
	FFunEffectSec.Diam4 = Turbine->getDiametroTurbinaIn();
	FFunEffectSec.CalculaGR = true;
	FFunEffectSec.Area = Area;
	FFunEffectSec.n_limit = n_limit;
	FFunEffectSec.Gamma = 1.4;
	FFunEffectSec.R = __R::Air;
	FFunEffectSec.BladeH = Turbine->getBladeHeight();
	FFunEffectSec.CritAngle = Turbine->getCriticalAngle();

	if (FFixedTurbine) {
		if (Turbine->IsVaneless()) {
			BladeMechanism = new stBMVaneless(Turbine->getDiametroRodete() / 2, Turbine->getBladeHeight());
		}
		else {
			BladeMechanism = new stBMVaned(0, Turbine->getVgTtoAlpha2(), Turbine->getDiametroRodete(), 0, Turbine->getFz0(),
				Turbine->getBladeHeight());
		}
	}
	else {
		BladeMechanism = new stBMVaned(Turbine->getVgTtoAlpha1(), Turbine->getVgTtoAlpha2(), Turbine->getFr2geom(),
			Turbine->getFlte(), Turbine->getFz0(), Turbine->getBladeHeight());
	}


	if (FExtrapolate) {
		FILE *fich = fopen("TurbineExtrp.txt", "w");
		fclose(fich);

		ExtrapolateTurbineMap(Turbine);

		fich = fopen("TurbineExtrp.txt", "a");
		fprintf(fich, "a = %lf\t", AeffFun->get_a());
		fprintf(fich, "b1 = %lf\t", AeffFun->get_b1());
		fprintf(fich, "b2 = %lf\t", AeffFun->get_b2());
		fprintf(fich, "c1 = %lf\t", AeffFun->get_c1());
		fprintf(fich, "c2 = %lf\t", AeffFun->get_c2());
		fprintf(fich, "d1 = %lf\t", AeffFun->get_d1());
		fprintf(fich, "d2 = %lf\n", AeffFun->get_d2());

		for (int i = 0; i < 6; i++) {
			fprintf(fich, "z%d = %lf\t", i, EffFun->get_z(i));
		}

		fprintf(fich, "\n");

		for (int i = 0; i < FNumPositions; ++i) {

			FTurbPosition[i].ExtrapolateTurbinePosition(fich, Turbine->getDiametroRodete(), Turbine->getBladeHeight(),
				Turbine->getDiametroRodeteOut(), Turbine->getDiametroTuerca(),
				Turbine->getDiametroTurbinaIn(), Turbine->getBeta());

		}
		fclose(fich);

		PrintExtrapolatedMap(Turbine->getDiametroRodete());
	}
	else{
		for (int i = 0; i < FNumPositions; ++i) {
			FTurbPosition[i].putAngle(BladeMechanism->ToAlpha(FTurbPosition[i].Rack()));
			if (FTurbPosition[i].getAngle() > Turbine->getCriticalAngle())
				CalculaGR = false;
			else
				CalculaGR = true;

			FTurbPosition[i].EffectiveArea(Area, CalculaGR, Turbine->getDiametroRodete(), Turbine->getDiametroRodeteOut(),
				Turbine->getDiametroTuerca(), n_limit);

			FTurbPosition[i].SearchMapLimits();
		}
	}


	SearchMaximumEfficiency();
}

void TTurbineMap::PrintExtrapolatedMap(double Dw) {
	double n;
	double gm = 1.4;
	double g1 = 2 / (gm - 1);
	double g2 = (gm - 1) / gm;
	double Cp = __R::Air / g2;
	dVector er;
	double rack, mred, eff;

	er.push_back(1.1);
	while(er.back() <= 4) {
		er.push_back(er.back() + 0.1);
	}
	FILE *fich = fopen("TurbinePlot.txt", "w");

	double K3;
	double bsr;

	if(FTurbPosition.front().Rack() > 0 && FTurbPosition.size() > 1) {
		rack = 0.;

		for(int j = 0; j < FTurbPosition.front().size(); j++) {
			n = FTurbPosition.front().getSpeedLine(j).Speed();
			K3 = pow2(n * __cons::Pi * Dw) / 2 / Cp;
			for(int k = 0; k < er.size(); k++) {
				bsr = sqrt(K3 / (1 - pow(1 / er[k], g2)));
				if(bsr < 1.5) {
					InterpolateTurbineMap(n, er[k], bsr, rack, mred, eff);

					if(eff > 0.05) {
						fprintf(fich, "%g\t", rack);
						fprintf(fich, "%g\t", n);
						fprintf(fich, "%g\t", er[k]);
						fprintf(fich, "%g\t", mred);
						fprintf(fich, "%g\t", eff);
						fprintf(fich, "%g\n", bsr);
					}
				}

			}
			fprintf(fich, "\n");
		}
	}

	for(int i = 0; i < FNumPositions; ++i) {
		rack = FTurbPosition[i].Rack();
		for(int j = 0; j < FTurbPosition[i].size(); j++) {
			n = FTurbPosition[i].getSpeedLine(j).Speed();
			K3 = pow2(n * __cons::Pi * Dw) / 2 / Cp;
			for(int k = 0; k < er.size(); k++) {
				bsr = sqrt(K3 / (1 - pow(1 / er[k], g2)));
				if(bsr < 1.5) {
					InterpolateTurbineMap(n, er[k], bsr, rack, mred, eff);

					if(eff > 0.05) {
						fprintf(fich, "%g\t", rack);
						fprintf(fich, "%g\t", n);
						fprintf(fich, "%g\t", er[k]);
						fprintf(fich, "%g\t", mred);
						fprintf(fich, "%g\t", eff);
						fprintf(fich, "%g\n", bsr);
					}
				}

			}
			fprintf(fich, "\n");
		}
	}

	if(FTurbPosition.back().Rack() < 100 && FTurbPosition.size() > 1) {
		rack = 100.;

		for(int j = 0; j < FTurbPosition.back().size(); j++) {
			n = FTurbPosition.back().getSpeedLine(j).Speed();
			K3 = pow2(n * __cons::Pi * Dw) / 2 / Cp;
			for(int k = 0; k < er.size(); k++) {
				bsr = sqrt(K3 / (1 - pow(1 / er[k], g2)));
				if(bsr < 1.5) {
					InterpolateTurbineMap(n, er[k], bsr, rack, mred, eff);

					if(eff > 0.05) {
						fprintf(fich, "%g\t", rack);
						fprintf(fich, "%g\t", n);
						fprintf(fich, "%g\t", er[k]);
						fprintf(fich, "%g\t", mred);
						fprintf(fich, "%g\t", eff);
						fprintf(fich, "%g\n", bsr);
					}
				}

			}
			fprintf(fich, "\n");
		}
	}
	fclose(fich);

}

void TTurbineMap::InterpExtendedMap(double n, double er, double bsr, double rack) {

	double Mass;

	InterpolateTurbineMap(n, er, bsr, rack, Mass, FEffTurb);

	FFunEffectSec.Angle = BladeMechanism->ToAlpha(rack);
	if(FFunEffectSec.Angle > FFunEffectSec.CritAngle)
		FFunEffectSec.CalculaGR = false;
	else
		FFunEffectSec.CalculaGR = true;
	FFunEffectSec.fun(FStatorES, FRotorES, Mass, er, FEffTurb, n);

}

void TTurbineMap::CorrectRotorEffArea(double T10T00){
	FRotorES = FRotorES * sqrt(T10T00);
}

void TTurbineMap::InterpolateTurbineMap(double n, double er, double bsr, double rack, double &mred, double &eff) {

	double Aeff, Aeff1, Aeff2, Eff1 = eff, Eff2 = eff;
	int i = 0;
	double nmin;
	double nmax;
	double p;

	double b = bsr;
	if (b > 1.5)
		b = 1.5;
	else if (b < 0.01)
		b = 0.01;

	ExtrapoltePoint(n, er, b, rack, Aeff, eff);

/*	if(FFixedTurbine) {
		nmin = FTurbPosition[0].getSpeedMin();
		nmax = FTurbPosition[0].getSpeedMax();
		if(n < nmin) {
			ExtrapoltePoint(n, er, bsr, rack, Aeff1, Eff1);
			FTurbPosition[0].MapInterpolation(nmin, bsr, er, Aeff2, Eff2);
			p = 1 / exp(20 * (nmin - n) / n);
			Aeff = Aeff1 * (1 - p) + Aeff2 * p;
			eff = Eff1 * (1 - p) + Eff2 * p;
		} else if(n > nmax) {
			ExtrapoltePoint(n, er, bsr, rack, Aeff1, Eff1);
			FTurbPosition[0].MapInterpolation(nmax, bsr, er, Aeff2, Eff2);
			p = 1 / exp(20 * (n - nmax) / n);
			Aeff = Aeff1 * (1 - p) + Aeff2 * p;
			eff = Eff1 * (1 - p) + Eff2 * p;
		} else {
			FTurbPosition[0].MapInterpolation(n, bsr, er, Aeff, eff);
		}
	} else {
		double n2 = n;

		if(rack < FMinOpening) {
			ExtrapoltePoint(n, er, bsr, rack, Aeff1, Eff1);
			FTurbPosition[0].MapInterpolation(n, bsr, er, Aeff2, Eff2);
			p = 1 / exp(20 * (FMinOpening - rack) / FMinOpening);
			Aeff = Aeff1 * (1 - p) + Aeff2 * p;
			eff = Eff1 * (1 - p) + Eff2 * p;
		} else if(rack > FMaxOpening) {
			ExtrapoltePoint(n, er, bsr, rack, Aeff1, Eff1);
			FTurbPosition.back().MapInterpolation(n, bsr, er, Aeff2, Eff2);
			p = 1 / exp(20 * (rack - FMaxOpening) / FMaxOpening);
			Aeff = Aeff1 * (1 - p) + Aeff2 * p;
			eff = Eff1 * (1 - p) + Eff2 * p;
		} else {
			int i = 1;
			while(rack > FTurbPosition[i].Rack() && i < FNumPositions) {
				++i;
			}
			double DeltaRack = (rack - FTurbPosition[i - 1].Rack()) / (FTurbPosition[i].Rack() - FTurbPosition[i - 1].Rack());
			nmin = Max(FTurbPosition[i].getSpeedMin(), FTurbPosition[i - 1].getSpeedMin());
			nmax = Min(FTurbPosition[i].getSpeedMax(), FTurbPosition[i - 1].getSpeedMax());
			if(n < nmin) {
				ExtrapoltePoint(n, er, bsr, rack, Aeff, eff);
				FTurbPosition[i].MapInterpolation(n, bsr, er, Aeff1, Eff1);
				FTurbPosition[i - 1].MapInterpolation(n, bsr, er, Aeff2, Eff2);
				p = 1 / exp(20 * (nmin - n) / n);
				Aeff = Aeff * (1 - p) + (Aeff2 * (1 - DeltaRack) + Aeff1 * DeltaRack) * p;
				eff = eff * (1 - p) + (Eff2 * (1 - DeltaRack) + Eff1 * DeltaRack) * p;

			} else if(n > nmax) {
				ExtrapoltePoint(n, er, bsr, rack, Aeff, eff);
				FTurbPosition[i].MapInterpolation(n, bsr, er, Aeff1, Eff1);
				FTurbPosition[i - 1].MapInterpolation(n, bsr, er, Aeff2, Eff2);
				p = 1 / exp(20 * (n - nmax) / n);
				Aeff = Aeff * (1 - p) + (Aeff2 * (1 - DeltaRack) + Aeff1 * DeltaRack) * p;
				eff = eff * (1 - p) + (Eff2 * (1 - DeltaRack) + Eff1 * DeltaRack) * p;
			} else {
				FTurbPosition[i].MapInterpolation(n, bsr, er, Aeff1, Eff1);
				FTurbPosition[i - 1].MapInterpolation(n, bsr, er, Aeff2, Eff2);

				Aeff = Aeff2 * (1 - DeltaRack) + Aeff1 * DeltaRack;
				eff = Eff2 * (1 - DeltaRack) + Eff1 * DeltaRack;
			}
		}

	}*/
	if(er < 1)
		er = 1 + 1e-3;

	mred = Aeff * __Gamma::G * 1e6 / sqrt(__Gamma::gxR) * pow(er, -1 / __Gamma::G) * sqrt(2 / __Gamma::G_1 * (1 - pow(er,
			-__Gamma::G_8)));

}

void TTurbineMap::ExtrapoltePoint(double n, double er, double bsr, double rack, double& aeff, double& eff) {

	double eff0 = 0.;
	double aeff0 = 0.;

	do {
		eff0 = eff;
		aeff0 = aeff;

		aeff = AeffFun->val(bsr, er, eff, rack);
		eff = EffFun->val(bsr, aeff, n, rack);

	} while(fabs((eff0 - eff) / eff0) > 1e-6 && fabs((aeff0 - aeff) / aeff0) > 1e-6);
}

void TTurbineMap::SearchMaximumEfficiency() {

	for(int i = 0; i < FNumPositions; i++) {
		FOpening.push_back(FTurbPosition[i].Rack());
		FMaxEff.push_back(FTurbPosition[i].MaximumEfficiency());
	}

	if(FNumPositions > 1)
		iMaxEff = new Hermite_interp(FOpening, FMaxEff);
}

double TTurbineMap::MaxEfficiency(double rack) {
	if(FNumPositions == 1)
		return FMaxEff[0];
	else
		return iMaxEff->interp(rack);
}

void TTurbineMap::ExtrapolateTurbineMap(TTurbina *Turbina) {

	dVector MapK3;
	dVector MapBSR;
	dVector MapAeff;

	double gm = 1.4;

	double g1 = 2 / (gm - 1);
	double g2 = (gm - 1) / gm;
	double Cp = __R::Air / g2;
	double gR = sqrt(gm * __R::Air);

	double K3 = 0.;

	for(int i = 0; i < FMapData.Mass.size(); ++i) {
		K3 = pow2(FMapData.Speed[i] * __cons::Pi * Turbina->getDiametroRodete()) / 2 / Cp;
		MapK3.push_back(K3);
		MapBSR.push_back(sqrt(K3 / (1 - pow(1 / FMapData.ERatio[i], g2))));

		MapAeff.push_back(1e-6 * FMapData.Mass[i] * gR / (gm * pow(1 / FMapData.ERatio[i],
						  1 / gm) * sqrt(g1 * (1 - pow(1 / FMapData.ERatio[i], g2)))));
	}

	double AR = __geom::Ring_area(Turbina->getDiametroTuerca(), Turbina->getDiametroRodeteOut());

	double diam = pow2((Turbina->getDiametroRodeteOut() + Turbina->getDiametroTuerca()) / 2 / Turbina->getDiametroRodete())
				  - 1;
	double RendMean = 0.8;

	AeffFun = new TAeffFunction_v2(AR, diam, gm, RendMean, BladeMechanism);

	AeffFun->fit(MapBSR, FMapData.Position, FMapData.Efficiency, FMapData.ERatio, MapAeff);

	double K1 = 2 * (diam + 1);
	double A0 = __geom::Circle_area(Turbina->getDiametroTurbinaIn());
	double A3 = __geom::Ring_area(Turbina->getDiametroTuerca(), Turbina->getDiametroRodeteOut());

	EffFun = new TEffTFunction_v2(K1, A0, A3, Turbina->getBeta(), BladeMechanism, AeffFun->get_c1(), AeffFun->get_c2(), gm,
								  Turbina->getDiametroRodete());

	EffFun->fit(MapBSR, FMapData.Speed, FMapData.Efficiency, FMapData.Position, MapAeff);

}

#pragma package(smart_init)
