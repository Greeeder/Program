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

#include "TTurbPosition.h"
#include "Globales.h"

// ---------------------------------------------------------------------------

TTurbPosition::TTurbPosition() {

}

TTurbPosition::~TTurbPosition() {
	FSpeedLine.clear();

}

void TTurbPosition::ReadTurbinPosition(FILE *Input, int rows, double pos, double ang) {
	double SP = 0., ER = 0., MF = 0., EF = 0., SP0 = 0., ER0 = 0.;

	FPosition = pos;
	FAngle = ang;

	for(int i = 0; i < rows; i++) {
		fscanf(Input, "%lf %lf %lf %lf", &SP, &ER, &MF, &EF);
		if(i == 0) {
			FLines = 1;
			FSpeedLine.resize(1);
			FSpeedLine[0].AsignaValores(ER, MF, EF);
			FSpeedLine[0].PutSpeed(SP);
		} else {
			if(SP == SP0) {
				if(ER > ER0) {
					FSpeedLine[FLines - 1].AsignaValores(ER, MF, EF);
				} else {
					cout << "ERROR: The turbine map is not correctely ordered" << endl;
					// error
				}
			} else if(SP > SP0) {
				FLines++;
				FSpeedLine.resize(FLines);
				FSpeedLine[FLines - 1].AsignaValores(ER, MF, EF);
				FSpeedLine[FLines - 1].PutSpeed(SP);
			}
		}
		SP0 = SP;
		ER0 = ER;
	}
	FSpeedMin = FSpeedLine[0].Speed();
	FSpeedMax = FSpeedLine[FLines - 1].Speed();

}

void TTurbPosition::LoadTurbinPosition(dVector SP, dVector ER, dVector MF, dVector EF, double pos) {

	FPosition = pos;

	for(int i = 0; i < SP.size(); i++) {
		if(i == 0) {
			FLines = 1;
			FSpeedLine.resize(1);
			FSpeedLine[0].AsignaValores(ER[i], MF[i], EF[i]);
			FSpeedLine[0].PutSpeed(SP[i]);
		} else {
			if(SP[i] == SP[i - 1]) {
				if(ER[i] > ER[i - 1]) {
					FSpeedLine[FLines - 1].AsignaValores(ER[i], MF[i], EF[i]);
				} else {
					cout << "ERROR: The turbine map is not correctely ordered" << endl;
					// error
				}
			} else if(SP[i] > SP[i - 1]) {
				FLines++;
				FSpeedLine.resize(FLines);
				FSpeedLine[FLines - 1].AsignaValores(ER[i], MF[i], EF[i]);
				FSpeedLine[FLines - 1].PutSpeed(SP[i]);
			}
		}
	}
	FSpeedMin = FSpeedLine[0].Speed();
	FSpeedMax = FSpeedLine[FLines - 1].Speed();

}

void TTurbPosition::ReadTurbinPositionXML(xml_node node_map) {
	double SP, ER, MF, EF, SP0, ER0;

	//rows = CountNodes(node_map,"TurbineMapPoint");
	FPosition = GetAttributeAsInt(node_map, "RackPos");
	//FAngle = GetXMLAngle(node_map, "BladeAngle");

	std::string RedSpeedUnit = "";
	std::string RedMassUnit = "";
	if(node_map.child("Units")) {
		xml_node node_unit = GetNodeChild(node_map, "Units");
		RedSpeedUnit = node_unit.attribute("ReducedSpeed").value();
		RedMassUnit = node_unit.attribute("ReducedMassFlow").value();
	}

	int i = 0;

	for(xml_node node_mpt = GetNodeChild(node_map, "Tmp_TurbineMapPoint"); node_mpt;
		node_mpt = node_mpt.next_sibling("Tmp_TurbineMapPoint")) {

		SP = GetXMLReducedSpeed(node_mpt, "Speed", RedSpeedUnit);
		MF = GetXMLReducedMassFlow(node_mpt, "Mass", RedMassUnit);
		ER = GetAttributeAsDouble(node_mpt, "ExpRatio");
		EF = GetAttributeAsDouble(node_mpt, "Efficiency");

		if(i == 0) {
			FLines = 1;
			FSpeedLine.resize(1);
			FSpeedLine[0].AsignaValores(ER, MF, EF);
			FSpeedLine[0].PutSpeed(SP);
		} else {
			if(SP == SP0) {
				if(ER > ER0) {
					FSpeedLine[FLines - 1].AsignaValores(ER, MF, EF);
				} else {
					// error
				}
			} else if(SP > SP0) {
				FLines++;
				FSpeedLine.resize(FLines);
				FSpeedLine[FLines - 1].AsignaValores(ER, MF, EF);
				FSpeedLine[FLines - 1].PutSpeed(SP);
			}
		}
		SP0 = SP;
		ER0 = ER;
		i++;
	}
	FSpeedMin = FSpeedLine[0].Speed();
	FSpeedMax = FSpeedLine[FLines - 1].Speed();

}

void TTurbPosition::EffectiveArea(double Area, bool CalculaGR, double Diam1, double Diam2, double Diam3,
								  double n_limit) {
	for(int i = 0; i < FLines; i++) {
		FSpeedLine[i].EffectiveSection(Area, CalculaGR, FAngle, Diam1, Diam2, Diam3, n_limit);
		FSpeedLine[i].Adimensionaliza();
	}
}

void TTurbPosition::CalculatePower(double Tin) {
	for(int i = 0; i < FLines; i++) {
		FSpeedLine[i].CalculatePower(Tin);
	}
}

void TTurbPosition::PutPosition(double Pos) {
	FPosition = Pos;
}

void TTurbPosition::SearchMapLimits() {
	for(int i = 0; i < FLines; i++) {
		FSpeed.push_back(FSpeedLine[i].Speed());
		FERMax.push_back(FSpeedLine[i].ERMax());
		FERMin.push_back(FSpeedLine[i].ERMin());
	}
	h_FERMax.resize(FLines);
	h_FERMin.resize(FLines);

	Hermite(FLines, FSpeed, FERMax, &h_FERMax);
	Hermite(FLines, FSpeed, FERMin, &h_FERMin);

}

void TTurbPosition::InterpolaPosicion(double n, double er) {
	double CDStator0, CDStator1, CDRotor0, CDRotor1, ERAdim, DeltaN, ERMax, ERMin, Eff0, Eff1;

	if(n <= FSpeedMin) {
		if(er < FSpeedLine[0].ERMin()) {
			ERAdim = 0;
		} else if(er > FSpeedLine[0].ERMax()) {
			ERAdim = 1;
		} else {
			ERAdim = (er - FSpeedLine[0].ERMin()) / (FSpeedLine[0].ERMax() - FSpeedLine[0].ERMin());
		}

		FStatorSec = FSpeedLine[0].Stator(ERAdim);
		FRotorSec = FSpeedLine[0].Rotor(ERAdim);
		FEfficiency = FSpeedLine[0].Efficiency(ERAdim);
	} else if(n >= FSpeedMax) {
		if(er < FSpeedLine[FLines - 1].ERMin()) {
			ERAdim = 0;
		} else if(er > FSpeedLine[FLines - 1].ERMax()) {
			ERAdim = 1;
		} else {
			ERAdim = (er - FSpeedLine[FLines - 1].ERMin()) / (FSpeedLine[FLines - 1].ERMax() - FSpeedLine[FLines - 1].ERMin());
		}

		FStatorSec = FSpeedLine[FLines - 1].Stator(ERAdim);
		FRotorSec = FSpeedLine[FLines - 1].Rotor(ERAdim);
		FEfficiency = FSpeedLine[FLines - 1].Efficiency(ERAdim);
	} else {
		int i = 0;
		while(FSpeedLine[i].Speed() < n) {
			++i;
		}

		ERMax = EvaluaHermite(n, FLines, FSpeed, FERMax, h_FERMax);
		ERMin = EvaluaHermite(n, FLines, FSpeed, FERMin, h_FERMin);
		ERAdim = (er - ERMin) / (ERMax - ERMin);
		CDStator0 = FSpeedLine[i - 1].Stator(ERAdim);
		CDStator1 = FSpeedLine[i].Stator(ERAdim);
		CDRotor0 = FSpeedLine[i - 1].Rotor(ERAdim);
		CDRotor1 = FSpeedLine[i].Rotor(ERAdim);
		Eff0 = FSpeedLine[i - 1].Efficiency(ERAdim);
		Eff1 = FSpeedLine[i].Efficiency(ERAdim);
		DeltaN = (n - FSpeed[i - 1]) / (FSpeed[i] - FSpeed[i - 1]);
		FStatorSec = CDStator0 * (1 - DeltaN) + CDStator1 * DeltaN;
		FRotorSec = CDRotor0 * (1 - DeltaN) + CDRotor1 * DeltaN;
		FEfficiency = Eff0 * (1 - DeltaN) + Eff1 * DeltaN;
	}
}

void TTurbPosition::PrintTurbinePosition(FILE *fich) {
	for(int i = 0; i < FLines; i++) {
		FSpeedLine[i].PrintEffectiveSection(fich);
	}
}

double TTurbPosition::StatorSec() {
	return FStatorSec;
}

double TTurbPosition::RotorSec() {
	return FRotorSec;
}

double TTurbPosition::Rack() {
	return FPosition;
}

double TTurbPosition::Efficiency() {
	return FEfficiency;
}

double TTurbPosition::MaxPowerLimit(double rtc) {

	int i = 0;
	while(FSpeedLine[i].Speed() < rtc) {
		++i;
	}
	double DeltaN = (rtc - FSpeed[i - 1]) / (FSpeed[i] - FSpeed[i - 1]);
	double PowerMax0 = FSpeedLine[i - 1].PowerMax();
	double PowerMax1 = FSpeedLine[i].PowerMax();
	return PowerMax0 * (1 - DeltaN) + PowerMax1 * DeltaN;

}

double TTurbPosition::MinPowerLimit(double rtc) {

	int i = 0;
	while(FSpeedLine[i].Speed() < rtc) {
		++i;
	}
	double DeltaN = (rtc - FSpeed[i - 1]) / (FSpeed[i] - FSpeed[i - 1]);
	double PowerMin0 = FSpeedLine[i - 1].PowerMin();
	double PowerMin1 = FSpeedLine[i].PowerMin();
	return PowerMin0 * (1 - DeltaN) + PowerMin1 * DeltaN;

}

//void TTurbPosition::AdiabaticEfficiency(TTC_HTM *HTM, double TinT, double TinC) {
//
//	double MaxEff = MaximumEfficiency();
//
//	for(int i = 0; i < FLines; ++i) {
//		FSpeedLine[i].GetAdiabaticEfficiency(HTM, TinT, TinC, MaxEff);
//	}
//}

void TTurbPosition::AdiabaticEfficiency(TTC_HTM2 *HTM, double TinT, double TinC) {

	double MaxEff = MaximumEfficiency();

	for(int i = 0; i < FLines; ++i) {
		FSpeedLine[i].GetAdiabaticEfficiency(HTM, TinT, TinC, MaxEff);
	}
}

double TTurbPosition::MaximumEfficiency() {

	double MaxEff = FSpeedLine[0].MaximumEfficiency();

	for(int i = 1; i < FSpeedLine.size(); i++) {
		if(FSpeedLine[i].MaximumEfficiency() > MaxEff) {
			MaxEff = FSpeedLine[i].MaximumEfficiency();
		}
	}
	return MaxEff;
}

void TTurbPosition::ExtrapolateTurbinePosition(FILE *fich, double Dwheel, double BladeH, double Dout, double Dnut,
		double Din, double beta) {

	double gm = 1.4;
	// double R = 287;

	double g1 = 2 / (gm - 1);
	double g2 = (gm - 1) / gm;
	double Cp = __R::Air / g2;
	double gR = sqrt(gm * __R::Air);
	double K3 = 0.;

	fprintf(fich, "Speed\tExpansion_Ratio\tReducedMass\tEfficiency\n");
	for(int i = 0; i < FLines; ++i) {
		K3 = pow2(FSpeedLine[i].MapSpeed() * __cons::Pi * Dwheel) / 2 / Cp;
		MapK3.push_back(K3);
		for(int j = 0; j < FSpeedLine[i].size(); ++j) {
			MapSpeed.push_back(FSpeedLine[i].MapSpeed());
			MapER.push_back(FSpeedLine[i].MapER(j));
			MapEff.push_back(FSpeedLine[i].MapEff(j));
			MapMass.push_back(FSpeedLine[i].MapMass(j));
			MapBSR.push_back(sqrt(K3 / (1 - pow(1 / FSpeedLine[i].MapER(j), g2))));

			MapAeff.push_back(1e-6 * FSpeedLine[i].MapMass(j) * gR / (gm * pow(1 / FSpeedLine[i].MapER(j),
							  1 / gm) * sqrt(g1 * (1 - pow(1 / FSpeedLine[i].MapER(j), g2)))));
			fprintf(fich, "%lf\t%lf\t%lf\t%lf\n", FSpeedLine[i].MapSpeed(), FSpeedLine[i].MapER(j), FSpeedLine[i].MapMass(j),
					FSpeedLine[i].MapEff(j));
		}
	}
	double AR = __geom::Ring_area(Dnut, Dout);
	double AN = __cons::Pi_x_2 * Dwheel * BladeH * (0.2 + FPosition / 100 * 0.8);

	double diam = pow2((Dout + Dnut) / 2 / Dwheel) - 1;
	double K1 = 2 * (diam + 1);
	double RendMean = Mean(MapEff);
	RendMean = 0.8;
	double A0 = __geom::Circle_area(Din);

	AeffFun = new TAeffFunction(AR, AN, diam, gm, __R::Air, RendMean, Dwheel);

	fprintf(fich, "TurbinePosition = %lf\n", FPosition);

	AeffFun->fit(MapSpeed, MapER, MapEff, MapMass);

	fprintf(fich, "AR = %lf\tAN = %lf\tdiam = %lf\tgamma = %lf\tR = %lf\tRendMean = %lf\tDwheel = %lf\n", AR, AN, diam, gm, __R::Air, RendMean, Dwheel);

	fprintf(fich, "a = %lf\t b = %lf\tc = %lf\td = %lf\n", AeffFun->get_a(), AeffFun->get_b(), AeffFun->get_c(), AeffFun->get_d());

	z_Eff.resize(2);
	for(int j = 0; j < 2; ++j)
		z_Eff[j].resize(FLines);
	dVector z1;
	dVector z2;
	dVector z3;

	fprintf(fich, "K1 = %lf\tA0 = %lf\tbeta = %lf\n", K1, A0, beta);

	for(int i = 0; i < FLines; ++i) {
		FSpeedLine[i].AsignEfficiencyData(K1, A0, beta, MapK3[i], gm, Dwheel);

		FSpeedLine[i].Fit_Efficiency(0, AeffFun);

		for(int j = 0; j < 2; ++j) {
			z_Eff[j][i] = FSpeedLine[i].GetEffCoef(j);
			// fprintf(fich, "%lf ", FSpeedLine[i].GetEffCoef(j));
		}

		z1.push_back(FSpeedLine[i].GetEffCoef(0));
		z2.push_back(FSpeedLine[i].GetEffCoef(1));
		z3.push_back(FSpeedLine[i].GetEffCoef(2));
		// fprintf(fich, "%lf\n", FSpeedLine[i].GetEffCoef(2));

	}
	double z1_mean = Mean(z1);
	double z2_mean = Mean(z2);
	z3_mean = Mean(z3);

	for(int i = 0; i < FLines; ++i) {

		FSpeed.push_back(FSpeedLine[i].Speed());

		FSpeedLine[i].Put_z(z1_mean, z2_mean, z3_mean);

		FSpeedLine[i].Fix_z3(z3_mean);

		FSpeedLine[i].Fit_Efficiency(1, AeffFun);

		for(int j = 0; j < 2; ++j) {
			z_Eff[j][i] = FSpeedLine[i].GetEffCoef(j);
			fprintf(fich, "z%d = %lf ", j, FSpeedLine[i].GetEffCoef(j));
		}
		// fprintf(fich, "%lf\n", z3_mean);
		fprintf(fich, "z2 = %lf\n", FSpeedLine[i].GetEffCoef(2));
	}
	z1_interp = new Hermite_interp(FSpeed, z_Eff[0]);
	z2_interp = new Hermite_interp(FSpeed, z_Eff[1]);

	EffTFun = new TEffTFunction(K1, A0, beta, K3, gm, Dwheel);

	EffTFun->put_z(2, z3_mean);

}

void TTurbPosition::MapInterpolation(double Speed, double BSR, double ER, double& Aeff_interp, double &Eff_interp) {

	double Eff_interp0 = 0.;
	double Aeff_interp0 = 0.;

	double z1 = z1_interp->interp(Speed);
	EffTFun->put_z(0, z1);
	double z2 = z2_interp->interp(Speed);
	EffTFun->put_z(1, z2);

	double ER_mod = 2 / (ER + 1);

	// Eff_interp = 0.7;
	Aeff_interp = AeffFun->val(BSR, ER_mod, Eff_interp);

	do {
		Eff_interp0 = Eff_interp;

		Aeff_interp0 = Aeff_interp;

		EffTFun->put_z(0, z1);
		EffTFun->put_z(1, z2);

		Eff_interp = EffTFun->val(BSR, Aeff_interp, ER);

		if(Eff_interp < 0.05)
			Eff_interp = 0.05;
		else if(Eff_interp > 0.99)
			Eff_interp = 0.99;

		Aeff_interp = AeffFun->val(BSR, ER_mod, Eff_interp);

	} while(fabs((Eff_interp0 - Eff_interp) / Eff_interp0) > 1e-6
			&& fabs((Aeff_interp0 - Aeff_interp) / Aeff_interp0) > 1e-6);

}

#pragma package(smart_init)
