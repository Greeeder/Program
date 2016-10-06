/**
* @file TSAEMap.cpp
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @date 19 de mar. de 2016
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
* Compressor map in SAE format.
*/
#include <Fit_MRQ.h>
#include <Math_wam.h>
#include <cmath>

#pragma hdrstop

#include "TSAEMap.h"

// ---------------------------------------------------------------------------

TSAEMap::TSAEMap(int i) :
	TCompressorMap() {

	FNumeroCompresor = i;

	FMassMAX_int = NULL;
	FPresMAX_int = NULL;
	FEffMAX_int = NULL;

	FIsAdiabatic = true;
}

TSAEMap::~TSAEMap() {
}

void TSAEMap::ReadSAECompressorMap(FILE *fich) {

	double speed = 0., mass = 0., pres = 0., eff = 0.;
	int i = 0; // Curva de isoregimen
	int j = 0; // Puntos de la curva
	int k = 0;
	int puntos = 0;
	double speedmax = 0, massmax = 0, presmax = 1, effmax = 0;
	int points = 0;
	double dwheel = 0.038;
	//TODO Read wheel diameter

	fscanf(fich, "%lf %lf %lf ", &FMassMultiplier, &FCRMultiplier, &FEffMultiplier);
	fscanf(fich, "%d", &points);
	FSpeed.resize(i + 1);
	FMass.resize(i + 1);
	FPres.resize(i + 1);
	FEff.resize(i + 1);

	while(k < points) {
		fscanf(fich, "%lf %lf %lf %lf", &speed, &mass, &pres, &eff);
		mass *= FMassMultiplier;
		pres = (pres - 1.) * FCRMultiplier + 1.;
		eff *= FEffMultiplier;
		if(k == 0)
			FSpeedVec.push_back(speed);
		if(!feof(fich)) {
			if(j > 0) {
				if(speed != FSpeed[i][j - 1]) {
					i++;
					j = 0;
					FSpeed.resize(i + 1);
					FMass.resize(i + 1);
					FPres.resize(i + 1);
					FEff.resize(i + 1);
					FSpeedVec.push_back(speed);
				}
			}
			FSpeed[i].push_back(speed);
			FMass[i].push_back(mass);
			FPres[i].push_back(pres);
			FEff[i].push_back(eff);
			j++;
		}
		k++;
	}

	FNumLines = FSpeed.size();

}

void TSAEMap::ReadSAECompressorMapXML(xml_node node_map) {

	double speed, mass, pres, eff;
	int i = 0; // Curva de isoregimen
	int j = 0; // Puntos de la curva
	int k = 0;
	int puntos;

	FMassMultiplier = GetAttributeAsDouble(node_map, "MassMultiplier");
	FCRMultiplier = GetAttributeAsDouble(node_map, "CRMultiplier");
	FEffMultiplier = GetAttributeAsDouble(node_map, "EffMultiplier");

	int points = CountNodes(node_map, "Cmp_CompMapPoint");
	FSpeed.resize(i + 1);
	FMass.resize(i + 1);
	FPres.resize(i + 1);
	FEff.resize(i + 1);
	xml_node unit_node = node_map.child("Units");
	std::string unitspeed = unit_node.attribute("RotationalSpeed").value();
	std::string unitmass = unit_node.attribute("MassFlow").value();

	for(xml_node node_mappoint = GetNodeChild(node_map, "Cmp_CompMapPoint"); node_mappoint;
		node_mappoint = node_mappoint.next_sibling("Cmp_CompMapPoint")) {
		speed = GetXMLRotationalSpeed(node_mappoint, "CorSpeed", unitspeed);
		mass = GetXMLMassFlow(node_mappoint, "CorMassFlow", unitmass) * FMassMultiplier;
		pres = (GetAttributeAsDouble(node_mappoint, "CompRatio") - 1) * FCRMultiplier + 1;
		eff = GetAttributeAsDouble(node_mappoint, "Efficiency") * FEffMultiplier;
		if(k == 0)
			FSpeedVec.push_back(speed);
		k++;
		if(j > 0) {
			if(speed != FSpeed[i][j - 1]) {
				i++;
				j = 0;
				FSpeed.resize(i + 1);
				FMass.resize(i + 1);
				FPres.resize(i + 1);
				FEff.resize(i + 1);
				FSpeedVec.push_back(speed);
			}
		}
		FSpeed[i].push_back(speed);
		FMass[i].push_back(mass);
		FPres[i].push_back(pres);
		FEff[i].push_back(eff);
		j++;

	}

	FNumLines = FSpeed.size();

}

void TSAEMap::ExtrapolateMap(double Dwheel) {

	int max_j;
	double max_p;

	FPres_ex = FPres;
	FSpeed_ex = FSpeed;
	FEff_ex = FEff;
	FMass_ex = FMass;

	//Look for the index of the maximum compression ratio
	for(int i = 0; i < FNumLines; i++) {
		max_j = 0;
		max_p = FPres[i][0];
		for(int j = 1; j < FPres[i].size(); j++) {
			if(FPres[i][j] > max_p) {
				max_j = j;
			}
		}
		FMaxIndex.push_back(max_j);
	}
	FPres_ns.resize(FNumLines);
	FSpeed_ns.resize(FNumLines);
	FEff_ns.resize(FNumLines);
	FMass_ns.resize(FNumLines);

	//Take only the map where the slope of the iso-speed lines is negative.
	for(int i = 0; i < FNumLines; i++) {
		for(int j = FMaxIndex[i]; j < FPres[i].size(); j++) {
			FPres_ns[i].push_back(FPres[i][j]);
			FSpeed_ns[i].push_back(FSpeed[i][j]);
			FEff_ns[i].push_back(FEff[i][j]);
			FMass_ns[i].push_back(FMass[i][j]);
		}
	}

	//Fit elipse for each iso-speed line.
	FitElipse();

	FitLeufvenModel();

	ExtrapolateLPR();

	ExtrapolateLSpeed(Dwheel);

	ExtrapolateHSpeed();

	MaximumEfficiencyCurve();

	ExtrapolateEfficiency();

	BuidExtrapolatedMap();

	AdimensionalizeMap();

}

void TSAEMap::FitElipse() {

	stCorrelation *Elipse;
	Elipse = new stElipse();

	dVector x;
	dVector y;
	dVector sig;
	dVector a;
	a.resize(3, 1);

	Fitmrq *MSR;
	MSR = new Fitmrq(x, y, sig, a, Elipse, 1e-5);

	// Maximum mass flow for each iso-speed line
	//FWmax.resize(FNumLines);
	// Mass flow for maximum compressor ratio
	//FWzs1.resize(FNumLines);
	// Maximum compressor ratio
	//FRCzs1.resize(FNumLines);

	// Fit elipse for each iso-speed line
	cout << "FITTING ELLIPSES" << endl;
	cout << "Line\tWmax\tWzsl\tRCzsl" << endl;

	for(int i = 0; i < FNumLines; i++) {
		if(FPres_ns[i].size() > 3) {
			x.resize(FPres_ns[i].size());
			y.resize(FPres_ns[i].size());
			sig.resize(FPres_ns[i].size(), 1);
			for(int j = 0; j < FPres_ns[i].size(); j++) {
				x[j] = FMass_ns[i][j];
				y[j] = FPres_ns[i][j];
			}
			a[0] = FMass_ns[i].front();
			MSR->max_limit(0, FMass_ns[i].front());
			MSR->min_limit(0, 0.002);
			a[1] = FMass_ns[i].back();
			MSR->min_limit(1, FMass_ns[i][FMass_ns[i].size() - 2]);
			a[2] = FPres_ns[i].front();
			MSR->min_limit(2, FPres_ns[i][1]);

			MSR->AsignData(x, y, sig, a, Elipse, 1e-5);

			MSR->fit();

			cout << i << "\t" << MSR->a[1] << "\t" << MSR->a[0] << "\t" << MSR->a[2] << endl;

			FSpeedEli.push_back(FSpeed[i][0]);
			FWmax.push_back(MSR->a[1]);
			FWzs1.push_back(MSR->a[0]);
			FRCzs1.push_back(MSR->a[2]);
			FWrate.push_back(MSR->a[0] / MSR->a[1]);
		}
	}

}

void TSAEMap::FitLeufvenModel() {

	dVector aa;
	dVector sig;
	sig.resize(FSpeedEli.size(), 1);

	stCorrelation *Wmax;
	Wmax = new stWmax();

	Fitmrq *fitWmax;
	aa.resize(2);
	aa[0] = 1.;
	aa[1] = 1.;
	fitWmax = new Fitmrq(FSpeedEli, FWmax, sig, aa, Wmax, 1e-5);
	fitWmax->fit();
	Wmax->p = fitWmax->a;
	cout << "FITTING WMAX" << endl;
	cout << Wmax->p[0] << "\t" << Wmax->p[1] << endl;

	stCorrelation *Wzs1;
	Wzs1 = new stWzs1();

	Fitmrq *fitWzs1;
	aa.resize(2);
	aa[0] = 1.;
	aa[1] = 1.;
	fitWzs1 = new Fitmrq(FSpeedEli, FWzs1, sig, aa, Wzs1, 1e-5);
	fitWzs1->fit();
	Wzs1->p = fitWzs1->a;
	cout << "FITTING WZSL" << endl;
	cout << Wzs1->p[0] << "\t" << Wzs1->p[1] << endl;

	stCorrelation *Wrate;
	Wrate = new stFparm();

	Fitmrq *fitWrate;
	aa.resize(3);
	aa[0] = 0.6;
	aa[1] = FSpeedEli[int(FSpeedEli.size() / 2)] * 1e-5;
	aa[2] = 2.0;
	fitWrate = new Fitmrq(FSpeedEli, FWrate, sig, aa, Wrate, 1e-3);
	fitWrate->min_limit(2, 0.0);
	fitWrate->min_limit(1, FSpeedEli.front() * 1e-5);
	fitWrate->fit();
	Wrate->p = fitWrate->a;
	cout << "FITTING WRATE" << endl;
	cout << Wrate->p[0] << "\t" << Wrate->p[1] << "\t" << Wrate->p[2] << endl;

	stCorrelation *RCzs1;
	RCzs1 = new stRCzs1();

	Fitmrq *fitRCzs1;
	aa.resize(2);
	aa[0] = 1.;
	aa[1] = 1.;
	fitRCzs1 = new Fitmrq(FSpeedEli, FRCzs1, sig, aa, RCzs1, 1e-5);
	fitRCzs1->fit();
	RCzs1->p = fitRCzs1->a;

	cout << "FITTING RCZSL" << endl;
	cout << RCzs1->p[0] << "\t" << RCzs1->p[1] << endl;

	//FLeufven = new stLeufven();

	//FitmrqMult *fitLeufven;
	//aa.resize(11);
	//aa[0] = 2.0;
	//aa[1] = 0.0;
	//aa[2] = 2.0;
	//aa[3] = 1e-5;
	//aa[4] = 1.0;
	//aa[5] = fitWmax->a[0];
	//aa[6] = fitWmax->a[1];
	//aa[7] = fitWzs1->a[0];
	//aa[8] = fitWzs1->a[1];
	//aa[9] = fitRCzs1->a[0];
	//aa[10] = fitRCzs1->a[1];

	//dMatrix x;
	//x.resize(0);
	//dVector xx;
	//dVector y;
	//xx.resize(2);
	//sig.resize(0);

	//for(int i = 0;i < FNumLines; i++){
	//	for(int j = 0;j < FSpeed_ns[i].size(); j++){
	//		xx[0] = FSpeed_ns[i][j];
	//		xx[1] = FMass_ns[i][j];
	//		x.push_back(xx);
	//		y.push_back(FPres_ns[i][j]);
	//		sig.push_back(1 / FPres_ns[i][j]);
	//	}
	//}
	////sig.resize(y.size(),1);

	//fitLeufven = new FitmrqMult(x, y, sig, aa, FLeufven, 1e-8);

	//fitLeufven->min_limit(0, 1.);
	//fitLeufven->min_limit(1, 0.);
	//fitLeufven->min_limit(2, 1.);
	//fitLeufven->min_limit(3, 1.e-5);
	//fitLeufven->min_limit(4, 0.);
	//fitLeufven->min_limit(5, 0.);
	//fitLeufven->min_limit(6, 0.);
	//fitLeufven->min_limit(7, 1.e-5);
	//fitLeufven->min_limit(8, 0.);
	//fitLeufven->min_limit(9, 1.e-5);
	//fitLeufven->min_limit(10, 0.);

	//fitLeufven->fit();

	//FLeufven->p = fitLeufven->a;

	//cout << "FITTING LEUFVEN" << endl;
	//cout << "C10:\t" << FLeufven->p[0] << endl;
	//cout << "C11:\t" << FLeufven->p[1] << endl;
	//cout << "C20:\t" << FLeufven->p[2] << endl;
	//cout << "C21:\t" << FLeufven->p[3] << endl;
	//cout << "C22:\t" << FLeufven->p[4] << endl;
	//cout << "C30:\t" << FLeufven->p[5] << endl;
	//cout << "C31:\t" << FLeufven->p[6] << endl;
	//cout << "C41:\t" << FLeufven->p[7] << endl;
	//cout << "C42:\t" << FLeufven->p[8] << endl;
	//cout << "C51:\t" << FLeufven->p[9] << endl;
	//cout << "C52:\t" << FLeufven->p[10] << endl;

	FLeufven = new stLeufvenMod();

	FitmrqMult *fitLeufven;
	aa.resize(12);
	aa[0] = 2.0;
	aa[1] = 0.0;
	aa[2] = 2.0;
	aa[3] = 1e-5;
	aa[4] = 1.0;
	aa[5] = fitWmax->a[0];
	aa[6] = fitWmax->a[1];
	aa[7] = fitWrate->a[0];
	aa[8] = fitWrate->a[1];
	aa[9] = fitWrate->a[2];
	aa[10] = fitRCzs1->a[0];
	aa[11] = fitRCzs1->a[1];

	dMatrix x;
	x.resize(0);
	dVector xx;
	dVector y;
	xx.resize(2);
	sig.resize(0);

	for(int i = 0; i < FNumLines; i++) {
		for(int j = 0; j < FSpeed_ns[i].size(); j++) {
			xx[0] = FSpeed_ns[i][j];
			xx[1] = FMass_ns[i][j];
			x.push_back(xx);
			y.push_back(FPres_ns[i][j]);
			sig.push_back(1 / FPres_ns[i][j]);
		}
	}
	//sig.resize(y.size(),1);

	fitLeufven = new FitmrqMult(x, y, sig, aa, FLeufven, 1e-5);

	fitLeufven->min_limit(0, 1.);
	fitLeufven->min_limit(1, 0.);
	fitLeufven->min_limit(2, 1.);
	fitLeufven->min_limit(3, 1.e-5);
	fitLeufven->min_limit(4, 0.);
	fitLeufven->min_limit(5, 0.);
	fitLeufven->min_limit(6, 0.);
	fitLeufven->min_limit(7, 1.e-5);
	fitLeufven->min_limit(8, FSpeedEli[0] * 1e-5);
	fitLeufven->min_limit(9, 1.);
	fitLeufven->min_limit(10, 1.e-5);
	fitLeufven->min_limit(11, 0.);

	fitLeufven->fit();

	FLeufven->p = fitLeufven->a;

	cout << "FITTING LEUFVEN" << endl;
	cout << "C10:\t" << FLeufven->p[0] << endl;
	cout << "C11:\t" << FLeufven->p[1] << endl;
	cout << "C20:\t" << FLeufven->p[2] << endl;
	cout << "C21:\t" << FLeufven->p[3] << endl;
	cout << "C22:\t" << FLeufven->p[4] << endl;
	cout << "C30:\t" << FLeufven->p[5] << endl;
	cout << "C31:\t" << FLeufven->p[6] << endl;
	cout << "C41:\t" << FLeufven->p[7] << endl;
	cout << "C42:\t" << FLeufven->p[8] << endl;
	cout << "C51:\t" << FLeufven->p[9] << endl;
	cout << "C52:\t" << FLeufven->p[10] << endl;

}

void TSAEMap::ExtrapolateLPR() {

	FPres_LPR.resize(FNumLines);
	FSpeed_LPR.resize(FNumLines);
	FEff_LPR.resize(FNumLines);
	FMass_LPR.resize(FNumLines);

	double DeltaPmax = 0.05;
	double DeltaP = DeltaPmax;
	int npoints;

	Hermite_interp *fun_Pre;

	dVector Inputs;
	Inputs.resize(2);
	double mass = 0.;
	double s = 0.;
	double drc = 0.;
	double rclow = 0.;

	dVector M;
	dVector RC;
	M.resize(3);
	RC.resize(3);
	std::string Option = "PlusOffset";

	for(int i = 0; i < FNumLines; i++) {

		drc = 0.;
		if(FPres_ns[i].back() > 1) {

			DeltaP = Min((FPres_ns[i].front() - FPres_ns[i].back()) / FPres_ns[i].size(), DeltaPmax);
			npoints = floor((FPres_ns[i].back() - 1) / DeltaP);
			// Extrapolation at low pressure ratio using Leufven
			Inputs[0] = FSpeed_ns[i][0];
			Inputs[1] = 1.;
			mass = dynamic_cast<stLeufvenMod *>(FLeufven)->val_m(Inputs);
			if(Option == "Hermite" || Option == "Hermite+Leufven" || Option == "Leufven") {
				if(mass <= FMass[i].back()) {
					Inputs[1] = FMass[i].back();
					double rc = dynamic_cast<stLeufvenMod *>(FLeufven)->val_rc(Inputs);
					drc = FPres[i].back() - rc;
				}
			} else {
				Inputs[1] = FMass[i].back();
				double rc = dynamic_cast<stLeufvenMod *>(FLeufven)->val_rc(Inputs);
				if(Option == "PlusOffset") {
					drc = FPres[i].back() - rc;
				} else {
					drc = (FPres[i].back() - 1) / (rc - 1);
				}
			}
			for(int j = 0; j < npoints; j++) {
				rclow = j * DeltaP + 1.;
				if(Option == "MultOffset") {
					Inputs[1] = rclow;
				} else {
					Inputs[1] = rclow - drc;
				}
				mass = dynamic_cast<stLeufvenMod *>(FLeufven)->val_m(Inputs);
				if(mass > FMass[i].back()) {
					FSpeed_LPR[i].push_back(FSpeed_ns[i][0]);
					if(Option == "MultOffset") {
						FPres_LPR[i].insert(FPres_LPR[i].begin(), (rclow - 1) * drc + 1);
					} else {
						FPres_LPR[i].insert(FPres_LPR[i].begin(), rclow);
					}
					FMass_LPR[i].insert(FMass_LPR[i].begin(), mass);
					FEff_LPR[i].push_back(0.5);
				} else {
					npoints = j;
					break;
				}
			}

			M[0] = FMass_ns[i].front();
			RC[0] = FPres_ns[i].front();
			M[1] = FMass_ns[i].back();
			RC[1] = FPres_ns[i].back();
			M[2] = FMass_LPR[i].back();
			RC[2] = 1;

			fun_Pre = new Hermite_interp(M, RC);
			// Correction at low pressure ratio to avoid discontinuities
			if(Option == "Hermite") {
				s = 0;
			} else {
				s = 1;
			}
			for(int j = 0; j < npoints; j++) {
				if(Option == "Hermite+Leufven") {
					s = (FPres_LPR[i][j] - FPres_ns[i].back()) / (1 - FPres_ns[i].back());
				}
				FPres_LPR[i][j] = (1 - s) * fun_Pre->interp(FMass_LPR[i][j]) + s * FPres_LPR[i][j];
			}

			FSpeed_ex[i].insert(FSpeed_ex[i].end(), FSpeed_LPR[i].begin(), FSpeed_LPR[i].end());
			FMass_ex[i].insert(FMass_ex[i].end(), FMass_LPR[i].begin(), FMass_LPR[i].end());
			FPres_ex[i].insert(FPres_ex[i].end(), FPres_LPR[i].begin(), FPres_LPR[i].end());
			FEff_ex[i].insert(FEff_ex[i].end(), FEff_LPR[i].begin(), FEff_LPR[i].end());
		}
	}
}

void TSAEMap::ExtrapolateLSpeed(double Dwheel) {

	dVector N_Surge;
	dVector RC_Surge;
	dVector sig;
	dVector aa;

	N_Surge.resize(FNumLines);

	for(int i = 0; i < FNumLines; i++) {
		N_Surge[i] = FSpeed_ex[i][0];
		RC_Surge.push_back(FPres_ex[i].front());
	}

	stCorrelation *Surge;
	FSurge = new stRCzs1();

	Fitmrq *fitSurge;
	aa.resize(2);
	aa[0] = 1.;
	aa[1] = 1.;

	sig.resize(RC_Surge.size(), 1);
	fitSurge = new Fitmrq(N_Surge, RC_Surge, sig, aa, FSurge, 1e-5);
	fitSurge->fit();

	FSurge->p = fitSurge->a;

	int MaxLowSpeedLines = 3;
	double deltaN = FSpeed[1][0] - FSpeed[0][0];
	MaxLowSpeedLines = Min(Max((int) floor(FSpeed[0][0] / deltaN), 0), MaxLowSpeedLines);

	dVector LS_ene;

	int newmaxn = MaxLowSpeedLines;
	double n;

	for(int i = 0; i < MaxLowSpeedLines; i++) {
		n = FSpeed_ex[0][0] - (MaxLowSpeedLines - i) * deltaN;
		if(n < 5000) {
			newmaxn = MaxLowSpeedLines - i;
		}
	}
	MaxLowSpeedLines = newmaxn;
	for(int i = 0; i < MaxLowSpeedLines; i++) {
		n = FSpeed_ex[0][0] - (MaxLowSpeedLines - i) * deltaN;
		if(n < 5000) {
			LS_ene.push_back(5000);
		} else {
			LS_ene.push_back(n);
		}
	}

	if(MaxLowSpeedLines > 0) {
		dVector Psi;
		dVector Phi;

		double Rho1 = FPresionRef / (__R::Air * FTempRef);
		double uc = Dwheel / 2 * __units::RPMToRad_s(FSpeed[0][0]);
		double cp = __Gamma::Cp;
		for(int j = 0; j < FSpeed_ex[0].size(); j++) {
			Phi.push_back(FMass_ex[0][j] / (Rho1 * __geom::Circle_area(Dwheel) * uc));
			Psi.push_back(cp * FTempRef * (pow(FPres_ex[0][j], __Gamma::G_8) - 1) * 2 / pow2(uc));
		}

//		stCorrelation *Martin;
//		Martin = new stMartin();
//
//		Fitmrq *fitMartin;
//
//		aa.resize(3);
//		aa[0] = 0.5;
//		aa[1] = -1.0;
//		aa[2] = 0.5;
//
//
//		sig.resize(FSpeed_ex[0].size(), 1.);
//		fitMartin = new Fitmrq(Phi, Psi, sig, aa, Martin, 1e-5);
//
//		fitMartin->fit();

		FSpeed_LS.resize(MaxLowSpeedLines);
		FMass_LS.resize(MaxLowSpeedLines);
		FPres_LS.resize(MaxLowSpeedLines);
		FEff_LS.resize(MaxLowSpeedLines);

		double deltaM;
		double phi_LS;
		double psi_LS;

		for(int i = 0; i < MaxLowSpeedLines; i++) {

			FSpeed_LS[i].resize(FMass_ex[0].size(), LS_ene[i]);
			FMass_LS[i].resize(FMass_ex[0].size(), 0);
			FPres_LS[i].resize(FMass_ex[0].size(), 1);
			FEff_LS[i].resize(FMass_ex[0].size(), 0.5);

			FPres_LS[i][0] = dynamic_cast<stRCzs1 *>(FSurge)->val(LS_ene[i]);
			FMass_LS[i][0] = FMass[0][0] + (FMass[0][0] - 0) / (FPres[0][0] - 1) * (FPres_LS[i][0] - FPres[0][0]);
			uc = Dwheel / 2 * __units::RPMToRad_s(LS_ene[i]);
//			psi_LS = 0.;
//			phi_LS = dynamic_cast<stMartin *>(Martin)->val_phi(psi_LS,
//					fitMartin->a);
			FMass_LS[i].back() = Phi.back() * (Rho1 * __geom::Circle_area(Dwheel) * uc);
			deltaM = (FMass_LS[i].back() - FMass_LS[i].front()) / (FMass_LS[0].size() - 1);
			for(int j = 0; j < FMass_LS[i].size() - 1; j++) {
				FMass_LS[i][j] = Phi[j] * (Rho1 * __geom::Circle_area(Dwheel) * uc);
//				FMass_LS[i][j] = FMass_LS[i][0] + deltaM * j;
//				phi_LS = FMass_LS[i][j] / (Rho1 * Pi * pow2(Dwheel) * uc / 4.);
//				psi_LS = dynamic_cast<stMartin *>(Martin)->val_psi(phi_LS,
//						fitMartin->a);
				FPres_LS[i][j] = pow(Psi[j] * pow2(uc) / 2 / cp / FTempRef + 1, __Gamma::G_9);
			}
		}

		FSpeed_ex.insert(FSpeed_ex.begin(), FSpeed_LS.begin(), FSpeed_LS.end());
		FMass_ex.insert(FMass_ex.begin(), FMass_LS.begin(), FMass_LS.end());
		FPres_ex.insert(FPres_ex.begin(), FPres_LS.begin(), FPres_LS.end());
	}

}

void TSAEMap::ExtrapolateHSpeed() {

	int nHS = 1; // Number of high speed lines extrapolated
	int npoints = 40;

	double deltaN = FSpeed[FSpeed.size() - 1][0] - FSpeed[FSpeed.size() - 2][0];
	double RCzs1 = 1.;
	double deltaRC = 0.;
	double M_surge;
	double RC_surge;

	dVector x;
	FSpeed_HS.resize(nHS);
	FPres_HS.resize(nHS);
	FMass_HS.resize(nHS);
	FEff_HS.resize(nHS);

	x.resize(2);

	for(int i = 0; i < nHS; i++) {
		FSpeed_HS[i].resize(npoints, FSpeed_ex[FSpeed_ex.size() - 1][0] + deltaN * (i + 1));
		RCzs1 = 1 + FLeufven->p[10] * pow(FSpeed_HS[i][0] * 1e-5, FLeufven->p[11]);
		deltaRC = (RCzs1 - 1) / (npoints - 1);
		x[0] = FSpeed_HS[i][0];
		for(int j = 0; j < npoints; j++) {
			FPres_HS[i].push_back(RCzs1 - deltaRC * j);
			x[1] = FPres_HS[i][j];
			FMass_HS[i].push_back(dynamic_cast<stLeufvenMod *>(FLeufven)->val_m(x));
		}
		FEff_HS[i].resize(npoints, 0.5);

		RC_surge = dynamic_cast<stRCzs1 *>(FSurge)->val(x[0]);
		M_surge = FMass_ex[FMass_ex.size() - 1][0] + (FMass_ex[FMass_ex.size() - 1][0] - 0) /
				  (FPres_ex[FPres_ex.size() - 1][0] - 1) * (RC_surge - FPres_ex[FPres_ex.size() - 1][0]);
		if(M_surge < FMass_HS[i][0]) {
			double A = -(-RC_surge + FPres_HS[i][0]) / pow2(M_surge - FMass_HS[i][0]);
			double B = -2 * FMass_HS[i][0] * (RC_surge - FPres_HS[i][0]) / pow2(M_surge - FMass_HS[i][0]);
			double C = FPres_HS[i][0] - pow2(FMass_HS[i][0]) * (-RC_surge + FPres_HS[i][0]) / pow2(M_surge - FMass_HS[i][0]);
			int np = (int) ceil((M_surge - FMass_HS[i][0]) / (FMass_HS[i][0] - FMass_HS[i][1]));
			double deltaM = -(M_surge - FMass_HS[i][0]) / np;
			double M;
			double P;

			for(int j = np - 1; j >= 0; j--) {
				M = M_surge + j * deltaM;
				FMass_HS[i].insert(FMass_HS[i].begin(), M);
				P = A * pow2(M) + B * M + C;
				FPres_HS[i].insert(FPres_HS[i].begin(), P);
				FSpeed_HS[i].insert(FSpeed_HS[i].begin(), x[0]);
				FEff_HS[i].insert(FEff_HS[i].begin(), 0.5);
			}

		} else {
			int np = FMass_HS[i].size();
			for(int j = 0; j < np; j++) {
				RC_surge = (FMass_HS[i].front() - FMass_ex[FMass_ex.size() - 1][0]) * (FPres_ex[FPres_ex.size() - 1][0] - 1) /
						   (FMass_ex[FMass_ex.size() - 1][0] - 0)
						   + FPres_ex[FPres_ex.size() - 1][0];
				if(RC_surge < FPres_HS[i].front()) {
					FMass_HS[i].erase(FMass_HS[i].begin());
					FSpeed_HS[i].erase(FSpeed_HS[i].begin());
					FPres_HS[i].erase(FPres_HS[i].begin());
					FEff_HS[i].erase(FEff_HS[i].begin());
				}

			}
			//while (M_surge > FMass_HS[i][0]) {
			//	FMass_HS[i].erase(FMass_HS[i].begin());
			//	FSpeed_HS[i].erase(FSpeed_HS[i].begin());
			//	FPres_HS[i].erase(FPres_HS[i].begin());
			//	FEff_HS[i].erase(FEff_HS[i].begin());
			//}
			//FMass_HS[i].insert(FMass_HS[i].begin(), M_surge);
			//FPres_HS[i].insert(FPres_HS[i].begin(), RC_surge);
			//FSpeed_HS[i].insert(FSpeed_HS[i].begin(), x[0]);
			//FEff_HS[i].insert(FEff_HS[i].begin(), 0.5);
		}
	}
	FMass_ex.insert(FMass_ex.end(), FMass_HS.begin(), FMass_HS.end());
	FSpeed_ex.insert(FSpeed_ex.end(), FSpeed_HS.begin(), FSpeed_HS.end());
	FPres_ex.insert(FPres_ex.end(), FPres_HS.begin(), FPres_HS.end());
	FEff_ex.insert(FEff_ex.end(), FEff_HS.begin(), FEff_HS.end());

}

void TSAEMap::MaximumEfficiencyCurve() {

	dVector SortedMass;
	dVector SortedEff;
	dVector SpeedMax;
	vector<size_t> SortedInd;

	dVector a;
	a.resize(3, 1);

	for(int i = 0; i < FNumLines; i++) {
		SortedInd = sort_indexes_down(FEff[i]);
		for(int j = 0; j < 3; j++) {
			SortedEff.push_back(FEff[i][SortedInd[j]]);
			SortedMass.push_back(FMass[i][SortedInd[j]]);
		}

		if(SortedInd[0] < SortedInd[1] && SortedInd[0] < SortedInd[2]) {
			if(SortedInd[0] > 0) {
				SortedEff[2] = FEff[i][SortedInd[0] - 1];
				SortedMass[2] = FMass[i][SortedInd[0] - 1];
			}
		} else if(SortedInd[0] > SortedInd[1] && SortedInd[0] > SortedInd[2]) {
			if(SortedInd[0] < FEff[i].size() - 1) {
				SortedEff[2] = FEff[i][SortedInd[0] + 1];
				SortedMass[2] = FMass[i][SortedInd[0] + 1];
			}
		}

		Poly2(SortedMass, SortedEff, a);

		FMassEffMax.push_back(-a[1] / 2 / a[0]);
		double max = MaxComponent(SortedMass);
		if(FMassEffMax[i] > MaxComponent(SortedMass) || FMassEffMax[i] < MinComponent(SortedMass)) {
			FMassEffMax[i] = FMass[i][SortedInd[0]];
			FEffMax.push_back(FEff[i][SortedInd[0]]);
		} else {
			FEffMax.push_back(-pow2(a[1]) / 4 / a[0] + a[2]);
		}
		SpeedMax.push_back(FSpeed[i][0] * 1e-5);

		SortedInd.clear();
		SortedMass.clear();
		SortedEff.clear();
	}

	FMaxEff = new stParabola();

	Fitmrq *fitMaxEff;
	dVector aa, sig;
	aa.resize(3);
	aa[0] = -0.1;
	aa[1] = 0.1;
	aa[2] = 0.6;
	sig.resize(SpeedMax.size(), 1.);
	fitMaxEff = new Fitmrq(SpeedMax, FEffMax, sig, aa, FMaxEff, 1e-10);
	fitMaxEff->fit();
	FMaxEff->p = fitMaxEff->a;

}

void TSAEMap::ExtrapolateEfficiency() {

	double kEff;
	dVector X_max;
	dMatrix Mass_norm;
	dMatrix Eff_norm;
	dVector m, n;
	Mass_norm.resize(FEff.size());
	Eff_norm.resize(FEff.size());
	dVector a, sig;
	a.resize(2, 1.);
	stCorrelation *EffCurve;
	EffCurve = new stEffCurve();
	double ms;

	Fitmrq *fitEffCurve;
	double SpeedRatio = FSpeed_HS.back()[0] / FSpeed.back()[0];
	double mmax = log(0.5 * (pow2(SpeedRatio) - 1) / pow2(SpeedRatio)) / log(0.5) - 0.05;

	//LOW PRESSURE RATIO EXTRAPOLATION
	for(int i = 0; i < FNumLines; i++) {
		if(FMass_LPR[i].size() > 0) {
			X_max.push_back(FMass_LPR[i].back());
			for(int j = 0; j < FMass[i].size(); j++) {
				Mass_norm[i].push_back(FMass[i][j] / X_max[i]);
				Eff_norm[i].push_back(FEff[i][j] / FEffMax[i]);
			}
		} else {
			X_max.push_back(FMass[i].back());
			for(int j = 0; j < FMass[i].size() - 1; j++) {
				Mass_norm[i].push_back(FMass[i][j] / X_max[i]);
				Eff_norm[i].push_back(FEff[i][j] / FEffMax[i]);
			}
		}

		a[0] = log(1 - FMassEffMax[i] / X_max[i]) / log(0.5);
		a[1] = 2;
		sig.resize(Mass_norm[i].size(), 1);
		fitEffCurve = new Fitmrq(Mass_norm[i], Eff_norm[i], sig, a, EffCurve, 1e-10);
		fitEffCurve->max_limit(0, mmax);
		fitEffCurve->fit();
		m.push_back(fitEffCurve->a[0]);
		n.push_back(fitEffCurve->a[1]);
		ms = FMass[i].back() / X_max[i];
		kEff = dynamic_cast<stEffCurve *>(EffCurve)->val(ms, fitEffCurve->a) * FEffMax[i] / FEff[i].back();
		for(int j = 0; j < FMass_LPR[i].size(); j++) {
			ms = FMass_LPR[i][j] / X_max[i];
			FEff_LPR[i][j] = dynamic_cast<stEffCurve *>(EffCurve)->val(ms, fitEffCurve->a) / kEff * FEffMax[i];
			if(FEff_LPR[i][j] < 0.001)
				FEff_LPR[i][j] = 0.001;
		}
		delete fitEffCurve;
	}
	// LOW SPEED EXTRAPOLATION

	double MaxEff;
	double M_norm;
	for(int i = 0; i < FSpeed_LS.size(); i++) {
		SpeedRatio = FSpeed_LS[i][0] / FSpeed[0][0];
		a[0] = log(1 - (((1 - pow(0.5, m[0])) - .5) * pow2(SpeedRatio) + .5)) / log(0.5);
		if(a[0] > 4.5)
			a[0] = 4.5;
		a[1] = n[0];
		MaxEff = dynamic_cast<stParabola *>(FMaxEff)->val(FSpeed_LS[i][0] * 1e-5);
		for(int j = 0; j < FMass_LS[i].size(); j++) {
			M_norm = FMass_LS[i][j] / FMass_LS[i].back();
			FEff_LS[i][j] = dynamic_cast<stEffCurve *>(EffCurve)->val(M_norm, a) * MaxEff;
			if(FEff_LS[i][j] < 0.001)
				FEff_LS[i][j] = 0.001;
		}
	}
	// HIGH SPEED EXTRAPOLATION
	for(int i = 0; i < FSpeed_HS.size(); i++) {
		SpeedRatio = FSpeed_HS[i][0] / FSpeed.back()[0];
		a[0] = log(1 - (((1 - pow(0.5, m.back())) - .5) * pow2(SpeedRatio) + .5)) / log(0.5);
		if(a[0] > 4.5)
			a[0] = 4.5;
		a[1] = n.back();
		MaxEff = dynamic_cast<stParabola *>(FMaxEff)->val(FSpeed_HS[i][0] * 1e-5);
		for(int j = 0; j < FMass_HS[i].size(); j++) {
			M_norm = FMass_HS[i][j] / FMass_HS[i].back();
			FEff_HS[i][j] = dynamic_cast<stEffCurve *>(EffCurve)->val(M_norm, a) * MaxEff;
			if(FEff_HS[i][j] < 0.001)
				FEff_HS[i][j] = 0.001;
		}
	}

}

void TSAEMap::BuidExtrapolatedMap() {

	for(int i = 0; i < FNumLines; i++) {

		FSpeed[i].insert(FSpeed[i].end(), FSpeed_LPR[i].begin(), FSpeed_LPR[i].end());
		FMass[i].insert(FMass[i].end(), FMass_LPR[i].begin(), FMass_LPR[i].end());
		FPres[i].insert(FPres[i].end(), FPres_LPR[i].begin(), FPres_LPR[i].end());
		FEff[i].insert(FEff[i].end(), FEff_LPR[i].begin(), FEff_LPR[i].end());
	}

	FSpeed.insert(FSpeed.begin(), FSpeed_LS.begin(), FSpeed_LS.end());
	FMass.insert(FMass.begin(), FMass_LS.begin(), FMass_LS.end());
	FPres.insert(FPres.begin(), FPres_LS.begin(), FPres_LS.end());
	FEff.insert(FEff.begin(), FEff_LS.begin(), FEff_LS.end());

	FMass.insert(FMass.end(), FMass_HS.begin(), FMass_HS.end());
	FSpeed.insert(FSpeed.end(), FSpeed_HS.begin(), FSpeed_HS.end());
	FPres.insert(FPres.end(), FPres_HS.begin(), FPres_HS.end());
	FEff.insert(FEff.end(), FEff_HS.begin(), FEff_HS.end());

	FNumLines = FSpeed.size();

	FSpeedVec.clear();
	for(int i = 0; i < FNumLines; i++) {
		FSpeedVec.push_back(FSpeed[i][0]);
	}

	FILE *fich = fopen("CompressorExtrp.txt", "w");
	for(int i = 0; i < FNumLines; i++) {
		for(int j = 0; j < FSpeed[i].size(); j++) {
			fprintf(fich, "%lf\t%lf\t%lf\t%lf\n", FSpeed[i][j], FMass[i][j], FPres[i][j], FEff[i][j]);
		}
		fprintf(fich, "\n");
	}

	fclose(fich);
}

void TSAEMap::AdimensionalizeMap() {

	double tmp = 0.;

	FMassMAX.resize(FNumLines);
	FPresMAX.resize(FNumLines);
	FEffMAX.resize(FNumLines);
	FSpeedMAX = FSpeed.back()[0];

	double pmax = 0., emax = 0.;
	FMassMAXMAX = 0;
	FPresMAXMAX = 1;
	FEffMAXMAX = 0;

	for(int i = 0; i < FNumLines; i++) {
		FMassMAX[i] = FMass[i].back();
		for(int j = 0; j < FPres[i].size(); j++) {
			if(FPres[i][j] > pmax)
				pmax = FPres[i][j];
			if(FEff[i][j] > emax)
				emax = FEff[i][j];
		}
		FPresMAX[i] = pmax;
		FEffMAX[i] = emax;

		pmax = 1.;
		emax = 0.;

		if(FMassMAX[i] > FMassMAXMAX)
			FMassMAXMAX = FMassMAX[i];
		if(FPresMAX[i] > FPresMAXMAX)
			FPresMAXMAX = FPresMAX[i];
		if(FEffMAX[i] > FEffMAXMAX)
			FEffMAXMAX = FEffMAX[i];
	}

	// FSpeedAdim.resize(FNumLines);
	FMassAdim.resize(FNumLines);
	FPresAdim.resize(FNumLines);
	FEffAdim.resize(FNumLines);

	FPre_MassCurve.resize(FNumLines);
	FEff_MassCurve.resize(FNumLines);

	for(int i = 0; i < FNumLines; i++) {
		FSpeedAdim.push_back(FSpeedVec[i] / FSpeedMAX);
		// printf("%lf %lf %lf\n",FSpeedVec[i] ,FSpeedAdim[i],FSpeedVec[i] / FSpeedMAX);
		FMassMAXAdim.push_back(FMassMAX[i] / FMassMAXMAX);
		FPresMAXAdim.push_back((FPresMAX[i] - 1) / (FPresMAXMAX - 1));
		FEffMAXAdim.push_back(FEffMAX[i] / FEffMAXMAX);
		for(unsigned int j = 0; j < FSpeed[i].size(); j++) {

			FMassAdim[i].push_back(FMass[i][j] / FMassMAX[i]);
			FPresAdim[i].push_back((FPres[i][j] - 1) / (FPresMAX[i] - 1));
			FEffAdim[i].push_back(FEff[i][j] / FEffMAX[i]);
		}
		FMassAdim[i].insert(FMassAdim[i].begin(), 0);
		FPresAdim[i].insert(FPresAdim[i].begin(), FPresAdim[i][0]);
		FEffAdim[i].insert(FEffAdim[i].begin(), FEffAdim[i][0]);

		FPre_MassCurve[i] = new Hermite_interp(FMassAdim[i], FPresAdim[i]);
		FEff_MassCurve[i] = new Hermite_interp(FMassAdim[i], FEffAdim[i]);
	}

	FSpeedAdim.insert(FSpeedAdim.begin(), 0.);
	FMassMAXAdim.insert(FMassMAXAdim.begin(), 0.);
	FPresMAXAdim.insert(FPresMAXAdim.begin(), 0.);
	FEffMAXAdim.insert(FEffMAXAdim.begin(), 0.);

	// printf("%d\n", FPresMAXAdim.size());

	FMassMAX_int = new Hermite_interp(FSpeedAdim, FMassMAXAdim);
	FPresMAX_int = new Hermite_interp(FSpeedAdim, FPresMAXAdim);
	FEffMAX_int = new Hermite_interp(FSpeedAdim, FEffMAXAdim);

}

void TSAEMap::InterpolateMAP(double RTC) {

	FRTCAdim = RTC / FSpeedMAX;
	FCurrentIND = FMassMAX_int->locate(FRTCAdim);

	// std::cout << FMassMAXAdim[4] << std::endl;

	if(FRTCAdim <= 1) {
		FCurrentIND = FMassMAX_int->locate(FRTCAdim);
		FCurrentMassMAX = FMassMAX_int->interp(FRTCAdim) * FMassMAXMAX;
		FCurrentPresMAX = FPresMAX_int->interp(FRTCAdim) * (FPresMAXMAX - 1.) + 1.;
		FCurrentEffMAX = FEffMAX_int->interp(FRTCAdim) * FEffMAXMAX;
		FDeltaLow = (FRTCAdim - FSpeedAdim[FCurrentIND]) / (FSpeedAdim[FCurrentIND + 1] - FSpeedAdim[FCurrentIND]);

	} else {
		FCurrentIND = FNumLines;
		FCurrentMassMAX = pow(FRTCAdim, 1.2) * FMassMAXMAX;
		FCurrentPresMAX = pow(FRTCAdim, 1.2) * (FPresMAXMAX - 1.) + 1.;
		FCurrentEffMAX = exp(1 - FRTCAdim) * FEffMAX_int->interp(FRTCAdim) * FEffMAXMAX;
	}

	// for (int i = 0; i <= FNumLines; ++i) {
	// printf("%lf %lf\n", FSpeedAdim[i], FSpeedAdim[i] * FSpeedMAX);
	// }

}

double TSAEMap::GetCurrentRC(double Mass) {

	double massadim = Mass / FCurrentMassMAX;
	double CurrentRC = 0.;

	if(FCurrentIND == 0) {
		CurrentRC = (FDeltaLow * FPre_MassCurve[FCurrentIND]->interp(massadim)) * (FCurrentPresMAX - 1) + 1;
	} else if(FCurrentIND == FNumLines) {
		CurrentRC = FPre_MassCurve[FCurrentIND - 1]->interp(massadim) * (FCurrentPresMAX - 1) + 1;
	} else {
		double pres_lo = FPre_MassCurve[FCurrentIND - 1]->interp(massadim);
		double pres_hi = FPre_MassCurve[FCurrentIND]->interp(massadim);

		CurrentRC = ((1 - FDeltaLow) * pres_lo + FDeltaLow * pres_hi) * (FCurrentPresMAX - 1) + 1;
	}
	return CurrentRC;

}

double TSAEMap::GetCurrentEff(double Mass) {

	double massadim = Mass / FCurrentMassMAX;
	double CurrentEff = 0.;

	if(FCurrentIND == 0) {
		CurrentEff = (FDeltaLow * FEff_MassCurve[FCurrentIND]->interp(massadim)) * FCurrentEffMAX;
	} else if(FCurrentIND == FNumLines) {
		CurrentEff = FEff_MassCurve[FCurrentIND - 1]->interp(massadim) * FCurrentEffMAX;
	} else {
		double pres_lo = FEff_MassCurve[FCurrentIND - 1]->interp(massadim);
		double pres_hi = FEff_MassCurve[FCurrentIND]->interp(massadim);

		CurrentEff = ((1 - FDeltaLow) * pres_lo + FDeltaLow * pres_hi) * FCurrentEffMAX;
	}
	return CurrentEff;

}

void TSAEMap::LeeMapa(FILE *fich) {

	int Adiab = 0;

	fscanf(fich, "%lf %lf ", &FPresionRef, &FTempRef);
	FTempRef = __units::degCToK(FTempRef);
	FPresionRef = __units::BarToPa(FPresionRef);

	fscanf(fich, "%lf %lf %lf ", &FMassMultiplier, &FCRMultiplier, &FEffMultiplier);
	fscanf(fich, "%d ", &Adiab);
	if(Adiab == 0) {
		FIsAdiabatic = false;
		fscanf(fich, "%lf ", &FTempMeasure);
	}

	ReadSAECompressorMap(fich);
}

void TSAEMap::LeeMapaXML(xml_node node_compressor) {

	xml_node node_map = GetNodeChild(node_compressor, "Com_SAE_Map");
	FPresionRef = GetXMLPressure(node_map, "PressureRef");
	FTempRef = GetXMLTemperature(node_map, "TemperatureRef");

	//FTempRef = __units::degCToK(FTempRef);
	//FPresionRef = __units::BarToPa(FPresionRef);

	FIsAdiabatic = GetAttributeAsBool(node_map, "IsAdiabatic");
	if(!FIsAdiabatic) {
		FTempMeasure = GetXMLTemperature(node_map, "MapTemperature");
	}

	ReadSAECompressorMapXML(node_map);
}

double TSAEMap::EvaluaRCHermite(double mass) {
	return GetCurrentRC(mass);
}

double TSAEMap::EvaluaRendSplines(double mass) {
	return GetCurrentEff(mass);
}

double TSAEMap::getGastoRelComp1() {
	return FCurrentMassMAX;
}

double TSAEMap::getRelCompBombeo() {
	if(FCurrentIND == 0) {
		return (1 - FDeltaLow) + FDeltaLow * FPres[FCurrentIND][0];
	}
	if(FCurrentIND == FNumLines) {
		return (FPres[FCurrentIND - 1][0] - 1) * pow(FRTCAdim, 1.2) + 1;
	} else {
		return (1 - FDeltaLow) * FPres[FCurrentIND - 1][0] + FDeltaLow * FPres[FCurrentIND][0];
	}
}

double TSAEMap::getGastoBombeo() {
	if(FCurrentIND == FNumLines) {
		return FMass[FCurrentIND - 1][0] * pow(FRTCAdim, 1.2);
	} else {
		return (1 - FDeltaLow) * FMass[FCurrentIND - 1][0] + FDeltaLow * FMass[FCurrentIND][0];
	}
}

void TSAEMap::InterpolaMapa(double rtc, double T10) {
	double rtc_cor = rtc * Sqrt(FTempRef / T10);

	InterpolateMAP(rtc_cor);
}

double TSAEMap::getRegimenCorregido() {
	return FRTCAdim * FSpeedMAX;
}

void TSAEMap::PutReference(double pref, double tref) {
	FTempRef = tref;
	FPresionRef = pref;
}

//void TSAEMap::CalculateAdiabaticEfficiency(TTC_HTM *HTM, double TinT) {
//
//	double m = 0., eff = 0., Rtc = 0.;
//
//	if(!FIsAdiabatic) {
//
//		FILE *fres = fopen("CompRes.dat", "w");
//
//		fprintf(fres, "Turbocharger_Speed(rmp)\t");
//		fprintf(fres, "Compressor_Mass(kg/s)\t");
//		fprintf(fres, "Compressor_Ratio(-)\t");
//		fprintf(fres, "Adiabatic_efficiency(-)\t");
//		fprintf(fres, "Efficiency(-)\t");
//		fprintf(fres, "Turbine_Heat_Flow(W)\t");
//		fprintf(fres, "Compressor_Heat_Flow(W)\t");
//		fprintf(fres, "Oil_Heat_Flow(W)\t");
//		fprintf(fres, "Water_Heat_Flow(W)\t");
//		fprintf(fres, "Ambient_Heat_Flow(W)\t");
//		fprintf(fres, "Temperature_G(K)\t");
//		fprintf(fres, "Temperature_AK)\t");
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
//		for(int i = 0; i < FNumLines; i++) {
//			for(unsigned int j = 0; j < FSpeed[i].size(); j++) {
//				m = FMass[i][j] / __units::PaToBar(FPresionRef) / sqrt(FTempMeasure / FTempRef);
//				Rtc = FSpeed[i][j] / sqrt(FTempRef / FTempMeasure);
//				if(FPres[i][j] > 1)
//					FEff[i][j] = HTM->CorrectCompressorMap(m, FPres[i][j], FEff[i][j], FTempMeasure, TinT, Rtc);
//			}
//		}
//	}
//
//}

void TSAEMap::CalculateAdiabaticEfficiency(TTC_HTM2 *HTM, double TinT) {

	double m = 0., eff = 0., Rtc = 0.;

	if(!FIsAdiabatic) {

		FILE *fres = fopen("CompRes.dat", "w");

		fprintf(fres, "Turbocharger_Speed(rmp)\t");
		fprintf(fres, "Compressor_Mass(kg/s)\t");
		fprintf(fres, "Compressor_Ratio(-)\t");
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
		fprintf(fres, "Temp_turbine_out(K)\t");
		fprintf(fres, "Temp_water_out(K)\t");
		fprintf(fres, "Compressor_Power(W)\t");
		fprintf(fres, "Turbine_Power(W)\t");
		fprintf(fres, "Mechanical_Efficiency(W)\t");
		fprintf(fres, "Turbine_In_Enthalpy(W)\n");

		fclose(fres);

		TurboMachine T, C;

		T.it = TinT;
		C.it = FTempMeasure;

		for(int i = 0; i < FNumLines; i++) {
			for(unsigned int j = 0; j < FSpeed[i].size(); j++) {
				C.m = FMass[i][j] / __units::PaToBar(FPresionRef) / sqrt(FTempMeasure / FTempRef);
				C.rtc = FSpeed[i][j] / sqrt(FTempRef / FTempMeasure);
				C.eff = FEff[i][j];
				C.pr = FPres[i][j];
				if(FPres[i][j] > 1)
					FEff[i][j] = HTM->getAdiabaticEfficiency("Compressor", T, C);
			}
		}
	}

}

#pragma package(smart_init)
