/* --------------------------------------------------------------------------------*\
|==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 | License
 |
 |	This file is part of OpenWAM.
 |
 |	OpenWAM is free software: you can redistribute it and/or modify
 |	it under the terms of the GNU General Public License as published by
 |	the Free Software Foundation, either version 3 of the License, or
 |	(at your option) any later version.
 |
 |	OpenWAM is distributed in the hope that it will be useful,
 |	but WITHOUT ANY WARRANTY; without even the implied warranty of
 |	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 |	GNU General Public License for more details.
 |
 |	You should have received a copy of the GNU General Public License
 |	along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 |
 \*-------------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------
#pragma hdrstop

#include "TTC_HTM.h"

// ---------------------------------------------------------------------------

TTC_HTM::TTC_HTM(stHTMoil *Oil) {

	FNumberNodes = 14;

	FNode.Gas = 0;
	FNode.Air = 1;
	FNode.Water = 2;
	FNode.T = 3;
	FNode.H1 = 4;
	FNode.H2 = 5;
	FNode.H3 = 6;
	FNode.C = 7;
	FNode.AirIn = 8;
	FNode.AirOut = 9;
	FNode.Amb = 10;
	FNode.O_H1 = 11;
	FNode.O_H3 = 12;
	FNode.O_H2 = 13;

	FNodeLabel.push_back("G");
	FNodeLabel.push_back("A");
	FNodeLabel.push_back("W");
	FNodeLabel.push_back("T");
	FNodeLabel.push_back("H1");
	FNodeLabel.push_back("H2");
	FNodeLabel.push_back("H3");
	FNodeLabel.push_back("C");
	FNodeLabel.push_back("Ai");
	FNodeLabel.push_back("Ao");
	FNodeLabel.push_back("Amb");
	FNodeLabel.push_back("OH1");
	FNodeLabel.push_back("OH3");
	FNodeLabel.push_back("OH2");

	FTime0 = 0;

	FExternalRadiation = false;
	FExternalConvection = false;
	FIsWaterCooled = false;

	FRePrmu.resize(39, 0);

	FMatrix_Kinv.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_Kinv[i].resize(FNumberNodes, 0.0);
		FMatrix_Kinv[i][i] = 1.0;
	}

	FMatrix_KS.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_KS[i].resize(FNumberNodes, 0.0);
		FMatrix_KS[i][i] = 1.0;
	}
	FMatrix_Kh.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_Kh[i].resize(FNumberNodes, 0.0);
		FMatrix_Kh[i][i] = 0.0;
	}
	FMatrix_Kr.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_Kr[i].resize(FNumberNodes, 0.0);
		FMatrix_Kr[i][i] = 0.0;
	}
	FMatrix_C.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_C[i].resize(FNumberNodes, 0.0);
	}
	FMatrix_KC.resize(2 * FNumberNodes);
	for(int i = 0; i < 2 * FNumberNodes; i++) {
		FMatrix_KC[i].resize(2 * FNumberNodes, 0.0);
		FMatrix_KC[i][i] = 1.0;
	}
	FMatrix_dT.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_dT[i].resize(FNumberNodes, 0.0);
	}
	FMatrix_HF.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_HF[i].resize(FNumberNodes, 0.0);
	}

	FMatrix_Conn.resize(FNumberNodes);
	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_Conn[i].resize(FNumberNodes, false);
	}

	FNode_Temp_K.resize(FNumberNodes);
	FNode_Temp_K_Tr.resize(2 * FNumberNodes);

	FC.Fluid = new stHTMair();
	FT.Fluid = new stHTMair();
	FOil = Oil;
	FWater = new stHTMwater();
	FAmbient = new stHTMair();

	FConverge = 1.;

	FMatrix_Conn[FNode.T][FNode.Gas] = true;
	FMatrix_Conn[FNode.T][FNode.H1] = true;

	FMatrix_Conn[FNode.H1][FNode.O_H1] = true;
	FMatrix_Conn[FNode.H1][FNode.T] = true;
	FMatrix_Conn[FNode.H1][FNode.H2] = true;
	// FMatrix_KS[FNode.H1][FNode.Amb] = true;

	// FMatrix_Conn[FNode.H2][FNode.Water] = true;
	FMatrix_Conn[FNode.H2][FNode.H1] = true;
	FMatrix_Conn[FNode.H2][FNode.H3] = true;
	FMatrix_Conn[FNode.H2][FNode.O_H2] = true;

	FMatrix_Conn[FNode.H3][FNode.H2] = true;
	FMatrix_Conn[FNode.H3][FNode.C] = true;
	// FMatrix_Conn[FNode.H3][FNode.AirIn] = true;
	// FMatrix_Conn[FNode.H3][FNode.O_H3] = true;

	FMatrix_Conn[FNode.C][FNode.Air] = true;
	FMatrix_Conn[FNode.C][FNode.H3] = true;

	FFitTurbExtHeat = 1.7;

}

TTC_HTM::~TTC_HTM() {

	FMechLosses = NULL;

}

void TTC_HTM::InputData(double T_AF, double T_Humidity, double T_MassFlow, double T_IT_C, double T_IP, double T_PR,
						double C_Humidity, double C_MassFlow, double C_IT_C, double C_IP, double C_PR, double O_MassFlow, double O_IT_C,
						double O_IP, double RTC) {

	FT.AF = T_AF; // Turbine A/F
	FT.Humidity = T_Humidity; // Turbine Humidity
	FT.MassFlow = T_MassFlow; // Turbine mass flow (kg/s)
	FT.IT_C = T_IT_C; // Turbine Inlet Temperature (�C)
	FT.IT_K = FT.IT_C + 273.; // Turbine Inlet Temperature (K)
	FT.IP = T_IP; // Turbine Inlet Pressure (bar)
	FT.PR = T_PR; // Turbine expansion ratio (-)

	// Input compressor data
	FC.Humidity = C_Humidity; // Compressor Humidity;
	FC.MassFlow = C_MassFlow; // Compressor mass flow (kg/s)
	FC.IT_C = C_IT_C; // Compressor inlet temperature (�C)
	FC.IT_K = FC.IT_C + 273.; // Compressor inlet temperature (K)
	FC.IP = C_IP; // Compressor inelt pressure (bar)
	FC.PR = C_PR; // Compressor compression ratio (-)

	// Input oil data
	FO_MassFlow = O_MassFlow; // Oil mass flow (kg/s)
	FO_IT_C = O_IT_C; // Oil inlet temperature (�C)
	FO_IT_K = FO_IT_C + 273.; // Oil inlet temperature (K)
	FO_IP = O_IP; // Oil inlet pressure (bar)

	// Input data turbocharger
	FRTC = RTC; // Turbocharger Speed

	FT.MassFlowC = FT.funCORRMass();
	FT.RTC_C = FT.funCORRRTC(RTC);
	FT.TISO_K = FT.funTiso(-1.0);

	FC.MassFlowC = FC.funCORRMass();
	FC.RTC_C = FC.funCORRRTC(RTC);
	FC.TISO_K = FC.funTiso(1.0);

}

void TTC_HTM::TurbochargerData(double DShaft, double HD, double Doil, double DWater, double DT, double LT, double DC,
							   double LC, double DH, double LH) {

	FDShaft = DShaft;
	FHD = HD;
	FDOil = Doil;
	FDWater = DWater;

	FDiamExtTurb = DT;
	FLengExtTurb = LT;
	FAreaExtTurb = Pi * (pow2(DT * FFitTurbExtHeat) - pow2(DH)) / 4.;
	FDiamExtHous = DH;
	FLengExtHous = LH;
	FAreaExtHous = Pi * DH * LH;
	FDiamExtComp = DC;
	FLengExtComp = LC;
	FAreaExtComp = Pi * (pow2(DC) - pow2(DH)) / 4.;
	FAlpha1 = 0.58;
	FAlpha2 = 0.21;
	FAlpha3 = 0.21;

	if(FExternalRadiation) {
		View_Factors(DT * FFitTurbExtHeat, DC, DH, LH);
	}

}

void TTC_HTM::TurbochargerWorkingPoint(double RTC, double MechEff, double O_MassFlow, double O_IT, double O_IP,
									   double W_IT, double W_MassFlow) {
	FRTC = RTC;
	FMechEff = MechEff;
	FO_MassFlow = O_MassFlow;
	FO_IT_K = O_IT;
	FO_IT_C = FO_IT_K - 273;
	FO_IP = O_IP;
	FW_IT_K = W_IT;
	FW_IT_C = W_IT - 273;
	FW_IP = 1;
	FW_MassFlow = W_MassFlow;

	FOilTemps.OI = FO_IT_K;
	FOilTemps.OIH1 = FO_IT_K;
	FOilTemps.OOH2 = FO_IT_K;
	FOilTemps.OO = FO_IT_K;

	FOil->CalcProperties(FO_IP, FO_IT_K);
	if(FIsWaterCooled)
		FWater->CalcProperties(FW_IP, FW_IT_K);
}

void TTC_HTM::CompressorData(double PREF, double TREF, double TMAP_K, double Din, double Dout) {
	FC.PREF = PREF;
	FC.TREF_K = TREF;
	FC.TMAP_K = TMAP_K;
	FC.DIn = Din;
	FC.DOut = Dout;
	FC.SecIn = Din * Din * Pi / 4;
	FC.SecOut = Dout * Dout * Pi / 4;
}

void TTC_HTM::CompressorWorkingPoint(double C_Humidity, double C_MassFlow, double C_IT_K, double C_IP, double C_PR,
									 double EFF) {

	// Input compressor data
	FC.Humidity = C_Humidity; // Compressor Humidity;
	FC.MassFlow = C_MassFlow; // Compressor mass flow (kg/s)
	FC.IT_K = C_IT_K; // Compressor inlet temperature (K)
	FC.IT_C = FC.IT_K - 273.; // Compressor inlet temperature (�C)
	FC.IP0 = C_IP; // Compressor inelt pressure (bar)
	FC.PR = C_PR; // Compressor compression ratio (-)
	FC.Fluid->Hum = C_Humidity;
	FC.EFF = EFF;

	// double a=FC.Fluid->fun_mu(200);
	FC.Fluid->CalcProperties(FC.IP0, FC.IT_K);

	// FC.MassFlowC = FC.funCORRMass();
	// FC.RTC_C = FC.funCORRRTC(FRTC);

	int iter = 0;
	double error = 1;
	FC.IP = FC.IP0;
	do {
		FC.IT0_K = FC.funT0_in();
		error = FC.funP_toStatic(FC.IP0, FC.IT_K, FC.IT0_K) - FC.IP;
		FC.IP = FC.IP + error;
	} while(iter < 100 && fabs(error) > 0.001);

	FC.TISO_K = FC.funTiso(1.0);
	FC.OP0 = FC.PR * FC.IP0;
}

void TTC_HTM::TurbineData(double PREF, double TREF, double TMAP_K, double Din, double Dout) {
	FT.PREF = PREF;
	FT.TREF_K = TREF;
	FT.TMAP_K = TMAP_K;
	FT.DIn = Din;
	FT.DOut = Dout;
	FT.SecIn = Din * Din * Pi / 4;
}

void TTC_HTM::TurbineWorkingPoint(double T_AF, double T_Humidity, double T_MassFlow, double T_IT_K, double T_IP,
								  double T_PR, double EFF) {
	FT.AF = T_AF; // Turbine A/F
	FT.Humidity = T_Humidity; // Turbine Humidity
	FT.MassFlow = T_MassFlow; // Turbine mass flow (kg/s)
	FT.IT_K = T_IT_K; // Turbine Inlet Temperature (K)
	FT.IT_C = FT.IT_K - 273.; // Turbine Inlet Temperature (�C)
	FT.IP0 = T_IP; // Turbine Inlet Pressure (bar)
	FT.PR = T_PR; // Turbine expansion ratio (-)
	FT.EFF = EFF;

	FT.Fluid->CalcProperties(FT.IP0, FT.IT_K);

	// FT.MassFlowC = FT.funCORRMass();
	// FT.RTC_C = FT.funCORRRTC(FRTC);
	FT.IP = FT.IP0;
	int iter = 0;
	double error = 1;

	do {
		FT.IT0_K = FT.funT0_in();
		error = FT.funP_toStatic(FT.IP0, FT.IT_K, FT.IT0_K) - FT.IP;
		FT.IP = FT.IP + error;
	} while(iter < 100 && fabs(error) > 0.001);
	FT.TISO_K = FT.funTiso(-1.0);
}

void TTC_HTM::BuildMatrix() {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			if(i != j)
				FMatrix_KS[i][j] = 0.;
		}
	}

	FT.Fluid->CalcProperties(FT.IP, FT.IT_K);
	FC.Fluid->CalcProperties(FC.IP, FC.IT_K);

	FOil->CalcProperties(1, FO_IT_K);
	FWater->CalcProperties(FW_IP, FW_IT_K);

	// **** Viscosity relations ****
	FRePrmu[nmmuO_H1] = FOil->fun_mu(FNode_Temp_K[FNode.O_H1]) / FOil->fun_mu(FNode_Temp_K[FNode.H1]);
	FRePrmu[nmmuO_H2] = FOil->fun_mu(FNode_Temp_K[FNode.O_H2]) / FOil->fun_mu(FNode_Temp_K[FNode.H2]);
	FRePrmu[nmmuG_T] = FT.Fluid->fun_mu(FNode_Temp_K[FNode.Gas]) / FT.Fluid->fun_mu(FNode_Temp_K[FNode.T]);

	// **** Reynolds numbers ****
	FRePrmu[nmReMassTur] = 4 * fabs(FT.MassFlow) / (Pi * FT.Fluid->fun_mu(FNode_Temp_K[FNode.Gas]) * FT.DIn);
	FRePrmu[nmReMassCom] = 4 * fabs(FC.MassFlow) / (Pi * FC.Fluid->fun_mu(FNode_Temp_K[FNode.Air]) * FC.DOut);
	FRePrmu[nmReOilH1] = 4 * FO_MassFlow / (Pi * FOil->fun_mu(FNode_Temp_K[FNode.O_H1]) * FDOil);
	FRePrmu[nmReOilH2] = 4 * FO_MassFlow / (Pi * FOil->fun_mu(FNode_Temp_K[FNode.O_H2]) * FDOil);

	if(FIsWaterCooled) {
		FRePrmu[nmReWater] = 4 * FW_MassFlow / (Pi * FWater->mu * FDWater);
		FRePrmu[nmReMassCom2] = 4 * fabs(FC.MassFlow) / (Pi * FC.Fluid->fun_mu(FNode_Temp_K[FNode.AirIn]) * FC.DOut);
	} else {
		FRePrmu[nmReMassCom2] = 4 * fabs(FC.MassFlow) / (Pi * FC.Fluid->fun_mu(FNode_Temp_K[FNode.AirOut]) * FC.DIn);
	}
	FRePrmu[nmReShaft] = FOil->fun_rho(FNode_Temp_K[FNode.O_H2]) * (Pi * FRTC / 60 * FDShaft) * (FHD) / FOil->fun_mu(
							 FNode_Temp_K[FNode.O_H2]);

	// **** Prandt numbers ****
	FRePrmu[nmPrTur] = FT.Fluid->fun_Pr(FNode_Temp_K[FNode.Gas]);
	if(FIsWaterCooled)
		FRePrmu[nmPrWater] = FWater->Pr;
	FRePrmu[nmPrOilH1] = FOil->fun_Pr(FNode_Temp_K[FNode.O_H1]);
	FRePrmu[nmPrOilH2] = FOil->fun_Pr(FNode_Temp_K[FNode.O_H2]);
	FRePrmu[nmPrCom] = FC.Fluid->fun_Pr(FNode_Temp_K[FNode.Air]);
	if(FIsWaterCooled)
		FRePrmu[nmPrCom2] = FC.Fluid->fun_Pr(FNode_Temp_K[FNode.AirIn]);
	else
		FRePrmu[nmPrCom2] = FC.Fluid->fun_Pr(FNode_Temp_K[FNode.AirOut]);

	// **** Conductances ****
	// --- [ T ] --- //
	FMatrix_KS[FNode.T][FNode.Gas] = FConv.GAS_T.Value(FRePrmu) * FT.Fluid->k * pow(FT.EFFMAX,
									 1.467) * Pi * pow2(FDiamExtTurb) / 4 / FT.DIn;

	//if (std::isnan(FMatrix_KS[FNode.T][FNode.Gas]))
	//	cout << FMatrix_KS[FNode.T][FNode.Gas] << endl;
	FMatrix_KS[FNode.T][FNode.H1] = K_T_H1();

	// --- [ H1 ] --- //
	FRePrmu[nmReShaft] = FOil->fun_rho(FNode_Temp_K[FNode.O_H1]) * (Pi * FRTC / 60 * FDShaft) * (FHD) / FOil->fun_mu(
							 FNode_Temp_K[FNode.O_H1]);
	FMatrix_KS[FNode.H1][FNode.O_H1] = FConv.OIL_H1.Value(FRePrmu) * FOil->fun_k(FNode_Temp_K[FNode.O_H1]) * Pi *
									   FLengExtHous;
	FMatrix_KS[FNode.H1][FNode.T] = K_T_H1();
	FMatrix_KS[FNode.H1][FNode.H2] = K_H1_H2();
	// FMatrix_KS[FNode.H1][FNode.Amb] = 0;

	// --- [ H2 ] --- //
	FRePrmu[nmReShaft] = FOil->fun_rho(FNode_Temp_K[FNode.O_H2]) * (Pi * FRTC / 60 * FDShaft) * (FHD) / FOil->fun_mu(
							 FNode_Temp_K[FNode.O_H2]);
	if(FIsWaterCooled)
		FMatrix_KS[FNode.H2][FNode.Water] = FConv.WAT_H2.Value(FRePrmu) * FWater->k * Pi * FLengExtHous;
	FMatrix_KS[FNode.H2][FNode.H1] = K_H1_H2();
	FMatrix_KS[FNode.H2][FNode.H3] = K_H2_H3();
	if(FNode_Temp_K[FNode.H2] < FNode_Temp_K[FNode.O_H2])
		FMatrix_KS[FNode.H2][FNode.O_H2] = FConv.OIL_H2_C.Value(FRePrmu) * FOil->fun_k(FNode_Temp_K[FNode.O_H2]) * Pi *
										   FLengExtHous;
	else
		FMatrix_KS[FNode.H2][FNode.O_H2] = FConv.OIL_H2_H.Value(FRePrmu) * FOil->fun_k(FNode_Temp_K[FNode.O_H2]) * Pi *
										   FLengExtHous;

	// --- [ H3 ] --- //
	FMatrix_KS[FNode.H3][FNode.H2] = K_H2_H3();
	FMatrix_KS[FNode.H3][FNode.C] = K_H3_C();
	if(!FIsWaterCooled) {
		if(FNode_Temp_K[FNode.H3] < FNode_Temp_K[FNode.AirOut])
			FMatrix_KS[FNode.H3][FNode.AirOut] = FConv.AIR_H3_C.Value(FRePrmu) * FC.Fluid->fun_k(
					FNode_Temp_K[FNode.AirOut]) * Pi * FDiamExtComp;
		else
			FMatrix_KS[FNode.H3][FNode.AirOut] = FConv.AIR_H3_H.Value(FRePrmu) * FC.Fluid->fun_k(
					FNode_Temp_K[FNode.AirOut]) * Pi * FDiamExtComp;
	}
	// FMatrix_KS[FNode.H3][FNode.O_H3] = FConv.OIL_H3.Value(FRePrmu) * FOil->fun_k
	// (FNode_Temp_K[FNode.O_H3]);

	// --- [ C ] --- //
	if(FNode_Temp_K[FNode.C] < FNode_Temp_K[FNode.Air])
		FMatrix_KS[FNode.C][FNode.Air] = FConv.AIR_C_C.Value(FRePrmu) * FC.Fluid->k * Pi * FDiamExtComp;
	else
		FMatrix_KS[FNode.C][FNode.Air] = FConv.AIR_C_H.Value(FRePrmu) * FC.Fluid->k * Pi * FDiamExtComp;

	FMatrix_KS[FNode.C][FNode.H3] = K_H3_C();

	FMatrix_KS[FNode.T][FNode.T] = 0;
	FMatrix_KS[FNode.H1][FNode.H1] = 0;
	FMatrix_KS[FNode.H2][FNode.H2] = 0;
	FMatrix_KS[FNode.H3][FNode.H3] = 0;
	FMatrix_KS[FNode.C][FNode.C] = 0;

	if(FExternalRadiation) {

		FMatrix_Kr[FNode.T][FNode.H1] = KRadiationSurfaces(FAreaExtTurb, FAreaExtHous * FAlpha1, FEmisExtTurb, FEmisExtHous,
										FVF.T_H1, FNode_Temp_K[FNode.T], FNode_Temp_K[FNode.H1]) * FFitRadiation;
		FMatrix_Kr[FNode.T][FNode.H2] = KRadiationSurfaces(FAreaExtTurb, FAreaExtHous * FAlpha2, FEmisExtTurb, FEmisExtHous,
										FVF.T_H2, FNode_Temp_K[FNode.T], FNode_Temp_K[FNode.H2]) * FFitRadiation;
		FMatrix_Kr[FNode.T][FNode.H3] = KRadiationSurfaces(FAreaExtTurb, FAreaExtHous * FAlpha3, FEmisExtTurb, FEmisExtHous,
										FVF.T_H3, FNode_Temp_K[FNode.T], FNode_Temp_K[FNode.H3]) * FFitRadiation;
		FMatrix_Kr[FNode.T][FNode.C] = KRadiationSurfaces(FAreaExtTurb, FAreaExtComp, FEmisExtTurb, FEmisExtComp, FVF.T_C,
									   FNode_Temp_K[FNode.T], FNode_Temp_K[FNode.C]) * FFitRadiation;
		FMatrix_Kr[FNode.T][FNode.Amb] = KRadiationExternal(FLengExtTurb, FDiamExtTurb * FFitTurbExtHeat, FEmisExtTurb,
										 FVF.T_Amb, FNode_Temp_K[FNode.T], FNode_Temp_K[FNode.Amb]) * FFitRadiation;

		FMatrix_Kr[FNode.C][FNode.T] = KRadiationSurfaces(FAreaExtComp, FAreaExtTurb, FEmisExtComp, FEmisExtTurb, FVF.C_T,
									   FNode_Temp_K[FNode.C], FNode_Temp_K[FNode.T]) * FFitRadiation;
		FMatrix_Kr[FNode.C][FNode.H1] = KRadiationSurfaces(FAreaExtComp, FAreaExtHous * FAlpha1, FEmisExtComp, FEmisExtHous,
										FVF.C_H1, FNode_Temp_K[FNode.C], FNode_Temp_K[FNode.H1]) * FFitRadiation;
		FMatrix_Kr[FNode.C][FNode.H2] = KRadiationSurfaces(FAreaExtComp, FAreaExtHous * FAlpha2, FEmisExtComp, FEmisExtHous,
										FVF.C_H2, FNode_Temp_K[FNode.C], FNode_Temp_K[FNode.H2]) * FFitRadiation;
		FMatrix_Kr[FNode.C][FNode.H3] = KRadiationSurfaces(FAreaExtComp, FAreaExtHous * FAlpha3, FEmisExtComp, FEmisExtHous,
										FVF.C_H3, FNode_Temp_K[FNode.C], FNode_Temp_K[FNode.H3]) * FFitRadiation;
		FMatrix_Kr[FNode.C][FNode.Amb] = KRadiationExternal(FLengExtComp, FDiamExtComp, FEmisExtComp, FVF.C_Amb,
										 FNode_Temp_K[FNode.C], FNode_Temp_K[FNode.Amb]) * FFitRadiation;

		FMatrix_Kr[FNode.H1][FNode.T] = KRadiationSurfaces(FAreaExtHous * FAlpha1, FAreaExtTurb, FEmisExtHous, FEmisExtTurb,
										FVF.H1_T, FNode_Temp_K[FNode.H1], FNode_Temp_K[FNode.T]) * FFitRadiation;
		FMatrix_Kr[FNode.H1][FNode.C] = KRadiationSurfaces(FAreaExtHous * FAlpha1, FAreaExtComp, FEmisExtHous, FEmisExtComp,
										FVF.H1_C, FNode_Temp_K[FNode.H1], FNode_Temp_K[FNode.C]) * FFitRadiation;
		FMatrix_Kr[FNode.H1][FNode.Amb] = KRadiationExternalHous(FAreaExtHous * FAlpha1, FEmisExtHous, FVF.H1_Amb,
										  FNode_Temp_K[FNode.H1], FNode_Temp_K[FNode.Amb]) * FFitRadiation;

		FMatrix_Kr[FNode.H2][FNode.T] = KRadiationSurfaces(FAreaExtHous * FAlpha2, FAreaExtTurb, FEmisExtHous, FEmisExtTurb,
										FVF.H2_T, FNode_Temp_K[FNode.H2], FNode_Temp_K[FNode.T]) * FFitRadiation;
		FMatrix_Kr[FNode.H2][FNode.C] = KRadiationSurfaces(FAreaExtHous * FAlpha2, FAreaExtComp, FEmisExtHous, FEmisExtComp,
										FVF.H2_C, FNode_Temp_K[FNode.H2], FNode_Temp_K[FNode.C]) * FFitRadiation;
		FMatrix_Kr[FNode.H2][FNode.Amb] = KRadiationExternalHous(FAreaExtHous * FAlpha2, FEmisExtHous, FVF.H1_Amb,
										  FNode_Temp_K[FNode.H2], FNode_Temp_K[FNode.Amb]) * FFitRadiation;

		FMatrix_Kr[FNode.H3][FNode.T] = KRadiationSurfaces(FAreaExtHous * FAlpha3, FAreaExtTurb, FEmisExtHous, FEmisExtTurb,
										FVF.H3_T, FNode_Temp_K[FNode.H3], FNode_Temp_K[FNode.T]) * FFitRadiation;
		FMatrix_Kr[FNode.H3][FNode.C] = KRadiationSurfaces(FAreaExtHous * FAlpha2, FAreaExtComp, FEmisExtHous, FEmisExtComp,
										FVF.H3_C, FNode_Temp_K[FNode.H3], FNode_Temp_K[FNode.C]) * FFitRadiation;
		FMatrix_Kr[FNode.H3][FNode.Amb] = KRadiationExternalHous(FAreaExtHous * FAlpha3, FEmisExtHous, FVF.H1_Amb,
										  FNode_Temp_K[FNode.H3], FNode_Temp_K[FNode.Amb]) * FFitRadiation;

		AddKConvRad(FMatrix_KS, FMatrix_Kr);
	}
	if(FExternalConvection) {
		double Tprop = 0., Beta = 0., rho = 0., mu = 0., visc = 0., Gr = 0., Re = 0., Pr = 0., Nu = 0.;
		double g = 9.8;
		double pamb = 1.;

		// Turbine to ambient

		Tprop = (FNode_Temp_K[FNode.T] + FNode_Temp_K[FNode.Amb]) / 2;

		Beta = 1 / Tprop;
		rho = pamb * 1e5 / FAmbient->fun_R() / Tprop;
		mu = FAmbient->fun_mu(Tprop);
		visc = mu / rho;
		Pr = FAmbient->fun_Pr(Tprop);

		Gr = g * Beta * (FNode_Temp_K[FNode.T] - FNode_Temp_K[FNode.Amb]) * pow3(FDiamExtTurb * FFitTurbExtHeat) / pow2(visc);

		Re = rho * FExtVelocity * FFitTurbExtHeat * FDiamExtTurb / mu;

		if(Gr > 10 * pow2(Re)) {
			Nu = NusseltFreeConv(Gr, Pr);
		} else if(Gr < 0.1 * pow2(Re)) {
			Nu = NusseltForcConv(Re, Pr);
		} else {
			Nu = pow(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)), 1 / 3);
		}
		FMatrix_Kh[FNode.T][FNode.Amb] = Pi * FLengExtTurb * Nu * FAmbient->fun_k(Tprop) * FFitConvection;

		// H1 to ambient

		Tprop = (FNode_Temp_K[FNode.H1] + FNode_Temp_K[FNode.Amb]) / 2;

		Beta = 1 / Tprop;
		rho = pamb * 1e5 / FAmbient->fun_R() / Tprop;
		mu = FAmbient->fun_mu(Tprop);
		visc = mu / rho;
		Pr = FAmbient->fun_Pr(Tprop);

		Gr = g * Beta * (FNode_Temp_K[FNode.H1] - FNode_Temp_K[FNode.Amb]) * pow3(FDiamExtHous) / pow2(visc);

		Re = rho * FExtVelocity * FDiamExtHous / mu;

		if(Gr > 10 * pow2(Re)) {
			Nu = NusseltFreeConv(Gr, Pr);
		} else if(Gr < 0.1 * pow2(Re)) {
			Nu = NusseltForcConv(Re, Pr);
		} else {
			Nu = pow(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)), 1 / 3);
		}
		FMatrix_Kh[FNode.H1][FNode.Amb] = Pi * FLengExtHous * FAlpha1 * Nu * FAmbient->fun_k(Tprop) * FFitConvection;

		// H2 to ambient

		Tprop = (FNode_Temp_K[FNode.H2] + FNode_Temp_K[FNode.Amb]) / 2;

		Beta = 1 / Tprop;
		rho = pamb * 1e5 / FAmbient->fun_R() / Tprop;
		mu = FAmbient->fun_mu(Tprop);
		visc = mu / rho;
		Pr = FAmbient->fun_Pr(Tprop);

		Gr = g * Beta * (FNode_Temp_K[FNode.H2] - FNode_Temp_K[FNode.Amb]) * pow3(FDiamExtHous) / pow2(visc);

		Re = rho * FExtVelocity * FDiamExtHous / mu;

		if(Gr > 10 * pow2(Re)) {
			Nu = NusseltFreeConv(Gr, Pr);
		} else if(Gr < 0.1 * pow2(Re)) {
			Nu = NusseltForcConv(Re, Pr);
		} else {
			Nu = pow(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)), 1 / 3);
		}
		FMatrix_Kh[FNode.H2][FNode.Amb] = Pi * FLengExtHous * FAlpha2 * Nu * FAmbient->fun_k(Tprop) * FFitConvection;

		// H3 to ambient

		Tprop = (FNode_Temp_K[FNode.H3] + FNode_Temp_K[FNode.Amb]) / 2;

		Beta = 1 / Tprop;
		rho = pamb * 1e5 / FAmbient->fun_R() / Tprop;
		mu = FAmbient->fun_mu(Tprop);
		visc = mu / rho;
		Pr = FAmbient->fun_Pr(Tprop);

		Gr = g * Beta * (FNode_Temp_K[FNode.H3] - FNode_Temp_K[FNode.Amb]) * pow3(FDiamExtHous) / pow2(visc);

		Re = rho * FExtVelocity * FDiamExtHous / mu;

		if(Gr > 10 * pow2(Re)) {
			Nu = NusseltFreeConv(Gr, Pr);
		} else if(Gr < 0.1 * pow2(Re)) {
			Nu = NusseltForcConv(Re, Pr);
		} else {
			Nu = pow(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)), 1 / 3);
		}
		FMatrix_Kh[FNode.H3][FNode.Amb] = Pi * FLengExtHous * FAlpha3 * Nu * FAmbient->fun_k(Tprop) * FFitConvection;

		// Compressor to ambient

		Tprop = (FNode_Temp_K[FNode.C] + FNode_Temp_K[FNode.Amb]) / 2;

		Beta = 1 / Tprop;
		rho = pamb * 1e5 / FAmbient->fun_R() / Tprop;
		mu = FAmbient->fun_mu(Tprop);
		visc = mu / rho;
		Pr = FAmbient->fun_Pr(Tprop);

		Gr = g * Beta * (FNode_Temp_K[FNode.C] - FNode_Temp_K[FNode.Amb]) * pow3(FDiamExtComp) / pow2(visc);

		Re = rho * FExtVelocity * FDiamExtComp / mu;

		if(Gr > 10 * pow2(Re)) {
			Nu = NusseltFreeConv(Gr, Pr);
		} else if(Gr < 0.1 * pow2(Re)) {
			Nu = NusseltForcConv(Re, Pr);
		} else {
			Nu = pow(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)), 1 / 3);
		}
		FMatrix_Kh[FNode.C][FNode.Amb] = Pi * FLengExtComp * Nu * FAmbient->fun_k(Tprop) * FFitConvection;

		AddKConvRad(FMatrix_KS, FMatrix_Kh);
	}

	for(int i = 0; i < FNumberNodes; i++) {

		if(i != FNode.T)
			FMatrix_KS[FNode.T][FNode.T] -= FMatrix_KS[FNode.T][i];

		if(i != FNode.H1)
			FMatrix_KS[FNode.H1][FNode.H1] -= FMatrix_KS[FNode.H1][i];

		if(i != FNode.H2)
			FMatrix_KS[FNode.H2][FNode.H2] -= FMatrix_KS[FNode.H2][i];

		if(i != FNode.H3)
			FMatrix_KS[FNode.H3][FNode.H3] -= FMatrix_KS[FNode.H3][i];

		if(i != FNode.C)
			FMatrix_KS[FNode.C][FNode.C] -= FMatrix_KS[FNode.C][i];

	}
}

void TTC_HTM::SolveNodeTemperatures(double TET, double TSC, double TEC, double TOIL, double MassOil, double MechLosses,
									double TW, double TAMB) {

	FOilTemps.OIH1 = FOilTemps.OI - FMatrix_HF[FNode.H1][FNode.O_H1] / FO_MassFlow / FOil->fun_Cp(FNode_Temp_K[FNode.O_H1]);
	if(FMatrix_HF[FNode.H1][FNode.O_H1] < 0 && FOilTemps.OIH1 > FNode_Temp_K[FNode.H1])
		FOilTemps.OIH1 = FNode_Temp_K[FNode.H1];
	if(FMatrix_HF[FNode.H1][FNode.O_H1] > 0 && FOilTemps.OIH1 < FNode_Temp_K[FNode.H1])
		FOilTemps.OIH1 = FNode_Temp_K[FNode.H1];

	FOilTemps.OOH2 = FOilTemps.OIH1 + MechLosses / FO_MassFlow / FOil->fun_Cp((FOilTemps.OIH1 + FOilTemps.OOH2) / 2);

	FOilTemps.OO = FOilTemps.OOH2 - FMatrix_HF[FNode.H2][FNode.O_H2] / FO_MassFlow / FOil->fun_Cp(FNode_Temp_K[FNode.O_H2]);
	if(FMatrix_HF[FNode.H2][FNode.O_H2] < 0 && FOilTemps.OO > FNode_Temp_K[FNode.H2])
		FOilTemps.OO = FNode_Temp_K[FNode.H2];
	if(FMatrix_HF[FNode.H2][FNode.O_H2] > 0 && FOilTemps.OO < FNode_Temp_K[FNode.H2])
		FOilTemps.OO = FNode_Temp_K[FNode.H2];

	FNode_Temp_K[FNode.O_H1] = (FOilTemps.OI + FOilTemps.OIH1) / 2;
	FNode_Temp_K[FNode.O_H2] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;

	BuildMatrix();

	LUdcmp TempSolver(FMatrix_KS);

	// TempSolver.inverse(FMatrix_Kinv);

	dVector InputTemp;
	InputTemp.resize(FNumberNodes, 0.0);

	InputTemp[FNode.Gas] = TET;
	InputTemp[FNode.Air] = TSC;
	if(FC.MassFlow > 1e-6) {
		double DeltaT = (FMatrix_HF[FNode.H3][FNode.AirOut] + FMatrix_HF[FNode.C][FNode.Air] / 2.) / FC.MassFlow /
						FC.Fluid->fun_Cp(FNode_Temp_K[FNode.AirOut]);

		InputTemp[FNode.Air] -= DeltaT;
		if(DeltaT > 0 && InputTemp[FNode.Air] < FOilTemps.OI && TSC > FOilTemps.OI)
			InputTemp[FNode.Air] = FOilTemps.OI;
		if(DeltaT < 0 && InputTemp[FNode.Air] > TET)
			InputTemp[FNode.Air] = TET;
	}
	InputTemp[FNode.O_H1] = (FOilTemps.OI + FOilTemps.OIH1) / 2;
	InputTemp[FNode.O_H3] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;
	InputTemp[FNode.O_H2] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;
	InputTemp[FNode.AirIn] = TEC;
	InputTemp[FNode.AirOut] = TSC;
	if(FC.MassFlow > 1e-6) {
		double DeltaT = FMatrix_HF[FNode.H3][FNode.AirOut] / 2 / FC.MassFlow / FC.Fluid->fun_Cp(FNode_Temp_K[FNode.AirOut]);

		InputTemp[FNode.AirOut] -= DeltaT;
		if(DeltaT > 0 && InputTemp[FNode.Air] < FOilTemps.OI)
			InputTemp[FNode.AirOut] = FOilTemps.OI;
		if(DeltaT < 0 && InputTemp[FNode.Air] > TET)
			InputTemp[FNode.AirOut] = TET;
	}
	InputTemp[FNode.Water] = TW;
	InputTemp[FNode.Amb] = TAMB;

//	FILE *fheat=fopen("TempsIn.txt","w");
//	for (int i = 0; i < FNumberNodes; i++) {
//		fprintf(fheat,"%g\t",InputTemp[i]);
//		fprintf(fheat,"\n");
//	}
//	fclose(fheat);

	TempSolver.solve(InputTemp, FNode_Temp_K);

	for(int i = 0; i < FNumberNodes; i++) {
		FNode_Temp_K_Tr[i] = (FNode_Temp_K[i] + FNode_Temp_K_Tr[i]) / 2.;
	}
}

void TTC_HTM::SolveDeltaTemp() {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			FMatrix_dT[i][j] = FNode_Temp_K[j] - FNode_Temp_K[i];
		}
	}

}

void TTC_HTM::SolveDeltaTempTr() {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			FMatrix_dT[i][j] = FNode_Temp_K_Tr[j] - FNode_Temp_K_Tr[i];
		}
	}

}

void TTC_HTM::SolveHeatFlowMatix() {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			FMatrix_HF[i][j] = FMatrix_KS[i][j] * FMatrix_dT[i][j];
		}
	}
}

double TTC_HTM::Oil_Heat_Flow() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.O_H1] + FMatrix_HF[i][FNode.O_H2] + FMatrix_HF[i][FNode.O_H3];
	}
	return heat;
}

double TTC_HTM::Ambient_Heat_Flow() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.Amb];
	}
	return heat;
}

double TTC_HTM::Water_Heat_Flow() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.Water];
	}
	return heat;
}

double TTC_HTM::Turb_Heat_Flow() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.Gas];
	}
	return heat;
}

double TTC_HTM::Comp_Heat_Flow() {

	return Comp_Heat_Flow_1() + Comp_Heat_Flow_2();
}

double TTC_HTM::Comp_Heat_Flow_1() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.Air];
	}
	return heat;
}

double TTC_HTM::Comp_Heat_Flow_2() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.AirOut];
	}
	return heat;
}

double TTC_HTM::Comp_Heat_Flow_In() {

	double heat = 0;

	for(int i = 0; i < FNumberNodes; i++) {
		heat += FMatrix_HF[i][FNode.AirIn];
	}
	return heat;
}

double TTC_HTM::AdiabaticEff(nmSide Case) {

	double HeatFlow_C = 0., HeatFlow_T = 0., HeatFlow_O = 0., HeatFlow_CI = 0.;
	double TOT0 = 0., OOT0 = 0.;
	double C_RealPower = 0., C_IsoPower = 0., T_IsoPower = 0., T_RealPower = 0., O_RealPower = 0.;
	double DOT = 0., DOT0 = 0., DOTs = 0., TIT = 0., TOT = 0., OOT = 0., ONT = 0.;
	bool error = true;
	double DT_max = 0.05;
	double AdiEff = 0.;
	int itera = 0;
	double CIT = FC.IT_K;
	double TAMB = 298.; //  = 0., PASAR EL VALOR DE TEMPERATURA AMBIENTE
	double DeltaTMech = 0.;

	C_IsoPower = FC.funIsoPower(1.0);

	FC.OT0_K = FC.funT0_in() + C_IsoPower / FC.MassFlow / FC.Fluid->Cp / FC.EFF;
	FC.OT_K = FC.funT_toStatic(FC.OT0_K, FC.OP0, FC.MassFlow);
	DOTs = FC.TISO_K;
	DOTs = FC.funT0_in() + C_IsoPower / FC.MassFlow / FC.Fluid->Cp;

	if(FC.OT_K < FO_IT_K) {
		DOT = FC.OT_K - 10;
		if(DOT < DOTs)
			DOT = DOTs;
	} else
		DOT = FC.OT_K + 10;

	TOT = FT.funTiso(-1.0);

	O_RealPower = FMechLosses->P_oil(FO_IT_K, Pi * FRTC / 30, 1., FC.PR, FT.IP, 1., FO_MassFlow);

	O_RealPower = FT.funIsoPower(-1) * FT.EFF * (1 - FMechEff);
	OOT = FO_IT_K + O_RealPower / FO_MassFlow / FOil->Cp;
	;

	DOT0 = DOT;
	TOT0 = TOT;
	OOT0 = OOT;

	// Initial Temperatures;

	FNode_Temp_K[FNode.Gas] = FT.TMAP_K;
	FNode_Temp_K[FNode.Air] = DOT;
	FNode_Temp_K[FNode.Water] = FW_IT_K;
	FNode_Temp_K[FNode.T] = 0.80 * (FT.TMAP_K - FC.TMAP_K) + FC.TMAP_K;
	FNode_Temp_K[FNode.H1] = 0.50 * (FT.TMAP_K - FC.TMAP_K) + FC.TMAP_K;
	FNode_Temp_K[FNode.H2] = 0.15 * (FT.TMAP_K - FC.TMAP_K) + FC.TMAP_K;
	FNode_Temp_K[FNode.H3] = 0.10 * (FT.TMAP_K - FC.TMAP_K) + FC.TMAP_K;
	FNode_Temp_K[FNode.C] = 0.05 * (FT.TMAP_K - FC.TMAP_K) + FC.TMAP_K;
	FNode_Temp_K[FNode.O_H1] = FO_IT_K;
	FNode_Temp_K[FNode.O_H2] = OOT;
	FNode_Temp_K[FNode.O_H3] = OOT;
	FNode_Temp_K[FNode.AirIn] = CIT;
	FNode_Temp_K[FNode.AirOut] = DOT;
	FNode_Temp_K[FNode.Amb] = TAMB;

	InitializeHeatFlow();

	while(error) {

		SolveNodeTemperatures(FT.TMAP_K, DOT, CIT, FO_IT_K, FO_MassFlow, O_RealPower, FW_IT_K, TAMB);

		SolveDeltaTemp();

		SolveHeatFlowMatix();

		HeatFlow_T = Turb_Heat_Flow();
		FT.IT_K = FT.TMAP_K - HeatFlow_T / FT.MassFlow / FT.Fluid->Cp;
		FT.TISO_K = FT.funTiso(-1.0);
		T_IsoPower = FT.MassFlow * FT.Fluid->Cp * (FT.funT0_in() - FT.TISO_K);

		HeatFlow_C = Comp_Heat_Flow();
		// HeatFlow_CI = Comp_Heat_Flow_In();

		// CIT = FC.IT_K + HeatFlow_CI / FC.MassFlow / FC.Fluid->Cp;
		DOT = FC.OT_K + HeatFlow_C / FC.MassFlow / FC.Fluid->Cp;
		if(DOT < DOTs)
			DOT = DOTs;

		FC.OT0_K = DOT + pow2(FC.MassFlow * 287 * DOT / (FC.OP * 1e5 * FC.SecOut)) / (2 * FC.Fluid->Cp);
		FC.IT0_K = CIT + pow2(FC.MassFlow * 287 * CIT / (FC.IP * 1e5 * FC.SecIn)) / (2 * FC.Fluid->Cp);
		C_RealPower = FC.MassFlow * FC.Fluid->Cp * (FC.OT0_K - FC.IT0_K);

		//DOTs = CIT * pow(FC.PR, (FC.Fluid->g - 1) / FC.Fluid->g);
		//C_IsoPower = FC.MassFlow * FC.Fluid->Cp * (DOTs - CIT);
		// if (C_RealPower < C_IsoPower) {
		// C_RealPower = FT.EFF * T_IsoPower;
		// DOT = CIT + C_RealPower / FC.MassFlow / FC.Fluid->Cp;
		// }

		ONT = FO_IT_K - FMatrix_HF[FNode.H1][FNode.O_H1] / FO_MassFlow / FOil->fun_Cp(FNode_Temp_K[FNode.O_H1]);

		O_RealPower = FMechLosses->P_oil(ONT, Pi * FRTC / 30, 1., FC.PR, FT.IP, 1., FO_MassFlow);

		T_RealPower = C_RealPower + O_RealPower;

		FT.IT0_K = FT.IT_K + pow2(FT.MassFlow * 287 * FT.IT_K / (FT.IP * 1e5 * FC.SecIn)) / (2 * FT.Fluid->Cp);
		TOT = FT.IT0_K - T_RealPower / FT.MassFlow / FT.Fluid->Cp;

		// O_RealPower = T_RealPower - C_RealPower;

		DeltaTMech = O_RealPower / FO_MassFlow / FOil->fun_Cp((FNode_Temp_K[FNode.O_H1] + FNode_Temp_K[FNode.O_H2]) / 2.);

		HeatFlow_O = Oil_Heat_Flow();
		OOT = FO_IT_K - HeatFlow_O / FO_MassFlow / FOil->fun_Cp((FNode_Temp_K[FNode.O_H1] + FNode_Temp_K[FNode.O_H3]) / 2.) +
			  DeltaTMech;

		if((fabs(DOT - DOT0) > DT_max || fabs(TOT - TOT0) > DT_max || fabs(OOT - OOT0) > DT_max) && itera < 100) {
			DOT0 = DOT;
			TOT0 = TOT;
			OOT0 = OOT;
			itera += 1;
		} else {
			error = false;
			if(itera == 100) {
				printf("Process does not converge. Error: Turbine %lf K, Compressor %lf K, Oil %lf K\n", TOT - TOT0, DOT - DOT0,
					   OOT - OOT0);
			}
		}

	}

	FC.Power = C_RealPower;
	FT.Power = T_RealPower;
	FMechEff = C_RealPower / T_RealPower;

	if(Case == nmCompressor) {
		AdiEff = C_IsoPower / C_RealPower;
	} else if(Case == nmTurbine) {
		AdiEff = T_RealPower / T_IsoPower;
	}
	return AdiEff;
}

void TTC_HTM::BuildCMatrix(double dt) {

	FMatrix_C[FNode.T][FNode.T] = FCap.T * FConverge / dt;
	FMatrix_C[FNode.H1][FNode.H1] = FCap.H1 * FConverge / dt;
	FMatrix_C[FNode.H2][FNode.H2] = FCap.H2 * FConverge / dt;
	FMatrix_C[FNode.H3][FNode.H3] = FCap.H3 * FConverge / dt;
	FMatrix_C[FNode.C][FNode.C] = FCap.C * FConverge / dt;
}

void TTC_HTM::BuildKCMatrix(double dt) {

	BuildMatrix();

	BuildCMatrix(dt);

	for(int i = 0; i < FNumberNodes; i++) {
		FMatrix_KC[FNode.T][i] = -FMatrix_KS[FNode.T][i] + FMatrix_C[FNode.T][i];
		FMatrix_KC[FNode.H1][i] = -FMatrix_KS[FNode.H1][i] + FMatrix_C[FNode.H1][i];
		FMatrix_KC[FNode.H2][i] = -FMatrix_KS[FNode.H2][i] + FMatrix_C[FNode.H2][i];
		FMatrix_KC[FNode.H3][i] = -FMatrix_KS[FNode.H3][i] + FMatrix_C[FNode.H3][i];
		FMatrix_KC[FNode.C][i] = -FMatrix_KS[FNode.C][i] + FMatrix_C[FNode.C][i];
	}
	FMatrix_KC[FNode.T][FNode.T + FNumberNodes] = -FMatrix_C[FNode.T][FNode.T];
	FMatrix_KC[FNode.H1][FNode.H1 + FNumberNodes] = -FMatrix_C[FNode.H1][FNode.H1];
	FMatrix_KC[FNode.H2][FNode.H2 + FNumberNodes] = -FMatrix_C[FNode.H2][FNode.H2];
	FMatrix_KC[FNode.H3][FNode.H3 + FNumberNodes] = -FMatrix_C[FNode.H3][FNode.H3];
	FMatrix_KC[FNode.C][FNode.C + FNumberNodes] = -FMatrix_C[FNode.C][FNode.C];

}

void TTC_HTM::SolveNodeTemperaturesTransient(double TET, double TSC, double TEC, double TOIL, double MassOil,
		double MechLosses, double TW, double TAMB, double time) {

	double dt = time - FTime0;
	FTime0 = time;
	double tmp = 0.;

	FOilTemps.OIH1 = FOilTemps.OI - FMatrix_HF[FNode.H1][FNode.O_H1] / FO_MassFlow / FOil->fun_Cp(FNode_Temp_K[FNode.O_H1]);
	if(FMatrix_HF[FNode.H1][FNode.O_H1] < 0 && FOilTemps.OIH1 > FNode_Temp_K[FNode.H1])
		FOilTemps.OIH1 = FNode_Temp_K[FNode.H1];
	if(FMatrix_HF[FNode.H1][FNode.O_H1] > 0 && FOilTemps.OIH1 < FNode_Temp_K[FNode.H1])
		FOilTemps.OIH1 = FNode_Temp_K[FNode.H1];

	FOilTemps.OOH2 = FOilTemps.OIH1 + MechLosses / FO_MassFlow / FOil->fun_Cp((FOilTemps.OIH1 + FOilTemps.OOH2) / 2);

	FOilTemps.OO = FOilTemps.OOH2 - FMatrix_HF[FNode.H2][FNode.O_H2] / FO_MassFlow / FOil->fun_Cp(FNode_Temp_K[FNode.O_H2]);
	if(FMatrix_HF[FNode.H2][FNode.O_H2] < 0 && FOilTemps.OO > FNode_Temp_K[FNode.H2])
		FOilTemps.OO = FNode_Temp_K[FNode.H2];
	if(FMatrix_HF[FNode.H2][FNode.O_H2] > 0 && FOilTemps.OO < FNode_Temp_K[FNode.H2])
		FOilTemps.OO = FNode_Temp_K[FNode.H2];

	FNode_Temp_K[FNode.O_H1] = (FOilTemps.OI + FOilTemps.OIH1) / 2;
	FNode_Temp_K[FNode.O_H2] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;

	// tmp = 20 * pow4(time);
	// if (tmp < 100)
	// FConverge = 1 - exp(-tmp);
	// else
	FConverge = 1;

	BuildKCMatrix(dt);

	dVector InputTemp;
	InputTemp.resize(2 * FNumberNodes, 0.0);

	InputTemp[FNode.Gas] = TET;
	InputTemp[FNode.Air] = TSC;
	if(FC.MassFlow > 1e-6) {
		double DeltaT = (FMatrix_HF[FNode.H3][FNode.AirOut] + FMatrix_HF[FNode.C][FNode.Air] / 2.) / FC.MassFlow /
						FC.Fluid->fun_Cp(FNode_Temp_K[FNode.AirOut]);

		InputTemp[FNode.Air] -= DeltaT;
		if(DeltaT > 0 && InputTemp[FNode.Air] < TEC)
			InputTemp[FNode.Air] = TEC;
		if(DeltaT < 0 && InputTemp[FNode.Air] > TET)
			InputTemp[FNode.Air] = TET;
	}
	InputTemp[FNode.O_H1] = (FOilTemps.OI + FOilTemps.OIH1) / 2;
	InputTemp[FNode.O_H3] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;
	InputTemp[FNode.O_H2] = (FOilTemps.OOH2 + FOilTemps.OO) / 2;
	InputTemp[FNode.AirIn] = TEC;
	InputTemp[FNode.AirOut] = TSC;
	if(FC.MassFlow > 1e-6) {
		double DeltaT = FMatrix_HF[FNode.H3][FNode.AirOut] / 2 / FC.MassFlow / FC.Fluid->fun_Cp(FNode_Temp_K[FNode.AirOut]);

		InputTemp[FNode.AirOut] -= DeltaT;
		if(DeltaT > 0 && InputTemp[FNode.Air] < TEC)
			InputTemp[FNode.AirOut] = TEC;
		if(DeltaT < 0 && InputTemp[FNode.Air] > TET)
			InputTemp[FNode.AirOut] = TET;
	}
	InputTemp[FNode.Water] = TW;
	InputTemp[FNode.Amb] = TAMB;

	for(int i = 0; i < FNumberNodes; i++) {
		InputTemp[FNumberNodes + i] = FNode_Temp_K_Tr[i];
	}

	LUdcmp TempSolver(FMatrix_KC);

	TempSolver.solve(InputTemp, FNode_Temp_K_Tr);

	for(int i = 0; i < FNumberNodes; i++) {
		FNode_Temp_K[i] = FNode_Temp_K_Tr[i];
	}

}

void TTC_HTM::Read_HTM(FILE *fich) {

	double tmp = 0.;
	dVector vtmp;
	int np = 0;
	int water = 0;

	fscanf(fich, "%d ", &water);
	if(water == 1) {
		FIsWaterCooled = true;
		FMatrix_Conn[FNode.H2][FNode.Water] = true;
	} else {
		FMatrix_Conn[FNode.H3][FNode.AirOut] = true;
	}

	fscanf(fich, "%lf ", &FCond.T_H1);
	fscanf(fich, "%lf ", &FCond.H1_H2);
	fscanf(fich, "%lf ", &FCond.H2_H3);
	fscanf(fich, "%lf ", &FCond.H3_C);

	fscanf(fich, "%lf ", &FCap.T);
	fscanf(fich, "%lf ", &FCap.H1);
	fscanf(fich, "%lf ", &FCap.H2);
	fscanf(fich, "%lf ", &FCap.H3);
	fscanf(fich, "%lf ", &FCap.C);

	/** Correlation GAS - T * */
	np = 4;
	FConv.GAS_T.Coef.resize(np);
	FConv.GAS_T.Ind.resize(np - 1);
	FConv.GAS_T.Ind[0] = nmReMassTur;
	FConv.GAS_T.Ind[1] = nmPrTur;
	FConv.GAS_T.Ind[2] = nmmuG_T;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.GAS_T.Coef[i]);
	}

	if(FIsWaterCooled) {
		/** Correlation Water - H2 * */
		np = 3;
		FConv.WAT_H2.Coef.resize(np);
		FConv.WAT_H2.Ind.resize(np - 1);
		FConv.WAT_H2.Ind[0] = nmReWater;
		FConv.WAT_H2.Ind[1] = nmPrWater;

		for(int i = 0; i < np; i++) {
			fscanf(fich, "%lf ", &FConv.WAT_H2.Coef[i]);
		}
	}

	/** Correlation Oil - H1 * */
	np = 5;
	FConv.OIL_H1.Coef.resize(np);
	FConv.OIL_H1.Ind.resize(np - 1);
	FConv.OIL_H1.Ind[0] = nmReOilH1;
	FConv.OIL_H1.Ind[1] = nmPrOilH1;
	FConv.OIL_H1.Ind[2] = nmmuO_H1;
	FConv.OIL_H1.Ind[3] = nmReShaft;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.OIL_H1.Coef[i]);
	}

	/** Corelation Oil - H2 * */
	np = 4;
	FConv.OIL_H2_C.Coef.resize(np);
	FConv.OIL_H2_C.Ind.resize(np - 1);
	FConv.OIL_H2_C.Ind[0] = nmReOilH2;
	FConv.OIL_H2_C.Ind[1] = nmPrOilH2;
	FConv.OIL_H2_C.Ind[2] = nmmuO_H2;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.OIL_H2_C.Coef[i]);
	}

	np = 4;
	FConv.OIL_H2_H.Coef.resize(np);
	FConv.OIL_H2_H.Ind.resize(np - 1);
	FConv.OIL_H2_H.Ind[0] = nmReOilH2;
	FConv.OIL_H2_H.Ind[1] = nmPrOilH2;
	FConv.OIL_H2_H.Ind[2] = nmmuO_H2;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.OIL_H2_H.Coef[i]);
	}

	/** Correlation Oil - H3 * */
	// np = 3;
	// FConv.OIL_H3.Coef.resize(np);
	// FConv.OIL_H3.Ind.resize(np - 1);
	// FConv.OIL_H3.Ind[0] = nmReOilH3;
	// FConv.OIL_H3.Ind[1] = nmPrOilH3;
	//
	// for (int i = 0; i < np; i++) {
	// fscanf(fich, "%lf ", &FConv.OIL_H3.Coef[i]);
	// }
	if(!FIsWaterCooled) {

		/** Correlation Air Out - H3 * */
		np = 3;
		FConv.AIR_H3_C.Coef.resize(np);
		FConv.AIR_H3_C.Ind.resize(np - 1);
		FConv.AIR_H3_C.Ind[0] = nmReMassCom2;
		FConv.AIR_H3_C.Ind[1] = nmPrCom2;

		for(int i = 0; i < np; i++) {
			fscanf(fich, "%lf ", &FConv.AIR_H3_C.Coef[i]);
		}

		np = 3;
		FConv.AIR_H3_H.Coef.resize(np);
		FConv.AIR_H3_H.Ind.resize(np - 1);
		FConv.AIR_H3_H.Ind[0] = nmReMassCom2;
		FConv.AIR_H3_H.Ind[1] = nmPrCom2;

		for(int i = 0; i < np; i++) {
			fscanf(fich, "%lf ", &FConv.AIR_H3_H.Coef[i]);
		}
	}

	np = 3;
	FConv.AIR_C_C.Coef.resize(np);
	FConv.AIR_C_C.Ind.resize(np - 1);
	FConv.AIR_C_C.Ind[0] = nmReMassCom;
	FConv.AIR_C_C.Ind[1] = nmPrCom;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.AIR_C_C.Coef[i]);
	}

	np = 3;
	FConv.AIR_C_H.Coef.resize(np);
	FConv.AIR_C_H.Ind.resize(np - 1);
	FConv.AIR_C_H.Ind[0] = nmReMassCom;
	FConv.AIR_C_H.Ind[1] = nmPrCom;

	for(int i = 0; i < np; i++) {
		fscanf(fich, "%lf ", &FConv.AIR_C_H.Coef[i]);
	}

	fscanf(fich, "%lf ", &FFitRadiation);
	if(FFitRadiation == 0.) {
		FExternalRadiation = false;
	} else {
		FExternalRadiation = true;
		fscanf(fich, "%lf %lf %lf", &FEmisExtTurb, &FEmisExtHous, &FEmisExtComp);
		int Shield = 0;
		fscanf(fich, "%d ", &Shield);
		if(Shield == 0) {
			FTurbineShielded = false;
		} else {
			FTurbineShielded = true;
			fscanf(fich, "%lf ", &FEmisExtShie);
		}
	}

	fscanf(fich, "%lf ", &FFitConvection);
	if(FFitConvection == 0.) {
		FExternalConvection = false;
	} else {
		FExternalConvection = true;
		fscanf(fich, "%lf ", &FExtVelocity);
	}

}

void TTC_HTM::Read_HTMXML(xml_node node_ht) {

	double tmp = 0.;
	dVector vtmp;
	int np = 0;
	int water = 0;

	xml_node node_set = GetNodeChild(node_ht, "Htr_Settings");

	FIsWaterCooled = GetAttributeAsBool(node_set, "WaterCooled");

	if(FIsWaterCooled) {
		FMatrix_Conn[FNode.H2][FNode.Water] = true;
	} else {
		FMatrix_Conn[FNode.H3][FNode.AirOut] = true;
	}

	xml_node node_cond = GetNodeChild(node_ht, "Htr_Conduction");

	FCond.T_H1 = GetAttributeAsDouble(node_cond, "T_H1");
	FCond.H1_H2 = GetAttributeAsDouble(node_cond, "H1_H2");
	FCond.H2_H3 = GetAttributeAsDouble(node_cond, "H2_H3");
	FCond.H3_C = GetAttributeAsDouble(node_cond, "H3_C");

	xml_node node_cap = GetNodeChild(node_ht, "Htr_Capacity");

	FCap.T = GetAttributeAsDouble(node_cap, "T");
	FCap.H1 = GetAttributeAsDouble(node_cap, "H1");
	FCap.H2 = GetAttributeAsDouble(node_cap, "H2");
	FCap.H3 = GetAttributeAsDouble(node_cap, "H3");
	FCap.C = GetAttributeAsDouble(node_cap, "C");

	xml_node node_iconv = GetNodeChild(node_ht, "Htr_Convection");

	/** Correlation GAS - T * */
	np = 4;
	FConv.GAS_T.Coef.resize(np);
	FConv.GAS_T.Ind.resize(np - 1);
	FConv.GAS_T.Ind[0] = nmReMassTur;
	FConv.GAS_T.Ind[1] = nmPrTur;
	FConv.GAS_T.Ind[2] = nmmuG_T;

	xml_node node_gast = GetNodeChild(node_iconv, "Cnv_Gas_T");

	FConv.GAS_T.Coef[0] = GetAttributeAsDouble(node_gast, "k1");
	FConv.GAS_T.Coef[1] = GetAttributeAsDouble(node_gast, "k2");
	FConv.GAS_T.Coef[2] = GetAttributeAsDouble(node_gast, "k3");
	FConv.GAS_T.Coef[3] = GetAttributeAsDouble(node_gast, "k4");

	if(FIsWaterCooled) {
		/** Correlation Water - H2 * */
		np = 3;
		FConv.WAT_H2.Coef.resize(np);
		FConv.WAT_H2.Ind.resize(np - 1);
		FConv.WAT_H2.Ind[0] = nmReWater;
		FConv.WAT_H2.Ind[1] = nmPrWater;

		xml_node node_wath2 = GetNodeChild(node_iconv, "Cnv_Wat_H2");

		FConv.WAT_H2.Coef[0] = GetAttributeAsDouble(node_wath2, "k1");
		FConv.WAT_H2.Coef[1] = GetAttributeAsDouble(node_wath2, "k2");
		FConv.WAT_H2.Coef[2] = GetAttributeAsDouble(node_wath2, "k3");
	}

	/** Correlation Oil - H1 * */
	np = 5;
	FConv.OIL_H1.Coef.resize(np);
	FConv.OIL_H1.Ind.resize(np - 1);
	FConv.OIL_H1.Ind[0] = nmReOilH1;
	FConv.OIL_H1.Ind[1] = nmPrOilH1;
	FConv.OIL_H1.Ind[2] = nmmuO_H1;
	FConv.OIL_H1.Ind[3] = nmReShaft;

	xml_node node_oilh1 = GetNodeChild(node_iconv, "Cnv_Oil_H1");

	FConv.OIL_H1.Coef[0] = GetAttributeAsDouble(node_oilh1, "k1");
	FConv.OIL_H1.Coef[1] = GetAttributeAsDouble(node_oilh1, "k2");
	FConv.OIL_H1.Coef[2] = GetAttributeAsDouble(node_oilh1, "k3");
	FConv.OIL_H1.Coef[3] = GetAttributeAsDouble(node_oilh1, "k4");
	FConv.OIL_H1.Coef[4] = GetAttributeAsDouble(node_oilh1, "k5");

	/** Corelation Oil - H2 * */
	np = 4;
	FConv.OIL_H2_C.Coef.resize(np);
	FConv.OIL_H2_C.Ind.resize(np - 1);
	FConv.OIL_H2_C.Ind[0] = nmReOilH2;
	FConv.OIL_H2_C.Ind[1] = nmPrOilH2;
	FConv.OIL_H2_C.Ind[2] = nmmuO_H2;

	xml_node node_oilh2c = GetNodeChild(node_iconv, "Cnv_Oil_H2_C");

	FConv.OIL_H2_C.Coef[0] = GetAttributeAsDouble(node_oilh2c, "k1");
	FConv.OIL_H2_C.Coef[1] = GetAttributeAsDouble(node_oilh2c, "k2");
	FConv.OIL_H2_C.Coef[2] = GetAttributeAsDouble(node_oilh2c, "k3");
	FConv.OIL_H2_C.Coef[3] = GetAttributeAsDouble(node_oilh2c, "k4");

	np = 4;
	FConv.OIL_H2_H.Coef.resize(np);
	FConv.OIL_H2_H.Ind.resize(np - 1);
	FConv.OIL_H2_H.Ind[0] = nmReOilH2;
	FConv.OIL_H2_H.Ind[1] = nmPrOilH2;
	FConv.OIL_H2_H.Ind[2] = nmmuO_H2;

	xml_node node_oilh2h = GetNodeChild(node_iconv, "Cnv_Oil_H2_H");

	FConv.OIL_H2_H.Coef[0] = GetAttributeAsDouble(node_oilh2h, "k1");
	FConv.OIL_H2_H.Coef[1] = GetAttributeAsDouble(node_oilh2h, "k2");
	FConv.OIL_H2_H.Coef[2] = GetAttributeAsDouble(node_oilh2h, "k3");
	FConv.OIL_H2_H.Coef[3] = GetAttributeAsDouble(node_oilh2h, "k4");

	if(!FIsWaterCooled) {

		/** Correlation Air Out - H3 * */
		np = 3;
		FConv.AIR_H3_C.Coef.resize(np);
		FConv.AIR_H3_C.Ind.resize(np - 1);
		FConv.AIR_H3_C.Ind[0] = nmReMassCom2;
		FConv.AIR_H3_C.Ind[1] = nmPrCom2;

		xml_node node_airh3c = GetNodeChild(node_iconv, "Cnv_Air_H3_C");

		FConv.AIR_H3_C.Coef[0] = GetAttributeAsDouble(node_airh3c, "k1");
		FConv.AIR_H3_C.Coef[1] = GetAttributeAsDouble(node_airh3c, "k2");
		FConv.AIR_H3_C.Coef[2] = GetAttributeAsDouble(node_airh3c, "k3");

		np = 3;
		FConv.AIR_H3_H.Coef.resize(np);
		FConv.AIR_H3_H.Ind.resize(np - 1);
		FConv.AIR_H3_H.Ind[0] = nmReMassCom2;
		FConv.AIR_H3_H.Ind[1] = nmPrCom2;

		xml_node node_airh3h = GetNodeChild(node_iconv, "Cnv_Air_H3_H");

		FConv.AIR_H3_H.Coef[0] = GetAttributeAsDouble(node_airh3h, "k1");
		FConv.AIR_H3_H.Coef[1] = GetAttributeAsDouble(node_airh3h, "k2");
		FConv.AIR_H3_H.Coef[2] = GetAttributeAsDouble(node_airh3h, "k3");

	}

	np = 3;
	FConv.AIR_C_C.Coef.resize(np);
	FConv.AIR_C_C.Ind.resize(np - 1);
	FConv.AIR_C_C.Ind[0] = nmReMassCom;
	FConv.AIR_C_C.Ind[1] = nmPrCom;

	xml_node node_aircc = GetNodeChild(node_iconv, "Cnv_Air_C_C");

	FConv.AIR_C_C.Coef[0] = GetAttributeAsDouble(node_aircc, "k1");
	FConv.AIR_C_C.Coef[1] = GetAttributeAsDouble(node_aircc, "k2");
	FConv.AIR_C_C.Coef[2] = GetAttributeAsDouble(node_aircc, "k3");

	np = 3;
	FConv.AIR_C_H.Coef.resize(np);
	FConv.AIR_C_H.Ind.resize(np - 1);
	FConv.AIR_C_H.Ind[0] = nmReMassCom;
	FConv.AIR_C_H.Ind[1] = nmPrCom;

	xml_node node_airch = GetNodeChild(node_iconv, "Cnv_Air_C_H");

	FConv.AIR_C_H.Coef[0] = GetAttributeAsDouble(node_airch, "k1");
	FConv.AIR_C_H.Coef[1] = GetAttributeAsDouble(node_airch, "k2");
	FConv.AIR_C_H.Coef[2] = GetAttributeAsDouble(node_airch, "k3");

	xml_node node_rad = GetNodeChild(node_ht, "Htr_Radiation");

	FExternalRadiation = GetAttributeAsBool(node_rad, "Radiation");

	if(FExternalRadiation) {
		FEmisExtTurb = GetAttributeAsDouble(node_rad, "TurbineEmissivity");
		FEmisExtHous = GetAttributeAsDouble(node_rad, "HousingEmissivity");
		FEmisExtComp = GetAttributeAsDouble(node_rad, "CompressorEmissivity");

		FTurbineShielded = GetAttributeAsBool(node_rad, "Shield");
		if(FTurbineShielded)
			FEmisExtShie = GetAttributeAsDouble(node_rad, "ShieldEmissivity");
	}

	xml_node node_econv = GetNodeChild(node_ht, "Htr_ConvectionExt");

	FFitConvection = GetAttributeAsBool(node_econv, "ExtConvection");

	if(FFitConvection)
		FExtVelocity = GetAttributeAsDouble(node_econv, "ExtVelocity");

}

double TTC_HTM::CorrectCompressorMap(double m, double cr, double eff, double TinC, double TinT, double Rtc, double mt) {

	double EffT = 0.6;
	if(mt == 0)
		mt = m;
	// double OilPower = FMechLosses->Value();

	CompressorWorkingPoint(0, m, TinC, 1.0, cr, eff);

	double PowerC = FC.funIsoPower(1.0) / eff;

	TurbineWorkingPoint(0, 0, mt, TinT, 1.0, 1.0, EffT);

	double Cpgas = FT.Fluid->Cp;
	double g = FT.Fluid->g;
	double gam = (g - 1) / g;

	double er = pow(1 - PowerC / (EffT * mt * Cpgas * (TinT + 273)), -1 / gam);
	double Lim1 = er - 0.1;
	if(Lim1 < 1)
		Lim1 = 1;
	double Lim2 = er + 0.1;

	stTurbinePower TurbinePower(PowerC, mt, EffT, TinT + 273, Cpgas, gam);

	if(zbrac(TurbinePower, Lim1, Lim2)) {
		er = FindRoot(TurbinePower, Lim1, Lim2);
	}

	TurbineWorkingPoint(0, 0, mt, TinT, er, er, EffT);

	double OilPower = FMechLosses->P_oil(FO_IT_K, Pi * Rtc / 30., 1.0, cr, er, 1.0, FO_MassFlow);
	double EffM = PowerC / (PowerC + OilPower);

	TurbochargerWorkingPoint(Rtc, EffM, FO_MassFlow, FO_IT_K, FO_IP, FW_IT_K, FW_MassFlow);

	double effadi = AdiabaticEff(nmCompressor);

	//printf("%lf %lf %lf %lf %lf %lf\n", Rtc, m, cr, eff, effadi,
	//	Comp_Heat_Flow());

	FILE *fres = fopen("CompRes.dat", "a");

	fprintf(fres, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", Rtc, m, cr, effadi, eff, Turb_Heat_Flow(),
			Comp_Heat_Flow(), Oil_Heat_Flow(), Water_Heat_Flow(), Ambient_Heat_Flow());

	for(int i = 0; i < FNumberNodes; i++) {
		fprintf(fres, "\t%lf", FNode_Temp_K[i]);
	}

	fprintf(fres, "\t%lf\t%lf\t%lf\t%lf", FC.Power, FT.Power, FMechEff, mt * FT.Fluid->fun_Cp(TinT) * TinT);

	fprintf(fres, "\n");
	fclose(fres);

	return effadi;

}

double TTC_HTM::CorrectTurbineMap(double m, double er, double eff, double TinC, double TinT, double Rtc, double mc) {

	double EffC = 0.6;
	// double EffM = FMechLosses->Value();
	if(mc == 0)
		mc = m;

	TurbineWorkingPoint(0, 0, m, TinT, er, er, eff);

	double PowerT = FT.funIsoPower(-1) * eff;

	CompressorWorkingPoint(0, mc, TinC, 1.0, 1.0, EffC);

	double Cpgas = FC.Fluid->Cp;
	double g = FC.Fluid->g;
	double gam = (g - 1) / g;

	double cr = pow(PowerT * EffC / (mc * Cpgas * (TinC)) + 1, 1 / gam);

	double Lim1 = cr - 0.1;
	if(Lim1 < 1)
		Lim1 = 1;
	double Lim2 = cr + 0.1;

	stCompressorPower CompressorPower(PowerT, mc, EffC, TinC, Cpgas, gam);

	if(zbrac(CompressorPower, Lim1, Lim2)) {
		cr = FindRoot(CompressorPower, Lim1, Lim2);
	}

	CompressorWorkingPoint(0, mc, TinC, 1.0, cr, EffC);

	double OilPower = FMechLosses->P_oil(FO_IT_K, Pi * Rtc / 30., 1.0, cr, er, 1.0, FO_MassFlow);
	double EffM = 1 - OilPower / PowerT;

	if(EffM < 0.1) {
		EffM = 0.1;
	}

	TurbochargerWorkingPoint(Rtc, EffM, FO_MassFlow, FO_IT_K, FO_IP, FW_IT_K, FW_MassFlow);

	double effadi = AdiabaticEff(nmTurbine);

	//printf("%lf %lf %lf %lf %lf\n", Rtc, m, er, eff, effadi);

	FILE *fres = fopen("TurbRes.dat", "a");

	fprintf(fres, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", Rtc, m, er, effadi, eff, Turb_Heat_Flow(),
			Comp_Heat_Flow(), Oil_Heat_Flow(), Water_Heat_Flow(), Ambient_Heat_Flow());

	for(int i = 0; i < FNumberNodes; i++) {
		fprintf(fres, "\t%lf", FNode_Temp_K[i]);
	}

	fprintf(fres, "\t%lf\t%lf\t%lf\t%lf", FC.Power, FT.Power, FMechEff, m * FT.Fluid->fun_Cp(TinT) * TinT);

	fprintf(fres, "\n");
	fclose(fres);

	return effadi;

}

void TTC_HTM::InitializeTemp(double TIT, double COT, double CIT, double OIT, double WIT, double TAMB) {
	FNode_Temp_K[FNode.Gas] = TIT;
	FNode_Temp_K[FNode.Air] = COT;
	FNode_Temp_K[FNode.Water] = WIT;
	FNode_Temp_K[FNode.T] = 0.80 * (TIT - COT) + COT;
	FNode_Temp_K[FNode.H1] = 0.50 * (TIT - COT) + COT;
	FNode_Temp_K[FNode.H2] = 0.15 * (TIT - COT) + COT;
	FNode_Temp_K[FNode.H3] = 0.10 * (TIT - COT) + COT;
	FNode_Temp_K[FNode.C] = 0.05 * (TIT - COT) + COT;
	FNode_Temp_K[FNode.O_H1] = OIT;
	FNode_Temp_K[FNode.O_H2] = OIT;
	FNode_Temp_K[FNode.O_H3] = OIT;
	FNode_Temp_K[FNode.AirIn] = CIT;
	FNode_Temp_K[FNode.AirOut] = COT;
	FNode_Temp_K[FNode.Amb] = TAMB;

	FNode_Temp_K_Tr[FNode.Gas] = TIT;
	FNode_Temp_K_Tr[FNode.Air] = COT;
	FNode_Temp_K_Tr[FNode.Water] = WIT;
	FNode_Temp_K_Tr[FNode.T] = 0.80 * (TIT - COT) + COT;
	FNode_Temp_K_Tr[FNode.H1] = 0.50 * (TIT - COT) + COT;
	FNode_Temp_K_Tr[FNode.H2] = 0.15 * (TIT - COT) + COT;
	FNode_Temp_K_Tr[FNode.H3] = 0.10 * (TIT - COT) + COT;
	FNode_Temp_K_Tr[FNode.C] = 0.05 * (TIT - COT) + COT;
	FNode_Temp_K_Tr[FNode.O_H1] = OIT;
	FNode_Temp_K_Tr[FNode.O_H2] = OIT;
	FNode_Temp_K_Tr[FNode.O_H3] = OIT;
	FNode_Temp_K_Tr[FNode.AirIn] = CIT;
	FNode_Temp_K_Tr[FNode.AirOut] = COT;
	FNode_Temp_K_Tr[FNode.Amb] = TAMB;

}

void TTC_HTM::PrintInsTemperatures(stringstream & insoutput) {

	for(int i = 0; i < FNumberNodes; i++) {
		insoutput << "\t" << FNode_Temp_K_Tr[i] + unKToC;
	}

}

void TTC_HTM::HeaderInsTemperatures(stringstream & insoutput, int num) {

	std::string Label;

	for(int i = 0; i < FNumberNodes; i++) {
		Label = "\t" + PutLabel(5007) + "/" + std::to_string(num) + "/" + PutLabel(4005) + PutLabel(910) + "/" + FNodeLabel[i];
		insoutput << Label.c_str();
	}

}

void TTC_HTM::PrintInsHeatFlow(stringstream & insoutput) {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			if(FMatrix_Conn[i][j])
				insoutput << "\t" << FMatrix_HF[i][j];
		}
	}

}

void TTC_HTM::HeaderInsHeatFlow(stringstream & insoutput, int num) {

	std::string Label;

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			if(FMatrix_Conn[i][j]) {
				Label = "\t" + PutLabel(5007) + "/" + std::to_string(num) + "/" + PutLabel(4010) + PutLabel(
							903) + "/" + "From_" + FNodeLabel[j] + "_To_" + FNodeLabel[i];
				insoutput << Label.c_str();
			}
		}
	}

}

void TTC_HTM::View_Factors(double Dt, double Dc, double Dh, double L) {

	double A1 = FAlpha1;
	double A2 = FAlpha2;
	double A3 = FAlpha3;

	double Dt2 = pow2(Dt);
	double Dc2 = pow2(Dc);
	double Dh2 = pow2(Dh);

	double DtDh = Dt2 - Dh2;
	double DcDh = Dc2 - Dh2;
	double LDh4 = L * Dh * 4;

	FVF.C_T = ViewFactorDiskDisk(Dc / 2, Dt / 2, Dh / 2, L);
	FVF.T_C = FVF.C_T * DcDh / DtDh;
	FVF.H1_T = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L * A1);
	FVF.T_H1 = FVF.H1_T * LDh4 * A1 / DtDh;
	FVF.T_H2 = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L * (A1 + A2)) * (A1 + A2) * LDh4 / DtDh - FVF.T_H1;
	FVF.H2_T = FVF.T_H2 * DtDh / A2 / LDh4;
	FVF.T_H3 = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L) * LDh4 / DtDh - FVF.T_H1 - FVF.T_H2;
	FVF.H3_T = FVF.T_H3 * DtDh / LDh4 / A3;
	FVF.T_Amb = 1 - (FVF.T_C + FVF.T_H1 + FVF.T_H2 + FVF.T_H3);
	FVF.H3_C = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L * A3);
	FVF.C_H3 = FVF.H3_C * LDh4 * A3 / DcDh;
	FVF.C_H2 = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L * (A3 + A2)) * (A3 + A2) * LDh4 / DcDh - FVF.C_H3;
	FVF.H2_C = FVF.C_H2 * DcDh / LDh4 / A2;
	FVF.C_H1 = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L) * LDh4 / DcDh - FVF.C_H3 - FVF.C_H2;
	FVF.H1_C = FVF.C_H1 * DcDh / LDh4 / A1;
	FVF.C_Amb = 1 - (FVF.C_H1 + FVF.C_H2 + FVF.C_H3 + FVF.C_T);
	FVF.H1_Amb = 1 - (FVF.H1_T + FVF.H1_C);
	FVF.H2_Amb = 1 - (FVF.H2_T + FVF.H2_C);
	FVF.H3_Amb = 1 - (FVF.H3_T + FVF.H3_C);

	if(FTurbineShielded) {
		FVF.T_C = FVF.T_C * FEmisExtShie / (2 * (FEmisExtShie + FVF.T_C * (1 - FEmisExtShie)));
		FVF.C_T = FVF.T_C * DtDh / DcDh;
		FVF.T_H1 = FVF.T_H1 * FEmisExtShie / (2 * (FEmisExtShie + FVF.T_H1 * (1 - FEmisExtShie)));
		FVF.H1_T = FVF.T_H1 / (LDh4 * A1 / DtDh);
		FVF.T_H2 = FVF.T_H2 * FEmisExtShie / (2 * (FEmisExtShie + FVF.T_H2 * (1 - FEmisExtShie)));
		FVF.H2_T = FVF.T_H2 / (LDh4 * A2 / DtDh);
		FVF.T_H3 = FVF.T_H3 * FEmisExtShie / (2 * (FEmisExtShie + FVF.T_H3 * (1 - FEmisExtShie)));
		FVF.H3_T = FVF.T_H3 / (LDh4 * A3 / DtDh);
		FVF.T_Amb = FVF.T_Amb * FEmisExtShie / (FEmisExtShie + 2 * FVF.T_H3 * (1 - FEmisExtShie));
		;

	}

}

double TTC_HTM::ViewFactorDiskDisk(double r1, double r2, double rc, double h) {

	double R1 = r1 / h;
	double R2 = r2 / h;
	double Rc = rc / h;
	double A = pow2(R1) - pow2(Rc);
	double B = pow2(R2) - pow2(Rc);
	double C = R1 + R2;
	double D = R2 - R1;
	double Y = sqrt(A) + sqrt(B);

	double VF1 = A / 2 * acos(Rc / R2) + B / 2 * acos(Rc / R1) + 2 * Rc * (atan(Y) - atan(sqrt(A)) - atan(sqrt(B)));
	double VF2 = sqrt((1 + pow2(C)) * (1 + pow2(D))) * atan(sqrt((1 + pow2(C)) * (pow2(Y) - pow2(D)) / (1 + pow2(D)) /
				 (pow2(C) - pow2(Y))));
	double VF3 = sqrt((1 + pow2(R1 + Rc)) * (1 + pow2(R1 - Rc))) * atan(sqrt((1 + pow2(R1 + Rc)) * (R1 - Rc) / (1 + pow2(
					 R1 - Rc)) / (R1 + Rc)));
	double VF4 = sqrt((1 + pow2(R2 + Rc)) * (1 + pow2(R2 - Rc))) * atan(sqrt((1 + pow2(R2 + Rc)) * (R2 - Rc) / (1 + pow2(
					 R2 - Rc)) / (R2 + Rc)));

	return (VF1 - VF2 + VF3 + VF4) / A / Pi;

}

double TTC_HTM::ViewFactorDiskCylinder(double r1, double r2, double h) {

	double R = r1 / r2;
	double H = h / r2;
	double A = pow2(H) + pow2(R) - 1;
	double B = pow2(H) - pow2(R) + 1;

	return B / 8 / R / H + 1 / (2 * Pi) * (acos(A / B) - 1 / (2 * H) * sqrt(pow2(A + 2) / pow2(R) - 4) * acos(
			A * R / B) - A / 2 / R / H * asin(R));

}

double TTC_HTM::KRadiationSurfaces(double A1, double A2, double e1, double e2, double F12, double T1, double T2) {

	double R1 = (1 - e1) / A1 / e1;
	double R2 = (1 - e2) / A2 / e2;
	double R12 = 1 / A1 / F12;

	double K = StephanBoltzmann * (pow2(T1) + pow2(T2)) * (T1 + T2) / (R1 + R12 + R2);

	return K;

}

double TTC_HTM::KRadiationExternal(double L, double D, double e, double F, double T, double Tamb) {

	return (D * (1 + F / (e * (1 - F) + F)) + 4 * L) * Pi / 4 * D * e * StephanBoltzmann * (pow2(T) + pow2(Tamb)) *
		   (T + Tamb);

}

double TTC_HTM::KRadiationExternalHous(double A, double e, double F, double T, double Tamb) {

	return e * F / (e * (1 - F) + F) * A * StephanBoltzmann * (pow2(T) + pow2(Tamb)) * (T + Tamb);

}

double TTC_HTM::NusseltFreeConv(double Gr, double Pr) {

	double Ra = Gr * Pr;

	return pow2(0.6 + 0.387 * pow(Ra, 1 / 6) / pow(1 + 0.559 / pow(Pr, 9 / 16), 8 / 27));

}

double TTC_HTM::NusseltForcConv(double Re, double Pr) {

	return 0.3 + 0.62 * sqrt(Re) * pow(Pr, 1 / 3) / pow(1 + 0.4 / pow(Pr, 2 / 3), 1 / 4) * pow(1 + pow(Re / 282000, 5 / 8),
			4 / 5);

}

void TTC_HTM::AddKConvRad(dMatrix& K, dMatrix K2) {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			K[i][j] += K2[i][j];
		}
	}
}

void TTC_HTM::PutHeatC(double HeatC1, double HeatC2) {
	FMatrix_HF[FNode.H3][FNode.AirOut] = HeatC2;
	FMatrix_HF[FNode.C][FNode.Air] = HeatC1;
}

void TTC_HTM::PutEffTurbMax(double EffMax) {
	FT.EFFMAX = EffMax;
}

void TTC_HTM::PrintTemperatures() {
	cout << "%%%% TEMPERATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
	for(int i = 0; i < FNumberNodes; i++) {
		cout << "Temperature node " << FNodeLabel[i] << ":" << "\t" << FNode_Temp_K_Tr[i] << "C" << endl;
	}
	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
}

void TTC_HTM::InitializeHeatFlow() {

	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			FMatrix_HF[i][j] = 0.;
		}
	}
}

void TTC_HTM::PrintHeatFlows() {
	FILE *fheat = fopen("heat.txt", "w");
	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			fprintf(fheat, "%g\t", FMatrix_HF[i][j]);
		}
		fprintf(fheat, "\n");
	}
	fclose(fheat);
}

void TTC_HTM::PrintKs() {
	FILE *fheat = fopen("ks.txt", "w");
	for(int i = 0; i < FNumberNodes; i++) {
		for(int j = 0; j < FNumberNodes; j++) {
			fprintf(fheat, "%g\t", FMatrix_KS[i][j]);
		}
		fprintf(fheat, "\n");
	}
	fclose(fheat);
}

void TTC_HTM::PrintTsTr() {
	FILE *fheat = fopen("TempsTr.txt", "w");
	for(int i = 0; i < FNumberNodes; i++) {
		fprintf(fheat, "%g\t", FNode_Temp_K_Tr[i]);
		fprintf(fheat, "\n");
	}
	fclose(fheat);
}

void TTC_HTM::PrintTs() {
	FILE *fheat = fopen("Temps.txt", "w");
	for(int i = 0; i < FNumberNodes; i++) {
		fprintf(fheat, "%g\t", FNode_Temp_K[i]);
		fprintf(fheat, "\n");
	}
	fclose(fheat);
}

#pragma package(smart_init)
