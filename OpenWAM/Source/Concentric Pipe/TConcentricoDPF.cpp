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
#pragma hdrstop

#include "TConcentricoDPF.h"
#include "TTubo.h"
#include "TBloqueMotor.h"
#include "TCondicionContorno.h"
#include "TCCUnionEntreTubos.h"
#include "TDPF.h"

#ifdef ParticulateFilter
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentricoDPF::TConcentricoDPF(int NumeroConcentrico) :
	TConcentrico() {

	FNumeroConcentrico = NumeroConcentrico;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentricoDPF::~TConcentricoDPF() {

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoDPF::LeeDatosTuboConcentrico(const char *FileWAM, fpos_t &filepos, TTubo **Tubo, TDPF **DPF) {
	try {

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d %d ", &FNumDPFInterna, &FNumTuboExterno);

		FHayDPF = true;
		FNumTuboInterno = 0;

		fgetpos(fich, &filepos);
		fclose(fich);

		FDPF = new TDPF*[1];
		FDPF[0] = DPF[FNumDPFInterna - 1];
		FDPF[0]->PutConcentric(1);
		FNumeroTubos = 1;

		FTg = new double[FNumeroTubos];
		FTpantpos = new double[FNumeroTubos];
		FTpant0 = new double[FNumeroTubos];
		FTpant1 = new double[FNumeroTubos];
		FTpant2 = new double[FNumeroTubos];
		FTpantant = new double[FNumeroTubos];

		FTubo = new TTubo*[FNumeroTubos];

		FTubo[0] = Tubo[FNumTuboExterno - 1];
		FTubo[0]->PutConcentric(1);

		for(int j = 0; j < FNumeroTubos; j++) {
			FTg[j] = FTubo[j]->getTemperaturaInicial();
			FTpantant[j] = FTubo[j]->getTempWallIni(); // Valor inicial de temperatura de pared tubo j
			FTpantpos[j] = __units::degCToK(FTubo[j]->getTempWallIni()); // Valor inicial de temperatura de pared tubo j
			FTpant0[j] = __units::degCToK(FTubo[j]->getTempWallIni());
			FTpant1[j] = __units::degCToK(FTubo[j]->getTempWallIni());
			FTpant2[j] = __units::degCToK(FTubo[j]->getTempWallIni());
		}

		FTParedAnt = new double**[FNumeroTubos];
		FTPared = new double**[FNumeroTubos];
		for(int j = 0; j < FNumeroTubos; j++) {
			FTParedAnt[j] = new double*[FNumeroNodosTemp];
			FTPared[j] = new double*[FNumeroNodosTemp];
			for(int i = 0; i < FNumeroNodosTemp; i++) {
				FTParedAnt[j][i] = new double[FTubo[j]->getNin()];
				FTPared[j][i] = new double[FTubo[j]->getNin()];
			}
		}

		for(int j = 0; j < FNumeroTubos; j++) {
			for(int i = 0; i < FNumeroNodosTemp; i++) {
				for(int k = 0; k < FTubo[j]->getNin(); k++) {
					FTParedAnt[j][i][k] = __units::degCToK(FTubo[j]->getTempWallIni());
					FTPared[j][i][k] = __units::degCToK(FTubo[j]->getTempWallIni());
				}
			}
		}

		FSUMTPTuboExtPro = new double**[2];
		for(int i = 0; i < 2; i++) {
			FSUMTPTuboExtPro[i] = new double*[2];
		}
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
				FSUMTPTuboExtPro[i][j] = new double[FTubo[0]->getNin()];
			}
		}

		for(int k = 0; k < FTubo[0]->getNin(); k++) {
			FSUMTPTuboExtPro[0][0][k] = 0.;
			FSUMTPTuboExtPro[0][1][k] = 0.;
			FSUMTPTuboExtPro[1][0][k] = 0.;
			FSUMTPTuboExtPro[1][1][k] = 0.;
		}

		FRg_int = new double[FTubo[0]->getNin()];
		FRg_int_ext = new double[FTubo[0]->getNin()];
		FRg_ext_int = new double[FTubo[0]->getNin()];
		FR_ext = new double[FTubo[0]->getNin()];
		FR_int_radiacion = new double[FTubo[0]->getNin()];
		FR_int_RadExt = new double[FTubo[0]->getNin()];
		FR_int_AxiAnt = new double[FTubo[0]->getNin()];
		FR_int_AxiPos = new double[FTubo[0]->getNin()];
		FR_int_RadInt = new double[FTubo[0]->getNin()];
		FR_ext_RadExt = new double[FTubo[0]->getNin()];
		FR_ext_AxiAnt = new double[FTubo[0]->getNin()];
		FR_ext_AxiPos = new double[FTubo[0]->getNin()];
		FR_ext_RadInt = new double[FTubo[0]->getNin()];

		FCapIntExt = new double[FTubo[0]->getNin()];
		FCapIntMed = new double[FTubo[0]->getNin()];
		FCapIntInt = new double[FTubo[0]->getNin()];
		FCapExtExt = new double[FTubo[0]->getNin()];
		FCapExtMed = new double[FTubo[0]->getNin()];
		FCapExtInt = new double[FTubo[0]->getNin()];

		FDuracionCiclo = FTubo[0]->getDuracionCiclo();
		FNumCiclosSinInerciaTermica = FTubo[0]->getNumCiclosSinInerciaTermica();

	} catch(exception &N) {
		cout << "ERROR: TConcentricoDPF::LeeDatosTuboConcentrico en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoDPF::CalculaTemperaturaPared(TBloqueMotor **Motor, double theta, TCondicionContorno **CC) {
	double Tg = 0.;
	double zzz = 0., czz = 0., cz1 = 0., uq1 = 0.;
	double DeltaTTPared = 0.;
	double R_Equiv_1 = 0., R_Equiv = 0., Text = 0., Ri = 0., Re = 0., ErrorTp = 0., Gamma = 0., RMezcla = 0., Asonido = 0.;
	bool EsPrimeraVez;
	int extremo = 0, nodo = 0;

	try {

		DeltaTTPared = FTubo[0]->getTime1() - FTubo[0]->getTime0();

		for(int j = 0; j < FNumeroTubos; j++) {    // Recorre los dos tubos conc�ntricos
			for(int i = 0; i < FTubo[j]->getNin(); i++) {
				FTParedAnt[j][0][i] = FTPared[j][0][i];
				FTParedAnt[j][1][i] = FTPared[j][1][i];
				FTParedAnt[j][2][i] = FTPared[j][2][i];

				zzz = 0.013 / DeltaTTPared;
				czz = 2 / (zzz + 1);
				uq1 = fabs(FTubo[j]->GetVelocidad(i) * __cons::ARef);
				if(FTubo[j]->getTipoTransCal() == nmTuboEscape)
					cz1 = czz;
				else
					cz1 = 1;
				FTubo[j]->PutVelPro(i, cz1 * uq1 + (1 - cz1) * FTubo[j]->GetVelPro(i));
			}
		}

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante
		   && FTubo[0]->getCoefAjustTC() != 0) {   // Tiene que existir motor.

			if(theta > FTubo[0]->getAnguloTotalCiclo()) {
				FSUMTime += DeltaTTPared;
			}
			for(int i = 0; i < FTubo[0]->getNin(); i++) {
				//Establece la temperatura del fluido en el exterior del tubo conc�ntrico.
				if(FTubo[0]->getTipoTransCal() == nmTuboAdmision || FTubo[0]->getTipoTransCal() == nmTuboEscape) {
					Text = FTubo[0]->getTExt();
				} else {
					Text = __units::degCToK(Motor[0]->getTempRefrigerante());
				}

				//Establece la temperatura de la pared en el instante de c�lculo anterior...
				//...en los nodos anterior y posterior.
				for(int j = 0; j < FNumeroTubos; j++) {
					//Establece la temperatura del gas en el interior del conducto.
					FTg[j] = pow(FTubo[j]->GetAsonido(i) * __cons::ARef, 2.) / (FTubo[j]->GetGamma(i) * FTubo[j]->GetRMezcla(i));

					if(i == 0) {
						if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
							if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoIzq() - 1])->getConductividad() > 0) {
								if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
									extremo = 1;
									if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								} else {
									extremo = 0;
									if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								}
								FTpantant[j] = __units::degCToK(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
							}
						}
						FTpantpos[j] = FTParedAnt[j][1][i + 1];
					} else if(i == FTubo[j]->getNin() - 1) {
						FTpantant[j] = FTParedAnt[j][1][i - 1];
						if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
							if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoDer() - 1])->getConductividad() > 0) {
								if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
									extremo = 1;
									if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								} else {
									extremo = 0;
									if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								}
								FTpantpos[j] = __units::degCToK(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
							}
						}
					} else {
						FTpantant[j] = FTParedAnt[j][1][i - 1];
						FTpantpos[j] = FTParedAnt[j][1][i + 1];
					}
					//...en los nodos interior, medio y exterior.
					FTpant0[j] = FTParedAnt[j][0][i];   // Nodo Interior
					FTpant1[j] = FTParedAnt[j][1][i];   // Nodo Medio
					FTpant2[j] = FTParedAnt[j][2][i];   // Nodo Exterior

					//Resistencias t�rmicas de convecci�n.
					if(FTubo[j]->Gethi(i) == 0) {
						FRg_int_ext[i] = 100000000.;
						FRg_ext_int[i] = 100000000.;
						if(FTpant0[0] == FDPF[0]->GetTSuperficie(i, 2)) {
							FR_int_radiacion[i] = 100000000.;
						} else {
							FR_int_radiacion[i] = (1. / FDPF[0]->getEmisividad() + FDPF[0]->getDiametroExt() / (FDPF[0]->getDiametroExt() + 2 *
												   FEspesor_fluido) * (1 / FTubo[j]->getEmisividad() - 1))
												  / (__cons::Sigma * __cons::Pi * FDPF[0]->getDiametroExt() * FTubo[j]->getXRef()) * (FDPF[0]->GetTSuperficie(i,
														  2) - FTpant0[0])
												  / (pow(FDPF[0]->GetTSuperficie(i, 2), 4) - pow(FTpant0[0], 4));
						}
					} else {
						FRg_int_ext[i] = 1 / __cons::Pi / FDPF[0]->getDiametroExt() / FTubo[j]->Gethi(i) / FTubo[j]->getXRef();
						FRg_ext_int[i] = 1 / __cons::Pi / (FDPF[0]->getDiametroExt() + 2 * FEspesor_fluido) / FTubo[j]->Gethi(
											 i) / FTubo[j]->getXRef();
						if(FTpant0[0] == FDPF[0]->GetTSuperficie(i, 2)) {
							FR_int_radiacion[i] = 100000000.;
						} else {
							FR_int_radiacion[i] = (1. / FDPF[0]->getEmisividad() + FDPF[0]->getDiametroExt() / (FDPF[0]->getDiametroExt() + 2 *
												   FEspesor_fluido) * (1 / FTubo[j]->getEmisividad() - 1))
												  / (__cons::Sigma * __cons::Pi * FDPF[0]->getDiametroExt() * FTubo[j]->getXRef()) * (FDPF[0]->GetTSuperficie(i,
														  2) - FTpant0[0])
												  / (pow(FDPF[0]->GetTSuperficie(i, 2), 4) - pow(FTpant0[0], 4));
						}
					}
					if(FTubo[j]->Gethe(i) == 0) {
						FR_ext[i] = 100000000.;
					} else {
						FR_ext[i] = 1 / __cons::Pi / (FDPF[0]->getDiametroExt() + 2 * (FEspesor_fluido + FEspesorTotalExt)) / FTubo[j]->Gethe(
										i) / FTubo[j]->getXRef();
					}
				}

				// C�lculo de las temperaturas de pared.
				// DPF wall temperature calculated at TDPF.cpp
				// Tubo exterior
				FTPared[0][2][i] = DeltaTTPared / FCapExtExt[i] * (1 / FR_ext[i] * (Text - FTpant2[0]) + 1 / FR_ext_RadExt[i] *
								   (FTpant1[0] - FTpant2[0])) + FTpant2[0];
				if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[0] - FTpant1[0])
										  + 1 / FR_ext_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_ext_AxiAnt[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_ext_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i] * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 /
									   FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0])) + FTpant1[0];
				}
				FTPared[0][0][i] = DeltaTTPared / FCapExtInt[i]
								   * (1 / FRg_ext_int[i] * (FTg[0] - FTpant0[0]) + 1 / FR_int_radiacion[i] * (FDPF[0]->GetTSuperficie(i,
										   2) - FTpant0[1]) + 1 / FR_ext_RadInt[i] * (FTpant1[0] - FTpant0[0]))
								   + FTpant0[0];

				for(int k = 0; k < 3; k++) {
					FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
				}
				// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
				if(FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica
				   || theta / FTubo[0]->getAnguloTotalCiclo() <= Motor[0]->getNumCiclosSinInerciaTermica()) {
					if(theta > FTubo[0]->getAnguloTotalCiclo()) {
						R_Equiv_1 = (FRg_ext_int[i] + FRg_int_ext[i]) * FR_int_radiacion[i] / (FRg_ext_int[i] + FRg_int_ext[i] +
									FR_int_radiacion[i]);
						R_Equiv = FR_ext_RadInt[i] + R_Equiv_1;
						// Tubo exterior
						// Sumatorio de R*Twall_DPF*incrt (para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][0][i] += __units::degCToK(FDPF[0]->GetTSuperficie(i, 2)) * DeltaTTPared / R_Equiv_1;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][0][i] += DeltaTTPared / R_Equiv_1;

						// Sumatorio de R*Twall_DPF*incrt (para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][1][i] += DeltaTTPared * __units::degCToK(FDPF[0]->GetTSuperficie(i, 2)) / R_Equiv;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][1][i] += DeltaTTPared / R_Equiv;
					}
				}

			}

			// Si est� al final del ciclo...
			if(FCicloTubo != Motor[0]->getCiclo() && FSUMTime > 0.) {
				// ...si (el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando) y est� en el segundo ciclo, calcula la temperatura de convergencia
				if((FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica
					|| theta / FTubo[0]->getAnguloTotalCiclo() <= Motor[0]->getNumCiclosSinInerciaTermica())
				   && theta > FTubo[0]->getAnguloTotalCiclo() + 1) {
					ErrorTp = 1.;
					EsPrimeraVez = true;
					while(ErrorTp >=
						  1) {   //Itera hasta conseguir una diferencia entre las temperaturas de pared menor a 1�C entre pasos.
						ErrorTp = 0.;
						for(int i = 0; i < FTubo[0]->getNin(); i++) {
							for(int j = 0; j < FNumeroTubos; j++) {
								//Establece la temperatura de la pared en el instante de c�lculo anterior...
								//...en los nodos anterior y posterior.
								FTpantant[j] = FTPared[j][1][i];
								FTpantpos[j] = FTPared[j][1][i];
								if(i == 0) {
									if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
										if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoIzq() - 1])->getConductividad() > 0) {
											if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
												extremo = 1;
												if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											} else {
												extremo = 0;
												if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											}
											FTpantant[j] = __units::degCToK(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
										}
									}
									FTpantpos[j] = FTPared[j][1][i + 1];
								} else if(i == FTubo[j]->getNin() - 1) {
									FTpantant[j] = FTPared[j][1][i - 1];
									if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
										if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoDer() - 1])->getConductividad() > 0) {
											if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
												extremo = 1;
												if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											} else {
												extremo = 0;
												if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											}
											FTpantpos[j] = __units::degCToK(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
										}
									}
								} else {
									FTpantant[j] = FTPared[j][1][i - 1];
									FTpantpos[j] = FTPared[j][1][i + 1];
								}
								//...en los nodos interior, medio y exterior.
								FTpant0[j] = FTPared[j][0][i];
								FTpant1[j] = FTPared[j][1][i];
								FTpant2[j] = FTPared[j][2][i];
							}

							// Tubo exterior
							if(EsPrimeraVez) {
								FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
												   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
							} else {
								if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[0] / FR_ext_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_ext_AxiPos[i]
														+ FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime * (1 / FR_ext_AxiAnt[i] + 1 / FR_ext_AxiPos[i]) + FSUMTime /
														  (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiAnt[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[0] / FR_ext_AxiAnt[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiAnt[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantpos[0] / FR_ext_AxiPos[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiPos[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								}
							}
							FTPared[0][0][i] = (FSUMTime * FTpant1[0] / FR_ext_RadInt[i] + FSUMTPTuboExtPro[0][0][i]) /
											   (FSUMTime / FR_ext_RadInt[i] + FSUMTPTuboExtPro[1][0][i]);
							FTPared[0][2][i] = (FSUMTime * FTpant1[0] / FR_ext_RadExt[i] + FSUMTime * Text / FR_ext[i]) / (FSUMTime *
											   (1 / FR_ext_RadExt[i] + 1 / FR_ext[i]));

							if(ErrorTp < fabs(FTpant1[0] - FTPared[0][1][i])) {
								ErrorTp = fabs(FTpant1[0] - FTPared[0][1][i]);
							}
							for(int k = 0; k < 3; k++) {
								FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
							}
						}
						EsPrimeraVez = false;
					}
				}
				for(int i = 0; i < FTubo[0]->getNin(); i++) {
					for(int j = 0; j < 2; j++) {
						for(int k = 0; k < 2; k++) {
							FSUMTPTuboExtPro[j][k][i] = 0.;
						}
					}
				}
			}
		}

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante && FTubo[0]->getCoefAjustTC() != 0) {
			if(FCicloTubo != Motor[0]->getCiclo()) {
				FSUMTime = 0.;
				FCicloTubo = Motor[0]->getCiclo();
			}
		}

	} catch(exception &N) {
		cout << "ERROR: TConcentricoDPF::CalculaTemperaturaPared en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void TConcentricoDPF::CalculaTemperaturaParedSinMotor(TCondicionContorno **CC) {
	double Tg = 0.;
	double zzz = 0., czz = 0., cz1 = 0., uq1 = 0.;
	double DeltaTTPared = 0.;
	double R_Equiv_1 = 0., R_Equiv = 0., Text = 0., Ri = 0., Re = 0., ErrorTp = 0., Gamma = 0., RMezcla = 0., Asonido = 0.;
	bool EsPrimeraVez;
	int extremo = 0, nodo = 0;

	try {

		DeltaTTPared = FTubo[0]->getTime1() - FTubo[0]->getTime0();

		for(int j = 0; j < FNumeroTubos; j++) {    // Recorre los dos tubos conc�ntricos
			for(int i = 0; i < FTubo[j]->getNin(); i++) {
				FTParedAnt[j][0][i] = FTPared[j][0][i];
				FTParedAnt[j][1][i] = FTPared[j][1][i];
				FTParedAnt[j][2][i] = FTPared[j][2][i];

				zzz = 0.013 / DeltaTTPared;
				czz = 2 / (zzz + 1);
				uq1 = fabs(FTubo[j]->GetVelocidad(i) * __cons::ARef);
				if(FTubo[j]->getTipoTransCal() == nmTuboEscape)
					cz1 = czz;
				else
					cz1 = 1;
				FTubo[j]->PutVelPro(i, cz1 * uq1 + (1 - cz1) * FTubo[j]->GetVelPro(i));
			}
		}

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante
		   && FTubo[0]->getCoefAjustTC() != 0) {   // Tiene que existir motor.
			FCicloActual = FTubo[0]->getTime1() / FDuracionCiclo;
			if(FTubo[0]->getTime1() > FDuracionCiclo) {
				FSUMTime += DeltaTTPared;
			}
			for(int i = 0; i < FTubo[0]->getNin(); i++) {
				//Establece la temperatura del fluido en el exterior del tubo conc�ntrico.
				//Como no hay motor quito los if referentes a tipo de tubo.
				Text = FTubo[0]->getTExt();

				//Establece la temperatura de la pared en el instante de c�lculo anterior...
				//...en los nodos anterior y posterior.
				for(int j = 0; j < FNumeroTubos; j++) {
					//Establece la temperatura del gas en el interior del conducto.
					FTg[j] = pow(FTubo[j]->GetAsonido(i) * __cons::ARef, 2.) / (FTubo[j]->GetGamma(i) * FTubo[j]->GetRMezcla(i));

					if(i == 0) {
						if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
							if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoIzq() - 1])->getConductividad() > 0) {
								if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
									extremo = 1;
									if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								} else {
									extremo = 0;
									if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								}
								FTpantant[j] = __units::degCToK(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
							}
						}
						FTpantpos[j] = FTParedAnt[j][1][i + 1];
					} else if(i == FTubo[j]->getNin() - 1) {
						FTpantant[j] = FTParedAnt[j][1][i - 1];
						if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
							if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoDer() - 1])->getConductividad() > 0) {
								if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
									extremo = 1;
									if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								} else {
									extremo = 0;
									if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
										nodo = 0;
									} else {
										nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
									}
								}
								FTpantpos[j] = __units::degCToK(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
							}
						}
					} else {
						FTpantant[j] = FTParedAnt[j][1][i - 1];
						FTpantpos[j] = FTParedAnt[j][1][i + 1];
					}
					//...en los nodos interior, medio y exterior.
					FTpant0[j] = FTParedAnt[j][0][i];   // Nodo Interior
					FTpant1[j] = FTParedAnt[j][1][i];   // Nodo Medio
					FTpant2[j] = FTParedAnt[j][2][i];   // Nodo Exterior

					//Resistencias t�rmicas de convecci�n.
					if(FTubo[j]->Gethi(i) == 0) {
						FRg_int_ext[i] = 100000000.;
						FRg_ext_int[i] = 100000000.;
						if(FTpant0[0] == FDPF[0]->GetTSuperficie(i, 2)) {
							FR_int_radiacion[i] = 100000000.;
						} else {
							FR_int_radiacion[i] = (1. / FDPF[0]->getEmisividad() + FDPF[0]->getDiametroExt() / (FDPF[0]->getDiametroExt() + 2 *
												   FEspesor_fluido) * (1 / FTubo[j]->getEmisividad() - 1))
												  / (__cons::Sigma * __cons::Pi * FDPF[0]->getDiametroExt() * FTubo[j]->getXRef()) * (FDPF[0]->GetTSuperficie(i,
														  2) - FTpant0[0])
												  / (pow(FDPF[0]->GetTSuperficie(i, 2), 4) - pow(FTpant0[0], 4));
						}
					} else {
						FRg_int_ext[i] = 1 / __cons::Pi / FDPF[0]->getDiametroExt() / FTubo[j]->Gethi(i) / FTubo[j]->getXRef();
						FRg_ext_int[i] = 1 / __cons::Pi / (FDPF[0]->getDiametroExt() + 2 * FEspesor_fluido) / FTubo[j]->Gethi(
											 i) / FTubo[j]->getXRef();
						if(FTpant0[0] == FDPF[0]->GetTSuperficie(i, 2)) {
							FR_int_radiacion[i] = 100000000.;
						} else {
							FR_int_radiacion[i] = (1. / FDPF[0]->getEmisividad() + FDPF[0]->getDiametroExt() / (FDPF[0]->getDiametroExt() + 2 *
												   FEspesor_fluido) * (1 / FTubo[j]->getEmisividad() - 1))
												  / (__cons::Sigma * __cons::Pi * FDPF[0]->getDiametroExt() * FTubo[j]->getXRef()) * (FDPF[0]->GetTSuperficie(i,
														  2) - FTpant0[0])
												  / (pow(FDPF[0]->GetTSuperficie(i, 2), 4) - pow(FTpant0[0], 4));
						}
					}
					if(FTubo[j]->Gethe(i) == 0) {
						FR_ext[i] = 100000000.;
					} else {
						FR_ext[i] = 1 / __cons::Pi / (FDPF[0]->getDiametroExt() + 2 * (FEspesor_fluido + FEspesorTotalExt)) / FTubo[j]->Gethe(
										i) / FTubo[j]->getXRef();
					}
				}

				// C�lculo de las temperaturas de pared.
				// DPF wall temperature calculated at TDPF.cpp
				// Tubo exterior
				FTPared[0][2][i] = DeltaTTPared / FCapExtExt[i] * (1 / FR_ext[i] * (Text - FTpant2[0]) + 1 / FR_ext_RadExt[i] *
								   (FTpant1[0] - FTpant2[0])) + FTpant2[0];
				if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[0] - FTpant1[0])
										  + 1 / FR_ext_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_ext_AxiAnt[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_ext_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_ext_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else {
					FTPared[0][1][i] = DeltaTTPared / FCapExtMed[i] * (1 / FR_ext_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 /
									   FR_ext_RadExt[i] * (FTpant2[0] - FTpant1[0])) + FTpant1[0];
				}
				FTPared[0][0][i] = DeltaTTPared / FCapExtInt[i]
								   * (1 / FRg_ext_int[i] * (FTg[0] - FTpant0[0]) + 1 / FR_int_radiacion[i] * (FDPF[0]->GetTSuperficie(i,
										   2) - FTpant0[1]) + 1 / FR_ext_RadInt[i] * (FTpant1[0] - FTpant0[0]))
								   + FTpant0[0];

				for(int k = 0; k < 3; k++) {
					FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
				}
				// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
				if(FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica || FCicloActual <= FNumCiclosSinInerciaTermica) {
					if(FTubo[0]->getTime1() > FDuracionCiclo) {
						R_Equiv_1 = (FRg_ext_int[i] + FRg_int_ext[i]) * FR_int_radiacion[i] / (FRg_ext_int[i] + FRg_int_ext[i] +
									FR_int_radiacion[i]);
						R_Equiv = FR_ext_RadInt[i] + R_Equiv_1;
						// Tubo exterior
						// Sumatorio de R*Twall_DPF*incrt (para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][0][i] += __units::degCToK(FDPF[0]->GetTSuperficie(i, 2)) * DeltaTTPared / R_Equiv_1;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][0][i] += DeltaTTPared / R_Equiv_1;

						// Sumatorio de R*Twall_DPF*incrt (para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][1][i] += DeltaTTPared * __units::degCToK(FDPF[0]->GetTSuperficie(i, 2)) / R_Equiv;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][1][i] += DeltaTTPared / R_Equiv;
					}
				}

			}

			// Si est� al final del ciclo...
			if(FCicloTubo != FCicloActual && FSUMTime > 0. && FCicloActual > 1) {
				// ...si (el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando) y est� en el segundo ciclo, calcula la temperatura de convergencia
				if((FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica || FCicloActual <= FNumCiclosSinInerciaTermica)
				   && FTubo[0]->getTime1() > FDuracionCiclo) {
					ErrorTp = 1.;
					EsPrimeraVez = true;
					while(ErrorTp >=
						  1) {   //Itera hasta conseguir una diferencia entre las temperaturas de pared menor a 1�C entre pasos.
						ErrorTp = 0.;
						for(int i = 0; i < FTubo[0]->getNin(); i++) {
							for(int j = 0; j < FNumeroTubos; j++) {
								//Establece la temperatura de la pared en el instante de c�lculo anterior...
								//...en los nodos anterior y posterior.
								FTpantant[j] = FTPared[j][1][i];
								FTpantpos[j] = FTPared[j][1][i];
								if(i == 0) {
									if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
										if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoIzq() - 1])->getConductividad() > 0) {
											if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
												extremo = 1;
												if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											} else {
												extremo = 0;
												if(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											}
											FTpantant[j] = __units::degCToK(CC[FTubo[j]->getNodoIzq() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
										}
									}
									FTpantpos[j] = FTPared[j][1][i + 1];
								} else if(i == FTubo[j]->getNin() - 1) {
									FTpantant[j] = FTPared[j][1][i - 1];
									if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
										if(dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[j]->getNodoDer() - 1])->getConductividad() > 0) {
											if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(0).Pipe->getNumeroTubo() == FTubo[j]->getNumeroTubo()) {
												extremo = 1;
												if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											} else {
												extremo = 0;
												if(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).TipoExtremo == nmLeft) {
													nodo = 0;
												} else {
													nodo = CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->getNin() - 1;
												}
											}
											FTpantpos[j] = __units::degCToK(CC[FTubo[j]->getNodoDer() - 1]->GetTuboExtremo(extremo).Pipe->GetTPTuboAnt(1, nodo));
										}
									}
								} else {
									FTpantant[j] = FTPared[j][1][i - 1];
									FTpantpos[j] = FTPared[j][1][i + 1];
								}
								//...en los nodos interior, medio y exterior.
								FTpant0[j] = FTPared[j][0][i];
								FTpant1[j] = FTPared[j][1][i];
								FTpant2[j] = FTPared[j][2][i];
							}

							// Tubo exterior
							if(EsPrimeraVez) {
								FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
												   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
							} else {
								if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[0] / FR_ext_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_ext_AxiPos[i]
														+ FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime * (1 / FR_ext_AxiAnt[i] + 1 / FR_ext_AxiPos[i]) + FSUMTime /
														  (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiAnt[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[0] / FR_ext_AxiAnt[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiAnt[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantpos[0] / FR_ext_AxiPos[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiPos[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else {
									FTPared[0][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								}
							}
							FTPared[0][0][i] = (FSUMTime * FTpant1[0] / FR_ext_RadInt[i] + FSUMTPTuboExtPro[0][0][i]) /
											   (FSUMTime / FR_ext_RadInt[i] + FSUMTPTuboExtPro[1][0][i]);
							FTPared[0][2][i] = (FSUMTime * FTpant1[0] / FR_ext_RadExt[i] + FSUMTime * Text / FR_ext[i]) / (FSUMTime *
											   (1 / FR_ext_RadExt[i] + 1 / FR_ext[i]));

							if(ErrorTp < fabs(FTpant1[0] - FTPared[0][1][i])) {
								ErrorTp = fabs(FTpant1[0] - FTPared[0][1][i]);
							}
							for(int k = 0; k < 3; k++) {
								FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
							}
						}
						EsPrimeraVez = false;
					}
				}
				for(int i = 0; i < FTubo[0]->getNin(); i++) {
					for(int j = 0; j < 2; j++) {
						for(int k = 0; k < 2; k++) {
							FSUMTPTuboExtPro[j][k][i] = 0.;
						}
					}
				}
			}
		}

		if(FCicloTubo != FCicloActual) {
			FCicloTubo = FCicloActual;
			FSUMTime = 0.;
		}

	} catch(exception &N) {
		cout << "ERROR: TConcentricoDPF::CalculaTemperaturaParedSinMotor en el tubo concentrico: " << FNumeroConcentrico <<
			 endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoDPF::CalculaResistenciasdePared(TCondicionContorno **CC) {
	try {
		double Dext = 0., Dint = 0., DIntPrin = 0., Rcond = 0., UnionEspes = 0., UnionConduct = 0., Cap = 0.;
		bool EsInterior;

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante && FTubo[0]->getCoefAjustTC() != 0) {
			FEspesorTotalInt = 0.;
			for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
				FEspesorTotalExt += FTubo[0]->GetCapa(j).Espesor;
			}
			for(int i = 0; i < FTubo[0]->getNin(); i++) {
				// C�lculo en el Tubo exterior
				//C�lculo de las resistencias t�rmicas radiales.
				EsInterior = true;
				FDiametroExtGap = pow(pow(FTubo[0]->GetDiametro(i), 2.) + pow(FDPF[0]->getDiametroExt(), 2), 0.5);
				FEspesor_fluido = (FDiametroExtGap - FDPF[0]->getDiametroExt()) / 2.;
				Dint = FDiametroExtGap;
				FR_ext_RadInt[i] = 0.;
				FR_ext_RadExt[i] = 0.;
				for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[0]->GetCapa(j).Espesor;
					if(FTubo[0]->GetCapa(j).EsPrincipal == false) {
						Rcond = log(Dext / Dint) / __cons::Pi_x_2 / FTubo[0]->GetCapa(j).Conductividad / FTubo[0]->getXRef();
						if(EsInterior) {
							//C�lculo de la resistencia t�rmica radial interior.
							FR_ext_RadInt[i] += Rcond;
						} else {
							//C�lculo de la resistencia t�rmica radial exterior.
							FR_ext_RadExt[i] += Rcond;
						}
					} else {
						//C�lculo de la resistencia t�rmica radial exterior e interior de la capa principal.
						FR_ext_RadInt[i] += log((Dint + FTubo[0]->GetCapa(j).Espesor) / Dint) / __cons::Pi_x_2 / FTubo[0]->GetCapa(
												j).Conductividad / FTubo[0]->getXRef();
						FR_ext_RadExt[i] += log(Dext / (Dint + FTubo[0]->GetCapa(j).Espesor)) / __cons::Pi_x_2 / FTubo[0]->GetCapa(
												j).Conductividad / FTubo[0]->getXRef();
						EsInterior = false;
					}
					Dint = Dext;
				}

				//C�lculo de las resistencias t�rmicas axiales.
				FR_ext_AxiAnt[i] = 0.;
				FR_ext_AxiPos[i] = 0.;
				DIntPrin = FDiametroExtGap + 2 * FTubo[0]->getEspesorIntPrin();
				if(i == 0) {
					if(CC[FTubo[0]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoIzq() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoIzq() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_ext_AxiAnt[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[0]->getEspesorPrin()) *
											   FTubo[0]->getEspesorPrin());
						}
					}
					FR_ext_AxiPos[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
				} else if(i == FTubo[0]->getNin() - 1) {
					FR_ext_AxiAnt[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
					if(CC[FTubo[0]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoDer() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoDer() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_ext_AxiPos[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[0]->getEspesorPrin()) *
											   FTubo[0]->getEspesorPrin());
						}
					}
				} else {
					FR_ext_AxiAnt[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
					FR_ext_AxiPos[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
				}

				//C�lculo de las capacidades t�rmicas.
				EsInterior = true;
				Dint = FDiametroExtGap;
				FCapExtInt[i] = 0.;
				FCapExtMed[i] = 0.;
				FCapExtExt[i] = 0.;
				for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[0]->GetCapa(j).Espesor;
					if(FTubo[0]->GetCapa(j).EsPrincipal == false) {
						if(EsInterior) {
							//C�lculo de la capacidad t�rmica interior.
							Cap = FTubo[0]->GetCapa(j).Density * FTubo[0]->GetCapa(j).CalorEspecifico * __geom::Ring_area(Dint,
									Dext) * FTubo[0]->getXRef();
							FCapExtInt[i] += Cap;
						} else {
							//C�lculo de la capacidad t�rmica exterior.
							Cap = FTubo[0]->GetCapa(j).Density * FTubo[0]->GetCapa(j).CalorEspecifico * __geom::Ring_area(Dint,
									Dext) * FTubo[0]->getXRef();
							FCapExtExt[i] += Cap;
						}
					} else {
						//C�lculo de la capacidad t�rmica exterior, media e interior de la capa principal.
						FCapExtInt[i] += FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin() * __geom::Ring_area(DIntPrin,
										 DIntPrin + 0.5 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						FCapExtMed[i] = FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin()
										* __geom::Ring_area(DIntPrin + 0.5 * FTubo[0]->getEspesorPrin(),
															DIntPrin + 1.5 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						FCapExtExt[i] += FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin()
										 * __geom::Ring_area(DIntPrin + 1.5 * FTubo[0]->getEspesorPrin(),
															 DIntPrin + 2 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						EsInterior = false;
					}
					Dint = Dext;
				}

			}
		}

	} catch(exception &N) {
		cout << "ERROR: TConcentricoTubos::CalculaResistenciasdePared en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#endif
#pragma package(smart_init)
