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

#include "TConcentricoTubos.h"
#include "TTubo.h"
#include "TDPF.h"
#include "TBloqueMotor.h"
#include "TCondicionContorno.h"
#include "TCCUnionEntreTubos.h"
#include "TCCCilindro.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentricoTubos::TConcentricoTubos(int NumeroConcentrico) :
	TConcentrico() {

	FNumeroConcentrico = NumeroConcentrico + 1;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentricoTubos::~TConcentricoTubos() {

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoTubos::LeeDatosTuboConcentrico(const char *FileWAM, fpos_t &filepos, TTubo **Tubo, TDPF **DPF) {
	try {

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d %d ", &FNumTuboInterno, &FNumTuboExterno);

		FNumDPFInterna = 0;

		fgetpos(fich, &filepos);
		fclose(fich);

		FNumeroTubos = 2;

		FTg = new double[FNumeroTubos];
		FTpantpos = new double[FNumeroTubos];
		FTpant0 = new double[FNumeroTubos];
		FTpant1 = new double[FNumeroTubos];
		FTpant2 = new double[FNumeroTubos];
		FTpantant = new double[FNumeroTubos];

		FTubo = new TTubo*[FNumeroTubos];

		FTubo[0] = Tubo[FNumTuboInterno - 1];
		FTubo[1] = Tubo[FNumTuboExterno - 1];

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
		FSUMTPTuboIntPro = new double**[2];
		for(int i = 0; i < 2; i++) {
			FSUMTPTuboExtPro[i] = new double*[2];
			FSUMTPTuboIntPro[i] = new double*[3];
		}
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
				FSUMTPTuboExtPro[i][j] = new double[FTubo[0]->getNin()];
			}
			for(int j = 0; j < 3; j++) {
				FSUMTPTuboIntPro[i][j] = new double[FTubo[0]->getNin()];
			}
		}

		for(int k = 0; k < FTubo[0]->getNin(); k++) {
			FSUMTPTuboIntPro[0][0][k] = 0.;
			FSUMTPTuboIntPro[0][1][k] = 0.;
			FSUMTPTuboIntPro[0][2][k] = 0.;
			FSUMTPTuboIntPro[1][0][k] = 0.;
			FSUMTPTuboIntPro[1][1][k] = 0.;
			FSUMTPTuboIntPro[1][2][k] = 0.;
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
		cout << "ERROR: TConcentricoTubos::LeeDatosTuboConcentrico en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoTubos::CalculaTemperaturaPared(TBloqueMotor **Motor, double theta, TCondicionContorno **CC) {
	double Tg = 0.;
	double zzz = 0., czz = 0., cz1 = 0., uq1 = 0.;
	double DeltaTTPared = 0.;
	double R_Equiv_1, R_Equiv, R_Equiv2, Text, Ri, Re, ErrorTp, Gamma, RMezcla, Asonido;
	bool EsPrimeraVez;
	int extremo = 0, nodo = 0;

	try {

		DeltaTTPared = FTubo[1]->getTime1() - FTubo[1]->getTime0();

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
					Text = FTubo[1]->getTExt();
				} else {
					Text = __units::degCToK(Motor[0]->getTempRefrigerante());
				}

				//Establece la temperatura de la pared en el instante de c�lculo anterior...
				//...en los nodos anterior y posterior.
				//Tpantant=FTParedAnt[1][i]+273.;
				//Tpantpos=FTParedAnt[1][i]+273.;
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
						} else if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
							if(FTubo[j]->getHayDPFNodoIzq()) {
#ifdef ParticulateFilter
								FTpantant[j] = __units::degCToK((FTubo[j]->getDPFEntrada())->GetTSuperficie(FTubo[j]->getNodoDPFEntrada(), 2));
#endif
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
						} else if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
							if(FTubo[j]->getHayDPFNodoDer()) {
#ifdef ParticulateFilter
								FTpantpos[j] = __units::degCToK((FTubo[j]->getDPFSalida())->GetTSuperficie(FTubo[j]->getNodoDPFSalida(), 2));
#endif
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
						if(j == 0) {
							FRg_int[i] = 100000000.;
						} else {
							FRg_int_ext[i] = 100000000.;
							FRg_ext_int[i] = 100000000.;
							if(FTpant0[1] == FTpant2[0]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FTubo[0]->getEmisividad()
													   + (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / (FTubo[0]->GetDiametro(i) + 2 *
															   (FEspesorTotalInt + FEspesor_fluido)) * (1 / FTubo[1]->getEmisividad() - 1))
													  / (__cons::Sigma * __cons::Pi * (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) * FTubo[0]->getXRef()) *
													  (FTpant2[0] - FTpant0[1])
													  / (pow(FTpant2[0], 4) - pow(FTpant0[1], 4));
							}
						}
					} else {
						if(j == 0) {
							FRg_int[i] = 1 / __cons::Pi / FTubo[j]->GetDiametro(i) / FTubo[j]->Gethi(i) / FTubo[j]->getXRef();
						} else {
							FRg_int_ext[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / FTubo[j]->Gethi(
												 i) / FTubo[j]->getXRef();
							FRg_ext_int[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * (FEspesorTotalInt + FEspesor_fluido)) /
											 FTubo[j]->Gethi(i) / FTubo[j]->getXRef();
							if(FTpant0[1] == FTpant2[0]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FTubo[0]->getEmisividad()
													   + (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / (FTubo[0]->GetDiametro(i) + 2 *
															   (FEspesorTotalInt + FEspesor_fluido)) * (1 / FTubo[1]->getEmisividad() - 1))
													  / (__cons::Sigma * __cons::Pi * (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) * FTubo[0]->getXRef()) *
													  (FTpant2[0] - FTpant0[1])
													  / (pow(FTpant2[0], 4) - pow(FTpant0[1], 4));
							}
						}
					}
					if(j == 1) {
						if(FTubo[j]->Gethe(i) == 0) {
							FR_ext[i] = 100000000.;
						} else {
							FR_ext[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * (FEspesorTotalInt + FEspesor_fluido + FEspesorTotalExt)) /
										FTubo[j]->Gethe(i) / FTubo[j]->getXRef();
						}
					}
				}
				//C�lculo de las temperaturas de pared.
				// Tubo interior
				FTPared[0][2][i] = DeltaTTPared / FCapIntExt[i]
								   * (1 / FR_int_RadExt[i] * (FTpant1[0] - FTpant2[0]) + 1 / FRg_int_ext[i] * (FTg[1] - FTpant2[0]) + 1 /
									  FR_int_radiacion[i] * (FTpant0[1] - FTpant2[0])) + FTpant2[0];
				if(FR_int_AxiAnt[i] > 0. && FR_int_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiAnt[i] * (FTpantant[0] - FTpant1[0])
										  + 1 / FR_int_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_int_AxiAnt[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiAnt[i] * (FTpantant[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_int_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i] * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 /
									   FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0])) + FTpant1[0];
				}
				FTPared[0][0][i] = DeltaTTPared / FCapIntInt[i] * (1 / FRg_int[i] * (FTg[0] - FTpant0[0]) + 1 / FR_int_RadInt[i] *
								   (FTpant1[0] - FTpant0[0])) + FTpant0[0];

				for(int k = 0; k < 3; k++) {
					FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
				}

				// Tubo exterior
				FTPared[1][2][i] = DeltaTTPared / FCapExtExt[i] * (1 / FR_ext[i] * (Text - FTpant2[1]) + 1 / FR_ext_RadExt[i] *
								   (FTpant1[1] - FTpant2[1])) + FTpant2[1];
				if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[1] - FTpant1[1])
										  + 1 / FR_ext_AxiPos[i] * (FTpantpos[1] - FTpant1[1])) + FTpant1[1];
				} else if(FR_ext_AxiAnt[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[1] - FTpant1[1])) + FTpant1[1];
				} else if(FR_ext_AxiPos[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiPos[i] * (FTpantpos[1] - FTpant1[1])) + FTpant1[1];
				} else {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i] * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 /
									   FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1])) + FTpant1[1];
				}
				FTPared[1][0][i] = DeltaTTPared / FCapExtInt[i]
								   * (1 / FRg_ext_int[i] * (FTg[1] - FTpant0[1]) + 1 / FR_int_radiacion[i] * (FTpant2[0] - FTpant0[1]) + 1 /
									  FR_ext_RadInt[i] * (FTpant1[1] - FTpant0[1])) + FTpant0[1];

				for(int k = 0; k < 3; k++) {
					FTubo[1]->PutTPTubo(k, i, __units::KTodegC(FTPared[1][k][i]));
				}
				// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
				if(FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica
				   || theta / FTubo[0]->getAnguloTotalCiclo() <= Motor[0]->getNumCiclosSinInerciaTermica()) {
					if(theta > FTubo[0]->getAnguloTotalCiclo()) {
						// Tubo interior
						R_Equiv_1 = (FRg_ext_int[i] + FRg_int_ext[i]) * FR_int_radiacion[i] / (FRg_ext_int[i] + FRg_int_ext[i] +
									FR_int_radiacion[i]);
						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboIntPro[0][0][i] += FTg[0] * DeltaTTPared / FRg_int[i];
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][0][i] += DeltaTTPared / FRg_int[i];

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						R_Equiv = FR_int_RadExt[i] + R_Equiv_1;
						FSUMTPTuboIntPro[0][1][i] += (FTg[0] / (FR_int_RadInt[i] + FRg_int[i]) + FTPared[1][0][i] / R_Equiv) * DeltaTTPared;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][1][i] += (1 / (FR_int_RadInt[i] + FRg_int[i]) + 1 / R_Equiv) * DeltaTTPared;

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboIntPro[0][2][i] += DeltaTTPared * FTPared[1][0][i] / R_Equiv_1;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][2][i] += DeltaTTPared / R_Equiv_1;

						// Tubo exterior
						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][0][i] += FTPared[0][2][i] * DeltaTTPared / R_Equiv_1;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][0][i] += DeltaTTPared / R_Equiv_1;

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						R_Equiv2 = FR_ext_RadInt[i] + R_Equiv_1;
						FSUMTPTuboExtPro[0][1][i] += DeltaTTPared * FTPared[0][2][i] / R_Equiv2;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][1][i] += DeltaTTPared / R_Equiv2;

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
									} else if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
										if(FTubo[j]->getHayDPFNodoIzq()) {
#ifdef ParticulateFilter
											FTpantant[j] = __units::degCToK((FTubo[j]->getDPFEntrada())->GetTSuperficie(FTubo[j]->getNodoDPFEntrada(), 2));
#endif
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
									} else if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
										if(FTubo[j]->getHayDPFNodoDer()) {
#ifdef ParticulateFilter
											FTpantpos[j] = __units::degCToK((FTubo[j]->getDPFSalida())->GetTSuperficie(FTubo[j]->getNodoDPFSalida(), 2));
#endif
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
							// Tubo interior
							if(EsPrimeraVez) {
								FTPared[0][1][i] = FSUMTPTuboIntPro[0][1][i] / FSUMTPTuboIntPro[1][1][i];
							} else {
								if(FR_int_AxiAnt[i] > 0. && FR_int_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantant[0] / FR_int_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_int_AxiPos[i])
													   / (FSUMTPTuboIntPro[1][1][i] + FSUMTime * (1 / FR_int_AxiAnt[i] + 1 / FR_int_AxiPos[i]));
								} else if(FR_int_AxiAnt[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantant[0] / FR_int_AxiAnt[i]) /
													   (FSUMTPTuboIntPro[1][1][i] + FSUMTime / FR_int_AxiAnt[i]);
								} else if(FR_int_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantpos[0] / FR_int_AxiPos[i]) /
													   (FSUMTPTuboIntPro[1][1][i] + FSUMTime / FR_int_AxiPos[i]);
								} else {
									FTPared[0][1][i] = FSUMTPTuboIntPro[0][1][i] / FSUMTPTuboIntPro[1][1][i];
								}
							}
							FTPared[0][0][i] = (FSUMTime * FTpant1[0] / FR_int_RadInt[i] + FSUMTPTuboIntPro[0][0][i]) /
											   (FSUMTime / FR_int_RadInt[i] + FSUMTPTuboIntPro[1][0][i]);
							FTPared[0][2][i] = (FSUMTime * FTpant1[0] / FR_int_RadExt[i] + FSUMTPTuboIntPro[0][2][i]) /
											   (FSUMTime / FR_int_RadExt[i] + FSUMTPTuboIntPro[1][2][i]);

							if(ErrorTp < fabs(FTpant1[0] - FTPared[0][1][i])) {
								ErrorTp = fabs(FTpant1[0] - FTPared[0][1][i]);
							}

							// Tubo exterior
							if(EsPrimeraVez) {
								FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
												   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
							} else {
								if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[1] / FR_ext_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_ext_AxiPos[i]
														+ FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime * (1 / FR_ext_AxiAnt[i] + 1 / FR_ext_AxiPos[i]) + FSUMTime /
														  (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiAnt[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[1] / FR_ext_AxiAnt[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiAnt[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiPos[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantpos[1] / FR_ext_AxiPos[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiPos[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								}
							}
							FTPared[1][0][i] = (FSUMTime * FTpant1[1] / FR_ext_RadInt[i] + FSUMTPTuboExtPro[0][0][i]) /
											   (FSUMTime / FR_ext_RadInt[i] + FSUMTPTuboExtPro[1][0][i]);
							FTPared[1][2][i] = (FSUMTime * FTpant1[1] / FR_ext_RadExt[i] + FSUMTime * Text / FR_ext[i]) / (FSUMTime *
											   (1 / FR_ext_RadExt[i] + 1 / FR_ext[i]));

							if(ErrorTp < fabs(FTpant1[1] - FTPared[1][1][i])) {
								ErrorTp = fabs(FTpant1[1] - FTPared[1][1][i]);
							}
							for(int k = 0; k < 3; k++) {
								FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
								FTubo[1]->PutTPTubo(k, i, __units::KTodegC(FTPared[1][k][i]));
							}
						}
						EsPrimeraVez = false;
					}
				}
				for(int i = 0; i < FTubo[0]->getNin(); i++) {
					for(int j = 0; j < 2; j++) {
						for(int k = 0; k < 3; k++) {
							FSUMTPTuboIntPro[j][k][i] = 0.;
						}
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
		cout << "ERROR: TConcentricoTubos::CalculaTemperaturaPared en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoTubos::CalculaTemperaturaParedSinMotor(TCondicionContorno **CC) {
	double Tg = 0.;
	double zzz = 0., czz = 0., cz1 = 0., uq1 = 0.;
	double DeltaTTPared = 0.;
	double R_Equiv_1, R_Equiv, R_Equiv2, Text, Ri, Re, ErrorTp, Gamma, RMezcla, Asonido;
	bool EsPrimeraVez;
	int extremo = 0, nodo = 0;

	try {

		DeltaTTPared = FTubo[1]->getTime1() - FTubo[1]->getTime0();

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

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante && FTubo[0]->getCoefAjustTC() != 0) {
			FCicloActual = FTubo[0]->getTime1() / FDuracionCiclo;
			if(FTubo[0]->getTime1() > FDuracionCiclo) {
				FSUMTime += DeltaTTPared;
			}
			for(int i = 0; i < FTubo[0]->getNin(); i++) {
				//Establece la temperatura del fluido en el exterior del tubo conc�ntrico.
				//Como no hay motor quito los if referentes a tipo de tubo.
				Text = FTubo[1]->getTExt();

				//Establece la temperatura de la pared en el instante de c�lculo anterior...
				//...en los nodos anterior y posterior.
				//Tpantant=FTParedAnt[1][i]+273.;
				//Tpantpos=FTParedAnt[1][i]+273.;
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
						} else if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
							if(FTubo[j]->getHayDPFNodoIzq()) {
#ifdef ParticulateFilter
								FTpantant[j] = __units::degCToK((FTubo[j]->getDPFEntrada())->GetTSuperficie(FTubo[j]->getNodoDPFEntrada(), 2));
#endif
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
							} else if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
								if(FTubo[j]->getHayDPFNodoDer()) {
#ifdef ParticulateFilter
									FTpantpos[j] = __units::degCToK((FTubo[j]->getDPFSalida())->GetTSuperficie(FTubo[j]->getNodoDPFSalida(), 2));
#endif
								}
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
						if(j == 0) {
							FRg_int[i] = 100000000.;
						} else {
							FRg_int_ext[i] = 100000000.;
							FRg_ext_int[i] = 100000000.;
							if(FTpant0[1] == FTpant2[0]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FTubo[0]->getEmisividad()
													   + (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / (FTubo[0]->GetDiametro(i) + 2 *
															   (FEspesorTotalInt + FEspesor_fluido)) * (1 / FTubo[1]->getEmisividad() - 1))
													  / (__cons::Sigma * __cons::Pi * (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) * FTubo[0]->getXRef()) *
													  (FTpant2[0] - FTpant0[1])
													  / (pow(FTpant2[0], 4) - pow(FTpant0[1], 4));
							}
						}
					} else {
						if(j == 0) {
							FRg_int[i] = 1 / __cons::Pi / FTubo[j]->GetDiametro(i) / FTubo[j]->Gethi(i) / FTubo[j]->getXRef();
						} else {
							FRg_int_ext[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / FTubo[j]->Gethi(
												 i) / FTubo[j]->getXRef();
							FRg_ext_int[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * (FEspesorTotalInt + FEspesor_fluido)) /
											 FTubo[j]->GetDiametro(i) / FTubo[j]->getXRef();
							if(FTpant0[1] == FTpant2[0]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FTubo[0]->getEmisividad()
													   + (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) / (FTubo[0]->GetDiametro(i) + 2 *
															   (FEspesorTotalInt + FEspesor_fluido)) * (1 / FTubo[1]->getEmisividad() - 1))
													  / (__cons::Sigma * __cons::Pi * (FTubo[0]->GetDiametro(i) + 2 * FEspesorTotalInt) * FTubo[0]->getXRef()) *
													  (FTpant2[0] - FTpant0[1])
													  / (pow(FTpant2[0], 4) - pow(FTpant0[1], 4));
							}
						}
					}
					if(j == 1) {
						if(FTubo[j]->Gethe(i) == 0) {
							FR_ext[i] = 100000000.;
						} else {
							FR_ext[i] = 1 / __cons::Pi / (FTubo[0]->GetDiametro(i) + 2 * (FEspesorTotalInt + FEspesor_fluido + FEspesorTotalExt)) /
										FTubo[j]->Gethe(i) / FTubo[j]->getXRef();
						}
					}
				}
				//C�lculo de las temperaturas de pared.
				// Tubo interior
				FTPared[0][2][i] = DeltaTTPared / FCapIntExt[i]
								   * (1 / FR_int_RadExt[i] * (FTpant1[0] - FTpant2[0]) + 1 / FRg_int_ext[i] * (FTg[1] - FTpant2[0]) + 1 /
									  FR_int_radiacion[i] * (FTpant0[1] - FTpant2[0])) + FTpant2[0];
				if(FR_int_AxiAnt[i] > 0. && FR_int_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiAnt[i] * (FTpantant[0] - FTpant1[0])
										  + 1 / FR_int_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_int_AxiAnt[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiAnt[i] * (FTpantant[0] - FTpant1[0])) + FTpant1[0];
				} else if(FR_int_AxiPos[i] > 0.) {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i]
									   * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 / FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0]) + 1 /
										  FR_int_AxiPos[i] * (FTpantpos[0] - FTpant1[0])) + FTpant1[0];
				} else {
					FTPared[0][1][i] = DeltaTTPared / FCapIntMed[i] * (1 / FR_int_RadInt[i] * (FTpant0[0] - FTpant1[0]) + 1 /
									   FR_int_RadExt[i] * (FTpant2[0] - FTpant1[0])) + FTpant1[0];
				}
				FTPared[0][0][i] = DeltaTTPared / FCapIntInt[i] * (1 / FRg_int[i] * (FTg[0] - FTpant0[0]) + 1 / FR_int_RadInt[i] *
								   (FTpant1[0] - FTpant0[0])) + FTpant0[0];

				for(int k = 0; k < 3; k++) {
					FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
				}

				// Tubo exterior
				FTPared[1][2][i] = DeltaTTPared / FCapExtExt[i] * (1 / FR_ext[i] * (Text - FTpant2[1]) + 1 / FR_ext_RadExt[i] *
								   (FTpant1[1] - FTpant2[1])) + FTpant2[1];
				if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[1] - FTpant1[1])
										  + 1 / FR_ext_AxiPos[i] * (FTpantpos[1] - FTpant1[1])) + FTpant1[1];
				} else if(FR_ext_AxiAnt[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiAnt[i] * (FTpantant[1] - FTpant1[1])) + FTpant1[1];
				} else if(FR_ext_AxiPos[i] > 0.) {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i]
									   * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 / FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1]) + 1 /
										  FR_ext_AxiPos[i] * (FTpantpos[1] - FTpant1[1])) + FTpant1[1];
				} else {
					FTPared[1][1][i] = DeltaTTPared / FCapExtMed[i] * (1 / FR_ext_RadInt[i] * (FTpant0[1] - FTpant1[1]) + 1 /
									   FR_ext_RadExt[i] * (FTpant2[1] - FTpant1[1])) + FTpant1[1];
				}
				FTPared[1][0][i] = DeltaTTPared / FCapExtInt[i]
								   * (1 / FRg_ext_int[i] * (FTg[1] - FTpant0[1]) + 1 / FR_int_radiacion[i] * (FTpant2[0] - FTpant0[1]) + 1 /
									  FR_ext_RadInt[i] * (FTpant1[1] - FTpant0[1])) + FTpant0[1];

				for(int k = 0; k < 3; k++) {
					FTubo[1]->PutTPTubo(k, i, __units::KTodegC(FTPared[1][k][i]));
				}
				// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
				if(FTubo[0]->getTipoCalcTempPared() == nmVariableSinInerciaTermica || FCicloActual <= FNumCiclosSinInerciaTermica) {
					if(FTubo[0]->getTime1() > FDuracionCiclo) {
						// Tubo interior
						R_Equiv_1 = (FRg_ext_int[i] + FRg_int_ext[i]) * FR_int_radiacion[i] / (FRg_ext_int[i] + FRg_int_ext[i] +
									FR_int_radiacion[i]);
						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboIntPro[0][0][i] += FTg[0] * DeltaTTPared / FRg_int[i];
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][0][i] += DeltaTTPared / FRg_int[i];

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						R_Equiv = FR_int_RadExt[i] + R_Equiv_1;
						FSUMTPTuboIntPro[0][1][i] += (FTg[0] / (FR_int_RadInt[i] + FRg_int[i]) + FTPared[1][0][i] / R_Equiv) * DeltaTTPared;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][1][i] += (1 / (FR_int_RadInt[i] + FRg_int[i]) + 1 / R_Equiv) * DeltaTTPared;

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboIntPro[0][2][i] += DeltaTTPared * FTPared[1][0][i] / R_Equiv_1;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboIntPro[1][2][i] += DeltaTTPared / R_Equiv_1;

						// Tubo exterior
						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						FSUMTPTuboExtPro[0][0][i] += FTPared[0][2][i] * DeltaTTPared / R_Equiv_1;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][0][i] += DeltaTTPared / R_Equiv_1;

						// Sumatorio de h*Tg*incrt interior(para el c�lculo del numerador de la integral).
						R_Equiv2 = FR_ext_RadInt[i] + R_Equiv_1;
						FSUMTPTuboExtPro[0][1][i] += DeltaTTPared * FTPared[0][2][i] / R_Equiv2;
						// Sumatorio de h*incrt interior(para el c�lculo del denominador de la integral).
						FSUMTPTuboExtPro[1][1][i] += DeltaTTPared / R_Equiv2;

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
									} else if(CC[FTubo[j]->getNodoIzq() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
										if(FTubo[j]->getHayDPFNodoIzq()) {
#ifdef ParticulateFilter
											FTpantant[j] = __units::degCToK((FTubo[j]->getDPFEntrada())->GetTSuperficie(FTubo[j]->getNodoDPFEntrada(), 2));
#endif
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
									} else if(CC[FTubo[j]->getNodoDer() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
										if(FTubo[j]->getHayDPFNodoDer()) {
#ifdef ParticulateFilter
											FTpantpos[j] = __units::degCToK((FTubo[j]->getDPFSalida())->GetTSuperficie(FTubo[j]->getNodoDPFSalida(), 2));
#endif
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
							// Tubo interior
							if(EsPrimeraVez) {
								FTPared[0][1][i] = FSUMTPTuboIntPro[0][1][i] / FSUMTPTuboIntPro[1][1][i];
							} else {
								if(FR_int_AxiAnt[i] > 0. && FR_int_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantant[0] / FR_int_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_int_AxiPos[i])
													   / (FSUMTPTuboIntPro[1][1][i] + FSUMTime * (1 / FR_int_AxiAnt[i] + 1 / FR_int_AxiPos[i]));
								} else if(FR_int_AxiAnt[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantant[0] / FR_int_AxiAnt[i]) /
													   (FSUMTPTuboIntPro[1][1][i] + FSUMTime / FR_int_AxiAnt[i]);
								} else if(FR_int_AxiPos[i] > 0.) {
									FTPared[0][1][i] = (FSUMTPTuboIntPro[0][1][i] + FSUMTime * FTpantpos[0] / FR_int_AxiPos[i]) /
													   (FSUMTPTuboIntPro[1][1][i] + FSUMTime / FR_int_AxiPos[i]);
								} else {
									FTPared[0][1][i] = FSUMTPTuboIntPro[0][1][i] / FSUMTPTuboIntPro[1][1][i];
								}
							}
							FTPared[0][0][i] = (FSUMTime * FTpant1[0] / FR_int_RadInt[i] + FSUMTPTuboIntPro[0][0][i]) /
											   (FSUMTime / FR_int_RadInt[i] + FSUMTPTuboIntPro[1][0][i]);
							FTPared[0][2][i] = (FSUMTime * FTpant1[0] / FR_int_RadExt[i] + FSUMTPTuboIntPro[0][2][i]) /
											   (FSUMTime / FR_int_RadExt[i] + FSUMTPTuboIntPro[1][2][i]);

							if(ErrorTp < fabs(FTpant1[0] - FTPared[0][1][i])) {
								ErrorTp = fabs(FTpant1[0] - FTPared[0][1][i]);
							}

							// Tubo exterior
							if(EsPrimeraVez) {
								FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
												   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
							} else {
								if(FR_ext_AxiAnt[i] > 0. && FR_ext_AxiPos[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[1] / FR_ext_AxiAnt[i] + FSUMTime * FTpantpos[0] /
														FR_ext_AxiPos[i]
														+ FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime * (1 / FR_ext_AxiAnt[i] + 1 / FR_ext_AxiPos[i]) + FSUMTime /
														  (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiAnt[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantant[1] / FR_ext_AxiAnt[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiAnt[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else if(FR_ext_AxiPos[i] > 0.) {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * FTpantpos[1] / FR_ext_AxiPos[i] + FSUMTime * Text /
														(FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / FR_ext_AxiPos[i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								} else {
									FTPared[1][1][i] = (FSUMTPTuboExtPro[0][1][i] + FSUMTime * Text / (FR_ext_RadExt[i] + FR_ext[i]))
													   / (FSUMTPTuboExtPro[1][1][i] + FSUMTime / (FR_ext_RadExt[i] + FR_ext[i]));
								}
							}
							FTPared[1][0][i] = (FSUMTime * FTpant1[1] / FR_ext_RadInt[i] + FSUMTPTuboExtPro[0][0][i]) /
											   (FSUMTime / FR_ext_RadInt[i] + FSUMTPTuboExtPro[1][0][i]);
							FTPared[1][2][i] = (FSUMTime * FTpant1[1] / FR_ext_RadExt[i] + FSUMTime * Text / FR_ext[i]) / (FSUMTime *
											   (1 / FR_ext_RadExt[i] + 1 / FR_ext[i]));

							if(ErrorTp < fabs(FTpant1[1] - FTPared[1][1][i])) {
								ErrorTp = fabs(FTpant1[1] - FTPared[1][1][i]);
							}
							for(int k = 0; k < 3; k++) {
								FTubo[0]->PutTPTubo(k, i, __units::KTodegC(FTPared[0][k][i]));
								FTubo[1]->PutTPTubo(k, i, __units::KTodegC(FTPared[1][k][i]));
							}
						}
						EsPrimeraVez = false;
					}
				}
				for(int i = 0; i < FTubo[0]->getNin(); i++) {
					for(int j = 0; j < 2; j++) {
						for(int k = 0; k < 3; k++) {
							FSUMTPTuboIntPro[j][k][i] = 0.;
						}
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
		cout << "ERROR: TConcentricoTubos::CalculaTemperaturaParedSinMotor en el tubo concentrico: " << FNumeroConcentrico <<
			 endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentricoTubos::CalculaResistenciasdePared(TCondicionContorno **CC) {
	try {
		double Dext = 0., Dint = 0., DIntPrin = 0., Rcond = 0., Rrad = 0., UnionEspes = 0., UnionConduct = 0., Cap = 0.;
		bool EsInterior;

		if(FTubo[0]->getTipoCalcTempPared() != nmTempConstante && FTubo[0]->getCoefAjustTC() != 0) {
			FEspesorTotalInt = 0.;
			for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
				FEspesorTotalInt += FTubo[0]->GetCapa(j).Espesor;
			}
			for(int j = 0; j < FTubo[1]->getNumCapas(); j++) {
				FEspesorTotalExt += FTubo[1]->GetCapa(j).Espesor;
			}
			for(int i = 0; i < FTubo[0]->getNin(); i++) {
				// C�lculo en el Tubo interior
				//C�lculo de las resistencias t�rmicas radiales.
				EsInterior = true;
				Dint = FTubo[0]->GetDiametro(i);
				FR_int_RadInt[i] = 0.;

				for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[0]->GetCapa(j).Espesor;
					if(FTubo[0]->GetCapa(j).EsPrincipal == false) {
						Rcond = log(Dext / Dint) / __cons::Pi_x_2 / FTubo[0]->GetCapa(j).Conductividad / FTubo[0]->getXRef();
						if(EsInterior) {
							//C�lculo de la resistencia t�rmica radial interior.
							FR_int_RadInt[i] += Rcond;
						} else {
							//C�lculo de la resistencia t�rmica radial exterior.
							FR_int_RadExt[i] += Rcond;
						}
					} else {
						//C�lculo de la resistencia t�rmica radial exterior e interior de la capa principal.
						FR_int_RadInt[i] += log((Dint + FTubo[0]->GetCapa(j).Espesor) / Dint) / __cons::Pi_x_2 / FTubo[0]->GetCapa(
												j).Conductividad / FTubo[0]->getXRef();
						FR_int_RadExt[i] += log(Dext / (Dint + FTubo[0]->GetCapa(j).Espesor)) / __cons::Pi_x_2 / FTubo[0]->GetCapa(
												j).Conductividad / FTubo[0]->getXRef();
						EsInterior = false;
					}
					Dint = Dext;
				}

				//C�lculo de las resistencias t�rmicas axiales.
				FR_int_AxiAnt[i] = 0.;
				FR_int_AxiPos[i] = 0.;
				DIntPrin = FTubo[0]->GetDiametro(i) + 2 * FTubo[0]->getEspesorIntPrin();
				if(i == 0) {
					if(CC[FTubo[0]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoIzq() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoIzq() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_int_AxiAnt[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[0]->getEspesorPrin()) *
											   FTubo[0]->getEspesorPrin());
						}
					} else if(CC[FTubo[0]->getNodoIzq() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
						if(FTubo[0]->getHayDPFNodoIzq()) {
#ifdef ParticulateFilter
							if(FTubo[0]->GetTipoCanal(0) == 0) {    // Junction to a DPF inlet channel
								FR_int_AxiAnt[i] = (FTubo[0]->getDPFEntrada())->getAjustRAxAnt()
												   / (__cons::Pi * (FTubo[0]->getDPFEntrada())->getConductividadMetal() * (FTubo[0]->getDPFEntrada())->getEspesorMetal()
													  * ((FTubo[0]->getDPFEntrada())->getDiametroEfect() + 2 * (FTubo[0]->getDPFEntrada())->getEspesorAislante() + 2 *
														 (FTubo[0]->getDPFEntrada())->getEspesorAire()
														 + (FTubo[0]->getDPFEntrada())->getEspesorMetal()));
							} else {   // Junction to a DPF outlet channel
								FR_int_AxiAnt[i] = (FTubo[0]->getDPFEntrada())->getAjustRAxPos()
												   / (__cons::Pi * (FTubo[0]->getDPFEntrada())->getConductividadMetal() * (FTubo[0]->getDPFEntrada())->getEspesorMetal()
													  * ((FTubo[0]->getDPFEntrada())->getDiametroEfect() + 2 * (FTubo[0]->getDPFEntrada())->getEspesorAislante() + 2 *
														 (FTubo[0]->getDPFEntrada())->getEspesorAire()
														 + (FTubo[0]->getDPFEntrada())->getEspesorMetal()));
							}
#endif
						}
					}
					FR_int_AxiPos[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
				} else if(i == FTubo[0]->getNin() - 1) {
					FR_int_AxiAnt[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
					if(CC[FTubo[0]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoDer() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[0]->getNodoDer() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_int_AxiPos[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[0]->getEspesorPrin()) *
											   FTubo[0]->getEspesorPrin());
						}
					} else if(CC[FTubo[0]->getNodoDer() - 1]->getTipoCC() == nmPipeToPlenumConnection) {
						if(FTubo[0]->getHayDPFNodoDer()) {
#ifdef ParticulateFilter
							if(FTubo[0]->GetTipoCanal(1) == 0) {    // Junction to a DPF inlet channel
								FR_int_AxiPos[i] = (FTubo[0]->getDPFSalida())->getAjustRAxAnt()
												   / (__cons::Pi * (FTubo[0]->getDPFSalida())->getConductividadMetal() * (FTubo[0]->getDPFSalida())->getEspesorMetal()
													  * ((FTubo[0]->getDPFSalida())->getDiametroEfect() + 2 * (FTubo[0]->getDPFSalida())->getEspesorAislante() + 2 *
														 (FTubo[0]->getDPFSalida())->getEspesorAire()
														 + (FTubo[0]->getDPFSalida())->getEspesorMetal()));
							} else {   // Junction to a DPF outlet channel
								FR_int_AxiPos[i] = (FTubo[0]->getDPFSalida())->getAjustRAxPos()
												   / (__cons::Pi * (FTubo[0]->getDPFSalida())->getConductividadMetal() * (FTubo[0]->getDPFSalida())->getEspesorMetal()
													  * ((FTubo[0]->getDPFSalida())->getDiametroEfect() + 2 * (FTubo[0]->getDPFSalida())->getEspesorAislante() + 2 *
														 (FTubo[0]->getDPFSalida())->getEspesorAire()
														 + (FTubo[0]->getDPFSalida())->getEspesorMetal()));
							}
#endif
						}
					}
				} else {
					FR_int_AxiAnt[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
					FR_int_AxiPos[i] = FTubo[0]->getXRef() / FTubo[0]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[0]->getEspesorPrin()) * FTubo[0]->getEspesorPrin());
				}

				//C�lculo de las capacidades t�rmicas.
				EsInterior = true;
				Dint = FTubo[0]->GetDiametro(i);
				FCapIntInt[i] = 0.;
				FCapIntMed[i] = 0.;
				FCapIntExt[i] = 0.;
				for(int j = 0; j < FTubo[0]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[0]->GetCapa(j).Espesor;
					if(FTubo[0]->GetCapa(j).EsPrincipal == false) {
						if(EsInterior) {
							//C�lculo de la capacidad t�rmica interior.
							Cap = FTubo[0]->GetCapa(j).Density * FTubo[0]->GetCapa(j).CalorEspecifico * __geom::Ring_area(Dint,
									Dext) * FTubo[0]->getXRef();
							FCapIntInt[i] += Cap;
						} else {
							//C�lculo de la capacidad t�rmica exterior.
							Cap = FTubo[0]->GetCapa(j).Density * FTubo[0]->GetCapa(j).CalorEspecifico * __geom::Ring_area(Dint,
									Dext) * FTubo[0]->getXRef();
							FCapIntExt[i] += Cap;
						}
					} else {
						//C�lculo de la capacidad t�rmica exterior, media e interior de la capa principal.
						FCapIntInt[i] += FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin() * __geom::Ring_area(DIntPrin,
										 DIntPrin + 0.5 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						FCapIntMed[i] = FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin()
										* __geom::Ring_area(DIntPrin + 0.5 * FTubo[0]->getEspesorPrin(),
															DIntPrin + 1.5 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						FCapIntExt[i] += FTubo[0]->getDensidadPrin() * FTubo[0]->getCalEspPrin()
										 * __geom::Ring_area(DIntPrin + 1.5 * FTubo[0]->getEspesorPrin(),
															 DIntPrin + 2.0 * FTubo[0]->getEspesorPrin()) * FTubo[0]->getXRef();
						EsInterior = false;
					}
					Dint = Dext;
				}

				// C�lculo en el Tubo exterior
				//C�lculo de las resistencias t�rmicas radiales.
				EsInterior = true;
				FDiametroExtGap = sqrt(pow2(FTubo[1]->GetDiametro(i)) + pow2(FTubo[0]->GetDiametro(i) + 2. * FEspesorTotalInt));
				FEspesor_fluido = (FDiametroExtGap - (FTubo[0]->GetDiametro(i) + 2. * FEspesorTotalInt)) / 2.;
				Dint = FDiametroExtGap;
				FR_ext_RadInt[i] = 0.;
				FR_ext_RadExt[i] = 0.;
				for(int j = 0; j < FTubo[1]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[1]->GetCapa(j).Espesor;
					if(FTubo[1]->GetCapa(j).EsPrincipal == false) {
						Rcond = log(Dext / Dint) / __cons::Pi_x_2 / FTubo[1]->GetCapa(j).Conductividad / FTubo[1]->getXRef();
						if(EsInterior) {
							//C�lculo de la resistencia t�rmica radial interior.
							FR_ext_RadInt[i] += Rcond;
						} else {
							//C�lculo de la resistencia t�rmica radial exterior.
							FR_ext_RadExt[i] += Rcond;
						}
					} else {
						//C�lculo de la resistencia t�rmica radial exterior e interior de la capa principal.
						FR_ext_RadInt[i] += log((Dint + FTubo[1]->GetCapa(j).Espesor) / Dint) / __cons::Pi_x_2 / FTubo[1]->GetCapa(
												j).Conductividad / FTubo[1]->getXRef();
						FR_ext_RadExt[i] += log(Dext / (Dint + FTubo[1]->GetCapa(j).Espesor)) / __cons::Pi_x_2 / FTubo[1]->GetCapa(
												j).Conductividad / FTubo[1]->getXRef();
						EsInterior = false;
					}
					Dint = Dext;
				}

				//C�lculo de las resistencias t�rmicas axiales.
				FR_ext_AxiAnt[i] = 0.;
				FR_ext_AxiPos[i] = 0.;
				DIntPrin = FDiametroExtGap + 2 * FTubo[1]->getEspesorIntPrin();
				if(i == 0) {
					if(CC[FTubo[1]->getNodoIzq() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[1]->getNodoIzq() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[1]->getNodoIzq() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_ext_AxiAnt[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[1]->getEspesorPrin()) *
											   FTubo[1]->getEspesorPrin());
						}
					}
					FR_ext_AxiPos[i] = FTubo[1]->getXRef() / FTubo[1]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[1]->getEspesorPrin()) * FTubo[1]->getEspesorPrin());
				} else if(i == FTubo[1]->getNin() - 1) {
					FR_ext_AxiAnt[i] = FTubo[1]->getXRef() / FTubo[1]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[1]->getEspesorPrin()) * FTubo[1]->getEspesorPrin());
					if(CC[FTubo[1]->getNodoDer() - 1]->getTipoCC() == nmPipesConnection) {
						UnionEspes = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[1]->getNodoDer() - 1])->getEspesor();
						UnionConduct = dynamic_cast<TCCUnionEntreTubos*>(CC[FTubo[1]->getNodoDer() - 1])->getConductividad();
						if(UnionConduct > 0) {
							FR_ext_AxiPos[i] = UnionEspes / UnionConduct / (__cons::Pi * (DIntPrin + FTubo[1]->getEspesorPrin()) *
											   FTubo[1]->getEspesorPrin());
						}
					}
				} else {
					FR_ext_AxiAnt[i] = FTubo[1]->getXRef() / FTubo[1]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[1]->getEspesorPrin()) * FTubo[1]->getEspesorPrin());
					FR_ext_AxiPos[i] = FTubo[1]->getXRef() / FTubo[1]->getConductPrin() / (__cons::Pi *
									   (DIntPrin + FTubo[1]->getEspesorPrin()) * FTubo[1]->getEspesorPrin());
				}

				//C�lculo de las capacidades t�rmicas.
				EsInterior = true;
				Dint = FDiametroExtGap;
				FCapExtInt[i] = 0.;
				FCapExtMed[i] = 0.;
				FCapExtExt[i] = 0.;
				for(int j = 0; j < FTubo[1]->getNumCapas(); j++) {
					Dext = Dint + 2 * FTubo[1]->GetCapa(j).Espesor;
					if(FTubo[1]->GetCapa(j).EsPrincipal == false) {
						if(EsInterior) {
							//C�lculo de la capacidad t�rmica interior.
							Cap = FTubo[1]->GetCapa(j).Density * FTubo[1]->GetCapa(j).CalorEspecifico * FTubo[1]->getXRef();
							FCapExtInt[i] += Cap;
						} else {
							//C�lculo de la capacidad t�rmica exterior.
							Cap = FTubo[1]->GetCapa(j).Density * FTubo[1]->GetCapa(j).CalorEspecifico * __geom::Ring_area(Dint,
									Dext) * FTubo[1]->getXRef();
							FCapExtExt[i] += Cap;
						}
					} else {
						//C�lculo de la capacidad t�rmica exterior, media e interior de la capa principal.
						FCapExtInt[i] += FTubo[1]->getDensidadPrin() * FTubo[1]->getCalEspPrin() * __geom::Ring_area(DIntPrin,
										 DIntPrin + 0.5 * FTubo[1]->getEspesorPrin()) * FTubo[1]->getXRef();
						FCapExtMed[i] = FTubo[1]->getDensidadPrin() * FTubo[1]->getCalEspPrin()
										* __geom::Ring_area(DIntPrin + 0.5 * FTubo[1]->getEspesorPrin(),
															DIntPrin + 1.5 * FTubo[1]->getEspesorPrin()) * FTubo[1]->getXRef();
						FCapExtExt[i] += FTubo[1]->getDensidadPrin() * FTubo[1]->getCalEspPrin()
										 * __geom::Ring_area(DIntPrin + 1.5 * FTubo[1]->getEspesorPrin(),
															 DIntPrin + 2 * FTubo[1]->getEspesorPrin()) * FTubo[1]->getXRef();
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

#pragma package(smart_init)
