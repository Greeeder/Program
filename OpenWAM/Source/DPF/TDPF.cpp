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

#include <iostream>
#ifdef __BORLANDC__
#include <vcl.h>
#endif
#include <cmath>

#include "TDPF.h"
#include "TBloqueMotor.h"
#include "TTubo.h"
#include "TConcentrico.h"
#include "TCondicionContorno.h"
#include "TCCDeposito.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TDPF::TDPF(int numeroDPF, TBloqueMotor **Motor, int NumeroEspecies) {

	FNumeroDPF = numeroDPF;

	FNumeroEspecies = NumeroEspecies;

	FCanal = NULL;
	FTuboEntradaDPF = NULL;
	FTuboSalidaDPF = NULL;

	FFactorFrecuenciaO2 = 2.8e-2;
	FFactorFrecuenciaNO2 = 5e-1;
	FEnergiaActO2 = 150000;
	FEnergiaActNO2 = 73300;

	if(Motor != NULL) {
		FAnguloTotalCiclo = Motor[0]->getAngTotalCiclo();
	} else if(Motor == NULL) {
		FAnguloTotalCiclo = 720.;
	}

	FResMediosDPF = NULL;

	FEspesorSoot = NULL;
	FEspesorSootIn = NULL;
	FEspesorSootDep = NULL;
	FVelocidadPared = NULL;
	FVelocidadPared0 = NULL;
	FKwall = NULL;
	FKwallClean = NULL;
	FKwallLoaded = NULL;
	FKsoot = NULL;
	FKsootIn = NULL;
	FKsootDep = NULL;
	FPorosidad = NULL;
	FDiametroUnidadColectora = NULL;
	FShapeFactor = NULL;
	FWallSootMass = NULL;
	FParametroIntercepcion = NULL;
	FEficiencia = NULL;
	FEficienciaPL = NULL;
	FEficienciaBrown = NULL;
	FEficienciaInter = NULL;
	FEficienciaIner = NULL;
	FCoeficienteParticion = NULL;
	FNumeroUnidadesCelulares = NULL;
	FLongitudVC = NULL;
	FMasaSootUC = NULL;
	FMasaSootNodo = NULL;
	FKreg1 = NULL;
	FKreg2 = NULL;
	FR1 = NULL;
	FR2 = NULL;
	FR3 = NULL;
	FR4 = NULL;
	FR5 = NULL;
	FSupActCatalizador = NULL;
	FSupActSoot = NULL;
	FSupEspecifica = NULL;
	FTasaFraccionMasicaEspecie = NULL;
	FFraccionMasicaEspecieSalida = NULL;
	FFraccionMasicaEspecieEntrantePared = NULL;
	FFraccionMasicaNO2Salida = NULL;
	FFraccionMasicaNO2Entrada = NULL;
	FQreg = NULL;
	FQ1 = NULL;
	FQ2 = NULL;
	FAreaVCCanalEntrada = NULL;
	FAreaVCCanalSalida = NULL;
	FNodoIzq = NULL;
	FNodoDer = NULL;

	FResistRadial = NULL;
	FResistConv = NULL;
	FResistAxialAnt = NULL;
	FResistAxialPost = NULL;
	FResistEntreHacesAnt = NULL;
	FResistEntreHacesPost = NULL;
	FCapEntrada = NULL;
	FCapPared = NULL;
	FCapSalida = NULL;
	FDiametroHazInt = NULL;
	FDiametroHazExt = NULL;
	FDiametroInicialHaz = NULL;
	FTPared = NULL;
	FTParedAnt = NULL;
	FSUMTPPromedio = NULL;
	FSUMTPPromedioConc = NULL;
	FSUMTime = NULL;
	FResistConduccionAislante = NULL;
	FResistConduccionAire = NULL;
	FResistRadiacionAire = NULL;
	FResistConduccionMetal = NULL;
	FResistConveccionExt = NULL;
	FResistAxialMetalExtAnt = NULL;
	FResistAxialMetalExtPost = NULL;
	FCapNodoExteriorSuperficie = NULL;
	FCapNodoMedioSuperficie = NULL;
	FCapNodoInteriorSuperficie = NULL;
	FTSuperficie = NULL;
	FNumCanalesHaz = NULL;
	FSuperficieTCHazAnt = NULL;
	FSuperficieTCHazPost = NULL;
	FAreaInternaHaz = NULL;

	FRg_int_ext = NULL;
	FR_int_radiacion = NULL;
	FR_ext_RadInt = NULL;

	FVolumenTotal = NULL;
	FNumeroUnidadesCelularesHaz = NULL;
	FMasaSootSaturacionHaz = NULL;
	FMasaSootInicialHaz = NULL;
	FBeamSootMass = NULL;
	FLayerSootMass = NULL;
	FLayerSootMassDep = NULL;
	FWallSootMass = NULL;
	FCVSootMass = NULL;
	FMasaSootLayerTotal = NULL;
	FDiametroPoro = NULL;
	FPorcentajeSootNodo = NULL;
	FPorcentajeSootBeam = NULL;
	FPorosidadCapaParticulas_in = 0.;
	FDensidadSootCapa_in = 0.;
	FPorosidadCapaParticulas_dep = 0.;
	FDensidadSootCapa_dep = 0.;

	FTime0DPF = 0.;
	FTime1DPF = 0.;

	FFraccionNOxTabla = NULL;
	FTemperaturasTabla = NULL;
	FMapa_ProporcionNO2 = NULL;

	FCicloActual = 0;
	FEspesorAire = 0.002;

	FDPFSootMass = 0.;
	FCoefAjusTC = 0.02;
	FHayConcentrico = false;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TDPF::~TDPF() {

	for(int i = 0; i < FNumeroHacesCanales; i++) {
		for(int j = 0; j < FCanal[i][0]->getNin(); j++) {
			for(int k = 0; k < 3; k++) {
				delete[] FSUMTPPromedio[i][j][k];
			}
		}
	}

	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 3; j++) {
			delete[] FSUMTPPromedioConc[i][j];
		}
		delete[] FSUMTPPromedioConc[i];
	}

	for(int i = 0; i < FCanal[0][0]->getNin(); i++) {
		delete[] FTSuperficie[i];
	}

	for(int i = 0; i < FNumeroHacesCanales; i++) {
		for(int j = 0; j < 2; j++) {
			delete[] FVelocidadPared[i][j];
			delete[] FVelocidadPared0[i][j];

		}
		for(int j = 0; j < FCanal[i][0]->getNin(); j++) {
			delete[] FTasaFraccionMasicaEspecie[i][j];
			delete[] FFraccionMasicaEspecieSalida[i][j];
			delete[] FFraccionMasicaEspecieEntrantePared[i][j];
			delete[] FTPared[i][j];
			delete[] FTParedAnt[i][j];
			delete[] FResistRadial[i][j];
			delete[] FResistConv[i][j];
			delete[] FSUMTPPromedio[i][j];
		}
	}

	for(int i = 0; i < FNumeroHacesCanales; i++) {
		delete[] FCanal[i];
		delete[] FEspesorSoot[i];
		delete[] FEspesorSootIn[i];
		delete[] FEspesorSootDep[i];
		delete[] FKwall[i];
		delete[] FKwallClean[i];
		delete[] FKwallLoaded[i];
		delete[] FKsoot[i];
		delete[] FKsootIn[i];
		delete[] FKsootDep[i];
		delete[] FPorosidad[i];
		delete[] FDiametroUnidadColectora[i];
		delete[] FShapeFactor[i];
		delete[] FCVSootMass[i];
		delete[] FParametroIntercepcion[i];
		delete[] FEficiencia[i];
		delete[] FEficienciaPL[i];
		delete[] FEficienciaBrown[i];
		delete[] FEficienciaInter[i];
		delete[] FEficienciaIner[i];
		delete[] FCoeficienteParticion[i];
		delete[] FLayerSootMass[i];
		delete[] FLayerSootMassDep[i];
		delete[] FWallSootMass[i];
		delete[] FNumeroUnidadesCelulares[i];
		delete[] FDiametroPoro[i];
		delete[] FLongitudVC[i];
		delete[] FMasaSootUC[i];
		delete[] FMasaSootNodo[i];
		delete[] FKreg1[i];
		delete[] FKreg2[i];
		delete[] FR1[i];
		delete[] FR2[i];
		delete[] FR3[i];
		delete[] FR4[i];
		delete[] FR5[i];
		delete[] FSupActCatalizador[i];
		delete[] FSupActSoot[i];
		delete[] FSupEspecifica[i];
		delete[] FQreg[i];
		delete[] FQ1[i];
		delete[] FQ2[i];
		delete[] FAreaVCCanalEntrada[i];
		delete[] FAreaVCCanalSalida[i];
		delete[] FResistRadial[i];
		delete[] FResistConv[i];
		delete[] FResistAxialAnt[i];
		delete[] FResistAxialPost[i];
		delete[] FResistEntreHacesAnt[i];
		delete[] FResistEntreHacesPost[i];
		delete[] FCapEntrada[i];
		delete[] FCapPared[i];
		delete[] FCapSalida[i];
		delete[] FTPared[i];
		delete[] FTParedAnt[i];
		delete[] FSUMTPPromedio[i];
		delete[] FVelocidadPared[i];
		delete[] FVelocidadPared0[i];
		delete[] FTasaFraccionMasicaEspecie[i];
		delete[] FFraccionMasicaEspecieSalida[i];
		delete[] FFraccionMasicaEspecieEntrantePared[i];
		delete[] FFraccionMasicaNO2Salida[i];
		delete[] FFraccionMasicaNO2Entrada[i];
		if(FCalculoDistribucionSoot) {
			delete[] FPorcentajeSootNodo[i];
		}
		delete[] FMasaSootLimitInfUC[i];
	}
	delete[] FCanal;
	delete[] FEspesorSoot;
	delete[] FEspesorSootIn;
	delete[] FEspesorSootDep;
	delete[] FKwall;
	delete[] FKwallClean;
	delete[] FKwallLoaded;
	delete[] FKsoot;
	delete[] FKsootIn;
	delete[] FKsootDep;
	delete[] FPorosidad;
	delete[] FDiametroUnidadColectora;
	delete[] FCVSootMass;
	delete[] FShapeFactor;
	delete[] FParametroIntercepcion;
	delete[] FEficiencia;
	delete[] FEficienciaPL;
	delete[] FEficienciaBrown;
	delete[] FEficienciaInter;
	delete[] FEficienciaIner;
	delete[] FCoeficienteParticion;
	delete[] FLayerSootMass;
	delete[] FLayerSootMassDep;
	delete[] FWallSootMass;
	delete[] FNumeroUnidadesCelulares;
	delete[] FDiametroPoro;
	delete[] FLongitudVC;
	delete[] FMasaSootUC;
	delete[] FMasaSootNodo;
	delete[] FKreg1;
	delete[] FKreg2;
	delete[] FR1;
	delete[] FR2;
	delete[] FR3;
	delete[] FR4;
	delete[] FR5;
	delete[] FSupActCatalizador;
	delete[] FSupActSoot;
	delete[] FSupEspecifica;
	delete[] FQreg;
	delete[] FQ1;
	delete[] FQ2;
	delete[] FAreaVCCanalEntrada;
	delete[] FAreaVCCanalSalida;
	delete[] FResistAxialAnt;
	delete[] FResistAxialPost;
	delete[] FResistEntreHacesAnt;
	delete[] FResistEntreHacesPost;
	delete[] FCapEntrada;
	delete[] FCapPared;
	delete[] FCapSalida;
	delete[] FResistConduccionAislante;
	delete[] FResistConduccionAire;
	delete[] FResistRadiacionAire;
	delete[] FResistConduccionMetal;
	delete[] FResistConveccionExt;
	delete[] FResistAxialMetalExtAnt;
	delete[] FResistAxialMetalExtPost;
	delete[] FCapNodoExteriorSuperficie;
	delete[] FCapNodoMedioSuperficie;
	delete[] FCapNodoInteriorSuperficie;
	delete[] FNumCanalesHaz;
	delete[] FSuperficieTCHazAnt;
	delete[] FSuperficieTCHazPost;
	delete[] FAreaInternaHaz;
	delete[] FTSuperficie;
	if(FCalculoDistribucionSoot) {
		delete[] FPorcentajeSootBeam;
		delete[] FPorcentajeSootNodo;
	}
	delete[] FRg_int_ext;
	delete[] FR_int_radiacion;
	delete[] FR_ext_RadInt;

	delete[] FTPared;
	delete[] FTParedAnt;
	delete[] FVelocidadPared;
	delete[] FVelocidadPared0;
	delete[] FTasaFraccionMasicaEspecie;
	delete[] FFraccionMasicaEspecieSalida;
	delete[] FFraccionMasicaEspecieEntrantePared;
	delete[] FFraccionMasicaNO2Salida;
	delete[] FFraccionMasicaNO2Entrada;
	delete[] FResistRadial;
	delete[] FResistConv;

	delete[] FSUMTPPromedio;
	delete[] FSUMTPPromedioConc;

	delete[] FNodoIzq;
	delete[] FNodoDer;
	delete[] FSUMTime;
	delete[] FDiametroHazInt;
	delete[] FDiametroHazExt;
	delete[] FDiametroInicialHaz;

	delete[] FNumeroUnidadesCelularesHaz;
	delete[] FVolumenTotal;
	delete[] FMasaSootSaturacionHaz;
	delete[] FMasaSootInicialHaz;
	delete[] FBeamSootMass;
	delete[] FMasaSootLayerTotal;
	delete[] FMasaSootLimitInfUC;

	for(int i = 0; i < FNumeroDatos_Temperaturas; i++) {
		delete[] FMapa_ProporcionNO2[i];
	}
	delete[] FMapa_ProporcionNO2;
	delete[] FFraccionNOxTabla;
	delete[] FTemperaturasTabla;

	if(FResMediosDPF != NULL) {
		for(int i = 0; i < FNumResMedios; i++) {
			delete[] FResMediosDPF[i];
		}
		delete[] FResMediosDPF;
	}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::LeeDatosDPF(const char *FileWAM, fpos_t &filepos, nmTipoCalculoEspecies CalculoEspecies,
					   nmCalculoGamma CalculoGamma, bool HayEGR, TBloqueMotor **Motor) {
	try {

		double SupEspecifica, AreaFiltradoDPF, VolumenDPF, FraccionDistribucion = 0.;
		int refrigerante, aislante, DPFCatalizada, CalculoRegeneracion, CalculoFiltrado, CalculoDistribucionSoot,
			HayConcentrico, numbernodos;

		if(HayEGR)
			FIntEGR = 0;
		else
			FIntEGR = 1;

// Canal[0] representa el canal de entrada a la trampa de part�culas.
// Canal[1] representa el canal de salida de la trampa de part�culas.
		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d %d", &FNumeroHacesCanales, &FNumeroCanalesTotales);

		fscanf(fich, "%d %d %d %d", &DPFCatalizada, &CalculoRegeneracion, &CalculoFiltrado, &CalculoDistribucionSoot);

		DPFCatalizada == 0 ? FDPFCatalizada = false : FDPFCatalizada = true;
		CalculoRegeneracion == 0 ? FCalculoRegeneracion = false : FCalculoRegeneracion = true;
		CalculoFiltrado == 0 ? FCalculoFiltrado = false : FCalculoFiltrado = true;
		CalculoDistribucionSoot == 0 ? FCalculoDistribucionSoot = false : FCalculoDistribucionSoot = true;

		if(FCalculoRegeneracion) {
			fscanf(fich, "%lf %lf %lf %lf", &FFactorFrecuenciaO2, &FFactorFrecuenciaNO2, &FEnergiaActO2, &FEnergiaActNO2);
			fscanf(fich, "%lf %lf ", &FIndiceCompletitud1, &FIndiceCompletitud2);

			// Lectura del mapa para determinar la proporci�n de NO2 en los NOx
			fscanf(fich, "%d %d ", &FNumeroDatos_FraccionesNOx, &FNumeroDatos_Temperaturas);
			FFraccionNOxTabla = new double[FNumeroDatos_FraccionesNOx];
			FTemperaturasTabla = new double[FNumeroDatos_Temperaturas];
			FMapa_ProporcionNO2 = new double*[FNumeroDatos_Temperaturas];
			for(int i = 0; i < FNumeroDatos_Temperaturas; i++) {
				FMapa_ProporcionNO2[i] = new double[FNumeroDatos_FraccionesNOx];
			}
			for(int i = 0; i < FNumeroDatos_FraccionesNOx; i++) {
				fscanf(fich, "%lf ", &FFraccionNOxTabla[i]);
			}
			for(int i = 0; i < FNumeroDatos_Temperaturas; i++) {
				fscanf(fich, "%lf ", &FTemperaturasTabla[i]);
			}
			for(int i = 0; i < FNumeroDatos_Temperaturas; i++) {
				for(int j = 0; j < FNumeroDatos_FraccionesNOx; j++) {
					fscanf(fich, "%lf ", &FMapa_ProporcionNO2[i][j]);
				}
			}
		}

// Geometr�a
		fscanf(fich, "%lf %lf %lf ", &FDiametroFiltroEfect, &FLongitudEfec, &FLongitudTapon);

		FDiametroInicialHaz = new double[FNumeroHacesCanales];
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			fscanf(fich, "%lf ", &FDiametroInicialHaz[i]);
		}

		FDiametroHazExt = new double[FNumeroHacesCanales];
		FDiametroHazInt = new double[FNumeroHacesCanales];

// Datos adicionales
		fscanf(fich, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &FEspesorParedPorosa, &FDiametroMedioPoroso,
			   &FPorosidadLimpia, &FPorosidadCapaParticulas_in,
			   &FPorosidadCapaParticulas_dep, &FDensidadSootPared, &FDensidadSootCapa_in, &FDensidadSootCapa_dep,
			   &FDiametroPrimarioSoot, &FDiametroAgregadoSoot, &FMasaSootInicial, &FFactorPercolacion,
			   &FFraccionParedSaturada, &FMasaSootWallTotalMax);
		if(FCalculoFiltrado) {
			fscanf(fich, "%lf %lf %lf ", &FaSF, &FbSF, &FMinShapeFactor);
			fscanf(fich, "%lf", &FStickingCoef);
			fscanf(fich, "%lf", &FFiltRelaxCoef);
		} else {
			fscanf(fich, "%lf", &FShapeFactorDiscrete);
		}
		if(FCalculoDistribucionSoot) {
			FPorcentajeSootNodo = new double*[FNumeroHacesCanales];
			FPorcentajeSootBeam = new double[FNumeroHacesCanales];
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				fscanf(fich, "%lf ", &FPorcentajeSootBeam[j]);
				fscanf(fich, "%d ", &numbernodos);
				FPorcentajeSootNodo[j] = new double[numbernodos];

				for(int i = 0; i < numbernodos; i++) {
					fscanf(fich, "%lf", &FPorcentajeSootNodo[j][i]); //introduccion de distribucion en [%]
					FraccionDistribucion += FPorcentajeSootNodo[j][i];
				}
				if(FraccionDistribucion != 100.) {
					std::cout << "ERROR: La fraccion de distribucion total no puede ser distinta de 100.% Repasa la lectura en el haz  " <<
							  j + 1 << std::endl;
					throw Exception(" ");
				}
			}
		}
		fscanf(fich, "%lf %lf %d %lf", &FTIniPared, &FTIniParedExt, &FTctpt, &FCoefAjusTC);

		switch(FTctpt) {
		case 0:
			FTipoCalcTempPared = nmVariableConInerciaTermica;
			break;
		case 1:
			FTipoCalcTempPared = nmVariableSinInerciaTermica;
			break;
		case 2:
			FTipoCalcTempPared = nmTempConstante;
			break;
		}

		if(FTipoCalcTempPared != nmTempConstante) {
			if(Motor == NULL) {
				fscanf(fich, "%lf %d ", &FDuracionCiclo, &FNumCiclosSinInerciaTermica);
			}
			fscanf(fich, "%lf %lf ", &FAjustRAxAnt, &FAjustRAxPos);
			fscanf(fich, "%lf %lf %d ", &FCoefExt, &FEmisividad, &refrigerante);
			switch(refrigerante) {
			case 0:
				FTipoRefrig = nmAire;
				break;
			case 1:
				FTipoRefrig = nmAgua;
				break;
			}
			if(FTipoRefrig == nmAgua) {
				fscanf(fich, "%lf ", &FTRefrigerante); /* Esta en �C */
			}

			// Propiedades de los materiales en la trampa de part�culas
			fscanf(fich, "%lf %lf %lf %lf ", &FDensidadParedPorosa, &FCalorEspPared, &FConductividadPared,
				   &FConductividadParedRadial);
			fscanf(fich, "%lf %lf %lf %lf", &FEspesorAire, &FConductividadAire, &FEmisividadIntAire, &FEmisividadExtAire);
			fscanf(fich, "%lf %lf ", &FCalorEspSoot, &FConductividadSoot);
			fscanf(fich, "%lf %lf %lf %lf ", &FEspesorMetal, &FDensidadMetal, &FCalorEspMetal, &FConductividadMetal);

			fscanf(fich, "%d ", &aislante);
			aislante == 0 ? FAislante = false : FAislante = true;
			if(FAislante) {
				fscanf(fich, "%lf %lf %lf %lf ", &FEspesorAislante, &FDensidadAislante, &FCalorEspAislante, &FConductividadAislante);
			}
		}

		fgetpos(fich, &filepos);
		fclose(fich);

		FCanal = new TCanalDPF**[FNumeroHacesCanales];
		FNodoIzq = new int[FNumeroHacesCanales];
		FNodoDer = new int[FNumeroHacesCanales];

		for(int i = 0; i < FNumeroHacesCanales; i++) {
			FCanal[i] = new TCanalDPF*[2];
			FCanal[i][0] = new TCanalDPF(FNumeroEspecies, i, CalculoEspecies, CalculoGamma, HayEGR, this, 0, FNumeroDPF);
			FCanal[i][1] = new TCanalDPF(FNumeroEspecies, i, CalculoEspecies, CalculoGamma, HayEGR, this, 1, FNumeroDPF);
		}
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			fich = fopen(FileWAM, "r");
			fsetpos(fich, &filepos);

			fscanf(fich, "%d %d ", &FNodoIzq[i], &FNodoDer[i]);
			fgetpos(fich, &filepos);
			fclose(fich);

			FCanal[i][0]->LeeDatosGeneralesCanal(FileWAM, filepos, FNodoIzq[i], FNodoDer[i]);   // Canal de entrada
			FCanal[i][0]->LeeDatosGeometricosCanal(FileWAM, filepos); // Canal de entrada
			FCanal[i][1]->LeeDatosGeneralesCanal(FileWAM, filepos, FNodoIzq[i], FNodoDer[i]);   // Canal de salida
			FCanal[i][1]->LeeDatosGeometricosCanal(FileWAM, filepos); // Canal de salida
		}

		FVelocidadPared = new double**[FNumeroHacesCanales]; // Para cada haz de canales considerado
		FVelocidadPared0 = new double**[FNumeroHacesCanales];
		FTasaFraccionMasicaEspecie = new double**[FNumeroHacesCanales];
		FFraccionMasicaEspecieSalida = new double**[FNumeroHacesCanales];
		FFraccionMasicaEspecieEntrantePared = new double**[FNumeroHacesCanales];
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			FVelocidadPared[i] = new double*[2]; // Para el canal de entrada y el de salida
			FVelocidadPared0[i] = new double*[2];
			FTasaFraccionMasicaEspecie[i] = new double*[FCanal[i][0]->getNin()]; // Se define en cada node del canal de entrada
			FFraccionMasicaEspecieSalida[i] = new double*[FCanal[i][0]->getNin()]; // Se define en cada node del canal de entrada
			FFraccionMasicaEspecieEntrantePared[i] = new double*[FCanal[i][0]->getNin()];
		}
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			for(int j = 0; j < 2; j++) {
				FVelocidadPared[i][j] = new double[FCanal[i][0]->getNin()]; // Para cada nodo
				FVelocidadPared0[i][j] = new double[FCanal[i][0]->getNin()]; // Para cada nodo
			}
			for(int j = 0; j < FCanal[i][0]->getNin(); j++) {
				FTasaFraccionMasicaEspecie[i][j] = new double[FNumeroEspecies]; // Para cada especie qu�mica transportada
				FFraccionMasicaEspecieSalida[i][j] = new double[FNumeroEspecies]; // Para cada especie qu�mica transportada
				FFraccionMasicaEspecieEntrantePared[i][j] = new double[FNumeroEspecies];
			}
		}

		FNumeroUnidadesCelularesHaz = new long long int[FNumeroHacesCanales];
		FVolumenTotal = new double[FNumeroHacesCanales];
		FMasaSootSaturacionHaz = new double[FNumeroHacesCanales];
		FMasaSootInicialHaz = new double[FNumeroHacesCanales];
		FMasaSootLayerTotal = new double[FNumeroHacesCanales];
		FBeamSootMass = new double[FNumeroHacesCanales];
		FCVSootMass = new double*[FNumeroHacesCanales];
		FLayerSootMass = new double*[FNumeroHacesCanales];
		FLayerSootMassDep = new double*[FNumeroHacesCanales];
		FWallSootMass = new double*[FNumeroHacesCanales];

		FDiametroPoro = new double*[FNumeroHacesCanales];
		FKwall = new double*[FNumeroHacesCanales];
		FKwallClean = new double*[FNumeroHacesCanales];
		FKwallLoaded = new double*[FNumeroHacesCanales];
		FKsoot = new double*[FNumeroHacesCanales];
		FKsootIn = new double*[FNumeroHacesCanales];
		FKsootDep = new double*[FNumeroHacesCanales];
		FPorosidad = new double*[FNumeroHacesCanales];
		FEspesorSoot = new double*[FNumeroHacesCanales];
		FEspesorSootIn = new double*[FNumeroHacesCanales];
		FEspesorSootDep = new double*[FNumeroHacesCanales];
		FDiametroUnidadColectora = new double*[FNumeroHacesCanales];
		FShapeFactor = new double*[FNumeroHacesCanales];
		FParametroIntercepcion = new double*[FNumeroHacesCanales];
		FEficiencia = new double*[FNumeroHacesCanales];
		FEficienciaPL = new double*[FNumeroHacesCanales];
		FEficienciaBrown = new double*[FNumeroHacesCanales];
		FEficienciaInter = new double*[FNumeroHacesCanales];
		FEficienciaIner = new double*[FNumeroHacesCanales];
		FCoeficienteParticion = new double*[FNumeroHacesCanales];
		FNumeroUnidadesCelulares = new long long int*[FNumeroHacesCanales];
		FLongitudVC = new double*[FNumeroHacesCanales];
		FMasaSootUC = new double*[FNumeroHacesCanales];
		FMasaSootNodo = new double*[FNumeroHacesCanales];
		FKreg1 = new double*[FNumeroHacesCanales];
		FKreg2 = new double*[FNumeroHacesCanales];
		FR1 = new double*[FNumeroHacesCanales];
		FR2 = new double*[FNumeroHacesCanales];
		FR3 = new double*[FNumeroHacesCanales];
		FR4 = new double*[FNumeroHacesCanales];
		FR5 = new double*[FNumeroHacesCanales];
		FSupActCatalizador = new double*[FNumeroHacesCanales];
		FSupActSoot = new double*[FNumeroHacesCanales];
		FSupEspecifica = new double*[FNumeroHacesCanales];
		FQreg = new double*[FNumeroHacesCanales];
		FQ1 = new double*[FNumeroHacesCanales];
		FQ2 = new double*[FNumeroHacesCanales];
		FAreaVCCanalEntrada = new double*[FNumeroHacesCanales];
		FAreaVCCanalSalida = new double*[FNumeroHacesCanales];
		FFraccionMasicaNO2Salida = new double*[FNumeroHacesCanales];
		FFraccionMasicaNO2Entrada = new double*[FNumeroHacesCanales];
		FMasaSootLimitInfUC = new double*[FNumeroHacesCanales];
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			FKwall[i] = new double[FCanal[i][0]->getNin()];
			FKwallClean[i] = new double[FCanal[i][0]->getNin()];
			FKwallLoaded[i] = new double[FCanal[i][0]->getNin()];
			FDiametroPoro[i] = new double[FCanal[i][0]->getNin()];
			FKsoot[i] = new double[FCanal[i][0]->getNin()];
			FKsootIn[i] = new double[FCanal[i][0]->getNin()];
			FKsootDep[i] = new double[FCanal[i][0]->getNin()];
			FPorosidad[i] = new double[FCanal[i][0]->getNin()];
			FEspesorSoot[i] = new double[FCanal[i][0]->getNin()];
			FEspesorSootIn[i] = new double[FCanal[i][0]->getNin()];
			FEspesorSootDep[i] = new double[FCanal[i][0]->getNin()];
			FDiametroUnidadColectora[i] = new double[FCanal[i][0]->getNin()];
			FShapeFactor[i] = new double[FCanal[i][0]->getNin()];
			FCVSootMass[i] = new double[FCanal[i][0]->getNin()];
			FLayerSootMass[i] = new double[FCanal[i][0]->getNin()];
			FLayerSootMassDep[i] = new double[FCanal[i][0]->getNin()];
			FWallSootMass[i] = new double[FCanal[i][0]->getNin()];
			FParametroIntercepcion[i] = new double[FCanal[i][0]->getNin()];
			FEficiencia[i] = new double[FCanal[i][0]->getNin()];
			FEficienciaPL[i] = new double[FCanal[i][0]->getNin()];
			FEficienciaBrown[i] = new double[FCanal[i][0]->getNin()];
			FEficienciaInter[i] = new double[FCanal[i][0]->getNin()];
			FEficienciaIner[i] = new double[FCanal[i][0]->getNin()];
			FCoeficienteParticion[i] = new double[FCanal[i][0]->getNin()];
			FNumeroUnidadesCelulares[i] = new long long int[FCanal[i][0]->getNin()];
			FLongitudVC[i] = new double[FCanal[i][0]->getNin()];
			FMasaSootUC[i] = new double[FCanal[i][0]->getNin()];
			FMasaSootNodo[i] = new double[FCanal[i][0]->getNin()];
			FKreg1[i] = new double[FCanal[i][0]->getNin()];
			FKreg2[i] = new double[FCanal[i][0]->getNin()];
			FR1[i] = new double[FCanal[i][0]->getNin()];
			FR2[i] = new double[FCanal[i][0]->getNin()];
			FR3[i] = new double[FCanal[i][0]->getNin()];
			FR4[i] = new double[FCanal[i][0]->getNin()];
			FR5[i] = new double[FCanal[i][0]->getNin()];
			FSupActCatalizador[i] = new double[FCanal[i][0]->getNin()];
			FSupActSoot[i] = new double[FCanal[i][0]->getNin()];
			FSupEspecifica[i] = new double[FCanal[i][0]->getNin()];
			FQreg[i] = new double[FCanal[i][0]->getNin()];
			FQ1[i] = new double[FCanal[i][0]->getNin()];
			FQ2[i] = new double[FCanal[i][0]->getNin()];
			FAreaVCCanalEntrada[i] = new double[FCanal[i][0]->getNin()];
			FAreaVCCanalSalida[i] = new double[FCanal[i][0]->getNin()];
			FFraccionMasicaNO2Entrada[i] = new double[FCanal[i][0]->getNin()];
			FFraccionMasicaNO2Salida[i] = new double[FCanal[i][0]->getNin()];
			FMasaSootLimitInfUC[i] = new double[FCanal[i][0]->getNin()];
		}

		VolumenDPF = __geom::Cylinder_volume(FDiametroFiltroEfect, FLongitudEfec); // Volumen de la trampa de part�culas
		AreaFiltradoDPF = 2. * FNumeroCanalesTotales * FCanal[0][0]->GetLadoCanal(0) * FLongitudEfec;
		SupEspecifica = AreaFiltradoDPF / VolumenDPF;

		FDiametroPoroPL_in = 2 * FPorosidadCapaParticulas_in * FDiametroAgregadoSoot / 3. / (1. - FPorosidadCapaParticulas_in);
		FDiametroPoroPL_dep = 2 * FPorosidadCapaParticulas_dep * FDiametroAgregadoSoot / 3. /
							  (1. - FPorosidadCapaParticulas_dep);
		for(int i = 0; i < FNumeroHacesCanales; i++) {
			for(int j = 0; j < FCanal[i][0]->getNin(); j++) {
				FKsoot[i][j] = 2e-14;
				FSupEspecifica[i][j] = SupEspecifica;
				FEspesorSootIn[i][j] = 0.;
				FEspesorSootDep[i][j] =
					0.; // inicializo a 0 el espesor del soot de ambas las contribuciones (inicial y de deposicion sucesiva) para poder luego sumarlas
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::LeeDatos DPF en la trampa de part�culas: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoVelocidadPared() {
	try {
		double PresionCanalSalida, DensidadCanalSalida, VelocidadParedCanalEntrada, DensidadCanalEntrada, LadoCanalEntrada,
			   EspesorSoot, funcion_f;
		double Knudsen, CaminoLibreMedio, SCF, Temperatura, ViscosidadCinematica, PM, Knudsen_limpia, SCF_limpia;
		double beta = 0., a = 0., b = 0., c = 0.;

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				FVelocidadPared0[j][0][i] = FVelocidadPared[j][0][i];
				if(FCanal[j][0]->getNodoInicialFuente() > i) {
					FVelocidadPared[j][0][i] = 0.;
				} else {
					if(FCanal[j][0]->getNodoInicialFuente() == 0) {
						PresionCanalSalida = FCanal[j][1]->GetPresion(i);
						DensidadCanalSalida = FCanal[j][1]->GetDensidad(i);
					} else {
						if(i == FCanal[j][0]->getNin() - 1) {
							PresionCanalSalida = FCanal[j][1]->GetPresion(i - FCanal[j][0]->getNodoInicialFuente())
												 - (FCanal[j][1]->GetPresion(i - FCanal[j][0]->getNodoInicialFuente() - 1) - FCanal[j][1]->GetPresion(
														i - FCanal[j][0]->getNodoInicialFuente()))
												 / FCanal[j][1]->getXRef() * FCanal[j][1]->getDistanciaInterpolacion();
							DensidadCanalSalida = FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente())
												  - (FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente() - 1) - FCanal[j][1]->GetDensidad(
														 i - FCanal[j][0]->getNodoInicialFuente()))
												  / FCanal[j][1]->getXRef() * FCanal[j][1]->getDistanciaInterpolacion();
						} else {
							PresionCanalSalida = Interpola(FCanal[j][1]->GetPresion(i - FCanal[j][0]->getNodoInicialFuente()),
														   FCanal[j][1]->GetPresion(i - FCanal[j][0]->getNodoInicialFuente() + 1),
														   1., FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
							DensidadCanalSalida = Interpola(FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente()),
															FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente() + 1), 1.,
															FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
						}
					}
					ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
					Temperatura = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
									  2.) / FCanal[j][0]->GetGamma(i) / FCanal[j][0]->GetRMezcla(i);
					CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * Temperatura), 0.5);
					Knudsen = 2. * CaminoLibreMedio /
							  FDiametroPoro[j][i]; // Referido al di�metro medio de la pared del filtro de part�culas
					if(Knudsen < 0.0016) {
						SCF = 1 + Knudsen * 1.257;
					} else {
						SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
					}

					if(FWallSootMass != 0) {
						funcion_f = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidad[j][i],
															  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1 - FPorosidad[j][i]);
						Knudsen_limpia = 2. * CaminoLibreMedio /
										 FDiametroMedioPoroso; // Referido al diametro medio de la pared del filtro de part�culas
						if(Knudsen_limpia < 0.0016) {
							SCF_limpia = 1 + Knudsen_limpia * 1.257;
						} else {
							SCF_limpia = 1 + Knudsen_limpia * (1.257 + 0.4 * exp(-1.1 / Knudsen_limpia));
						}
						FKwallClean[j][i] = 0.09 * 2 * (2 - 1.8 * pow((1 - FPorosidadLimpia),
														1. / 3.) - FPorosidadLimpia - 0.2 * pow((1 - FPorosidadLimpia), 2.)) / (9 * (1 - FPorosidadLimpia))
											* pow(FDiametroUnidadColectoraLimpia, 2.) * SCF_limpia;
						FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
												 2.) * funcion_f / Ffuncion_f_Limpia * FKwallClean[j][i] * SCF / SCF_limpia;
						FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
									   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);
					} else {
						FKwallClean[j][i] = 0.09 * 2 * (2 - 1.8 * pow((1 - FPorosidad[j][i]),
														1. / 3.) - FPorosidad[j][i] - 0.2 * pow((1 - FPorosidad[j][i]), 2.)) / (9 * (1 - FPorosidad[j][i]))
											* pow(FDiametroUnidadColectora[j][i], 2.) * SCF;
						FKwall[j][i] = FKwallClean[j][i];
					}
					FVelocidadPared[j][0][i] = fabs(__units::BarToPa(FCanal[j][0]->GetPresion(i) - PresionCanalSalida))
											   / ((FCanal[j][0]->GetViscosidadDinamica(i) * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i])
												   * log(FCanal[j][0]->GetLadoCanal(i) / (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i])) / (2 * FKsoot[j][i]))
												  + (FCanal[j][0]->GetViscosidadDinamica(i) * FEspesorParedPorosa * FCanal[j][0]->GetDensidad(i) *
													 (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i])
													 / (FKwall[j][i] * DensidadCanalSalida * FCanal[j][1]->GetLadoCanal(i))));

					if(FCanal[j][0]->GetPresion(i) < PresionCanalSalida) {
						FVelocidadPared[j][0][i] = -FVelocidadPared[j][0][i];
					}
				}

			}
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				FVelocidadPared0[j][1][i] = FVelocidadPared[j][1][i];
				if(FCanal[j][1]->getNodoFinalFuente() >= i) {
					if(FCanal[j][0]->getNodoInicialFuente() == 0) {
						FVelocidadPared[j][1][i] = FVelocidadPared[j][0][i] * FCanal[j][0]->GetDensidad(i) * (FCanal[j][0]->GetLadoCanal(
													   i) - 2 * FEspesorSoot[j][i])
												   / (FCanal[j][1]->GetDensidad(i) * FCanal[j][1]->GetLadoCanal(i));
					} else {
						if(i == 0) {
							VelocidadParedCanalEntrada = FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente() + i]
														 - (FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente() + i + 1] -
															FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente() + i]) / FCanal[j][0]->getXRef()
														 * FCanal[j][0]->getDistanciaInterpolacion();
							DensidadCanalEntrada = FCanal[j][0]->GetDensidad(FCanal[j][0]->getNodoInicialFuente() + i)
												   - (FCanal[j][0]->GetDensidad(FCanal[j][0]->getNodoInicialFuente() + i + 1) - FCanal[j][0]->GetDensidad(
														  FCanal[j][0]->getNodoInicialFuente() + i))
												   / FCanal[j][0]->getXRef() * FCanal[j][0]->getDistanciaInterpolacion();
							LadoCanalEntrada = FCanal[j][0]->GetLadoCanal(FCanal[j][0]->getNodoInicialFuente() + i)
											   - (FCanal[j][0]->GetLadoCanal(FCanal[j][0]->getNodoInicialFuente() + i + 1) - FCanal[j][0]->GetLadoCanal(
													  FCanal[j][0]->getNodoInicialFuente() + i))
											   / FCanal[j][0]->getXRef() * FCanal[j][0]->getDistanciaInterpolacion();
							EspesorSoot = FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente() + i]
										  - (FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente() + i + 1] - FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente()
												  + i]) / FCanal[j][0]->getXRef()
										  * FCanal[j][0]->getDistanciaInterpolacion();
							FVelocidadPared[j][1][i] = VelocidadParedCanalEntrada * DensidadCanalEntrada * (LadoCanalEntrada - 2 * EspesorSoot)
													   / (FCanal[j][1]->GetDensidad(i) * FCanal[j][1]->GetLadoCanal(i));
						} else {
							VelocidadParedCanalEntrada = Interpola(FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente() - 1 + i],
																   FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente() + i],
																   1., (FCanal[j][0]->getXRef() - FCanal[j][0]->getDistanciaInterpolacion()) / FCanal[j][0]->getXRef());
							DensidadCanalEntrada = Interpola(FCanal[j][0]->GetDensidad(FCanal[j][0]->getNodoInicialFuente() - 1 + i),
															 FCanal[j][0]->GetDensidad(FCanal[j][0]->getNodoInicialFuente() + i), 1.,
															 (FCanal[j][0]->getXRef() - FCanal[j][0]->getDistanciaInterpolacion()) / FCanal[j][0]->getXRef());
							LadoCanalEntrada = Interpola(FCanal[j][0]->GetLadoCanal(FCanal[j][0]->getNodoInicialFuente() - 1 + i),
														 FCanal[j][0]->GetLadoCanal(FCanal[j][0]->getNodoInicialFuente() + i),
														 1., (FCanal[j][0]->getXRef() - FCanal[j][0]->getDistanciaInterpolacion()) / FCanal[j][0]->getXRef());
							EspesorSoot = Interpola(FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente() - 1 + i],
													FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente() + i], 1.,
													(FCanal[j][0]->getXRef() - FCanal[j][0]->getDistanciaInterpolacion()) / FCanal[j][0]->getXRef());
							FVelocidadPared[j][1][i] = VelocidadParedCanalEntrada * DensidadCanalEntrada * (LadoCanalEntrada - 2 * EspesorSoot)
													   / (FCanal[j][1]->GetDensidad(i) * FCanal[j][1]->GetLadoCanal(i));
						}
					}
				} else {
					FVelocidadPared[j][1][i] = 0.;
				}
			}
		}

		/*for(int j=0;j<FNumeroHacesCanales;j++){
		 for(int i=0;i<FCanal[j][0]->getNin();i++){
		 if(FCanal[j][0]->getNodoInicialFuente()>i){
		 FVelocidadPared[j][0][i]=0.;
		 }else{
		 if(FCanal[j][0]->getNodoInicialFuente()==0){
		 PresionCanalSalida=FCanal[j][1]->Presion[i];
		 DensidadCanalSalida=FCanal[j][1]->Densidad[i];
		 }else{
		 PresionCanalSalida=Interpola(FCanal[j][1]->Presion[i-FCanal[j][0]->getNodoInicialFuente()],
		 FCanal[j][1]->Presion[i-FCanal[j][0]->getNodoInicialFuente()+1],1.,FCanal[j][1]->getDistanciaInterpolacion());
		 DensidadCanalSalida=Interpola(FCanal[j][1]->Densidad[i-FCanal[j][0]->getNodoInicialFuente()],
		 FCanal[j][1]->Densidad[i-FCanal[j][0]->getNodoInicialFuente()+1],1.,FCanal[j][1]->getDistanciaInterpolacion());
		 }
		 ViscosidadCinematica=FCanal[j][0]->ViscosidadDinamica[i]/FCanal[j][0]->Densidad[i];
		 Temperatura=pow(FCanal[j][0]->Asonido[i]*__cons::ARef,2.)/FCanal[j][0]->Gamma[i]/FCanal[j][0]->R[i];
		 //CaminoLibreMedio=__cons::Kb*Temperatura/(pow(2,0.5)*Pi*pow(FDiametroMedioPoroso,2.)*FCanal[j][0]->Presion[i]*1e5);
		 //PM=__R::Universal/FCanal[j][0]->R[i];
		 CaminoLibreMedio=ViscosidadCinematica*pow(Pi/(2.*FCanal[j][0]->R[i]*Temperatura),0.5);
		 Knudsen=2.*CaminoLibreMedio/FDiametroMedioPoroso;  // Referido al di�metro medio de la pared del filtro de part�culas
		 if(Knudsen<0.0016){
		 SCF=1+Knudsen*1.257;
		 }else{
		 SCF=1+Knudsen*(1.257+0.4*exp(-1.1/Knudsen));
		 }
		 //FKwall[j][i]=pow(FPorosidad[j][i],5.5)/5.6*pow(FDiametroMedioPoroso,2.)*SCF;     // Rumpf and Gupte
		 FKwall[j][i]=0.09*2*(2-1.8*pow((1-FPorosidad[j][i]),1./3.)-FPorosidad[j][i]-0.2*pow((1-FPorosidad[j][i]),2.))/
		 (9*(1-FPorosidad[j][i]))*pow(FDiametroUnidadColectora[j][i],2.)*SCF;
		 beta=1.75*(1-FPorosidad[j][i])/FDiametroUnidadColectora[j][i]/pow(FPorosidad[j][i],3.);
		 a=beta*FCanal[j][0]->Densidad[i]*FEspesorParedPorosa*pow(FCanal[j][0]->Densidad[i]*(FCanal[j][0]->LadoCanal[i]-2*FEspesorSoot[j][i])/
		 (FCanal[j][1]->Densidad[i]*FCanal[j][0]->LadoCanal[i]),2.);
		 b=FCanal[j][0]->ViscosidadDinamica[i]*FEspesorParedPorosa*FCanal[j][0]->Densidad[i]*
		 (FCanal[j][0]->LadoCanal[i]-2*FEspesorSoot[j][i])/(FKwall[j][i]*
		 DensidadCanalSalida*FCanal[j][1]->LadoCanal[i]);
		 c=-fabs((FCanal[j][0]->Presion[i]-PresionCanalSalida)*1e5);
		 FVelocidadPared[j][0][i]=(-b+pow(pow(b,2.)-4.*a*c,0.5))/2./a;
		 if(FCanal[j][0]->Presion[i]<PresionCanalSalida){
		 FVelocidadPared[j][0][i]=-FVelocidadPared[j][0][i];
		 }
		 }
		 }
		 for(int i=0;i<FCanal[j][0]->getNin();i++){
		 if(FCanal[j][1]->NodoFinalFuente>=i){
		 if(FCanal[j][0]->getNodoInicialFuente()==0){
		 FVelocidadPared[j][1][i]=FVelocidadPared[j][0][i]*FCanal[j][0]->Densidad[i]*
		 (FCanal[j][0]->LadoCanal[i]-2*FEspesorSoot[j][i])/
		 (FCanal[j][1]->Densidad[i]*FCanal[j][1]->LadoCanal[i]);
		 }else{
		 VelocidadParedCanalEntrada=Interpola(FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente()-1+i],
		 FVelocidadPared[j][0][FCanal[j][0]->getNodoInicialFuente()+i],1.,
		 FCanal[j][0]->getXRef()-FCanal[j][0]->getDistanciaInterpolacion());
		 DensidadCanalEntrada=Interpola(FCanal[j][0]->Densidad[FCanal[j][0]->getNodoInicialFuente()-1+i],
		 FCanal[j][0]->Densidad[FCanal[j][0]->getNodoInicialFuente()+i],1.,
		 FCanal[j][0]->getXRef()-FCanal[j][0]->getDistanciaInterpolacion());
		 LadoCanalEntrada=Interpola(FCanal[j][0]->LadoCanal[FCanal[j][0]->getNodoInicialFuente()-1+i],
		 FCanal[j][0]->LadoCanal[FCanal[j][0]->getNodoInicialFuente()+i],1.,
		 FCanal[j][0]->getXRef()-FCanal[j][0]->getDistanciaInterpolacion());
		 EspesorSoot=Interpola(FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente()-1+i],
		 FEspesorSoot[j][FCanal[j][0]->getNodoInicialFuente()+i],1.,
		 FCanal[j][0]->getXRef()-FCanal[j][0]->getDistanciaInterpolacion());
		 FVelocidadPared[j][1][i]=VelocidadParedCanalEntrada*DensidadCanalEntrada*(LadoCanalEntrada-2*EspesorSoot)/
		 (FCanal[j][1]->Densidad[i]*FCanal[j][1]->LadoCanal[i]);
		 }
		 }else{
		 FVelocidadPared[j][1][i]=0.;
		 }
		 }
		 } */
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoVelocidadPared en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::InicializaDPF(int nconcentrico, TConcentrico **Concentrico) {
	try {
		double funcion_g, VelocidadPared, VelocidadIntersticial, ViscosidadCinematica, Temperatura, CaminoLibreMedio, Knudsen,
			   SCF, CoefDifusionParticulas, Peclet, EficienciaBrowniana, funcion_f,
			   EficienciaIntercepcion, Cs, Stokes, EficienciaInercial, EficienciaCombinada, VolumenControlPared, DensidadIni,
			   ViscosidadDinamica, SootFactor;

		double SCFPL, KnudsenPL, PecletPL, CoefDifusionParticulasPL, EficienciaBrownianaPL, StokesPL, EficienciaInercialPL,
			   EficienciaCombinadaPL, VelocidadIntersticialPL;

		long long int UnidadesCelularesCanal = 0, UnidadesCelulares = 0;

#ifdef ConcentricElement
		if(FHayConcentrico) {
			for(int i = 0; i < nconcentrico; i++) {
				if(Concentrico[i]->GetNumDPF() == FNumeroDPF) {
					FNumTuboExt = Concentrico[i]->GetNumTuboExterno();
					FNumConcentrico = i;
				}
			}
		}
#endif
// Previamiente han de inicializarse los canales.

// Inicializac��n variables de la regeneraci�n
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				if(FCanal[j][0]->getNodoInicialFuente() > i) {
					if(i == 0) {
						if(i + 1 == FCanal[j][0]->getNodoInicialFuente()) {
							FLongitudVC[j][i] = FCanal[j][0]->getXRef() - FCanal[j][0]->getDistanciaInterpolacion();
						} else {
							FLongitudVC[j][i] = FCanal[j][0]->getXRef() / 2.;
						}
					} else if(i + 1 == FCanal[j][0]->getNodoInicialFuente()) {
						FLongitudVC[j][i] = FCanal[j][0]->getXRef() + FCanal[j][0]->getXRef() / 2. - FCanal[j][0]->getDistanciaInterpolacion();
					} else {
						FLongitudVC[j][i] = FCanal[j][0]->getXRef();
					}
				} else if(i == FCanal[j][0]->getNodoInicialFuente()) {
					// C�lculo de la longitud del Volumen de Control de cada nodo del canal de entrada
					FLongitudVC[j][i] = FCanal[j][0]->getXRef() / 2. + FCanal[j][0]->getDistanciaInterpolacion();
				} else if(i == FCanal[j][0]->getNin() - 1) {
					FLongitudVC[j][i] = FCanal[j][0]->getXRef() / 2.;
				} else {
					FLongitudVC[j][i] = FCanal[j][0]->getXRef();
				}
				FR1[j][i] = 0.;
				FR2[j][i] = 0.;
				FR3[j][i] = 0.;
				FR4[j][i] = 0.;
				FR5[j][i] = 0.;
				FKreg1[j][i] = 0.;
				FKreg2[j][i] = 0.;
				FQreg[j][i] = 0.;
				FQ1[j][i] = 0.;
				FQ2[j][i] = 0.;
				FFraccionMasicaNO2Entrada[j][i] = 0.;
				FFraccionMasicaNO2Salida[j][i] = 0.;
				for(int k = 0; k < FNumeroEspecies; k++) {
					FTasaFraccionMasicaEspecie[j][i][k] = 0.;
					FFraccionMasicaEspecieSalida[j][i][k] = 0.;
					FFraccionMasicaEspecieEntrantePared[j][i][k] = 0.;
				}
			}
		}

// Entalp�a de reacci�n (J/mol)
		FIncrH_reg1 = -2 * (FIndiceCompletitud1 - 0.5) * __HFormacion::CO2 - 2 * (1 - FIndiceCompletitud1) * __HFormacion::CO;
		FIncrH_reg2 = -(1 - FIndiceCompletitud2) * __HFormacion::CO2 + (FIndiceCompletitud2 - 0 - 5) * __HFormacion::CO -
					  FIndiceCompletitud2 * __HFormacion::NO;
		FIncrH_R1 = -__HFormacion::CO2;
		FIncrH_R2 = -4 * __HFormacion::CO2 - 6.93 / 2 * __HFormacion::H2O;

// C�lculo de datos para el submodelo de filtrado
		Ffuncion_f_Limpia = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidadLimpia,
									   1. / 3.) - FPorosidadLimpia - 0.2 * pow(1. - FPorosidadLimpia, 2.)) / (1. - FPorosidadLimpia);
		FDiametroUnidadColectoraLimpia = 1.5 * (1 - FPorosidadLimpia) / FPorosidadLimpia * FDiametroMedioPoroso;
		FDiametroUnidadCelular = pow(pow(FDiametroUnidadColectoraLimpia, 3.) / (1 - FPorosidadLimpia), 1. / 3.);
		Ffuncion_gPL_in = pow(FPorosidadCapaParticulas_in / (2. - FPorosidadCapaParticulas_in - 1.8 * pow(
								  1 - FPorosidadCapaParticulas_in, 1. / 3.) - 0.2 * pow(1 - FPorosidadCapaParticulas_in, 2)),
							  1. / 3.);
		Ffuncion_gPL_dep = pow(
							   FPorosidadCapaParticulas_dep / (2. - FPorosidadCapaParticulas_dep - 1.8 * pow(1 - FPorosidadCapaParticulas_dep,
									   1. / 3.) - 0.2 * pow(1 - FPorosidadCapaParticulas_dep, 2)), 1. / 3.);

// Se inicializa la velocidad a trav�s de la pared a 0 y el n�mero de unidades celulares.
		FNumeroUnidadesCelularesFiltro = 0;
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FNumeroUnidadesCelularesHaz[j] = 0;
			FBeamSootMass[j] = 0.;
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				FVelocidadPared[j][0][i] = 0.;
				FVelocidadPared[j][1][i] = 0.;
				FLayerSootMass[j][i] = 0.;
				FWallSootMass[j][i] = 0.;
				FCVSootMass[j][i] = 0.;
				VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
										  0) * 4;
				FNumeroUnidadesCelulares[j][i] = (long long int)(floor(VolumenControlPared / (__cons::Pi * pow(FDiametroUnidadCelular,
												 3.) / 6.)));
				if(i >= FCanal[j][0]->getNodoInicialFuente()) {
					UnidadesCelularesCanal = UnidadesCelularesCanal + FNumeroUnidadesCelulares[j][i];
				}
			}
			FNumeroUnidadesCelularesHaz[j] = UnidadesCelularesCanal * FCanal[j][0]->getNumeroCanales();
			FNumeroUnidadesCelularesFiltro = FNumeroUnidadesCelularesFiltro + FNumeroUnidadesCelularesHaz[j];
			UnidadesCelularesCanal = 0;
		}

// Calculo masa soot de saturacion en una unidad celular. Se supone DiametroUnidadColectora=FactorPercolaci�n*DiametroUnidadCelular
		FMasaSootSaturacion = (pow(FFactorPercolacion * FDiametroUnidadCelular / 2.,
								   3.) - pow(FDiametroUnidadColectoraLimpia / 2., 3.)) * (4. * __cons::Pi * FDensidadSootPared) / 3.;
		FMasaSootWallMaxUC = FMasaSootWallTotalMax / FNumeroUnidadesCelularesFiltro;

//Calculo de la masa de soot acumulable dentro de toda la pared porosa. Se suponen todos los canales iguales, tambien los de los extremos.
		FMasaSootSaturacionTotal = FMasaSootSaturacion * FNumeroUnidadesCelularesFiltro;

		if(FMasaSootSaturacion < FMasaSootWallMaxUC) {
			FMasaSootLimitUC = FMasaSootSaturacion;
		} else {
			FMasaSootLimitUC = FMasaSootWallMaxUC;
		}

// Inicializaci�n Permeabilidad de la pared limpia (Kuwabara)
		ViscosidadCinematica = FCanal[0][0]->GetViscosidadDinamica(0) / FCanal[0][0]->GetDensidad(0);
		Temperatura = pow(FCanal[0][0]->GetAsonido(0) * __cons::ARef,
						  2.) / FCanal[0][0]->GetGamma(0) / FCanal[0][0]->GetRMezcla(0);
		CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[0][0]->GetRMezcla(0) * Temperatura), 0.5);
		Knudsen = 2. * CaminoLibreMedio /
				  FDiametroMedioPoroso; // Referido al di�metro medio de la pared del filtro de part�culas
		if(Knudsen < 0.0016) {
			SCF = 1 + Knudsen * 1.257;
		} else {
			SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
		}
		FKwallLimpia = 0.09 * 2 * (2 - 1.8 * pow((1 - FPorosidadLimpia),
								   1. / 3.) - FPorosidadLimpia - 0.2 * pow((1 - FPorosidadLimpia), 2.)) / (9 * (1 - FPorosidadLimpia))
					   * pow(FDiametroUnidadColectoraLimpia, 2.) * SCF;
		double masatotal = 0.;
		if(FMasaSootInicial != 0) {
			FDPFSootMass = FMasaSootInicial;
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				//C�lculo del n�mero de unidades celulares y la masa de soot acumulable dentro de toda la pared porosa en el haz. Se suponen todos los canales iguales.
				FVolumenTotal[j] = FFraccionParedSaturada * FLongitudEfec * FEspesorParedPorosa * FCanal[j][0]->GetLadoCanal(
									   0) * FCanal[j][0]->getNumeroCanales() * 4.; // Volumen de pared porosa en el haz.
				FMasaSootSaturacionHaz[j] = FMasaSootLimitUC * FNumeroUnidadesCelularesHaz[j];
				if(!FCalculoDistribucionSoot) {
					FMasaSootInicialHaz[j] = FMasaSootInicial / ((FNumeroCanalesTotales / 2)) * FCanal[j][0]->getNumeroCanales();
					FBeamSootMass[j] = FMasaSootInicialHaz[j];
				}

				if(FCalculoDistribucionSoot) {
					FMasaSootInicialHaz[j] = FMasaSootInicial * FPorcentajeSootBeam[j] / 100.;
					FBeamSootMass[j] = FMasaSootInicialHaz[j];
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						// C�lculo del n�mero de unidades celulares en cada volumen de control de los nodos del canal de entrada.
						VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
												  0) * 4;
						FMasaSootSaturacionNodo = FMasaSootLimitUC * FNumeroUnidadesCelulares[j][i];
						if(i > 0) {
							masatotal = masatotal + FMasaSootSaturacionNodo;
						}
						// A continuaci�n se inicializa el modelo de filtrado a las condiciones de saturaci�n.
						if(FCalculoFiltrado) {
							FMasaSootNodo[j][i] = (FPorcentajeSootNodo[j][i] / 100.) * FMasaSootInicialHaz[j] / FCanal[j][0]->getNumeroCanales();
							if(FMasaSootNodo[j][i] < FMasaSootSaturacionNodo) {   // Node not saturated  or below maximum allowed loading
								FMasaSootUC[j][i] = FMasaSootNodo[j][i] / FNumeroUnidadesCelulares[j][i];
								FWallSootMass[j][i] = FMasaSootNodo[j][i];
								if(FWallSootMass[j][i] != 0) {
									SootFactor = FDensidadSootPared * VolumenControlPared / FWallSootMass[j][i];
								}
								FShapeFactor[j][i] = FaSF * pow(SootFactor, FbSF) + FMinShapeFactor;
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactor[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
								FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
															3.) * (1 - FPorosidadLimpia);
								FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);
								FLayerSootMass[j][i] = 0.;
								FWallSootMass[j][i] = FMasaSootNodo[j][i];
								FEspesorSoot[j][i] = 0.;
								FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
								FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
								FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];
							} else { // Node saturated (total or up to maximum allowed)
								FMasaSootUC[j][i] = FMasaSootLimitUC;
								FWallSootMass[j][i] = FMasaSootSaturacionNodo;
								SootFactor = FDensidadSootPared * VolumenControlPared / FWallSootMass[j][i];
								FShapeFactor[j][i] = FaSF * pow(SootFactor, FbSF) + FMinShapeFactor;
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactor[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
								FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
															3.) * (1 - FPorosidadLimpia);
								FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);
								FLayerSootMass[j][i] = FMasaSootNodo[j][i] - FMasaSootSaturacionNodo;
								FWallSootMass[j][i] = FMasaSootSaturacionNodo;
								FEspesorSootIn[j][i] = (FCanal[j][0]->GetLadoCanal(i)
														- pow(pow(FCanal[j][0]->GetLadoCanal(i), 2.) - FLayerSootMass[j][i] / (FLongitudVC[j][i] * FDensidadSootCapa_in),
															  0.5)) / 2.;
								FEspesorSoot[j][i] = FEspesorSootIn[j][i];
								FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
								FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
								FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];
							}

							funcion_f = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1 - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);

							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));

							// C�lculo de Eficiencia de Deposici�n debido a Difusi�n Browniana
							funcion_g = pow(FPorosidad[j][i] / (2. - FPorosidad[j][i] - 1.8 * pow(1 - FPorosidad[j][i],
																1. / 3.) - 0.2 * pow(1 - FPorosidad[j][i], 2.)), 1. / 3.);
							if(i == FCanal[j][0]->getNin() - 1) {
								VelocidadPared = FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()]
												 - (FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() - 1] - FVelocidadPared[j][1][i -
														 FCanal[j][0]->getNodoInicialFuente()]) / FCanal[j][1]->getXRef()
												 * FCanal[j][1]->getDistanciaInterpolacion();
							} else {
								VelocidadPared = Interpola(FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()],
														   FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() + 1], 1.,
														   FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
							}

							VelocidadIntersticial = VelocidadPared / FPorosidad[j][i];
							ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
							CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * __units::degCToK(
												   FCanal[j][0]->getTemperaturaInicial())), 0.5);
							Knudsen = 2 * CaminoLibreMedio / FDiametroPoro[j][i]; // como en el submodelo de filtrado
							if(Knudsen < 0.0016) {
								SCF = 1 + Knudsen * 1.257;
							} else {
								SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
							}
							CoefDifusionParticulas = __cons::Kb * __units::degCToK(FCanal[j][0]->getTemperaturaInicial()) * SCF
													 / (3. * __cons::Pi * FCanal[j][0]->GetViscosidadDinamica(i) * FDiametroAgregadoSoot);
							Peclet = VelocidadIntersticial * FDiametroUnidadColectora[j][i] /
									 CoefDifusionParticulas; // como en el submodelo de filtrado

							if(Peclet == 0) {
								EficienciaBrowniana = 0.;
							} else {
								EficienciaBrowniana = 3.5 * funcion_g * pow(Peclet, -2. / 3.);
							}
							if(EficienciaBrowniana > 1)
								EficienciaBrowniana = 1.;

							// Calculo de Eficiencia de Deposici�n debido a Intercepcion
							FParametroIntercepcion[j][i] = FDiametroAgregadoSoot / FDiametroUnidadColectora[j][i];
							EficienciaIntercepcion = 1.5 * pow(FParametroIntercepcion[j][i], 2.) * pow(funcion_g, 3.)
													 / pow(1 + FParametroIntercepcion[j][i], (3. - 2. * FPorosidad[j][i]) / (3. * FPorosidad[j][i]));

							// Calculo de Eficiencia de Deposicion Inertial
							Stokes = (2. / 9. * SCF * FDensidadSootPared * VelocidadIntersticial * pow(FDiametroAgregadoSoot,
									  2.)) / (2. * ViscosidadCinematica * FDiametroUnidadColectora[j][i]);
							EficienciaInercial = pow(Stokes, 2) / pow((Stokes + 0.25), 2);

							// Calculo de la Eficiencia de Deposicion Combinada
							EficienciaCombinada = (EficienciaBrowniana + EficienciaIntercepcion + EficienciaInercial)
												  - (EficienciaBrowniana * EficienciaIntercepcion + EficienciaBrowniana * EficienciaInercial + EficienciaIntercepcion *
													 EficienciaInercial)
												  + (EficienciaBrowniana * EficienciaIntercepcion * EficienciaInercial);

							// Calculo de la Eficiencia de Deposicion Total
							FEficiencia[j][i] = 1.
												- exp(
													-3 * EficienciaCombinada * (1. - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							// Evaluacion de las diferentes contribuciones a la eficiencia global
							FEficienciaBrown[j][i] = 1.
													 - exp(
														 -3 * EficienciaBrowniana * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaInter[j][i] = 1.
													 - exp(
														 -3 * EficienciaIntercepcion * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaIner[j][i] = 1.
													- exp(
														-3 * EficienciaInercial * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							// Eficiencia de la capa de partículas
							if(FCoeficienteParticion[j][i] >= FFiltRelaxCoef && FFiltRelaxCoef < 1 && FMasaSootUC[j][i] < FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i] * (FCoeficienteParticion[j][i] - FFiltRelaxCoef) / (1 - FFiltRelaxCoef);
							} else if(FMasaSootUC[j][i] >= FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i];
							}
						} else {
							FMasaSootNodo[j][i] = (FPorcentajeSootNodo[j][i] / 100) * FMasaSootInicialHaz[j] / FCanal[j][0]->getNumeroCanales();
							if(FMasaSootNodo[j][i] < FMasaSootSaturacionNodo) {     // Node not saturated
								FMasaSootUC[j][i] = FMasaSootNodo[j][i] / FNumeroUnidadesCelulares[j][i];
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactorDiscrete / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
								FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
															3.) * (1 - FPorosidadLimpia);
								FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);
								FLayerSootMass[j][i] = 0.;
								FWallSootMass[j][i] = FMasaSootNodo[j][i];
								FEspesorSoot[j][i] = 0.;
								FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
								FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
								FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];
							} else { // Node saturated (total or up to maximum allowed)
								FMasaSootUC[j][i] = FMasaSootLimitUC;
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactorDiscrete / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
								FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
															3.) * (1 - FPorosidadLimpia);
								FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);
								FLayerSootMass[j][i] = FMasaSootNodo[j][i] - FMasaSootSaturacionNodo;
								FWallSootMass[j][i] = FMasaSootSaturacionNodo;
								FEspesorSootIn[j][i] = (FCanal[j][0]->GetLadoCanal(i)
														- pow(pow(FCanal[j][0]->GetLadoCanal(i), 2.) - FLayerSootMass[j][i] / (FLongitudVC[j][i] * FDensidadSootCapa_in),
															  0.5)) / 2.;
								FEspesorSoot[j][i] = FEspesorSootIn[j][i] + FEspesorSootDep[j][i];
								FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
								FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
								FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];
							}
							funcion_f = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1 - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);

							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));
							FParametroIntercepcion[j][i] = 0.;
							FEficiencia[j][i] = 0.;
							FEficienciaPL[j][i] = 0.;
							FEficienciaBrown[j][i] = 0.;
							FEficienciaInter[j][i] = 0.;
							FEficienciaIner[j][i] = 0.;
						}
					}
				} else if(FMasaSootSaturacionHaz[j] < FMasaSootInicialHaz[j]) {
					// Ya hay formada una capa de soot sobre las paredes porosas, que se encuentran saturadas. Soot homogeneamente distribuido.
					// Se calcula el n�mero de unidades celulares asociadas a cada uno de los volumenes de control
					// existentes en el canal de entrada donde existe velocidad a trav�s de la pared. (Los nodos extremos son diferentes)
					FMasaSootLayerTotal[j] = FMasaSootInicialHaz[j] - FMasaSootSaturacionHaz[j];
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						// C�lculo del n�mero de unidades celulares en cada volumen de control de los nodos del canal de entrada.
						VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
												  0) * 4;

						// A continuaci�n se inicializa el modelo de filtrado a las condiciones de saturaci�n.
						if(FCalculoFiltrado) {
							FMasaSootUC[j][i] = FMasaSootLimitUC;
							if(FCalculoFiltrado) {
								FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
								SootFactor = FDensidadSootPared * VolumenControlPared / FWallSootMass[j][i];
								FShapeFactor[j][i] = FaSF * pow(SootFactor, FbSF) + FMinShapeFactor;
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactor[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							} else {
								FDiametroUnidadColectora[j][i] = 2.
																 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactorDiscrete / (__cons::Pi * FDensidadSootPared) + pow(
																		 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							}
							FPorosidad[j][i] = 1 - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													   3.) * (1 - FPorosidadLimpia);
							FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);

							funcion_f = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1 - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);

							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));

							// C�lculo de Eficiencia de Deposici�n debido a Difusi�n Browniana
							funcion_g = pow(FPorosidad[j][i] / (2. - FPorosidad[j][i] - 1.8 * pow(1 - FPorosidad[j][i],
																1. / 3.) - 0.2 * pow(1 - FPorosidad[j][i], 2.)), 1. / 3.);
							if(i == FCanal[j][0]->getNin() - 1) {
								VelocidadPared = FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()]
												 - (FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() - 1] - FVelocidadPared[j][1][i -
														 FCanal[j][0]->getNodoInicialFuente()]) / FCanal[j][1]->getXRef()
												 * FCanal[j][1]->getDistanciaInterpolacion();
							} else {
								VelocidadPared = Interpola(FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()],
														   FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() + 1], 1.,
														   FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
							}

							VelocidadIntersticial = VelocidadPared / FPorosidad[j][i];
							ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
							CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * __units::degCToK(
												   FCanal[j][0]->getTemperaturaInicial())), 0.5);
							Knudsen = 2 * CaminoLibreMedio / FDiametroPoro[j][i]; // como en el submodelo de filtrado
							if(Knudsen < 0.0016) {
								SCF = 1 + Knudsen * 1.257;
							} else {
								SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
							}
							CoefDifusionParticulas = __cons::Kb * __units::degCToK(FCanal[j][0]->getTemperaturaInicial()) * SCF
													 / (3. * __cons::Pi * FCanal[j][0]->GetViscosidadDinamica(i) * FDiametroAgregadoSoot);
							Peclet = VelocidadIntersticial * FDiametroUnidadColectora[j][i] /
									 CoefDifusionParticulas; // como en el submodelo de filtrado

							if(Peclet == 0) {
								EficienciaBrowniana = 0.;
							} else {
								EficienciaBrowniana = 3.5 * funcion_g * pow(Peclet, -2. / 3.);
							}
							if(EficienciaBrowniana > 1)
								EficienciaBrowniana = 1.;

							// C�lculo de Eficiencia de Deposici�n debido a Intercepci�n
							FParametroIntercepcion[j][i] = FDiametroAgregadoSoot / FDiametroUnidadColectora[j][i];
							EficienciaIntercepcion = 1.5 * pow(FParametroIntercepcion[j][i], 2.) * pow(funcion_g, 3.)
													 / pow(1 + FParametroIntercepcion[j][i], (3. - 2. * FPorosidad[j][i]) / (3. * FPorosidad[j][i]));

							// Cï¿½lculo de Eficiencia de Deposicion Inertial
							// Cs=1+2/FDiametroAgregadoSoot*(1.23+0.41*exp(-0.88*(FDiametroAgregadoSoot/(2.*CaminoLibreMedio))));
							Stokes = (2. / 9. * SCF * FDensidadSootPared * VelocidadIntersticial * pow(FDiametroAgregadoSoot,
									  2.)) / (2. * ViscosidadCinematica * FDiametroUnidadColectora[j][i]);
							EficienciaInercial = pow(Stokes, 2) / pow((Stokes + 0.25), 2);

							// Cï¿½lculo de la Eficiencia  de Deposiciï¿½n Combinada
							EficienciaCombinada = (EficienciaBrowniana + EficienciaIntercepcion + EficienciaInercial)
												  - (EficienciaBrowniana * EficienciaIntercepcion + EficienciaBrowniana * EficienciaInercial + EficienciaIntercepcion *
													 EficienciaInercial)
												  + (EficienciaBrowniana * EficienciaIntercepcion * EficienciaInercial);

							// C�lculo de la Eficiencia de Deposici�n Total
							FEficiencia[j][i] = 1.
												- exp(
													-3 * EficienciaCombinada * (1. - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							// Evaluacion de las diferentes contribuciones a la eficiencia global
							FEficienciaBrown[j][i] = 1.
													 - exp(
														 -3 * EficienciaBrowniana * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaInter[j][i] = 1.
													 - exp(
														 -3 * EficienciaIntercepcion * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaIner[j][i] = 1.
													- exp(
														-3 * EficienciaInercial * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							FEspesorSootIn[j][i] = (FCanal[j][0]->GetLadoCanal(i)
													- pow(pow(FCanal[j][0]->GetLadoCanal(i),
															  2.) - FMasaSootLayerTotal[j] / (FCanal[j][0]->getNumeroCanales() * FLongitudEfec * FDensidadSootCapa_in), 0.5)) / 2.;
							FEspesorSoot[j][i] = FEspesorSootIn[j][i] + FEspesorSootDep[j][i];
							// C�lculo del espesor de soot inicial, por lo que se considera uniforme.
							FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
							FLayerSootMass[j][i] = FMasaSootLayerTotal[j] / FCanal[j][0]->getNumeroCanales() / (FCanal[j][0]->getNin() - 1);
							FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
							FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
							FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];

							// Eficiencia de la capa de partículas
							if(FCoeficienteParticion[j][i] >= FFiltRelaxCoef && FFiltRelaxCoef < 1 && FMasaSootUC[j][i] < FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i] * (FCoeficienteParticion[j][i] - FFiltRelaxCoef) / (1 - FFiltRelaxCoef);
							} else if(FMasaSootUC[j][i] >= FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i];
							}
						} else {
							FMasaSootUC[j][i] = FMasaSootLimitUC;
							FDiametroUnidadColectora[j][i] = 2.
															 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactorDiscrete / (__cons::Pi * FDensidadSootPared) + pow(
																	 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FPorosidad[j][i] = 1 - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													   3.) * (1 - FPorosidadLimpia);
							FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);

							funcion_f = 2. / 9. * (2. - 1.8 * pow(1 - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1 - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);

							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));

							// C�lculo del espesor de soot inicial, por lo que se considera uniforme.
							FEspesorSootIn[j][i] = (FCanal[j][0]->GetLadoCanal(i)
													- pow(pow(FCanal[j][0]->GetLadoCanal(i),
															  2.) - FMasaSootLayerTotal[j] / (FCanal[j][0]->getNumeroCanales() * FLongitudEfec * FDensidadSootCapa_in), 0.5)) / 2.;
							FEspesorSoot[j][i] = FEspesorSootIn[j][i] + FEspesorSootDep[j][i];
							FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
							FLayerSootMass[j][i] = FMasaSootLayerTotal[j] / FCanal[j][0]->getNumeroCanales() / (FCanal[j][0]->getNin() - 1);
							FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
							FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
							FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];

							FParametroIntercepcion[j][i] = 0.;
							FEficiencia[j][i] = 0.;
							FEficienciaPL[j][i] = 0.;
							FEficienciaBrown[j][i] = 0.;
							FEficienciaInter[j][i] = 0.;
							FEficienciaIner[j][i] = 0.;
						}
					}
				} else { // Solo existe soot en el interior de la capa porosa homogeneamente distribuido.
					// Se calcula el n�mero de unidades celulares asociadas a cada uno de los volumenes de control
					// existentes en el canal de entrada donde existe velocidad a trav�s de la pared. (Los nodos extremos son diferentes)
					FMasaSootLayerTotal[j] = 0.;
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						// C�lculo del n�mero de unidades celulares en cada volumen de control de los nodos del canal de entrada.
						VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
												  0) * 4;
						// A continuaci�n se inicializa el modelo de filtrado a las condiciones de ensuciamiento iniciales.
						if(FCalculoFiltrado) {
							FMasaSootUC[j][i] = FMasaSootInicialHaz[j] / FNumeroUnidadesCelularesHaz[j];

							FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
							SootFactor = FDensidadSootPared * VolumenControlPared / FWallSootMass[j][i];
							FShapeFactor[j][i] = FaSF * pow(SootFactor, FbSF) + FMinShapeFactor;
							FDiametroUnidadColectora[j][i] = 2.
															 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactor[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
																	 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);

							FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
														3.) * (1 - FPorosidadLimpia);
							FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);

							funcion_f = 2. / 9. * (2. - 1.8 * pow(1. - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1. - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);
							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));

							// Calculo de Eficiencia de Deposici�n debido a Difusi�n Browniana
							funcion_g = pow(FPorosidad[j][i] / (2. - FPorosidad[j][i] - 1.8 * pow(1 - FPorosidad[j][i],
																1. / 3.) - 0.2 * pow(1 - FPorosidad[j][i], 2.)), 1. / 3.);
							if(i == FCanal[j][0]->getNin() - 1) {
								VelocidadPared = FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()]
												 - (FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() - 1] - FVelocidadPared[j][1][i -
														 FCanal[j][0]->getNodoInicialFuente()]) / FCanal[j][1]->getXRef()
												 * FCanal[j][1]->getDistanciaInterpolacion();
							} else {
								VelocidadPared = Interpola(FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()],
														   FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() + 1], 1.,
														   FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
							}
							VelocidadIntersticial = VelocidadPared / FPorosidad[j][i];
							ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
							CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * __units::degCToK(
												   FCanal[j][0]->getTemperaturaInicial())), 0.5);
							Knudsen = 2 * CaminoLibreMedio / FDiametroPoro[j][i]; // como en el submodelo de filtrado

							if(Knudsen < 0.0016) {
								SCF = 1 + Knudsen * 1.257;
							} else {
								SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
							}
							CoefDifusionParticulas = __cons::Kb * (FCanal[j][0]->getTemperaturaInicial() + 273.) * SCF
													 / (3. * __cons::Pi * FCanal[j][0]->GetViscosidadDinamica(i) * FDiametroAgregadoSoot);
							Peclet = VelocidadIntersticial * FDiametroUnidadColectora[j][i] /
									 CoefDifusionParticulas; // como en el submodelo de filtrado

							if(Peclet == 0) {
								EficienciaBrowniana = 0.;
							} else {
								EficienciaBrowniana = 3.5 * funcion_g * pow(Peclet, -2. / 3.);
							}
							if(EficienciaBrowniana > 1)
								EficienciaBrowniana = 1.;

							// C�lculo de Eficiencia de Deposici�n debido a Intercepci�n
							FParametroIntercepcion[j][i] = FDiametroAgregadoSoot / FDiametroUnidadColectora[j][i];
							EficienciaIntercepcion = 1.5 * pow(FParametroIntercepcion[j][i], 2.) * pow(funcion_g, 3.)
													 / pow(1 + FParametroIntercepcion[j][i], (3 - 2 * FPorosidad[j][i]) / (3 * FPorosidad[j][i]));

							// Calculo de Eficiencia de Deposicion Inertial
							Stokes = (2. / 9. * SCF * FDensidadSootPared * VelocidadIntersticial * pow(FDiametroAgregadoSoot,
									  2.)) / (2. * ViscosidadCinematica * FDiametroUnidadColectora[j][i]);
							EficienciaInercial = pow(Stokes, 2) / pow((Stokes + 0.25), 2);

							// C�lculo de la Eficiencia  de Deposici�n Combinada
							EficienciaCombinada = (EficienciaBrowniana + EficienciaIntercepcion + EficienciaInercial)
												  - (EficienciaBrowniana * EficienciaIntercepcion + EficienciaBrowniana * EficienciaInercial + EficienciaIntercepcion *
													 EficienciaInercial)
												  + (EficienciaBrowniana * EficienciaIntercepcion * EficienciaInercial);

							// C�lculo de la Eficiencia de Deposici�n Total
							FEficiencia[j][i] = 1.
												- exp(
													-3 * EficienciaCombinada * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							// Evaluacion de las diferentes contribuciones a la eficiencia global
							FEficienciaBrown[j][i] = 1.
													 - exp(
														 -3 * EficienciaBrowniana * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaInter[j][i] = 1.
													 - exp(
														 -3 * EficienciaIntercepcion * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
							FEficienciaIner[j][i] = 1.
													- exp(
														-3 * EficienciaInercial * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
														/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

							// Calculo del espesor de soot inicial.
							FEspesorSoot[j][i] = 0.;
							FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
							FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
							FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
							FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];

							// Eficiencia de la capa de partículas
							if(FCoeficienteParticion[j][i] >= FFiltRelaxCoef && FFiltRelaxCoef < 1 && FMasaSootUC[j][i] < FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i] * (FCoeficienteParticion[j][i] - FFiltRelaxCoef) / (1 - FFiltRelaxCoef);
							} else if(FMasaSootUC[j][i] >= FMasaSootLimitUC) {
								FEficienciaPL[j][i] = FEficiencia[j][i];
							}
						} else {
							FMasaSootUC[j][i] = FMasaSootInicialHaz[j] / FNumeroUnidadesCelularesHaz[j];
							FDiametroUnidadColectora[j][i] = 2.
															 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactorDiscrete / (__cons::Pi * FDensidadSootPared) + pow(
																	 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
														3.) * (1 - FPorosidadLimpia);
							FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);

							funcion_f = 2. / 9. * (2. - 1.8 * pow(1. - FPorosidad[j][i],
																  1. / 3.) - FPorosidad[j][i] - 0.2 * pow(1. - FPorosidad[j][i], 2.)) / (1. - FPorosidad[j][i]);
							FKwallClean[j][i] = FKwallLimpia;
							FKwallLoaded[j][i] = pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													 2.) * funcion_f / Ffuncion_f_Limpia * FKwallLimpia;
							FKwall[j][i] = FKwallClean[j][i] * FKwallLoaded[j][i] / (FKwallClean[j][i] * FFraccionParedSaturada +
										   (1 - FFraccionParedSaturada) * FKwallLoaded[j][i]);
							double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
									FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
							FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
														  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));

							// Calculo del espesor de soot inicial.
							FEspesorSoot[j][i] = 0.;
							FWallSootMass[j][i] = FMasaSootUC[j][i] * FNumeroUnidadesCelulares[j][i];
							FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
							FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
							FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];

							FParametroIntercepcion[j][i] = 0.;
							FEficiencia[j][i] = 0.;
							FEficienciaPL[j][i] = 0.;
							FEficienciaBrown[j][i] = 0.;
							FEficienciaInter[j][i] = 0.;
							FEficienciaIner[j][i] = 0.;
						}
					}
				}
			}
		} else {
			// No hay que introducir el submodelo de distribucion de soot porque no hay soot inicialmente y el calculo es el mismo
			// No hay soot inicial, luego la trampa est� inicialmente completamente limpia.
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FVolumenTotal[j] = FFraccionParedSaturada * FLongitudEfec * FEspesorParedPorosa * FCanal[j][0]->GetLadoCanal(
									   0) * FCanal[j][0]->getNumeroCanales() * 4.; // Volumen de pared porosa en el haz.
				FMasaSootSaturacionHaz[j] = FMasaSootLimitUC * FNumeroUnidadesCelularesHaz[j];   // Para al 25% (0.1)

				FMasaSootInicialHaz[j] = 0.;
				FMasaSootLayerTotal[j] = 0.;

				for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
					// C�lculo del n�mero de unidades celulares en cada volumen de control de los nodos del canal de entrada.
					VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
											  0) * 4;

					// A continuaci�n se inicializa el modelo de filtrado a las condiciones de trampa limpia.
					/*if(FCanal[j][0]->getNodoInicialFuente()>i){
					 FEficiencia[j][i]=0.;
					 FEficienciaPL[j][i]=0.;
					 FEficienciaBrown[j][i]=0.;
					 FEficienciaInter[j][i]=0.;
					 FEficienciaIner[j][i]=0.;

					 }else{*/
					if(FCalculoFiltrado) {
						FMasaSootUC[j][i] = 0.;
						FDiametroUnidadColectora[j][i] = FDiametroUnidadColectoraLimpia;
						FPorosidad[j][i] = FPorosidadLimpia;
						FKwallClean[j][i] = FKwallLimpia;
						FKwallLoaded[j][i] = FKwallLimpia;
						FKwall[j][i] = FKwallLimpia;

						FCoeficienteParticion[j][i] = 0.;
						FDiametroPoro[j][i] = FDiametroMedioPoroso;

						// C�lculo de Eficiencia de Deposici�n debido a Difusi�n Browniana
						funcion_g = pow(FPorosidad[j][i] / (2. - FPorosidad[j][i] - 1.8 * pow(1 - FPorosidad[j][i],
															1. / 3.) - 0.2 * pow(1 - FPorosidad[j][i], 2.)), 1. / 3.);
						if(i == FCanal[j][0]->getNin() - 1) {
							VelocidadPared = FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()]
											 - (FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() - 1] - FVelocidadPared[j][1][i -
													 FCanal[j][0]->getNodoInicialFuente()]) / FCanal[j][1]->getXRef()
											 * FCanal[j][1]->getDistanciaInterpolacion();
						} else {
							VelocidadPared = Interpola(FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente()],
													   FVelocidadPared[j][1][i - FCanal[j][0]->getNodoInicialFuente() + 1], 1.,
													   FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
						}
						VelocidadIntersticial = FVelocidadPared[j][0][i] / FPorosidad[j][i];
						ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
						CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) *
										   (FCanal[j][0]->getTemperaturaInicial() + 273.)), 0.5);
						Knudsen = 2 * CaminoLibreMedio / FDiametroPoro[j][i]; // como en el submodelo de filtrado
						if(Knudsen < 0.0016) {
							SCF = 1 + Knudsen * 1.257;
						} else {
							SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
						}
						CoefDifusionParticulas = __cons::Kb * __units::degCToK(FCanal[j][0]->getTemperaturaInicial()) * SCF
												 / (3. * __cons::Pi * FCanal[j][0]->GetViscosidadDinamica(i) * FDiametroAgregadoSoot);
						Peclet = VelocidadIntersticial * FDiametroUnidadColectora[j][i] /
								 CoefDifusionParticulas; // como en el submodelo de filtrado
						if(Peclet == 0) {
							EficienciaBrowniana = 0.;
						} else {
							EficienciaBrowniana = 3.5 * funcion_g * pow(Peclet, -2. / 3.);
						}
						if(EficienciaBrowniana > 1)
							EficienciaBrowniana = 1.;

						// C�lculo de Eficiencia de Deposici�n debido a Intercepci�n
						FParametroIntercepcion[j][i] = FDiametroAgregadoSoot / FDiametroUnidadColectora[j][i];
						EficienciaIntercepcion = 1.5 * pow(FParametroIntercepcion[j][i], 2.) * pow(funcion_g, 3.)
												 / pow(1. + FParametroIntercepcion[j][i], (3. - 2. * FPorosidad[j][i]) / (3. * FPorosidad[j][i]));

						// Cï¿½lculo de Eficiencia de Deposicion Inertial
						// Cs=1+2/FDiametroAgregadoSoot*(1.23+0.41*exp(-0.88*(FDiametroAgregadoSoot/(2.*CaminoLibreMedio))));
						Stokes = (2. / 9. * SCF * FDensidadSootPared * VelocidadIntersticial * pow(FDiametroAgregadoSoot,
								  2.)) / (2. * ViscosidadCinematica * FDiametroUnidadColectora[j][i]);
						EficienciaInercial = pow(Stokes, 2) / pow((Stokes + 0.25), 2);

						// Cï¿½lculo de la Eficiencia  de Deposiciï¿½n Combinada
						EficienciaCombinada = (EficienciaBrowniana + EficienciaIntercepcion + EficienciaInercial)
											  - (EficienciaBrowniana * EficienciaIntercepcion + EficienciaBrowniana * EficienciaInercial + EficienciaIntercepcion *
												 EficienciaInercial)
											  + (EficienciaBrowniana * EficienciaIntercepcion * EficienciaInercial);

						// C�lculo de la Eficiencia de Deposici�n Total
						FEficiencia[j][i] = 1
											- exp(
												-3 * EficienciaCombinada * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
												/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

						// Evaluacion de las diferentes contribuciones a la eficiencia global
						FEficienciaBrown[j][i] = 1.
												 - exp(
													 -3 * EficienciaBrowniana * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
						FEficienciaInter[j][i] = 1.
												 - exp(
													 -3 * EficienciaIntercepcion * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
						FEficienciaIner[j][i] = 1.
												- exp(
													-3 * EficienciaInercial * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

						//Particulate layer efficiency
						FEficienciaPL[j][i] = 0.;

						// C�lculo del espesor de soot inicial.
						FEspesorSoot[j][i] = 0.;
						FLayerSootMass[j][i] = 0.;
						FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
						FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];
					} else {
						FMasaSootUC[j][i] = 0.;
						FDiametroUnidadColectora[j][i] = FDiametroUnidadColectoraLimpia;
						FPorosidad[j][i] = FPorosidadLimpia;
						FKwallClean[j][i] = FKwallLimpia;
						FKwallLoaded[j][i] = FKwallLimpia;
						FKwall[j][i] = FKwallLimpia;
						FCoeficienteParticion[j][i] = 0.;
						FDiametroPoro[j][i] = FDiametroMedioPoroso;

						// C�lculo del espesor de soot inicial.
						FEspesorSoot[j][i] = 0.;
						FLayerSootMass[j][i] = 0.;
						FAreaVCCanalEntrada[j][i] = 4. * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];
						FAreaVCCanalSalida[j][i] = 4. * FCanal[j][0]->GetLadoCanal(i) * FLongitudVC[j][i];

						FParametroIntercepcion[j][i] = 0.;
						FEficiencia[j][i] = 0.;
						FEficienciaPL[j][i] = 0.;
						FEficienciaBrown[j][i] = 0.;
						FEficienciaInter[j][i] = 0.;
						FEficienciaIner[j][i] = 0.;
					}
				}
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::InicializaDPF en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::SubmodeloFiltrado() {
	try {
		double funcion_g, funcion_gPL_in, VelocidadPared, VelocidadIntersticial, ViscosidadCinematica, Temperatura,
			   CaminoLibreMedio, Knudsen, SCF, SCF_in, CoefDifusionParticulas, Peclet,
			   EficienciaBrowniana, EficienciaIntercepcion, Cs, Stokes, EficienciaInercial, EficienciaCombinada, SootFactor,
			   VolumenControlPared;

// double VelocidadIntersticialPL=0., CoefDifusionParticulasPL=0., PecletPL=0., EficienciaBrownianaPL=0.,
//	   EficienciaIntercepcionPL=0., StokesPL=0., EficienciaInercialPL=0., EficienciaCombinadaPL=0.,

		double masasootatrapadaVC = 0., FLayerSootMassDEP = 0.;
		double densidad = 0.;

		FDPFSootMass = 0;

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FBeamSootMass[j] = 0.;
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				//FLayerSootMassDep[j][i]=0.;
				if(i >= FCanal[j][0]->getNodoInicialFuente()) {
					VolumenControlPared = FFraccionParedSaturada * FEspesorParedPorosa * FLongitudVC[j][i] * FCanal[j][0]->GetLadoCanal(
											  0) * 4;
					Temperatura = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
									  2.) / (FCanal[j][0]->GetGamma(i) * FCanal[j][0]->GetRMezcla(i));
					densidad = __units::BarToPa(FCanal[j][0]->GetPresion(i)) / (FCanal[j][0]->GetRMezcla(i) * Temperatura);

					if((FCoeficienteParticion[j][i] < 1 && FMasaSootUC[j][i] < FMasaSootLimitUC) || FCanal[j][0]->getTime1() < 0.001) {
						// Filtrado en la pared porosa
						// Calculo de Eficiencia de Deposicion debido a Difusion Browniana
						if(FVelocidadPared[j][0][i] < 0) {
							VelocidadIntersticial = 0.;
						} else {
							VelocidadIntersticial = FVelocidadPared[j][0][i] / FPorosidad[j][i];
						}

						if(FWallSootMass[j][i] <= 0.) {
							FDiametroUnidadColectora[j][i] = FDiametroUnidadColectoraLimpia;
							FShapeFactor[j][i] = 1;
						} else {
							SootFactor = FDensidadSootPared * VolumenControlPared / FWallSootMass[j][i];
							FShapeFactor[j][i] = FaSF * pow(SootFactor, FbSF) + FMinShapeFactor;
							FDiametroUnidadColectora[j][i] = 2.
															 * pow(0.75 * FMasaSootUC[j][i] / FShapeFactor[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
																	 FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
						}
						FPorosidad[j][i] = 1. - pow(FDiametroUnidadColectora[j][i] / FDiametroUnidadColectoraLimpia,
													3.) * (1 - FPorosidadLimpia);
						funcion_g = pow(FPorosidad[j][i] / (2. - FPorosidad[j][i] - 1.8 * pow(1 - FPorosidad[j][i],
															1. / 3.) - 0.2 * pow(1 - FPorosidad[j][i], 2)), 1. / 3.);
						double DiametroUnidadColectoraVirtual = 2. * pow(0.75 * FMasaSootUC[j][i] / (__cons::Pi * FDensidadSootPared) + pow(
								FDiametroUnidadColectoraLimpia / 2., 3.), 1. / 3.);
						FCoeficienteParticion[j][i] = (pow(DiametroUnidadColectoraVirtual, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.))
													  / (pow(FFactorPercolacion * FDiametroUnidadCelular, 3.) - pow(FDiametroUnidadColectoraLimpia, 3.));
						FDiametroPoro[j][i] = 2. / 3. * FDiametroUnidadColectora[j][i] * FPorosidad[j][i] / (1 - FPorosidad[j][i]);

						ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
						CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * Temperatura), 0.5);
						Knudsen = 2 * CaminoLibreMedio /
								  FDiametroPoro[j][i]; // Referido al diametro medio de la pared del filtro de part�culas
						if(Knudsen < 0.0016) {
							SCF = 1 + Knudsen * 1.257;
						} else {
							SCF = 1 + Knudsen * (1.257 + 0.4 * exp(-1.1 / Knudsen));
						}
						CoefDifusionParticulas = __cons::Kb * Temperatura * SCF / (3. * __cons::Pi * FCanal[j][0]->GetViscosidadDinamica(
													 i) * FDiametroAgregadoSoot);
						Peclet = VelocidadIntersticial * FDiametroUnidadColectora[j][i] / CoefDifusionParticulas;

						if(Peclet == 0 || Peclet < 0) {
							EficienciaBrowniana = 0.;
						} else {
							EficienciaBrowniana = 3.5 * funcion_g * pow(Peclet, -2. / 3.);
						}
						if(EficienciaBrowniana > 1)
							EficienciaBrowniana = 1.;

						// Calculo de Eficiencia de Deposicion debido a Intercepcion
						FParametroIntercepcion[j][i] = FDiametroAgregadoSoot / FDiametroUnidadColectora[j][i];
						EficienciaIntercepcion = 1.5 * pow(FParametroIntercepcion[j][i], 2.) * pow(funcion_g, 3.)
												 / pow(1. + FParametroIntercepcion[j][i], (3. - 2. * FPorosidad[j][i]) / (3. * FPorosidad[j][i]));

						// Calculo de Eficiencia de Deposicion Inertial
						// Cs=1+2/FDiametroAgregadoSoot*(1.23+0.41*exp(-0.88*(FDiametroAgregadoSoot/(2.*CaminoLibreMedio))));
						Stokes = (2. / 9. * SCF * FDensidadSootPared * VelocidadIntersticial * pow(FDiametroAgregadoSoot,
								  2.)) / (2. * ViscosidadCinematica * FDiametroUnidadColectora[j][i]);
						EficienciaInercial = pow(Stokes, 2) / pow((Stokes + 0.25), 2);

						// Calculo de la Eficiencia  de Deposicion Combinada
						EficienciaCombinada = (EficienciaBrowniana + EficienciaIntercepcion + EficienciaInercial)
											  - (EficienciaBrowniana * EficienciaIntercepcion + EficienciaBrowniana * EficienciaInercial + EficienciaIntercepcion *
												 EficienciaInercial)
											  + (EficienciaBrowniana * EficienciaIntercepcion * EficienciaInercial);

						// Calculo de la Eficiencia de Deposicion Total
						FEficiencia[j][i] = 1.
											- exp(
												-3 * EficienciaCombinada * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
												/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

						// Evaluacion de las diferentes contribuciones a la eficiencia global
						FEficienciaBrown[j][i] = 1.
												 - exp(
													 -3 * EficienciaBrowniana * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
						FEficienciaInter[j][i] = 1.
												 - exp(
													 -3 * EficienciaIntercepcion * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													 / (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));
						FEficienciaIner[j][i] = 1.
												- exp(
													-3 * EficienciaInercial * (1 - FPorosidad[j][i]) * FFraccionParedSaturada * FEspesorParedPorosa * FStickingCoef
													/ (2. * FPorosidad[j][i] * FDiametroUnidadColectora[j][i]));

					}
					// Influencia de la capa de particulas

					if(FCoeficienteParticion[j][i] >= FFiltRelaxCoef && FFiltRelaxCoef < 1 && FMasaSootUC[j][i] < FMasaSootLimitUC) {
						FEficienciaPL[j][i] = FEficiencia[j][i] * (FCoeficienteParticion[j][i] - FFiltRelaxCoef) / (1 - FFiltRelaxCoef);
					} else if(FMasaSootUC[j][i] >= FMasaSootLimitUC) {
						FEficienciaPL[j][i] = FEficiencia[j][i];
					}

					if(FCoeficienteParticion[j][i] < 1 && FCoeficienteParticion[j][i] <= FFiltRelaxCoef
					   && FMasaSootUC[j][i] < FMasaSootLimitUC) {
						// Balance de masas y especies quimicas
						// Pared porosa no saturada. Filtrado solo dentro de la pared

						masasootatrapadaVC = FEficiencia[j][i] * FCanal[j][0]->GetYespecie(i,
											 4) * FVelocidadPared[j][0][i] * densidad * FAreaVCCanalEntrada[j][i] * FDeltaTimeDPF;

						FMasaSootUC[j][i] = FMasaSootUC[j][i] + masasootatrapadaVC / FNumeroUnidadesCelulares[j][i];
						FWallSootMass[j][i] = FWallSootMass[j][i] + masasootatrapadaVC;
						FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
						FBeamSootMass[j] = FBeamSootMass[j] + FCVSootMass[j][i] * FCanal[j][0]->getNumeroCanales();

						FAreaVCCanalEntrada[j][i] = 4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];

						FFraccionMasicaEspecieSalida[j][i][0] = FCanal[j][0]->GetYespecie(i, 0)
																+ FCanal[j][0]->GetYespecie(i, 0) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][1] = FCanal[j][0]->GetYespecie(i, 1)
																+ FCanal[j][0]->GetYespecie(i, 1) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][2] = FCanal[j][0]->GetYespecie(i, 2)
																+ FCanal[j][0]->GetYespecie(i, 2) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][3] = FCanal[j][0]->GetYespecie(i, 3)
																+ FCanal[j][0]->GetYespecie(i, 3) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][4] = FCanal[j][0]->GetYespecie(i, 4) * (1 - FEficiencia[j][i]); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FCanal[j][0]->GetYespecie(i, 5)
																+ FCanal[j][0]->GetYespecie(i, 5) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][6] = FCanal[j][0]->GetYespecie(i, 6)
																+ FCanal[j][0]->GetYespecie(i, 6) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
																	+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
																	+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
																	+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
																	+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][9] = FCanal[j][0]->GetYespecie(i, 9)
																	+ FCanal[j][0]->GetYespecie(i, 9) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficiencia[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						}
						FMasaSootLimitInfUC[j][i] = FMasaSootUC[j][i];

					} else if(FCoeficienteParticion[j][i] < 1 && FCoeficienteParticion[j][i] > FFiltRelaxCoef
							  && FMasaSootUC[j][i] < FMasaSootLimitUC) {
						// Balance de masas y especies quimicas
						// Pared porosa no saturada que sigue cargándose, pero el filtrado produce también el aumento de la capa de part�culas cuyas propiedades se dan al inicio de la ejecuci�n
						// Hay que calcular la eficiencia de la capa

						/*	if(FVelocidadPared[j][0][i]<0){
						 VelocidadIntersticialPL=0.;
						 }else{
						 VelocidadIntersticialPL=FVelocidadPared[j][0][i]/FPorosidadCapaParticulas_in;
						 }
						 CoefDifusionParticulasPL=__cons::Kb*Temperatura*SCF_in/(3.*Pi*FCanal[j][0]->GetViscosidadDinamica(i)*FDiametroAgregadoSoot);
						 PecletPL=VelocidadIntersticialPL*FDiametroAgregadoSoot/CoefDifusionParticulasPL;

						 if(PecletPL==0 || PecletPL<0){
						 EficienciaBrownianaPL=0.;
						 }else{
						 EficienciaBrownianaPL=3.5*funcion_gPL_in*pow(PecletPL,-2./3.);
						 }
						 if(EficienciaBrownianaPL>1) EficienciaBrownianaPL=1.;

						 // Calculo de Eficiencia de Deposicion debido a Intercepcion
						 FParametroIntercepcionPL=FDiametroAgregadoSoot/FDiametroAgregadoSoot;
						 EficienciaIntercepcionPL=1.5*pow(FParametroIntercepcionPL,2.)*pow(funcion_gPL_in,3.)/pow(1.+FParametroIntercepcionPL,(3.-2.*FPorosidadCapaParticulas_in)/(3.*FPorosidadCapaParticulas_in));

						 // Calculo de Eficiencia de Deposicion Inertial
						 // Cs=1+2/FDiametroAgregadoSoot*(1.23+0.41*exp(-0.88*(FDiametroAgregadoSoot/(2.*CaminoLibreMedio))));
						 StokesPL=(2./9.*SCF_in*FDensidadSootCapa_in*VelocidadIntersticialPL*pow(FDiametroAgregadoSoot,2.))/(2.*ViscosidadCinematica*FDiametroAgregadoSoot);
						 EficienciaInercialPL=pow(StokesPL,2)/pow((StokesPL+0.25),2);

						 // Calculo de la Eficiencia  de Deposicion Combinada
						 EficienciaCombinadaPL=(EficienciaBrownianaPL+EficienciaIntercepcionPL+EficienciaInercialPL)-(EficienciaBrownianaPL*EficienciaIntercepcionPL+EficienciaBrownianaPL*EficienciaInercialPL+EficienciaIntercepcionPL*EficienciaInercialPL)+(EficienciaBrownianaPL*EficienciaIntercepcionPL*EficienciaInercialPL);

						 // Calculo de la Eficiencia de Deposicion Total
						 FEficienciaPL[j][i]=1.-exp(-3*EficienciaCombinadaPL*(1-FPorosidadCapaParticulas_in)*FEspesorSoot[j][i]/(2.*FPorosidadCapaParticulas_in*FDiametroAgregadoSoot));

						 */

						// FLayerSootMass[j][i]=FLayerSootMass[j][i]+FEficienciaPL[j][i]*FCanal[j][0]->GetYespecie(i,4)*
						// FVelocidadPared[j][0][i]*densidad*FAreaVCCanalEntrada[j][i]*FDeltaTimeDPF;
						FLayerSootMassDEP = FEficienciaPL[j][i] * FCanal[j][0]->GetYespecie(i,
											4) * FVelocidadPared[j][0][i] * densidad * FAreaVCCanalEntrada[j][i] * FDeltaTimeDPF;
						FLayerSootMassDep[j][i] = FLayerSootMassDep[j][i] + FLayerSootMassDEP;
						FLayerSootMass[j][i] = FLayerSootMass[j][i] + FLayerSootMassDEP;
						FFraccionMasicaEspecieEntrantePared[j][i][4] = FCanal[j][0]->GetYespecie(i, 4) * (1 - FEficienciaPL[j][i]); // Soot
						masasootatrapadaVC = FEficiencia[j][i] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FVelocidadPared[j][0][i] *
											 densidad * FAreaVCCanalEntrada[j][i] * FDeltaTimeDPF;

						FMasaSootUC[j][i] = FMasaSootUC[j][i] + masasootatrapadaVC / FNumeroUnidadesCelulares[j][i];
						FWallSootMass[j][i] = FWallSootMass[j][i] + masasootatrapadaVC;
						FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
						FBeamSootMass[j] = FBeamSootMass[j] + FCVSootMass[j][i] * FCanal[j][0]->getNumeroCanales();

						// FEspesorSoot[j][i]=(FCanal[j][0]->GetLadoCanal(i)-pow(pow(FCanal[j][0]->GetLadoCanal(i),2.)-FLayerSootMass[j][i]/(FLongitudVC[j][i]*FDensidadSootCapa),0.5))/2.;
						FEspesorSootDep[j][i] = ((FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSootIn[j][i]
												  - pow(pow(FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSootIn[j][i],
															2.) - FLayerSootMassDep[j][i] / (FLongitudVC[j][i] * FDensidadSootCapa_dep), 0.5)) / 2.)
												/ (FPorosidad[j][i] + (FMasaSootUC[j][i] - FMasaSootLimitInfUC[j][i]) * (1 - FPorosidad[j][i]) /
												   (FMasaSootLimitUC - FMasaSootLimitInfUC[j][i]));
						FEspesorSoot[j][i] = FEspesorSootIn[j][i] + FEspesorSootDep[j][i];
						FAreaVCCanalEntrada[j][i] = 4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];

						FFraccionMasicaEspecieEntrantePared[j][i][0] = FCanal[j][0]->GetYespecie(i, 0)
								+ FCanal[j][0]->GetYespecie(i, 0) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieEntrantePared[j][i][1] = FCanal[j][0]->GetYespecie(i, 1)
								+ FCanal[j][0]->GetYespecie(i, 1) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieEntrantePared[j][i][2] = FCanal[j][0]->GetYespecie(i, 2)
								+ FCanal[j][0]->GetYespecie(i, 2) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieEntrantePared[j][i][3] = FCanal[j][0]->GetYespecie(i, 3)
								+ FCanal[j][0]->GetYespecie(i, 3) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieEntrantePared[j][i][5] = FCanal[j][0]->GetYespecie(i, 5)
								+ FCanal[j][0]->GetYespecie(i, 5) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieEntrantePared[j][i][6] = FCanal[j][0]->GetYespecie(i, 6)
								+ FCanal[j][0]->GetYespecie(i, 6) * FCanal[j][0]->GetYespecie(i,
										4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieEntrantePared[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
									+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
											4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieEntrantePared[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
									+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
											4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieEntrantePared[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
									+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
											4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieEntrantePared[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
									+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
											4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieEntrantePared[j][i][9] = FCanal[j][0]->GetYespecie(i, 9)
									+ FCanal[j][0]->GetYespecie(i, 9) * FCanal[j][0]->GetYespecie(i,
											4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						}

						FFraccionMasicaEspecieSalida[j][i][0] = FFraccionMasicaEspecieEntrantePared[j][i][0]
																+ FFraccionMasicaEspecieEntrantePared[j][i][0] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						FFraccionMasicaEspecieSalida[j][i][1] = FFraccionMasicaEspecieEntrantePared[j][i][1]
																+ FFraccionMasicaEspecieEntrantePared[j][i][1] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						FFraccionMasicaEspecieSalida[j][i][2] = FFraccionMasicaEspecieEntrantePared[j][i][2]
																+ FFraccionMasicaEspecieEntrantePared[j][i][2] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						FFraccionMasicaEspecieSalida[j][i][3] = FFraccionMasicaEspecieEntrantePared[j][i][3]
																+ FFraccionMasicaEspecieEntrantePared[j][i][3] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						FFraccionMasicaEspecieSalida[j][i][4] = FFraccionMasicaEspecieEntrantePared[j][i][4] * (1 - FEficiencia[j][i]); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FFraccionMasicaEspecieEntrantePared[j][i][5]
																+ FFraccionMasicaEspecieEntrantePared[j][i][5] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						FFraccionMasicaEspecieSalida[j][i][6] = FFraccionMasicaEspecieEntrantePared[j][i][6]
																+ FFraccionMasicaEspecieEntrantePared[j][i][6] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FFraccionMasicaEspecieEntrantePared[j][i][7]
																	+ FFraccionMasicaEspecieEntrantePared[j][i][7] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																	(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
							FFraccionMasicaEspecieSalida[j][i][8] = FFraccionMasicaEspecieEntrantePared[j][i][8]
																	+ FFraccionMasicaEspecieEntrantePared[j][i][8] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																	(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FFraccionMasicaEspecieEntrantePared[j][i][7]
																	+ FFraccionMasicaEspecieEntrantePared[j][i][7] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																	(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
							FFraccionMasicaEspecieSalida[j][i][8] = FFraccionMasicaEspecieEntrantePared[j][i][8]
																	+ FFraccionMasicaEspecieEntrantePared[j][i][8] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																	(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
							FFraccionMasicaEspecieSalida[j][i][9] = FFraccionMasicaEspecieEntrantePared[j][i][9]
																	+ FFraccionMasicaEspecieEntrantePared[j][i][9] * FFraccionMasicaEspecieEntrantePared[j][i][4] * FEficiencia[j][i] /
																	(1 - FFraccionMasicaEspecieEntrantePared[j][i][4]);
						}

					} else {
						// Pared porosa saturada o con carga máxima; el filtrado produce el aumento de la capa de part�culas cuyas propiedades se dan al inicio de la ejecuci�n
						// Hay que calcular la eficiencia de la capa
						// FLayerSootMass[j][i]=FLayerSootMass[j][i]+FEficienciaPL[j][i]*FCanal[j][0]->GetYespecie(i,4)*
						// FVelocidadPared[j][0][i]*densidad*FAreaVCCanalEntrada[j][i]*FDeltaTimeDPF;

						/* if(FVelocidadPared[j][0][i]<0){
						 VelocidadIntersticialPL=0.;
						 }else{
						 VelocidadIntersticialPL=FVelocidadPared[j][0][i]/FPorosidadCapaParticulas_in;
						 }
						 CoefDifusionParticulasPL=__cons::Kb*Temperatura*SCF_in/(3.*Pi*FCanal[j][0]->GetViscosidadDinamica(i)*FDiametroAgregadoSoot);
						 PecletPL=VelocidadIntersticialPL*FDiametroAgregadoSoot/CoefDifusionParticulasPL;

						 if(PecletPL==0 || PecletPL<0){
						 EficienciaBrownianaPL=0.;
						 }else{
						 EficienciaBrownianaPL=3.5*funcion_gPL_in*pow(PecletPL,-2./3.);
						 }
						 if(EficienciaBrownianaPL>1) EficienciaBrownianaPL=1.;

						 // Calculo de Eficiencia de Deposicion debido a Intercepcion
						 FParametroIntercepcionPL=FDiametroAgregadoSoot/FDiametroAgregadoSoot;
						 EficienciaIntercepcionPL=1.5*pow(FParametroIntercepcionPL,2.)*pow(funcion_gPL_in,3.)/pow(1.+FParametroIntercepcionPL,(3.-2.*FPorosidadCapaParticulas_in)/(3.*FPorosidadCapaParticulas_in));

						 // Calculo de Eficiencia de Deposicion Inertial
						 // Cs=1+2/FDiametroAgregadoSoot*(1.23+0.41*exp(-0.88*(FDiametroAgregadoSoot/(2.*CaminoLibreMedio))));
						 StokesPL=(2./9.*SCF_in*FDensidadSootCapa_in*VelocidadIntersticialPL*pow(FDiametroAgregadoSoot,2.))/(2.*ViscosidadCinematica*FDiametroAgregadoSoot);
						 EficienciaInercialPL=pow(StokesPL,2)/pow((StokesPL+0.25),2);

						 // Calculo de la Eficiencia  de Deposicion Combinada
						 EficienciaCombinadaPL=(EficienciaBrownianaPL+EficienciaIntercepcionPL+EficienciaInercialPL)-(EficienciaBrownianaPL*EficienciaIntercepcionPL+EficienciaBrownianaPL*EficienciaInercialPL+EficienciaIntercepcionPL*EficienciaInercialPL)+(EficienciaBrownianaPL*EficienciaIntercepcionPL*EficienciaInercialPL);

						 // Calculo de la Eficiencia de Deposicion Total
						 FEficienciaPL[j][i]=1.-exp(-3*EficienciaCombinadaPL*(1-FPorosidadCapaParticulas_in)*FEspesorSoot[j][i]/(2.*FPorosidadCapaParticulas_in*FDiametroAgregadoSoot));

						 */

						FLayerSootMassDEP = FEficienciaPL[j][i] * FCanal[j][0]->GetYespecie(i,
											4) * FVelocidadPared[j][0][i] * densidad * FAreaVCCanalEntrada[j][i] * FDeltaTimeDPF;
						FLayerSootMassDep[j][i] = FLayerSootMassDep[j][i] + FLayerSootMassDEP;
						//-(FKreg1[j][i]+FKreg2[j][i])*__PM::C*FAreaVCCanalEntrada[j][i];
						FLayerSootMass[j][i] = FLayerSootMass[j][i] + FLayerSootMassDEP;
						FCVSootMass[j][i] = FWallSootMass[j][i] + FLayerSootMass[j][i];
						FBeamSootMass[j] = FBeamSootMass[j] + FCVSootMass[j][i] * FCanal[j][0]->getNumeroCanales();

						// FEspesorSoot[j][i]=(FCanal[j][0]->GetLadoCanal(i)-pow(pow(FCanal[j][0]->GetLadoCanal(i),2.)-FLayerSootMass[j][i]/(FLongitudVC[j][i]*FDensidadSootCapa),0.5))/2.;
						FEspesorSootDep[j][i] = ((FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSootIn[j][i]
												  - pow(pow(FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSootIn[j][i],
															2.) - FLayerSootMassDep[j][i] / (FLongitudVC[j][i] * FDensidadSootCapa_dep), 0.5)) / 2.)
												/ (FPorosidad[j][i] + (FMasaSootUC[j][i] - FMasaSootLimitInfUC[j][i]) * (1 - FPorosidad[j][i]) /
												   (FMasaSootLimitUC - FMasaSootLimitInfUC[j][i]));
						FEspesorSoot[j][i] = FEspesorSootIn[j][i] + FEspesorSootDep[j][i];
						FAreaVCCanalEntrada[j][i] = 4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FLongitudVC[j][i];

						FFraccionMasicaEspecieSalida[j][i][0] = FCanal[j][0]->GetYespecie(i, 0)
																+ FCanal[j][0]->GetYespecie(i, 0) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][1] = FCanal[j][0]->GetYespecie(i, 1)
																+ FCanal[j][0]->GetYespecie(i, 1) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][2] = FCanal[j][0]->GetYespecie(i, 2)
																+ FCanal[j][0]->GetYespecie(i, 2) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][3] = FCanal[j][0]->GetYespecie(i, 3)
																+ FCanal[j][0]->GetYespecie(i, 3) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][4] = FCanal[j][0]->GetYespecie(i, 4) * (1 - FEficienciaPL[j][i]); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FCanal[j][0]->GetYespecie(i, 5)
																+ FCanal[j][0]->GetYespecie(i, 5) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						FFraccionMasicaEspecieSalida[j][i][6] = FCanal[j][0]->GetYespecie(i, 6)
																+ FCanal[j][0]->GetYespecie(i, 6) * FCanal[j][0]->GetYespecie(i,
																		4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
																	+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
																	+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7)
																	+ FCanal[j][0]->GetYespecie(i, 7) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8)
																	+ FCanal[j][0]->GetYespecie(i, 8) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
							FFraccionMasicaEspecieSalida[j][i][9] = FCanal[j][0]->GetYespecie(i, 9)
																	+ FCanal[j][0]->GetYespecie(i, 9) * FCanal[j][0]->GetYespecie(i,
																			4) * FEficienciaPL[j][i] / (1 - FCanal[j][0]->GetYespecie(i, 4));
						}
					}
				}
			}
			FDPFSootMass = FDPFSootMass + FBeamSootMass[j];
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::SubmodeloFiltrado en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::SubmodeloRegeneracion() {
	try {

		double X_O2, X_HC, X_NO2, X_NO, X_CO, X_Combustible, Y_NO2, Ka1, Ka2, Ka3, Ka4, Ka5, Ka6, Ka7, Ka8, Ka9, Ka10, Ka11,
			   Ka12, Ka13, G1, G2, G3, G4, Kp3, Kp4, Eq3, Eq4, Kreg1, Kreg2;
		double Raire = 0., PMaire = 0., temperatura = 0.;

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[i][0]->getNin(); i++) {
				if(FCanal[j][0]->getNodoInicialFuente() <= i) {
					// C�lculo de las tasas de reacci�n debido a que el medio poroso posee catalizadores.
					if(FDPFCatalizada) {
						/*
						 // C�lculo de las fracciones molares
						 Raire=FCanal[j][0]->Yespecie[i][0]*__R::O2+FCanal[j][0]->Yespecie[i][1]*__R::CO2+FCanal[j][0]->Yespecie[i][2]*__R::H2O+
						 (1-0.01292*(1-FCanal[j][0]->Yespecie[i][2])-FCanal[j][0]->Yespecie[i][0]-FCanal[j][0]->Yespecie[i][1]-FCanal[j][0]->Yespecie[i][2])*__R::N2+
						 (0.01292*(1-FCanal[j][0]->Yespecie[i][2]))*__R::Ar;
						 PMaire=__R::Universal/Raire;
						 X_O2=FCanal[j][0]->Yespecie[i][0]*PMaire/__PM::O2;
						 X_HC=FCanal[j][0]->Yespecie[i][3]*PMaire/__PM::UHC;
						 Y_NO2=FRatioNO2_NOx*FCanal[j][0]->Yespecie[i][5];
						 X_NO2=FRatioNO2_NOx*FCanal[j][0]->Yespecie[i][5]*PMaire/__PM::NO2;
						 X_NO=(1-FRatioNO2_NOx)*FCanal[j][0]->Yespecie[i][5]*PMaire/__PM::NO;
						 X_CO=FCanal[j][0]->Yespecie[i][6]*PMaire/__PM::CO;

						 // C�lculo de los t�rminos de inhibici�n
						 Ka1=ka01*exp(-Ha1/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka2=ka02*exp(-Ha2/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka3=ka03*exp(-Ha3/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka4=ka04*exp(-Ha4/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka5=ka05*exp(-Ha5/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka6=ka06*exp(-Ha6/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka7=ka07*exp(-Ha7/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka8=ka08*exp(-Ha8/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka9=ka09*exp(-Ha9/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka10=ka010*exp(-Ha10/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka11=ka011*exp(-Ha11/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka12=ka012*exp(-Ha12/(__R::Universal*(FTPared[j][i][2]+273.)));
						 Ka13=ka013*exp(-Ha13/(__R::Universal*(FTPared[j][i][2]+273.)));

						 G1=(FTPared[j][i][2]+273.)*pow(1+Ka1*X_CO+Ka2*X_HC,2.)*(1+Ka3*pow(X_CO*X_HC,2.))*(1+Ka4*pow(X_NO,0.7));
						 G2=(FTPared[j][i][2]+273.)*pow(1+Ka5*X_CO+Ka6*X_HC,2.)*(1+Ka7*pow(X_CO*X_HC,2.))*(1+Ka8*pow(X_NO,0.7));
						 G3=(1+Ka9*pow(X_O2,1.5));
						 G4=(FTPared[j][i][2]+273.)*pow(1+Ka10*X_CO+Ka11*X_HC,2.)*(1+Ka12*pow(X_CO*X_HC,2.))*(1+Ka13*pow(X_NO,0.7));

						 // Factores que tienen en cuenta el equilibrio NO2-NO
						 Kp3=0; // �C�mo se calcula?
						 Kp4=0; // �C�mo se calcula?
						 Eq3=1-X_NO2/(X_NO*pow(X_O2,0.5)*Kp3);
						 Eq4=1-X_NO2*pow(X_O2,0.5)/(X_NO*Kp4);

						 // C�lculo de las tasas de reacci�n [mol/m3s]   A[K/m3s]
						 FR1[j][i]=A1*exp(-E1/(__R::Universal*(FTPared[j][i][2]+273.)))*X_CO*X_O2/G2; // Combusti�n de CO
						 FR2[j][i]=A2*exp(-E2/(__R::Universal*(FTPared[j][i][2]+273.)))*X_HC*X_O2/(G1*G3); // Combusti�n de HC
						 FR3[j][i]=A3*exp(-E3/(__R::Universal*(FTPared[j][i][2]+273.)))*X_NO*X_O2*Eq3/G4;
						 FR4[j][i]=A4*exp(-E4/(__R::Universal*(FTPared[j][i][2]+273.)))*X_NO2*Eq4/G4;
						 if(FNumeroEspecies==10){
						 X_Combustible=FCanal[j][0]->Yespecie[i][7]*PMaire/__PM::Diesel;
						 FR5[j][i]=A2*exp(-E2/(__R::Universal*(FTPared[j][i][2]+273.)))*X_Combustible*X_O2/(G1*G3);
						 }
						 */

						FKreg1[i][j] = FFactorFrecuenciaO2 * (__units::degCToK(FTPared[j][i][2])) * exp(-FEnergiaActO2 /
									   (__R::Universal * __units::degCToK(FTPared[j][i][2]))); // Reacci�n con O2
						FKreg2[i][j] = FFactorFrecuenciaNO2 * (__units::degCToK(FTPared[j][i][2])) * exp(-FEnergiaActNO2 /
									   (__R::Universal * __units::degCToK(FTPared[j][i][2]))); // Reacci�n con N2

						temperatura = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
										  2.) / (FCanal[j][0]->GetRMezcla(i) * FCanal[j][0]->GetGamma(i));
						FFraccionMasicaNO2Entrada[j][i] = Interpolacion_bidimensional(FCanal[j][0]->GetYespecie(i, 5), temperatura,
														  FTemperaturasTabla, FFraccionNOxTabla, FMapa_ProporcionNO2,
														  FNumeroDatos_FraccionesNOx, FNumeroDatos_Temperaturas);

						FFraccionMasicaEspecieSalida[j][i][0] = FCanal[j][0]->GetYespecie(i, 0)
																* (1 - exp(-FSupEspecifica[j][i] * FEspesorSoot[j][i] * FKreg1[i][j] * FIndiceCompletitud1 / FVelocidadPared[j][0][i]));
						FFraccionMasicaEspecieSalida[j][i][1] = FCanal[j][0]->GetYespecie(i, 1)
																* (1
																		+ exp(
																				FSupEspecifica[j][i] * FEspesorSoot[j][i] * (FKreg1[i][j] * (2 * (FIndiceCompletitud1 - 0.5)) + FKreg2[i][j] *
																						(FIndiceCompletitud2 - 1))
																				/ FVelocidadPared[j][0][i])); // CO2
						FFraccionMasicaEspecieSalida[j][i][2] = FCanal[j][0]->GetYespecie(i, 2); // H2O
						FFraccionMasicaEspecieSalida[j][i][3] = FCanal[j][0]->GetYespecie(i, 3);    // HC
						FFraccionMasicaEspecieSalida[j][i][4] = FCanal[j][0]->GetYespecie(i,
																4) * (1 - exp(-FEficiencia[j][i] * FEspesorSoot[j][i])); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FCanal[j][0]->GetYespecie(i, 5);   // NOx
						FFraccionMasicaEspecieSalida[j][i][6] =
							FCanal[j][0]->GetYespecie(i, 6)
							* (1
							   + exp(
								   FSupEspecifica[j][i] * FEspesorSoot[j][i] * (FKreg1[i][j] * (2 * (1 - FIndiceCompletitud1)) + FKreg2[i][j] *
										   (2 - FIndiceCompletitud2))
								   / FVelocidadPared[j][0][i])); // CO
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7);   // N2
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // EGR
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7); // Combustible
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // N2
							FFraccionMasicaEspecieSalida[j][i][9] = FCanal[j][0]->GetYespecie(i, 9);   // EGR
						}
						FFraccionMasicaNO2Salida[j][i] = FFraccionMasicaNO2Entrada[j][i]
														 * (1 - exp(-FSupEspecifica[j][i] * FEspesorSoot[j][i] * FKreg2[i][j] * FIndiceCompletitud2 /
																 FVelocidadPared[j][0][i])); // NO2 a la salida

						/*FTasaFraccionMasicaEspecie[j][i][0]=(-FIndiceCompletitud1*FKreg1[j][i]-0.5*FR1[j][i]-7.5*FR2[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][0]; // O2
						 FTasaFraccionMasicaEspecie[j][i][1]=(FR1[j][i]+5*FR2[j][i]+2*(FIndiceCompletitud1-0.5)*FKreg1[j][i]+(FIndiceCompletitud2-1)*FKreg2[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][1]; // CO2
						 FTasaFraccionMasicaEspecie[j][i][2]=2.5*FR2[j][i]*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][2]; // H2O
						 FTasaFraccionMasicaEspecie[j][i][3]=-FR2[j][i]*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][3];    // HC
						 FTasaFraccionMasicaEspecie[j][i][4]=-FEficiencia[j][i]*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][4]*FVelocidadPared[j][1][i]; // Soot
						 FTasaFraccionMasicaEspecie[j][i][5]=0.;   // NOx
						 FTasaFraccionMasicaEspecie[j][i][6]=(-FR1[j][i]+2*(1-FIndiceCompletitud1)*FKreg1[j][i]+(2-FIndiceCompletitud2)*FKreg2[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][6]; // CO
						 if(FNumeroEspecies==9){
						 FTasaFraccionMasicaEspecie[j][i][7]=0.;   // N2
						 FTasaFraccionMasicaEspecie[j][i][8]=0.;   // EGR
						 }else if(FNumeroEspecies==10){
						 FTasaFraccionMasicaEspecie[j][i][7]=-FR5[j][i]*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][5];   // Combustible
						 FTasaFraccionMasicaEspecie[j][i][8]=0.;   // N2
						 FTasaFraccionMasicaEspecie[j][i][9]=0.;   // EGR
						 }   */

						/* C�lculo del calor liberado durante el proceso de regeneraci�n */

						//FQreg[j][i]=(FIncrH_reg1/(__PM::O2*FIndiceCompletitud1/1000)*FCanal[j][0]->Yespecie[i][0]*FKreg1[j][i]+FIncrH_reg2/(__PM::NO2*FIndiceCompletitud2/1000)*Y_NO2*FKreg2[j][i])*FSupEspecifica[j][i];   // (J/kgs) Ha de estar en (W)
						//FQ1[j][i]=FIncrH_R1/(__PM::O2*0.5/1000)*FSupEspecifica[j][i]*FR1[j][i]*FCanal[j][0]->Yespecie[i][0];
						//FQ2[j][i]=FIncrH_R2/(__PM::O2*(4+6.93/4)/1000)*FSupEspecifica[j][i]*FR2[j][i]*FCanal[j][0]->Yespecie[i][0];
					} else {

						FKreg1[i][j] = FFactorFrecuenciaO2 * (__units::degCToK(FTPared[j][i][2])) * exp(-FEnergiaActO2 /
									   (__R::Universal * __units::degCToK(FTPared[j][i][2])));
						FKreg2[i][j] = FFactorFrecuenciaNO2 * (__units::degCToK(FTPared[j][i][2])) * exp(-FEnergiaActNO2 /
									   (__R::Universal * __units::degCToK(FTPared[j][i][2])));

						temperatura = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
										  2.) / (FCanal[j][0]->GetRMezcla(i) * FCanal[j][0]->GetGamma(i));
						FFraccionMasicaNO2Entrada[j][i] = Interpolacion_bidimensional(FCanal[j][0]->GetYespecie(i, 5), temperatura,
														  FTemperaturasTabla, FFraccionNOxTabla, FMapa_ProporcionNO2,
														  FNumeroDatos_FraccionesNOx, FNumeroDatos_Temperaturas);

						FFraccionMasicaEspecieSalida[j][i][0] = FCanal[j][0]->GetYespecie(i, 0)
																* (1 - exp(-FSupEspecifica[j][i] * FEspesorSoot[j][i] * FKreg1[i][j] * FIndiceCompletitud1 / FVelocidadPared[j][0][i]));
						FFraccionMasicaEspecieSalida[j][i][1] = FCanal[j][0]->GetYespecie(i, 1)
																* (1
																		+ exp(
																				FSupEspecifica[j][i] * FEspesorSoot[j][i] * (FKreg1[i][j] * (2 * (FIndiceCompletitud1 - 0.5)) + FKreg2[i][j] *
																						(FIndiceCompletitud2 - 1))
																				/ FVelocidadPared[j][0][i])); // CO2
						FFraccionMasicaEspecieSalida[j][i][2] = FCanal[j][0]->GetYespecie(i, 2); // H2O
						FFraccionMasicaEspecieSalida[j][i][3] = FCanal[j][0]->GetYespecie(i, 3);    // HC
						FFraccionMasicaEspecieSalida[j][i][4] = FCanal[j][0]->GetYespecie(i,
																4) * (1 - exp(-FEficiencia[j][i] * FEspesorSoot[j][i])); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FCanal[j][0]->GetYespecie(i, 5);   // NOx
						FFraccionMasicaEspecieSalida[j][i][6] =
							FCanal[j][0]->GetYespecie(i, 6)
							* (1
							   + exp(
								   FSupEspecifica[j][i] * FEspesorSoot[j][i] * (FKreg1[i][j] * (2 * (1 - FIndiceCompletitud1)) + FKreg2[i][j] *
										   (2 - FIndiceCompletitud2))
								   / FVelocidadPared[j][0][i])); // CO
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7);   // N2
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // EGR
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7); // Combustible
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // N2
							FFraccionMasicaEspecieSalida[j][i][9] = FCanal[j][0]->GetYespecie(i, 9);   // EGR
						}
						FFraccionMasicaNO2Salida[j][i] = FFraccionMasicaNO2Entrada[j][i]
														 * (1 - exp(-FSupEspecifica[j][i] * FEspesorSoot[j][i] * FKreg2[i][j] * FIndiceCompletitud2 /
																 FVelocidadPared[j][0][i])); // NO2 a la salida

						/* FTasaFraccionMasicaEspecie[j][i][0]=(-FIndiceCompletitud1*FKreg1[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][0]; // O2
						 FTasaFraccionMasicaEspecie[j][i][1]=(2*(FIndiceCompletitud1-0.5)*FKreg1[j][i]+(FIndiceCompletitud2-1)*FKreg2[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][1]; // CO2
						 FTasaFraccionMasicaEspecie[j][i][2]=0.;   // H2O
						 FTasaFraccionMasicaEspecie[j][i][3]=0.;   // HC
						 FTasaFraccionMasicaEspecie[j][i][4]=-FEficiencia[j][i]*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FVelocidadPared[j][1][i]*FCanal[j][0]->Yespecie[i][4]; // Soot
						 FTasaFraccionMasicaEspecie[j][i][5]=0.;   // NOx
						 FTasaFraccionMasicaEspecie[j][i][6]=(2*(1-FIndiceCompletitud1)*FKreg1[j][i]+(2-FIndiceCompletitud2)*FKreg2[j][i])*FSupEspecifica[j][i]*(4*FLongitudVC[j][i]*FCanal[j][0]->LadoCanal[i])*FCanal[j][0]->Yespecie[i][6]; // CO

						 if(FNumeroEspecies==9){
						 FTasaFraccionMasicaEspecie[j][i][7]=0.;   // N2
						 FTasaFraccionMasicaEspecie[j][i][8]=0.;   // EGR
						 }else if(FNumeroEspecies==10){
						 FTasaFraccionMasicaEspecie[j][i][7]=0.;   // Combustible
						 FTasaFraccionMasicaEspecie[j][i][8]=0.;   // N2
						 FTasaFraccionMasicaEspecie[j][i][9]=0.;   // EGR
						 } */

						/* C�lculo del calor liberado durante el proceso de regeneraci�n */
						//FQreg[j][i]=(FIncrH_reg1/(__PM::O2*FIndiceCompletitud1)*FCanal[j][0]->Yespecie[i][0]*FKreg1[j][i]+FIncrH_reg2/(__PM::NO2*FIndiceCompletitud2)*Y_NO2*FKreg2[j][i])*FSupEspecifica[j][i];   // (J/kgs)
						//FQ1[j][i]=FIncrH_R1/(__PM::O2*0.5)*FSupEspecifica[j][i]*FR1[j][i]*FCanal[j][0]->Yespecie[i][0];
						//FQ2[j][i]=FIncrH_R2/(__PM::O2*(4+6.93/4))*FSupEspecifica[j][i]*FR2[j][i]*FCanal[j][0]->Yespecie[i][0];
					}
					for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
						if(FFraccionMasicaEspecieSalida[j][i][k] < 0.) {
							FFraccionMasicaEspecieSalida[j][i][k] = 0.;
						}
					}

				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::SubmodeloRegeneracion en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoKsoot() {
	try {
		double Kuwabara_in, Kuwabara_dep, Knudsen_in, Knudsen_dep, ViscosidadCinematica, Temperatura, CaminoLibreMedio, SCF_in,
			   SCF_dep;
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				if(FLayerSootMass[j][i] > 0.) {
					Kuwabara_in = 2 * (2. - 1.8 * pow(1. - FPorosidadCapaParticulas_in,
													  1. / 3.) - FPorosidadCapaParticulas_in - 0.2 * pow(1 - FPorosidadCapaParticulas_in, 2.)) / 9
								  / (1. - FPorosidadCapaParticulas_in);
					Kuwabara_dep = 2 * (2. - 1.8 * pow(1. - FPorosidadCapaParticulas_dep,
													   1. / 3.) - FPorosidadCapaParticulas_dep - 0.2 * pow(1 - FPorosidadCapaParticulas_dep, 2.)) / 9
								   / (1. - FPorosidadCapaParticulas_dep);
					ViscosidadCinematica = FCanal[j][0]->GetViscosidadDinamica(i) / FCanal[j][0]->GetDensidad(i);
					Temperatura = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
									  2.) / (FCanal[j][0]->GetGamma(i) * FCanal[j][0]->GetRMezcla(i));
					CaminoLibreMedio = ViscosidadCinematica * pow(__cons::Pi / (2. * FCanal[j][0]->GetRMezcla(i) * Temperatura), 0.5);
					Knudsen_in = 2 * CaminoLibreMedio / FDiametroPoroPL_in; // Referido al di�metro de los poros de la capa
					Knudsen_dep = 2 * CaminoLibreMedio / FDiametroPoroPL_dep; // Referido al di�metro de los poros de la capa
					//(2*FPorosidadCapaParticulas*FDiametroAgregadoSoot/3./(1.-FPorosidadCapaParticulas));
					SCF_in = 1 + Knudsen_in * (1.257 + 0.4 * exp(-1.1 / Knudsen_in));
					SCF_dep = 1 + Knudsen_dep * (1.257 + 0.4 * exp(-1.1 / Knudsen_dep));
					FKsootIn[j][i] = Kuwabara_in * SCF_in * pow(FDiametroAgregadoSoot, 2.);
					FKsootDep[j][i] = Kuwabara_dep * SCF_dep * pow(FDiametroAgregadoSoot, 2.);
					// FKsoot[j][i]=(FKsoot_in[j][i]*FKsoot_dep[j][i])/(FKsoot_in[j][i]*FEspesorSoot_dep[j][i]/(FEspesorSoot_in[j][i]+FEspesorSoot_dep[j][i])+FKsoot_in[j][i]*FEspesorSoot_in[j][i]/(FEspesorSoot_in[j][i]+FEspesorSoot_dep[j][i]));
					if(FEspesorSootDep[j][i] == 0) {
						FKsoot[j][i] = FKsootIn[j][i];
					} else {
						FKsoot[j][i] = (FKsootIn[j][i] * FKsootDep[j][i])
									   / (FKsootIn[j][i] * FEspesorSootDep[j][i] / (FEspesorSootIn[j][i] + FEspesorSootDep[j][i])
										  + FKsootDep[j][i] * FEspesorSootIn[j][i] / (FEspesorSootIn[j][i] + FEspesorSootDep[j][i]));
					}

				}
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoKsoot en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

/*id TDPF::CalculoEspesorSoot()
 {
 try
 {
 double temp,densidad;

 for(int j=0;j<FNumeroHacesCanales;j++){
 for(int i=0;i<FCanal[j][0]->getNin();i++){
 temp=pow(FCanal[j][0]->Asonido[i]*__cons::ARef,2.)/(FCanal[j][0]->R[i]*FCanal[j][0]->Gamma[i]);
 densidad=FCanal[j][0]->Presion[i]*1e5/(FCanal[j][0]->R[i]*temp);

 FMasaSootLayer[j][i]=FMasaSootLayer[j][i]+FEficiencia[j][i]*FCanal[j][0]->Yespecie[i][4]*
 FVelocidadPared[j][0][i]*densidad*FAreaVCCanalEntrada[j][i]*FDeltaTimeDPF-
 (FKreg1[j][i]+FKreg2[j][i])*__PM::C*FAreaVCCanalEntrada[j][i];
 *//*(FKreg1[j][i]+FKreg2[j][i])*densidad*FAreaVCCanalEntrada[j][i]*FCanal[j][0]->getXRef()*/;

/*     FEspesorSoot[j][i]=(FCanal[j][0]->LadoCanal[i]-pow(pow(FCanal[j][0]->LadoCanal[i],2.)-
 FMasaSootLayer[j][i]/(FLongitudVC[j][i]*FDensidadSootCapa),0.5))/2.;

 FAreaVCCanalEntrada[j][i]=4*(FCanal[j][0]->LadoCanal[i]-FEspesorSoot[j][i])*FLongitudVC[j][i];

 }
 }

 }
 catch(Exception &N)
 {
 std::cout << "ERROR: TDPF::CalculoEspesorSoot en la DPF " << FNumeroDPF << std::endl;
 std::cout << "Tipo de error: " << N.what() << std::endl;
 throw Exception("");
 }
 } */

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
void TDPF::IniciaVariablesTransmisionCalor(double tamb) {
	try {

		int NumCanalesUltimoHaz = 0;

		FSUMTime = new double[FNumeroHacesCanales];
		FResistAxialAnt = new double*[FNumeroHacesCanales];
		FResistAxialPost = new double*[FNumeroHacesCanales];
		FResistEntreHacesAnt = new double*[FNumeroHacesCanales];
		FResistEntreHacesPost = new double*[FNumeroHacesCanales];
		FCapEntrada = new double*[FNumeroHacesCanales];
		FCapSalida = new double*[FNumeroHacesCanales];
		FCapPared = new double*[FNumeroHacesCanales];
		FResistRadial = new double**[FNumeroHacesCanales];
		FResistConv = new double**[FNumeroHacesCanales];
		FTPared = new double**[FNumeroHacesCanales];
		FTParedAnt = new double**[FNumeroHacesCanales];
		FSUMTPPromedio = new double***[FNumeroHacesCanales];
		FSUMTPPromedioConc = new double**[2];
		FDiametroHazExt = new double[FNumeroHacesCanales];
		FDiametroHazInt = new double[FNumeroHacesCanales];
		FResistConduccionAislante = new double[FCanal[0][0]->getNin()];
		FResistConduccionAire = new double[FCanal[0][0]->getNin()];
		FResistRadiacionAire = new double[FCanal[0][0]->getNin()];
		FResistConduccionMetal = new double[FCanal[0][0]->getNin()];
		FResistConveccionExt = new double[FCanal[0][0]->getNin()];
		FResistAxialMetalExtAnt = new double[FCanal[0][0]->getNin()];
		FResistAxialMetalExtPost = new double[FCanal[0][0]->getNin()];
		FCapNodoExteriorSuperficie = new double[FCanal[0][0]->getNin()];
		FCapNodoMedioSuperficie = new double[FCanal[0][0]->getNin()];
		FCapNodoInteriorSuperficie = new double[FCanal[0][0]->getNin()];
		FNumCanalesHaz = new double[FNumeroHacesCanales];
		FSuperficieTCHazAnt = new double[FNumeroHacesCanales];
		FSuperficieTCHazPost = new double[FNumeroHacesCanales];
		FAreaInternaHaz = new double[FNumeroHacesCanales];
		FTSuperficie = new double*[FCanal[0][0]->getNin()];

		FRg_int_ext = new double[FCanal[0][0]->getNin()];
		FR_int_radiacion = new double[FCanal[0][0]->getNin()];
		FR_ext_RadInt = new double[FCanal[0][0]->getNin()];

		for(int i = 0; i < 2; i++) {
			FSUMTPPromedioConc[i] = new double*[3];
		}
		for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 3; j++) {
				FSUMTPPromedioConc[i][j] = new double[FCanal[0][0]->getNin()];
			}
		}
		for(int k = 0; k < FCanal[0][0]->getNin(); k++) {
			FSUMTPPromedioConc[0][0][k] = 0.;
			FSUMTPPromedioConc[0][1][k] = 0.;
			FSUMTPPromedioConc[0][2][k] = 0.;
			FSUMTPPromedioConc[1][0][k] = 0.;
			FSUMTPPromedioConc[1][1][k] = 0.;
			FSUMTPPromedioConc[1][2][k] = 0.;
		}

		for(int i = 0; i < FCanal[0][0]->getNin(); i++) {
			FTSuperficie[i] = new double[3];
			FTSuperficie[i][0] = FTIniParedExt;
			FTSuperficie[i][1] = FTIniParedExt - 2.;
			FTSuperficie[i][2] = FTIniParedExt - 4.;
		}

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FSUMTPPromedio[j] = new double**[FCanal[j][0]->getNin()];
			FTPared[j] = new double*[FCanal[j][0]->getNin()];
			FTParedAnt[j] = new double*[FCanal[j][0]->getNin()];
			FResistAxialAnt[j] = new double[FCanal[j][0]->getNin()];
			FResistAxialPost[j] = new double[FCanal[j][0]->getNin()];
			FResistEntreHacesAnt[j] = new double[FCanal[j][0]->getNin()];
			FResistEntreHacesPost[j] = new double[FCanal[j][0]->getNin()];
			FCapEntrada[j] = new double[FCanal[j][0]->getNin()];
			FCapSalida[j] = new double[FCanal[j][0]->getNin()];
			FCapPared[j] = new double[FCanal[j][0]->getNin()];
			FResistRadial[j] = new double*[FCanal[j][0]->getNin()];
			FResistConv[j] = new double*[FCanal[j][0]->getNin()];
		}
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				FSUMTPPromedio[j][i] = new double*[3]; // Dimensionado de los nodos de pared
				FTPared[j][i] = new double[3]; // Dimensionado de los nodos de pared
				FTParedAnt[j][i] = new double[3]; // Dimensionado de los nodos de pared
				FResistRadial[j][i] = new double[2]; // Dimensionado para tipo de canal 0->Entrada 1->Salida
				FResistConv[j][i] = new double[2]; // Dimensionado para tipo de canal 0->Entrada 1->Salida
			}
		}
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				for(int k = 0; k < 3; k++) {
					FSUMTPPromedio[j][i][k] = new double[2];
				}
			}
		}

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FDiametroHazInt[j] = FDiametroInicialHaz[j];
			if(j + 1 == FNumeroHacesCanales) {
				FDiametroHazExt[j] = FDiametroFiltroEfect;
			} else {
				FDiametroHazExt[j] = FDiametroInicialHaz[j + 1];
			}
		}

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FSUMTime[j] = 0.;
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				for(int k = 0; k < 3; k++) {
					FTPared[j][i][k] = FTIniPared;
					FTParedAnt[j][i][k] = FTIniPared;
					for(int h = 0; h < 2; h++) {
						FSUMTPPromedio[j][i][k][h] = 0.;
					}
				}
			}
		}

		if(FTipoCalcTempPared != nmTempConstante) {
			if(FTipoRefrig == nmAgua)
				FTExt = __units::degCToK(FTRefrigerante);
			else
				FTExt = __units::degCToK(tamb);

		}

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FNumCanalesHaz[j] = floor(__cons::Pi * (pow(FDiametroHazExt[j], 2.) - pow(FDiametroHazInt[j],
													2.)) / 4. / pow(FCanal[j][0]->GetLadoCanal(0) + FEspesorParedPorosa, 2.));
			FSuperficieTCHazAnt[j] = __cons::Pi * FDiametroHazInt[j];
			FSuperficieTCHazPost[j] = __cons::Pi * FDiametroHazExt[j];
			FAreaInternaHaz[j] = 4 * FCanal[j][0]->GetLadoCanal(0) * FNumCanalesHaz[j];
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::IniciaVariablesTransmisionCalor en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoTransmisionCalor(TBloqueMotor **Motor, double theta, TTubo **Tubo, TConcentrico **Concentrico) {
	try {

		double GastoO2Entrada, GastoO2Salida, GastoNO2Entrada, GastoNO2Salida, O2Salida, N2Salida, DensidadSalida;

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			// Coeficientes de pel�cula del canal de entrada y del canal de salida
			FCanal[j][0]->CalculaCoeficientePeliculaInterior(); // Canal de entrada
			FCanal[j][1]->CalculaCoeficientePeliculaInterior(); // Canal de salida
		}

		if(FTipoCalcTempPared != nmTempConstante && FCoefAjusTC != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				if(FCalculoRegeneracion) {
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						if(FCanal[j][0]->getNodoInicialFuente() > i) {
							FQreg[j][i] = 0.;
						} else {
							if(FCanal[j][0]->getNodoInicialFuente() == 0) {
								DensidadSalida = FCanal[j][1]->GetDensidad(i);
							} else {
								if(i == FCanal[j][0]->getNin() - 1) {
									DensidadSalida = FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente())
													 - (FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente() - 1) - FCanal[j][1]->GetDensidad(
															i - FCanal[j][0]->getNodoInicialFuente()))
													 / FCanal[j][1]->getXRef() * FCanal[j][1]->getDistanciaInterpolacion();
								} else {
									DensidadSalida = Interpola(FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente()),
															   FCanal[j][1]->GetDensidad(i - FCanal[j][0]->getNodoInicialFuente() + 1), 1.,
															   FCanal[j][1]->getDistanciaInterpolacion() / FCanal[j][1]->getXRef());
								}
							}
							GastoO2Entrada = FAreaVCCanalEntrada[j][i] * FVelocidadPared[j][0][i] * FCanal[j][0]->GetDensidad(
												 i) * FCanal[j][0]->GetYespecie(i, 0);
							GastoO2Salida = FAreaVCCanalSalida[j][i] * FVelocidadPared[j][1][i] * DensidadSalida *
											FFraccionMasicaEspecieSalida[j][i][0];
							GastoNO2Entrada = FAreaVCCanalEntrada[j][i] * FVelocidadPared[j][0][i] * FCanal[j][0]->GetDensidad(
												  i) * FFraccionMasicaNO2Entrada[j][i];
							GastoNO2Salida = FAreaVCCanalSalida[j][i] * FVelocidadPared[j][1][i] * DensidadSalida * FFraccionMasicaNO2Salida[j][i];
							FQreg[j][i] = FIncrH_reg1 / (__PM::O2 * FIndiceCompletitud1) * (GastoO2Entrada - GastoO2Salida)
										  + FIncrH_reg2 / (__PM::NO2 * FIndiceCompletitud2) * (GastoNO2Entrada - GastoNO2Salida); // (J/s)

						}
					}
				}
				// C�lculo de la temperatura de pared
				if(Motor != NULL) {
#ifdef ConcentricElement
					CalculaTemperaturaPared(Motor, theta, j, Tubo, Concentrico);
#else
					CalculaTemperaturaPared(Motor, theta, j, Tubo, NULL);
#endif
				} else {
					FCicloActual = FTime1DPF / FDuracionCiclo;
#ifdef ConcentricElement
					CalculaTemperaturaParedSinMotor(j, Tubo, Concentrico);
#else
					CalculaTemperaturaParedSinMotor(j, Tubo, NULL);
#endif
				}
			}
			if(Motor != NULL) {
				if(FCicloDPF != Motor[0]->getCiclo()) {
					FCicloDPF = Motor[0]->getCiclo();
				}
			} else {
				if(FCicloDPF != FCicloActual) {
					FCicloDPF = FCicloActual;
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoTransmisionCalor en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoResistenciaTC(int j, TTubo **Tubo, TConcentrico **Concentrico) {
	try {

		double Rrad = 0.;

//for(int j=0;j<FNumeroHacesCanales;j++){
		if(FTipoCalcTempPared != nmTempConstante && FCoefAjusTC != 0) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				// C�lculo de las resistencias en nodos iniciales y nodo final
				if(i == 0) {      // XRef vale la mitad
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					if(FCanal[j][0]->Geth(i) == 0.) {
						FResistConv[j][i][0] = 100000000.;
					} else {
						FResistConv[j][i][0] = 1 / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
													FCanal[j][0]->Geth(i));
					}
					if(FCanal[j][1]->Geth(i) == 0.) {
						FResistConv[j][i][1] = 100000000.;
					} else {
						FResistConv[j][i][1] = 1 / (2 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));
					}

					// C�lculo de las resistencias exteriores
					if(j + 1 == FNumeroHacesCanales) {       // Haz exterior
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));

						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_2_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 2.
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
					}

					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef() / 2.
										+ 2 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();

				} else if(i + 1 == FCanal[j][0]->getNin()) {     // XRef vale la mitad
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					if(FCanal[j][0]->Geth(i) == 0.) {
						FResistConv[j][i][0] = 100000000.;
					} else {
						FResistConv[j][i][0] = 1 / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
													FCanal[j][0]->Geth(i));
					}
					if(FCanal[j][1]->Geth(i) == 0.) {
						FResistConv[j][i][1] = 100000000.;
					} else {
						FResistConv[j][i][1] = 1 / (2 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));
					}

					// C�lculo de las resistencias a haces de canales adyacentes
					if(j + 1 == FNumeroHacesCanales) {
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));

						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_2_Pi(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 2.
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
					}

					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef() / 2.
										+ 2 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();
				} else {
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (8 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					FResistConv[j][i][0] = 1 / (4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
												FCanal[j][0]->Geth(i));
					FResistConv[j][i][1] = 1 / (4 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));

					// C�lculo de las resistencias exteriores
					if(j + 1 == FNumeroHacesCanales) {       // Haz externp
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef()) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_1_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef()) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 1.
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
					}
					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef()
										+ 4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();
				}
			}
		}
//}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoResistenciaTC en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoResistenciaTC_First_Time(int j, TTubo **Tubo, TConcentrico **Concentrico) {
	try {

		double Rrad = 0.;

//for(int j=0;j<FNumeroHacesCanales;j++){
		if(FTipoCalcTempPared != nmTempConstante && FCoefAjusTC != 0) {
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				// C�lculo de las resistencias en nodos iniciales y nodo final
				if(i == 0) {      // XRef vale la mitad
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());
					FResistRadial[j][i][1] = FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared *
											 FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					if(FCanal[j][0]->Geth(i) == 0.) {
						FResistConv[j][i][0] = 100000000.;
					} else {
						FResistConv[j][i][0] = 1 / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
													FCanal[j][0]->Geth(i));
					}
					if(FCanal[j][1]->Geth(i) == 0.) {
						FResistConv[j][i][1] = 100000000.;
					} else {
						FResistConv[j][i][1] = 1 / (2 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));
					}

					// C�lculo de las resistencias axiales
					FResistAxialAnt[j][i] = 10000000.;
					FResistAxialPost[j][i] = FCanal[j][0]->getXRef() / (4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa *
											 FConductividadPared); // En este caso XRef no vale la mitad

					// C�lculo de las resistencias a haces de canales adyacentes
					if(j == 0) {     // Haz interior
						FResistEntreHacesAnt[j][i] = 10000000.;
						if(j + 1 == FNumeroHacesCanales) {
							FResistEntreHacesPost[j][i] = FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect / 0.0001)
														  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
							FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect) /
														   (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAislante);
							FResistRadiacionAire[i] = (1 / FEmisividadIntAire
													   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
													   (1 / FEmisividadExtAire - 1))
													  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
													  (FTSuperficie[i][1] - FTSuperficie[i][2])
													  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
							FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
														   (FDiametroFiltroEfect + 2 * FEspesorAislante))
													   / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAire);
							FResistConduccionMetal[i] = log(
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
														/ (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadMetal);

							if(FHayConcentrico) {
#ifdef ConcetricElement
								if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
									FRg_int_ext[i] = __cons::_2_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
													 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
								} else {
									FRg_int_ext[i] = 100000000.;
								}
								if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
									FR_int_radiacion[i] = 100000000.;
								} else {
									FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
														   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
																   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
														   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
																   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
														  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
														  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
								}
#endif
							} else {
								FResistConveccionExt[i] = 2.
														  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
															 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
							}
							//FResistAxialMetalExtAnt[i]=10000000.;
							FResistAxialMetalExtAnt[i] = FAjustRAxAnt
														 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															FEspesorMetal);
							FResistAxialMetalExtPost[i] = FCanal[j][0]->getXRef()
														  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															 FEspesorMetal);
						} else {
							FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
															  FDiametroHazExt[j + 1] / 0.0001)
														  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						}
					} else if(j + 1 == FNumeroHacesCanales) {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect /
													  FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
						FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect) /
													   (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAislante);
						FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
													   (FDiametroFiltroEfect + 2 * FEspesorAislante))
												   / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAire);
						FResistConduccionMetal[i] = log(
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
													/ (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadMetal);

						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_2_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 2.
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
						//FResistAxialMetalExtAnt[i]=10000000.;
						FResistAxialMetalExtAnt[i] = FAjustRAxAnt
													 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														FEspesorMetal);
						FResistAxialMetalExtPost[i] = FCanal[j][0]->getXRef()
													  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														 FEspesorMetal);

					} else if(j - 1 == 0) {
						FResistEntreHacesAnt[j][i] = 1 * (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / 0.0001)
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					} else {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					}

					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef() / 2.
										+ 2 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();
					FCapPared[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									  FCanal[j][0]->getXRef();
					FCapSalida[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									   FCanal[j][0]->getXRef() / 2.;
					if(j + 1 == FNumeroHacesCanales) {
						FCapNodoInteriorSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadAislante * FCalorEspAislante
														* (pow(FDiametroFiltroEfect + FEspesorAislante, 2.) - pow(FDiametroFiltroEfect, 2.));
						FCapNodoMedioSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadAislante * FCalorEspAislante
													 * (pow(FDiametroFiltroEfect + 2 * FEspesorAislante, 2.) - pow(FDiametroFiltroEfect + FEspesorAislante, 2.));
						FCapNodoExteriorSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadMetal * FCalorEspMetal
														* (pow(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal,
															   2.) - pow(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire, 2.));
					}
				} else if(i + 1 == FCanal[j][0]->getNin()) {     // XRef vale la mitad
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());
					FResistRadial[j][i][1] = FEspesorParedPorosa / (4 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared *
											 FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					if(FCanal[j][0]->Geth(i) == 0.) {
						FResistConv[j][i][0] = 100000000.;
					} else {
						FResistConv[j][i][0] = 1 / (2 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
													FCanal[j][0]->Geth(i));
					}
					if(FCanal[j][1]->Geth(i) == 0.) {
						FResistConv[j][i][1] = 100000000.;
					} else {
						FResistConv[j][i][1] = 1 / (2 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));
					}

					// C�lculo de las resistencias axiales
					FResistAxialAnt[j][i] = FCanal[j][0]->getXRef() / (4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa *
											FConductividadPared); // En este caso XRef no vale la mitad
					FResistAxialPost[j][i] = 10000000.;

					// C�lculo de las resistencias a haces de canales adyacentes
					if(j == 0) {     // Haz interior
						FResistEntreHacesAnt[j][i] = 10000000.;
						if(j + 1 == FNumeroHacesCanales) {
							FResistEntreHacesPost[j][i] = 1 * FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect / 0.0001)
														  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
							FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect) /
														   (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAislante);
							FResistRadiacionAire[i] = (1 / FEmisividadIntAire
													   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
													   (1 / FEmisividadExtAire - 1))
													  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
													  (FTSuperficie[i][1] - FTSuperficie[i][2])
													  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
							FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
														   (FDiametroFiltroEfect + 2 * FEspesorAislante))
													   / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAire);
							FResistConduccionMetal[i] = log(
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
														/ (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadMetal);

							if(FHayConcentrico) {
#ifdef ConcetricElement
								if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
									FRg_int_ext[i] = __cons::_2_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
													 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
								} else {
									FRg_int_ext[i] = 100000000.;
								}
								if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
									FR_int_radiacion[i] = 100000000.;
								} else {
									FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
														   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
																   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
														   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
																   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
														  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
														  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
								}
#endif
							} else {
								FResistConveccionExt[i] = 2.
														  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
															 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
							}
							FResistAxialMetalExtAnt[i] = FCanal[j][0]->getXRef()
														 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															FEspesorMetal);
							//FResistAxialMetalExtPost[i]=10000000.;
							FResistAxialMetalExtPost[i] = FAjustRAxPos
														  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															 FEspesorMetal);
						} else {
							FResistEntreHacesPost[j][i] = 1 * (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
															  FDiametroHazExt[j + 1] / 0.0001)
														  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						}
					} else if(j + 1 == FNumeroHacesCanales) {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect /
													  FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect) /
													   (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAislante);
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef() / 2.) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
						FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
													   (FDiametroFiltroEfect + 2 * FEspesorAislante))
												   / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadAire);
						FResistConduccionMetal[i] = log(
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
													/ (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadMetal);

						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_2_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef() / 2.) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 2.
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
						FResistAxialMetalExtAnt[i] = FCanal[j][0]->getXRef()
													 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														FEspesorMetal);
						FResistAxialMetalExtPost[i] = FAjustRAxPos
													  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														 FEspesorMetal);
					} else if(j - 1 == 0) {
						FResistEntreHacesAnt[j][i] = 1 * (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / 0.0001)
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					} else {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					}

					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef() / 2.
										+ 2 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();
					FCapPared[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									  FCanal[j][0]->getXRef();
					FCapSalida[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									   FCanal[j][0]->getXRef() / 2.;
					if(j + 1 == FNumeroHacesCanales) {
						FCapNodoInteriorSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadAislante * FCalorEspAislante
														* (pow(FDiametroFiltroEfect + FEspesorAislante, 2.) - pow(FDiametroFiltroEfect, 2.));
						FCapNodoMedioSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadAislante * FCalorEspAislante
													 * (pow(FDiametroFiltroEfect + 2 * FEspesorAislante, 2.) - pow(FDiametroFiltroEfect + FEspesorAislante, 2.));
						FCapNodoExteriorSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi / 8. * FDensidadMetal * FCalorEspMetal
														* (pow(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal,
															   2.) - pow(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire, 2.));
					}
				} else {
					// C�lculo de las resistencias t�rmicas radiales dentro de la pared porosa en nodos intermedios
					FResistRadial[j][i][0] = FEspesorSoot[j][i] / (4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) *
											 FConductividadSoot * FCanal[j][0]->getXRef())
											 + FEspesorParedPorosa / (8 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared * FCanal[j][0]->getXRef());
					FResistRadial[j][i][1] = FEspesorParedPorosa / (8 * FCanal[j][0]->GetLadoCanal(i) * FConductividadPared *
											 FCanal[j][0]->getXRef());

					// C�lculo de las resistencias t�rmicas radiales debidas a convecci�n
					FResistConv[j][i][0] = 1 / (4 * (FCanal[j][0]->GetLadoCanal(i) - 2 * FEspesorSoot[j][i]) * FCanal[j][0]->getXRef() *
												FCanal[j][0]->Geth(i));
					FResistConv[j][i][1] = 1 / (4 * FCanal[j][0]->GetLadoCanal(i) * FCanal[j][0]->getXRef() * FCanal[j][1]->Geth(i));

					// C�lculo de las resistencias axiales
					FResistAxialAnt[j][i] = FCanal[j][0]->getXRef() / (4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa *
											FConductividadPared);
					FResistAxialPost[j][i] = FCanal[j][0]->getXRef() / (4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa *
											 FConductividadPared);

					// C�lculo de las resistencias a haces de canales adyacentes
					if(j == 0) {     // Haz interior
						FResistEntreHacesAnt[j][i] = 10000000.;
						if(j + 1 == FNumeroHacesCanales) {
							FResistEntreHacesPost[j][i] = FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect / 0.0001)
														  / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
							FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect)
														   / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadAislante);
							FResistRadiacionAire[i] = (1 / FEmisividadIntAire
													   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
													   (1 / FEmisividadExtAire - 1))
													  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef()) *
													  (FTSuperficie[i][1] - FTSuperficie[i][2])
													  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
							FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
														   (FDiametroFiltroEfect + 2 * FEspesorAislante))
													   / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadAire);
							FResistConduccionMetal[i] = log(
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
															(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
														/ (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadMetal);

							if(FHayConcentrico) {
#ifdef ConcetricElement
								if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
									FRg_int_ext[i] = __cons::_1_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
													 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
								} else {
									FRg_int_ext[i] = 100000000.;
								}
								if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
									FR_int_radiacion[i] = 100000000.;
								} else {
									FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
														   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
																   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
														   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
																   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef()) *
														  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
														  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
								}
#endif
							} else {
								FResistConveccionExt[i] = 1
														  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
															 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
							}

							FResistAxialMetalExtAnt[i] = FCanal[j][0]->getXRef()
														 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															FEspesorMetal);
							FResistAxialMetalExtPost[i] = FCanal[j][0]->getXRef()
														  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
															 FEspesorMetal);
						} else {
							FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
															  FDiametroHazExt[j + 1] / 0.0001)
														  / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						}

					} else if(j + 1 == FNumeroHacesCanales) {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = FAreaInternaHaz[j] / FSuperficieTCHazPost[j] * log(FDiametroFiltroEfect /
													  FDiametroHazInt[j])
													  / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistConduccionAislante[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante) / FDiametroFiltroEfect) /
													   (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadAislante);
						FResistRadiacionAire[i] = (1 / FEmisividadIntAire
												   + (FDiametroFiltroEfect + 2 * FEspesorAislante) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) *
												   (1 / FEmisividadExtAire - 1))
												  / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante) * FCanal[j][0]->getXRef()) *
												  (FTSuperficie[i][1] - FTSuperficie[i][2])
												  / (pow(__units::degCToK(FTSuperficie[i][1]), 4) - pow(__units::degCToK(FTSuperficie[i][2]), 4));
						FResistConduccionAire[i] = log((FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire) /
													   (FDiametroFiltroEfect + 2 * FEspesorAislante))
												   / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadAire);
						FResistConduccionMetal[i] = log(
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
														(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire))
													/ (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadMetal);

						if(FHayConcentrico) {
#ifdef ConcetricElement
							if(Tubo[FNumTuboExt - 1]->Gethi(i) != 0.) {
								FRg_int_ext[i] = __cons::_1_Pi / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) /
												 Tubo[FNumTuboExt - 1]->Gethi(i) / Tubo[FNumTuboExt - 1]->getXRef();
							} else {
								FRg_int_ext[i] = 100000000.;
							}
							if(Concentrico[FNumConcentrico]->GetTPared(0, 0, i) == FTSuperficie[i][2]) {
								FR_int_radiacion[i] = 100000000.;
							} else {
								FR_int_radiacion[i] = (1. / FEmisividad + (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 *
													   FEspesorMetal) / (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal + 2 *
															   Concentrico[FNumConcentrico]->GetEspesorFluido()) *
													   (1 / Tubo[FNumTuboExt - 1]->getEmisividad() - 1)) / (__cons::Sigma * __cons::Pi * (FDiametroFiltroEfect + 2 *
															   FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) * Tubo[FNumTuboExt - 1]->getXRef()) *
													  (FTSuperficie[i][2] - Concentrico[FNumConcentrico]->GetTPared(0, 0, i)) /
													  (pow(FTSuperficie[i][2], 4) - pow(Concentrico[FNumConcentrico]->GetTPared(0, 0, i), 4));
							}
#endif
						} else {
							FResistConveccionExt[i] = 1
													  / (__cons::Pi * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal) *
														 FCanal[j][0]->getXRef() * FCanal[j][0]->Gethext(i));
						}
						FResistAxialMetalExtAnt[i] = FCanal[j][0]->getXRef()
													 / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														FEspesorMetal);
						FResistAxialMetalExtPost[i] = FCanal[j][0]->getXRef()
													  / (__cons::Pi * FConductividadMetal * (FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + FEspesorMetal) *
														 FEspesorMetal);

					} else if(j - 1 == 0) {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / 0.0001)
													 / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					} else {
						FResistEntreHacesAnt[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j - 1]) / FSuperficieTCHazAnt[j] * log(
														 FDiametroHazExt[j] / FDiametroHazInt[j - 1])
													 / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
						FResistEntreHacesPost[j][i] = (FAreaInternaHaz[j] + FAreaInternaHaz[j + 1]) / FSuperficieTCHazPost[j] * log(
														  FDiametroHazExt[j + 1] / FDiametroHazInt[j])
													  / (__cons::Pi_x_2 * FCanal[j][0]->getXRef() * FConductividadParedRadial);
					}

					// C�lculo de las capacitancias t�rmicas
					FCapEntrada[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
										FCanal[j][0]->getXRef()
										+ 4 * FCanal[j][0]->GetLadoCanal(i) * FEspesorSoot[j][i] * FDensidadSootCapa_in * FCalorEspSoot *
										FCanal[j][0]->getXRef();
					FCapPared[j][i] = 2 * FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									  FCanal[j][0]->getXRef();
					FCapSalida[j][i] = FCanal[j][0]->GetLadoCanal(i) * FEspesorParedPorosa * FDensidadParedPorosa * FCalorEspPared *
									   FCanal[j][0]->getXRef();
					if(j + 1 == FNumeroHacesCanales) {
						FCapNodoInteriorSuperficie[i] = FCanal[j][0]->getXRef() * FDensidadAislante * FCalorEspAislante
														* __geom::Ring_area(FDiametroFiltroEfect, FDiametroFiltroEfect + FEspesorAislante);
						FCapNodoMedioSuperficie[i] = FCanal[j][0]->getXRef() * FDensidadAislante * FCalorEspAislante
													 * __geom::Ring_area(FDiametroFiltroEfect + FEspesorAislante, FDiametroFiltroEfect + 2 * FEspesorAislante);
						FCapNodoExteriorSuperficie[i] = FCanal[j][0]->getXRef() * __cons::Pi_4 * FDensidadMetal * FCalorEspMetal
														* __geom::Ring_area(FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire,
																FDiametroFiltroEfect + 2 * FEspesorAislante + 2 * FEspesorAire + 2 * FEspesorMetal);
					}
				}
			}
		}
//}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoResistenciaTC_First_Time en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculaTemperaturaPared(TBloqueMotor **Motor, double theta, int j, TTubo **Tubo,
								   TConcentrico **Concentrico) {
	double TgEntrada = 0., TgSalida = 0.;
//double zzz,czz,uq1,uq2;
	double DeltaTTPared = 0.;
	double Tpant0 = 0., Tpant1 = 0., Tpant2 = 0., Tpantant = 0., Tpantpos = 0., ErrorTp = 0.;
	double Tphazant, Tphazpost, Tpsuperficie_int_ant, Tpsuperficie_ext_ant, Tpsuperficie_med_ant, Tpsup_ant_ant,
		   Tpsup_ant_post;
	double TgCon = 0., R_Equiv_Gap = 0.;
	bool EsPrimeraVez;

	try {

		DeltaTTPared = FTime1DPF - FTime0DPF;

		for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
			FTParedAnt[j][i][0] = FTPared[j][i][0];
			FTParedAnt[j][i][1] = FTPared[j][i][1];
			FTParedAnt[j][i][2] = FTPared[j][i][2];
		}

		if(theta > FAnguloTotalCiclo) {
			FSUMTime[j] += DeltaTTPared;
		}
		for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
			//Establece la temperatura del gas en el interior del conducto.
			TgEntrada = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
							2.) / (FCanal[j][0]->GetGamma(i) * FCanal[j][0]->GetRMezcla(i));
			TgSalida = pow(FCanal[j][1]->GetAsonido(i) * __cons::ARef,
						   2.) / (FCanal[j][1]->GetGamma(i) * FCanal[j][1]->GetRMezcla(i));

			//Establece la temperatura de la pared en el instante de c�lculo anterior...
			//...en los nodos anterior y posterior.
			if(i == 0) {
				Tpantant = __units::degCToK(FTParedAnt[j][i][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i + 1][1]);
			} else if(i == FCanal[j][0]->getNin() - 1) {
				Tpantant = __units::degCToK(FTParedAnt[j][i - 1][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i][1]);
			} else {
				Tpantant = __units::degCToK(FTParedAnt[j][i - 1][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i + 1][1]);
			}

			// Temperatura en haces posterior y anterior en el instante anterior
			if(j == 0 && j == FNumeroHacesCanales - 1) {
				Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
				Tphazpost = __units::degCToK(FTSuperficie[i][0]);
			} else if(j == 0) {
				Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
				Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
			} else if(j == FNumeroHacesCanales - 1) {
				Tphazant = __units::degCToK(FTParedAnt[j - 1][i][1]);
				Tphazpost = __units::degCToK(FTSuperficie[i][0]);
			} else {
				Tphazant = __units::degCToK(FTParedAnt[j - 1][i][1]);
				Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
			}

			//...en los nodos interior, medio y exterior.
			Tpant0 = __units::degCToK(FTParedAnt[j][i][0]);
			Tpant1 = __units::degCToK(FTParedAnt[j][i][1]);
			Tpant2 = __units::degCToK(FTParedAnt[j][i][2]);

			//C�lculo de las temperaturas de pared.
			// Verificar riesgos de resistencia igual a 0!!!!!!!!!!!!!!!!!!
			FTPared[j][i][2] = DeltaTTPared / FCapEntrada[j][i] * (1 / FResistRadial[j][i][0] *
							   (Tpant1 - Tpant2) + 1 / FResistConv[j][i][0] * (TgEntrada - Tpant2) + FQreg[j][i]) + Tpant2;
			FTPared[j][i][1] = DeltaTTPared / FCapPared[j][i]
							   * (1 / FResistRadial[j][i][1] * (Tpant0 - Tpant1) + 1 / FResistRadial[j][i][0] * (Tpant2 - Tpant1) + 1 /
								  FResistAxialAnt[j][i] * (Tpantant - Tpant1)
								  + 1 / FResistAxialPost[j][i] * (Tpantpos - Tpant1) + 1 / FResistEntreHacesAnt[j][i] *
								  (Tphazant - Tpant1) + 1 / FResistEntreHacesPost[j][i] * (Tphazpost - Tpant1)) + Tpant1;
			FTPared[j][i][0] = DeltaTTPared / FCapSalida[j][i] * (1 / FResistConv[j][i][1] * (TgSalida - Tpant0) + 1 /
							   FResistRadial[j][i][1] * (Tpant1 - Tpant0)) + Tpant0;

			if(j + 1 == FNumeroHacesCanales) {
				Tpsuperficie_int_ant = __units::degCToK(FTSuperficie[i][0]);
				Tpsuperficie_med_ant = __units::degCToK(FTSuperficie[i][1]);
				Tpsuperficie_ext_ant = __units::degCToK(FTSuperficie[i][2]);
				if(i == 0) {
					Tpsup_ant_ant = __units::degCToK(FTuboEntradaDPF->GetTPTubo(1, FNodoTuboEntrada));
					Tpsup_ant_post = __units::degCToK(FTSuperficie[i + 1][2]);
				} else if(i == FCanal[j][0]->getNin() - 1) {
					Tpsup_ant_ant = __units::degCToK(FTSuperficie[i - 1][2]);
					Tpsup_ant_post = __units::degCToK(FTuboSalidaDPF->GetTPTubo(1, FNodoTuboSalida));
				} else {
					Tpsup_ant_ant = __units::degCToK(FTSuperficie[i + 1][2]);
					Tpsup_ant_post = __units::degCToK(FTSuperficie[i - 1][2]);
				}
				FTSuperficie[i][0] = DeltaTTPared / FCapNodoInteriorSuperficie[i]
									 * (1 / FResistConduccionAislante[i] * (Tpsuperficie_med_ant - Tpsuperficie_int_ant) + 1 / FResistEntreHacesPost[j][i]
										* (Tpant1 - Tpsuperficie_int_ant)) + Tpsuperficie_int_ant;
				FTSuperficie[i][1] = DeltaTTPared / FCapNodoMedioSuperficie[i]
									 * (1 / FResistConduccionAislante[i] * (Tpsuperficie_int_ant - Tpsuperficie_med_ant)
										+ 1 / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i]) *
										(Tpsuperficie_ext_ant - Tpsuperficie_med_ant)) + Tpsuperficie_med_ant;

				if(FHayConcentrico) {
#ifdef ConcetricElement
					TgCon = pow(Tubo[FNumTuboExt - 1]->GetAsonido(i) * __cons::ARef,
								2.) / (Tubo[FNumTuboExt - 1]->GetGamma(i) * Tubo[FNumTuboExt - 1]->GetRMezcla(i));
					FTSuperficie[i][2] = DeltaTTPared / FCapNodoExteriorSuperficie[i] * (1 / (1 / (1 / FResistConduccionAire[i] + 1 /
										 FResistRadiacionAire[i]) + FResistConduccionMetal[i]) * (Tpsuperficie_med_ant - Tpsuperficie_ext_ant) +
										 1 / FRg_int_ext[i] * (TgCon - Tpsuperficie_ext_ant) + 1 / FR_int_radiacion[i] *
										 (Concentrico[FNumConcentrico]->GetTPared(0, 0, i) - Tpsuperficie_ext_ant) +
										 1 / FResistAxialMetalExtPost[i] * (Tpsup_ant_post - Tpsuperficie_ext_ant) +
										 1 / FResistAxialMetalExtAnt[i] * (Tpsup_ant_ant - Tpsuperficie_ext_ant)) + Tpsuperficie_ext_ant;
#endif
				} else {
					FTSuperficie[i][2] = DeltaTTPared / FCapNodoExteriorSuperficie[i]
										 * (1 / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i]) *
											(Tpsuperficie_med_ant - Tpsuperficie_ext_ant)
											+ 1 / FResistConveccionExt[i] * (FTExt - Tpsuperficie_ext_ant) + 1 / FResistAxialMetalExtPost[i] *
											(Tpsup_ant_post - Tpsuperficie_ext_ant)
											+ 1 / FResistAxialMetalExtAnt[i] * (Tpsup_ant_ant - Tpsuperficie_ext_ant)) + Tpsuperficie_ext_ant;
				}
				FTSuperficie[i][0] = __units::KTodegC(FTSuperficie[i][0]);
				FTSuperficie[i][1] = __units::KTodegC(FTSuperficie[i][1]);
				FTSuperficie[i][2] = __units::KTodegC(FTSuperficie[i][2]);
			}

			for(int k = 0; k < 3; k++) {
				FTPared[j][i][k] = __units::KTodegC(FTPared[j][i][k]);
			}

			// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
			if(FTipoCalcTempPared == nmVariableSinInerciaTermica
			   || theta / FAnguloTotalCiclo <= Motor[0]->getNumCiclosSinInerciaTermica()) {
				if(theta > FAnguloTotalCiclo) {
					if(FHayConcentrico) {
#ifdef ConcetricElement
						R_Equiv_Gap = (Concentrico[FNumConcentrico]->GetRg_ext_int(i) + FRg_int_ext[i]) * FR_int_radiacion[i] /
									  (Concentrico[FNumConcentrico]->GetRg_ext_int(i) + FRg_int_ext[i] + FR_int_radiacion[i]);
						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][2][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0, i) * DeltaTTPared / R_Equiv_Gap;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][2][i] += DeltaTTPared / R_Equiv_Gap;

						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][1][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0,
													   i) * DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
															   FResistConduccionMetal[i]);
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][1][i] += DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 /
													   FResistRadiacionAire[i]) + FResistConduccionMetal[i]);

						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][0][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0,
													   i) * DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
															   FResistConduccionMetal[i] + FResistConduccionAislante[i]);
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][0][i] += DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 /
													   FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistConduccionAislante[i]);
#endif
					}

					FSUMTPPromedio[j][i][2][0] += 1 / (FResistConv[j][i][0]) * TgEntrada * DeltaTTPared + FQreg[j][i] * DeltaTTPared;
					FSUMTPPromedio[j][i][2][1] += 1 / (FResistConv[j][i][0]) * DeltaTTPared;
					FSUMTPPromedio[j][i][1][0] += 1 / (FResistConv[j][i][0] + FResistRadial[j][i][0]) * TgEntrada * DeltaTTPared
												  + 1 / (FResistConv[j][i][1] + FResistRadial[j][i][1]) * TgSalida * DeltaTTPared;
					FSUMTPPromedio[j][i][1][1] += 1 / (FResistConv[j][i][0] + FResistRadial[j][i][0]) * DeltaTTPared + 1 /
												  (FResistConv[j][i][1] + FResistRadial[j][i][1]) * DeltaTTPared;
					FSUMTPPromedio[j][i][0][0] += 1 / (FResistConv[j][i][1]) * TgSalida * DeltaTTPared;
					FSUMTPPromedio[j][i][0][1] += 1 / (FResistConv[j][i][1]) * DeltaTTPared;
				}
			}

		}

		// Si est� al final del ciclo...
		if(FCicloDPF != Motor[0]->getCiclo() && FSUMTime[j] > 0.) {
			// ...si (el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando) y est� en el segundo ciclo, calcula la temperatura de convergencia
			if((FTipoCalcTempPared == nmVariableSinInerciaTermica
				|| theta / FAnguloTotalCiclo <= Motor[0]->getNumCiclosSinInerciaTermica()) && theta > FAnguloTotalCiclo + 1) {
				ErrorTp = 1;
				EsPrimeraVez = true;
				while(ErrorTp >=
					  1) {   //Itera hasta conseguir una diferencia entre las temperaturas de pared menor a 1�C entre pasos.
					ErrorTp = 0.;
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						//Establece la temperatura de la pared en el instante de c�lculo actual...
						//...en los nodos anterior y posterior.

						if(i == 0) {
							Tpantant = __units::degCToK(FTPared[j][i][1]);
							Tpantpos = __units::degCToK(FTPared[j][i + 1][1]);
						} else if(i == FCanal[j][0]->getNin() - 1) {
							Tpantant = __units::degCToK(FTPared[j][i - 1][1]);
							Tpantpos = __units::degCToK(FTPared[j][i][1]);
						} else {
							Tpantant = __units::degCToK(FTPared[j][i - 1][1]);
							Tpantpos = __units::degCToK(FTPared[j][i + 1][1]);
						}

						// Temperatura en haces posterior y anterior en el instante anterior
						if(j == 0 && j == FNumeroHacesCanales - 1) {
							Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
							Tphazpost = __units::degCToK(FTSuperficie[i][0]);
						} else if(j == 0) {
							Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
							Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
						} else if(j == FNumeroHacesCanales - 1) {
							Tphazant = __units::degCToK(FTPared[j - 1][i][1]);
							Tphazpost = __units::degCToK(FTSuperficie[i][0]);
						} else {
							Tphazant = __units::degCToK(FTPared[j - 1][i][1]);
							Tphazpost = __units::degCToK(FTPared[j + 1][i][1]);
						}

						//...en los nodos interior, medio y exterior.
						Tpant0 = __units::degCToK(FTPared[j][i][0]);
						Tpant1 = __units::degCToK(FTPared[j][i][1]);
						Tpant2 = __units::degCToK(FTPared[j][i][2]);

						if(EsPrimeraVez) {
							FTPared[j][i][1] = FSUMTPPromedio[j][i][1][0] / FSUMTPPromedio[j][i][1][1];

						} else {
							FTPared[j][i][1] = (FSUMTime[j] * Tpantant / FResistAxialAnt[j][i] + FSUMTime[j] * Tpantpos / FResistAxialPost[j][i] +
												FSUMTime[j] * Tphazant / FResistEntreHacesAnt[j][i]
												+ FSUMTime[j] * Tphazpost / FResistEntreHacesPost[j][i] + FSUMTPPromedio[j][i][1][0])
											   / (FSUMTime[j] / FResistAxialAnt[j][i] + FSUMTime[j] / FResistAxialPost[j][i] + FSUMTime[j] / FResistEntreHacesAnt[j][i]
												  + FSUMTime[j] / FResistEntreHacesPost[j][i]
												  + FSUMTPPromedio[j][i][1][1]);
						}
						FTPared[j][i][0] = (FSUMTime[j] * Tpant1 / FResistRadial[j][i][1] + FSUMTPPromedio[j][i][0][0]) /
										   (FSUMTime[j] / FResistRadial[j][i][1] + FSUMTPPromedio[j][i][0][1]);
						FTPared[j][i][2] = (FSUMTime[j] * Tpant1 / FResistRadial[j][i][0] + FSUMTPPromedio[j][i][2][0]) /
										   (FSUMTime[j] / FResistRadial[j][i][0] + FSUMTPPromedio[j][i][2][1]);

						if(ErrorTp < fabs(Tpant1 - FTPared[j][i][1])) {
							ErrorTp = fabs(Tpant1 - FTPared[j][i][1]);
						}

						if(j + 1 == FNumeroHacesCanales) {
							Tpsuperficie_int_ant = __units::degCToK(FTSuperficie[i][0]);
							Tpsuperficie_med_ant = __units::degCToK(FTSuperficie[i][1]);
							Tpsuperficie_ext_ant = __units::degCToK(FTSuperficie[i][2]);
							if(i == 0) {
								Tpsup_ant_ant = __units::degCToK(FTuboEntradaDPF->GetTPTubo(1, FNodoTuboEntrada));
								Tpsup_ant_post = __units::degCToK(FTSuperficie[i + 1][2]);
							} else if(i == FCanal[j][0]->getNin() - 1) {
								Tpsup_ant_ant = __units::degCToK(FTSuperficie[i - 1][2]);
								Tpsup_ant_post = __units::degCToK(FTuboSalidaDPF->GetTPTubo(1, FNodoTuboSalida));
							} else {
								Tpsup_ant_ant = __units::degCToK(FTSuperficie[i + 1][2]);
								Tpsup_ant_post = __units::degCToK(FTSuperficie[i - 1][2]);
							}

							if(FHayConcentrico) {
#ifdef ConcetricElement
								FTSuperficie[i][0] = (FSUMTime[j] * Tpant1 / FResistEntreHacesPost[j][i] + FSUMTPPromedioConc[0][0][i]) /
													 (FSUMTime[j] / FResistEntreHacesPost[j][i] + FSUMTPPromedioConc[1][0][i]);
								FTSuperficie[i][1] = (FSUMTime[j] * Tpant1 / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i]) +
													  FSUMTPPromedioConc[0][1][i]) /
													 (FSUMTime[j] / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i]) + FSUMTPPromedioConc[1][1][i]);
								FTSuperficie[i][2] = (FSUMTPPromedioConc[0][2][i] + FSUMTime[j] * Tpant1 / (FResistConduccionAislante[i] + 1 /
													  (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistEntreHacesPost[j][i])
													  +
													  FSUMTime[j] * Tpsup_ant_ant / FResistAxialMetalExtAnt[i] + FSUMTime[j] * Tpsup_ant_post / FResistAxialMetalExtPost[i]) /
													 (FSUMTPPromedioConc[1][2][i] + FSUMTime[j] / (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 /
															 FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistEntreHacesPost[j][i]) +
													  FSUMTime[j] / FResistAxialMetalExtAnt[i] + FSUMTime[j] / FResistAxialMetalExtPost[i]);
#endif
							} else {
								FTSuperficie[i][0] = (FSUMTime[j] * Tpant1 / FResistEntreHacesPost[j][i]
													  + FSUMTime[j] * FTExt
													  / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
														 FResistConveccionExt[i] + FResistConduccionAislante[i]))
													 / (FSUMTime[j] / FResistEntreHacesPost[j][i]
														+ FSUMTime[j]
														/ (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
														   FResistConveccionExt[i] + FResistConduccionAislante[i]));
								FTSuperficie[i][1] = (FSUMTime[j] * Tpant1 / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i])
													  + FSUMTime[j] * FTExt / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
															  FResistConveccionExt[i]))
													 / (FSUMTime[j] / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i])
														+ FSUMTime[j] / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
																FResistConveccionExt[i]));
								FTSuperficie[i][2] = (FSUMTime[j] * FTExt / FResistConveccionExt[i]
													  + FSUMTime[j] * Tpant1
													  / (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
														 FResistConduccionMetal[i] + FResistEntreHacesPost[j][i])
													  + FSUMTime[j] * Tpsup_ant_ant / FResistAxialMetalExtAnt[i] + FSUMTime[j] * Tpsup_ant_post / FResistAxialMetalExtPost[i])
													 / (FSUMTime[j] / FResistConveccionExt[i]
														+ FSUMTime[j]
														/ (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
														   FResistConduccionMetal[i] + FResistConduccionMetal[i]
														   + FResistEntreHacesPost[j][i]) + FSUMTime[j] / FResistAxialMetalExtAnt[i] + FSUMTime[j] / FResistAxialMetalExtPost[i]);
							}
							if(ErrorTp < fabs(Tpsuperficie_ext_ant - FTSuperficie[i][2])) {
								ErrorTp = fabs(Tpsuperficie_ext_ant - FTSuperficie[i][2]);
							}

							FTSuperficie[i][0] = __units::KTodegC(FTSuperficie[i][0]);
							FTSuperficie[i][1] = __units::KTodegC(FTSuperficie[i][1]);
							FTSuperficie[i][2] = __units::KTodegC(FTSuperficie[i][2]);
						}

						for(int k = 0; k < 3; k++) {
							FTPared[j][i][k] = __units::KTodegC(FTPared[j][i][k]);
						}

					}
					EsPrimeraVez = false;
				}
			}
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				for(int k = 0; k < 3; k++) {
					for(int h = 0; h < 2; h++) {
						FSUMTPPromedio[j][i][k][h] = 0.;
					}
				}
			}
			FSUMTime[j] = 0.;
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculaTemperaturaPared en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculaTemperaturaParedSinMotor(int j, TTubo **Tubo, TConcentrico **Concentrico) {
	double TgEntrada = 0., TgSalida = 0.;
	double DeltaTTPared = 0.;
	double Tpant0 = 0., Tpant1 = 0., Tpant2 = 0., Tpantant = 0., Tpantpos = 0., ErrorTp = 0.;
	double Tphazant, Tphazpost, Tpsuperficie_int_ant, Tpsuperficie_ext_ant, Tpsuperficie_med_ant, Tpsup_ant_ant,
		   Tpsup_ant_post;
	double TgCon = 0., R_Equiv_Gap = 0.;
	bool EsPrimeraVez;

	try {

		DeltaTTPared = FTime1DPF - FTime0DPF;

		for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
			FTParedAnt[j][i][0] = FTPared[j][i][0];
			FTParedAnt[j][i][1] = FTPared[j][i][1];
			FTParedAnt[j][i][2] = FTPared[j][i][2];
		}

		if(FTime1DPF > FDuracionCiclo) {
			FSUMTime[j] += DeltaTTPared;
		}
		for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
			//Establece la temperatura del gas en el interior del conducto.
			TgEntrada = pow(FCanal[j][0]->GetAsonido(i) * __cons::ARef,
							2.) / (FCanal[j][0]->GetGamma(i) * FCanal[j][0]->GetRMezcla(i));
			TgSalida = pow(FCanal[j][1]->GetAsonido(i) * __cons::ARef,
						   2.) / (FCanal[j][1]->GetGamma(i) * FCanal[j][1]->GetRMezcla(i));

			//Establece la temperatura de la pared en el instante de c�lculo anterior...
			//...en los nodos anterior y posterior.
			if(i == 0) {
				Tpantant = __units::degCToK(FTParedAnt[j][i][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i + 1][1]);
			} else if(i == FCanal[j][0]->getNin() - 1) {
				Tpantant = __units::degCToK(FTParedAnt[j][i - 1][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i][1]);
			} else {
				Tpantant = __units::degCToK(FTParedAnt[j][i - 1][1]);
				Tpantpos = __units::degCToK(FTParedAnt[j][i + 1][1]);
			}

			// Temperatura en haces posterior y anterior en el instante anterior
			if(j == 0 && j == FNumeroHacesCanales - 1) {
				Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
				Tphazpost = __units::degCToK(FTSuperficie[i][0]);
			} else if(j == 0) {
				Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
				Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
			} else if(j == FNumeroHacesCanales - 1) {
				Tphazant = __units::degCToK(FTParedAnt[j - 1][i][1]);
				Tphazpost = __units::degCToK(FTSuperficie[i][0]);
			} else {
				Tphazant = __units::degCToK(FTParedAnt[j - 1][i][1]);
				Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
			}

			//...en los nodos interior, medio y exterior.
			Tpant0 = __units::degCToK(FTParedAnt[j][i][0]);
			Tpant1 = __units::degCToK(FTParedAnt[j][i][1]);
			Tpant2 = __units::degCToK(FTParedAnt[j][i][2]);

			//C�lculo de las temperaturas de pared.
			// Verificar riesgos de resistencia igual a 0!!!!!!!!!!!!!!!!!!
			FTPared[j][i][2] = DeltaTTPared / FCapEntrada[j][i] * (1 / FResistRadial[j][i][0] *
							   (Tpant1 - Tpant2) + 1 / FResistConv[j][i][0] * (TgEntrada - Tpant2) + FQreg[j][i]) + Tpant2;
			FTPared[j][i][1] = DeltaTTPared / FCapPared[j][i]
							   * (1 / FResistRadial[j][i][1] * (Tpant0 - Tpant1) + 1 / FResistRadial[j][i][0] * (Tpant2 - Tpant1) + 1 /
								  FResistAxialAnt[j][i] * (Tpantant - Tpant1)
								  + 1 / FResistAxialPost[j][i] * (Tpantpos - Tpant1) + 1 / FResistEntreHacesAnt[j][i] *
								  (Tphazant - Tpant1) + 1 / FResistEntreHacesPost[j][i] * (Tphazpost - Tpant1)) + Tpant1;
			FTPared[j][i][0] = DeltaTTPared / FCapSalida[j][i] * (1 / FResistConv[j][i][1] * (TgSalida - Tpant0) + 1 /
							   FResistRadial[j][i][1] * (Tpant1 - Tpant0)) + Tpant0;

			if(j + 1 == FNumeroHacesCanales) {
				Tpsuperficie_int_ant = __units::degCToK(FTSuperficie[i][0]);
				Tpsuperficie_med_ant = __units::degCToK(FTSuperficie[i][1]);
				Tpsuperficie_ext_ant = __units::degCToK(FTSuperficie[i][2]);
				if(i == 0) {
					Tpsup_ant_ant = __units::degCToK(FTuboEntradaDPF->GetTPTubo(1, FNodoTuboEntrada));
					Tpsup_ant_post = __units::degCToK(FTSuperficie[i + 1][2]);
				} else if(i == FCanal[j][0]->getNin() - 1) {
					Tpsup_ant_ant = __units::degCToK(FTSuperficie[i - 1][2]);
					Tpsup_ant_post = __units::degCToK(FTuboSalidaDPF->GetTPTubo(1, FNodoTuboSalida));
				} else {
					Tpsup_ant_ant = __units::degCToK(FTSuperficie[i + 1][2]);
					Tpsup_ant_post = __units::degCToK(FTSuperficie[i - 1][2]);
				}
				FTSuperficie[i][0] = DeltaTTPared / FCapNodoInteriorSuperficie[i]
									 * (1 / FResistConduccionAislante[i] * (Tpsuperficie_med_ant - Tpsuperficie_int_ant) + 1 / FResistEntreHacesPost[j][i]
										* (Tpant1 - Tpsuperficie_int_ant)) + Tpsuperficie_int_ant;
				FTSuperficie[i][1] = DeltaTTPared / FCapNodoMedioSuperficie[i]
									 * (1 / FResistConduccionAislante[i] * (Tpsuperficie_int_ant - Tpsuperficie_med_ant)
										+ 1 / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i]) *
										(Tpsuperficie_ext_ant - Tpsuperficie_med_ant)) + Tpsuperficie_med_ant;
				if(FHayConcentrico) {
#ifdef ConcetricElement
					TgCon = pow(Tubo[FNumTuboExt - 1]->GetAsonido(i) * __cons::ARef,
								2.) / (Tubo[FNumTuboExt - 1]->GetGamma(i) * Tubo[FNumTuboExt - 1]->GetRMezcla(i));
					FTSuperficie[i][2] = DeltaTTPared / FCapNodoExteriorSuperficie[i] * (1 / (1 / (1 / FResistConduccionAire[i] + 1 /
										 FResistRadiacionAire[i]) + FResistConduccionMetal[i]) * (Tpsuperficie_med_ant - Tpsuperficie_ext_ant) +
										 1 / FRg_int_ext[i] * (TgCon - Tpsuperficie_ext_ant) + 1 / FR_int_radiacion[i] *
										 (Concentrico[FNumConcentrico]->GetTPared(0, 0, i) - Tpsuperficie_ext_ant) +
										 1 / FResistAxialMetalExtPost[i] * (Tpsup_ant_post - Tpsuperficie_ext_ant) +
										 1 / FResistAxialMetalExtAnt[i] * (Tpsup_ant_ant - Tpsuperficie_ext_ant)) + Tpsuperficie_ext_ant;
#endif
				} else {
					FTSuperficie[i][2] = DeltaTTPared / FCapNodoExteriorSuperficie[i]
										 * (1 / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i]) *
											(Tpsuperficie_med_ant - Tpsuperficie_ext_ant)
											+ 1 / FResistConveccionExt[i] * (FTExt - Tpsuperficie_ext_ant) + 1 / FResistAxialMetalExtPost[i] *
											(Tpsup_ant_post - Tpsuperficie_ext_ant)
											+ 1 / FResistAxialMetalExtAnt[i] * (Tpsup_ant_ant - Tpsuperficie_ext_ant)) + Tpsuperficie_ext_ant;
				}
				FTSuperficie[i][0] = __units::KTodegC(FTSuperficie[i][0]);
				FTSuperficie[i][1] = __units::KTodegC(FTSuperficie[i][1]);
				FTSuperficie[i][2] = __units::KTodegC(FTSuperficie[i][2]);
			}

			for(int k = 0; k < 3; k++) {
				FTPared[j][i][k] = __units::KTodegC(FTPared[j][i][k]);
			}

			// Si el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando...
			if(FTipoCalcTempPared == nmVariableSinInerciaTermica || FCicloActual <= FNumCiclosSinInerciaTermica) {
				if(FTime1DPF > FDuracionCiclo) {
					if(FHayConcentrico) {
#ifdef ConcetricElement
						R_Equiv_Gap = (Concentrico[FNumConcentrico]->GetRg_ext_int(i) + FRg_int_ext[i]) * FR_int_radiacion[i] /
									  (Concentrico[FNumConcentrico]->GetRg_ext_int(i) + FRg_int_ext[i] + FR_int_radiacion[i]);
						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][2][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0, i) * DeltaTTPared / R_Equiv_Gap;
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][2][i] += DeltaTTPared / R_Equiv_Gap;

						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][1][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0,
													   i) * DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
															   FResistConduccionMetal[i]);
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][1][i] += DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 /
													   FResistRadiacionAire[i]) + FResistConduccionMetal[i]);

						// Sumatorio de R*Twall_concentric*incrt (para el c�lculo del numerador de la integral).
						FSUMTPPromedioConc[0][0][i] += Concentrico[FNumConcentrico]->GetTPared(0, 0,
													   i) * DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
															   FResistConduccionMetal[i] + FResistConduccionAislante[i]);
						// Sumatorio de R*incrt (para el c�lculo del denominador de la integral).
						FSUMTPPromedioConc[1][0][i] += DeltaTTPared / (R_Equiv_Gap + 1 / (1 / FResistConduccionAire[i] + 1 /
													   FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistConduccionAislante[i]);
#endif
					}

					FSUMTPPromedio[j][i][2][0] += 1 / (FResistConv[j][i][0]) * TgEntrada * DeltaTTPared + FQreg[j][i] * DeltaTTPared;
					FSUMTPPromedio[j][i][2][1] += 1 / (FResistConv[j][i][0]) * DeltaTTPared;
					FSUMTPPromedio[j][i][1][0] += 1 / (FResistConv[j][i][0] + FResistRadial[j][i][0]) * TgEntrada * DeltaTTPared
												  + 1 / (FResistConv[j][i][1] + FResistRadial[j][i][1]) * TgSalida * DeltaTTPared;
					FSUMTPPromedio[j][i][1][1] += 1 / (FResistConv[j][i][0] + FResistRadial[j][i][0]) * DeltaTTPared + 1 /
												  (FResistConv[j][i][1] + FResistRadial[j][i][1]) * DeltaTTPared;
					FSUMTPPromedio[j][i][0][0] += 1 / (FResistConv[j][i][1]) * TgSalida * DeltaTTPared;
					FSUMTPPromedio[j][i][0][1] += 1 / (FResistConv[j][i][1]) * DeltaTTPared;
				}
			}

		}

		// Si est� al final del ciclo...
		if(FCicloDPF != FCicloActual && FSUMTime[j] > 0. && FCicloActual > 1) {
			// ...si (el tipo de calculo es sin inercia termica � lleva menos de "NumCiclosSinInerciaTermica" ciclos calculando) y est� en el segundo ciclo, calcula la temperatura de convergencia
			if((FTipoCalcTempPared == nmVariableSinInerciaTermica || FCicloActual <= FNumCiclosSinInerciaTermica)
			   && FTime1DPF > FDuracionCiclo) {
				ErrorTp = 1;
				EsPrimeraVez = true;
				while(ErrorTp >=
					  1) {   //Itera hasta conseguir una diferencia entre las temperaturas de pared menor a 1�C entre pasos.
					ErrorTp = 0.;
					for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
						//Establece la temperatura de la pared en el instante de c�lculo actual...
						//...en los nodos anterior y posterior.

						if(i == 0) {
							Tpantant = __units::degCToK(FTPared[j][i][1]);
							Tpantpos = __units::degCToK(FTPared[j][i + 1][1]);
						} else if(i == FCanal[j][0]->getNin() - 1) {
							Tpantant = __units::degCToK(FTPared[j][i - 1][1]);
							Tpantpos = __units::degCToK(FTPared[j][i][1]);
						} else {
							Tpantant = __units::degCToK(FTPared[j][i - 1][1]);
							Tpantpos = __units::degCToK(FTPared[j][i + 1][1]);
						}

						// Temperatura en haces posterior y anterior en el instante anterior
						if(j == 0 && j == FNumeroHacesCanales - 1) {
							Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
							Tphazpost = __units::degCToK(FTSuperficie[i][0]);
						} else if(j == 0) {
							Tphazant = __units::degCToK(FTParedAnt[j][i][1]); // As� no hay salto de temperaturas
							Tphazpost = __units::degCToK(FTParedAnt[j + 1][i][1]);
						} else if(j == FNumeroHacesCanales - 1) {
							Tphazant = __units::degCToK(FTPared[j - 1][i][1]);
							Tphazpost = __units::degCToK(FTSuperficie[i][0]);
						} else {
							Tphazant = __units::degCToK(FTPared[j - 1][i][1]);
							Tphazpost = __units::degCToK(FTPared[j + 1][i][1]);
						}

						//...en los nodos interior, medio y exterior.
						Tpant0 = __units::degCToK(FTPared[j][i][0]);
						Tpant1 = __units::degCToK(FTPared[j][i][1]);
						Tpant2 = __units::degCToK(FTPared[j][i][2]);

						if(EsPrimeraVez) {
							FTPared[j][i][1] = FSUMTPPromedio[j][i][1][0] / FSUMTPPromedio[j][i][1][1];

						} else {
							FTPared[j][i][1] = (FSUMTime[j] * Tpantant / FResistAxialAnt[j][i] + FSUMTime[j] * Tpantpos / FResistAxialPost[j][i] +
												FSUMTime[j] * Tphazant / FResistEntreHacesAnt[j][i]
												+ FSUMTime[j] * Tphazpost / FResistEntreHacesPost[j][i] + FSUMTPPromedio[j][i][1][0])
											   / (FSUMTime[j] / FResistAxialAnt[j][i] + FSUMTime[j] / FResistAxialPost[j][i] + FSUMTime[j] / FResistEntreHacesAnt[j][i]
												  + FSUMTime[j] / FResistEntreHacesPost[j][i]
												  + FSUMTPPromedio[j][i][1][1]);
						}
						FTPared[j][i][0] = (FSUMTime[j] * Tpant1 / FResistRadial[j][i][1] + FSUMTPPromedio[j][i][0][0]) /
										   (FSUMTime[j] / FResistRadial[j][i][1] + FSUMTPPromedio[j][i][0][1]);
						FTPared[j][i][2] = (FSUMTime[j] * Tpant1 / FResistRadial[j][i][0] + FSUMTPPromedio[j][i][2][0]) /
										   (FSUMTime[j] / FResistRadial[j][i][0] + FSUMTPPromedio[j][i][2][1]);

						if(ErrorTp < fabs(Tpant1 - FTPared[j][i][1])) {
							ErrorTp = fabs(Tpant1 - FTPared[j][i][1]);
						}

						if(j + 1 == FNumeroHacesCanales) {
							Tpsuperficie_int_ant = __units::degCToK(FTSuperficie[i][0]);
							Tpsuperficie_med_ant = __units::degCToK(FTSuperficie[i][1]);
							Tpsuperficie_ext_ant = __units::degCToK(FTSuperficie[i][2]);
							if(i == 0) {
								Tpsup_ant_ant = __units::degCToK(FTuboEntradaDPF->GetTPTubo(1, FNodoTuboEntrada));
								Tpsup_ant_post = __units::degCToK(FTSuperficie[i + 1][2]);
							} else if(i == FCanal[j][0]->getNin() - 1) {
								Tpsup_ant_ant = __units::degCToK(FTSuperficie[i - 1][2]);
								Tpsup_ant_post = __units::degCToK(FTuboSalidaDPF->GetTPTubo(1, FNodoTuboSalida));
							} else {
								Tpsup_ant_ant = __units::degCToK(FTSuperficie[i + 1][2]);
								Tpsup_ant_post = __units::degCToK(FTSuperficie[i - 1][2]);
							}

							if(FHayConcentrico) {
#ifdef ConcetricElement
								FTSuperficie[i][0] = (FSUMTime[j] * Tpant1 / FResistEntreHacesPost[j][i] + FSUMTPPromedioConc[0][0][i]) /
													 (FSUMTime[j] / FResistEntreHacesPost[j][i] + FSUMTPPromedioConc[1][0][i]);
								FTSuperficie[i][1] = (FSUMTime[j] * Tpant1 / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i]) +
													  FSUMTPPromedioConc[0][1][i]) /
													 (FSUMTime[j] / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i]) + FSUMTPPromedioConc[1][1][i]);
								FTSuperficie[i][2] = (FSUMTPPromedioConc[0][0][i] + FSUMTime[j] * Tpant1 / (FResistConduccionAislante[i] + 1 /
													  (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistEntreHacesPost[j][i])
													  +
													  FSUMTime[j] * Tpsup_ant_ant / FResistAxialMetalExtAnt[i] + FSUMTime[j] * Tpsup_ant_post / FResistAxialMetalExtPost[i]) /
													 (FSUMTPPromedioConc[1][0][i] + FSUMTime[j] / (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 /
															 FResistRadiacionAire[i]) + FResistConduccionMetal[i] + FResistEntreHacesPost[j][i]) +
													  FSUMTime[j] / FResistAxialMetalExtAnt[i] + FSUMTime[j] / FResistAxialMetalExtPost[i]);
#endif
							} else {
								FTSuperficie[i][0] = (FSUMTime[j] * Tpant1 / FResistEntreHacesPost[j][i]
													  + FSUMTime[j] * FTExt
													  / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
														 FResistConveccionExt[i] + FResistConduccionAislante[i]))
													 / (FSUMTime[j] / FResistEntreHacesPost[j][i]
														+ FSUMTime[j]
														/ (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
														   FResistConveccionExt[i] + FResistConduccionAislante[i]));
								FTSuperficie[i][1] = (FSUMTime[j] * Tpant1 / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i])
													  + FSUMTime[j] * FTExt / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
															  FResistConveccionExt[i]))
													 / (FSUMTime[j] / (FResistEntreHacesPost[j][i] + FResistConduccionAislante[i])
														+ FSUMTime[j] / (1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) + FResistConduccionMetal[i] +
																FResistConveccionExt[i]));
								FTSuperficie[i][2] =
									(FSUMTime[j] * FTExt / FResistConveccionExt[i]
									 + FSUMTime[j] * Tpant1
									 / (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
										FResistConduccionMetal[i] + FResistEntreHacesPost[j][i])
									 + FSUMTime[j] * Tpsup_ant_ant / FResistAxialMetalExtAnt[i] + FSUMTime[j] * Tpsup_ant_post / FResistAxialMetalExtPost[i])
									/ (FSUMTime[j] / FResistConveccionExt[i]
									   + FSUMTime[j]
									   / (FResistConduccionAislante[i] + 1 / (1 / FResistConduccionAire[i] + 1 / FResistRadiacionAire[i]) +
										  FResistConduccionMetal[i]
										  + FResistEntreHacesPost[j][i]) + FSUMTime[j] / FResistAxialMetalExtAnt[i] + FSUMTime[j] / FResistAxialMetalExtPost[i]);
							}
							if(ErrorTp < fabs(Tpsuperficie_ext_ant - FTSuperficie[i][2])) {
								ErrorTp = fabs(Tpsuperficie_ext_ant - FTSuperficie[i][2]);
							}

							FTSuperficie[i][0] = __units::KTodegC(FTSuperficie[i][0]);
							FTSuperficie[i][1] = __units::KTodegC(FTSuperficie[i][1]);
							FTSuperficie[i][2] = __units::KTodegC(FTSuperficie[i][2]);
						}

						for(int k = 0; k < 3; k++) {
							FTPared[j][i][k] = __units::KTodegC(FTPared[j][i][k]);
						}

					}
					EsPrimeraVez = false;
				}
			}
			for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
				for(int k = 0; k < 3; k++) {
					for(int h = 0; h < 2; h++) {
						FSUMTPPromedio[j][i][k][h] = 0.;
					}
				}
			}
			FSUMTime[j] = 0.;
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculaTemperaturaParedSinMotor en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::LeeResultadosMediosDPF(const char *FileWAM, fpos_t &filepos) {
	int NumVars = 0, TipoVar = 0;
	try {
		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d ", &FNumResMedios);
		FResMediosDPF = new stResMediosDPF*[FNumResMedios];
		for(int i = 0; i < FNumResMedios; i++) {
			FResMediosDPF[i] = new stResMediosDPF[FNumeroHacesCanales];
		}
		FTiempoMedSUM = 0.;
		FControlResMed = 1;

		for(int i = 0; i < FNumResMedios; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FResMediosDPF[i][j].VelocidadParedCanalEntrada = false;
				FResMediosDPF[i][j].VelocidadParedCanalEntradaSUM = 0.;
				FResMediosDPF[i][j].VelocidadParedCanalEntradaMED = 0;
				FResMediosDPF[i][j].VelocidadParedCanalSalida = false;
				FResMediosDPF[i][j].VelocidadParedCanalSalidaSUM = 0.;
				FResMediosDPF[i][j].VelocidadParedCanalSalidaMED = 0.;
				FResMediosDPF[i][j].DPFSootMass = false;
				FResMediosDPF[i][j].DPFSootMassSUM = 0.;
				FResMediosDPF[i][j].DPFSootMassMED = 0.;
				FResMediosDPF[i][j].BeamSootMass = false;
				FResMediosDPF[i][j].BeamSootMassSUM = 0.;
				FResMediosDPF[i][j].BeamSootMassMED = 0.;
				FResMediosDPF[i][j].CVSootMass = false;
				FResMediosDPF[i][j].CVSootMassSUM = 0.;
				FResMediosDPF[i][j].CVSootMassMED = 0.;
				FResMediosDPF[i][j].WallSootMass = false;
				FResMediosDPF[i][j].WallSootMassSUM = 0.;
				FResMediosDPF[i][j].WallSootMassMED = 0.;
				FResMediosDPF[i][j].LayerSootMass = false;
				FResMediosDPF[i][j].LayerSootMassSUM = 0.;
				FResMediosDPF[i][j].LayerSootMassMED = 0.;

				FResMediosDPF[i][j].LayerSootMassDep = false;
				FResMediosDPF[i][j].LayerSootMassDepSUM = 0.;
				FResMediosDPF[i][j].LayerSootMassDepMED = 0.;

				FResMediosDPF[i][j].EspesorSoot = false;
				FResMediosDPF[i][j].EspesorSootSUM = 0.;
				FResMediosDPF[i][j].EspesorSootMED = 0.;

				FResMediosDPF[i][j].EspesorSootIn = false;
				FResMediosDPF[i][j].EspesorSootInSUM = 0.;
				FResMediosDPF[i][j].EspesorSootInMED = 0.;

				FResMediosDPF[i][j].TemperaturaParedCE = false;
				FResMediosDPF[i][j].TemperaturaParedCESUM = 0.;
				FResMediosDPF[i][j].TemperaturaParedCEMED = 0.;
				FResMediosDPF[i][j].TemperaturaIntermediaPared = false;
				FResMediosDPF[i][j].TemperaturaIntermediaParedSUM = 0.;
				FResMediosDPF[i][j].TemperaturaIntermediaParedMED = 0.;
				FResMediosDPF[i][j].TemperaturaParedCS = false;
				FResMediosDPF[i][j].TemperaturaParedCSSUM = 0.;
				FResMediosDPF[i][j].TemperaturaParedCSMED = 0.;
				FResMediosDPF[i][j].TemperaturaInternaSuperficie = false;
				FResMediosDPF[i][j].TemperaturaInternaSuperficieSUM = 0.;
				FResMediosDPF[i][j].TemperaturaInternaSuperficieMED = 0.;
				FResMediosDPF[i][j].TemperaturaExternaSuperficie = false;
				FResMediosDPF[i][j].TemperaturaExternaSuperficieSUM = 0.;
				FResMediosDPF[i][j].TemperaturaExternaSuperficieMED = 0.;
				FResMediosDPF[i][j].TemperaturaMediaSuperficie = false;
				FResMediosDPF[i][j].TemperaturaMediaSuperficieSUM = 0.;
				FResMediosDPF[i][j].TemperaturaMediaSuperficieMED = 0.;
				FResMediosDPF[i][j].KwallClean = false;
				FResMediosDPF[i][j].KwallCleanSUM = 0;
				FResMediosDPF[i][j].KwallCleanMED = 0;
				FResMediosDPF[i][j].KwallLoaded = false;
				FResMediosDPF[i][j].KwallLoadedSUM = 0;
				FResMediosDPF[i][j].KwallLoadedMED = 0;
				FResMediosDPF[i][j].Kwall = false;
				FResMediosDPF[i][j].KwallSUM = 0;
				FResMediosDPF[i][j].KwallMED = 0;

				FResMediosDPF[i][j].KsootIn = false;
				FResMediosDPF[i][j].KsootInSUM = 0.;
				FResMediosDPF[i][j].KsootInMED = 0.;
				FResMediosDPF[i][j].KsootDep = false;
				FResMediosDPF[i][j].KsootDepSUM = 0.;
				FResMediosDPF[i][j].KsootDepMED = 0.;

				FResMediosDPF[i][j].Ksoot = false;
				FResMediosDPF[i][j].KsootSUM = 0.;
				FResMediosDPF[i][j].KsootMED = 0.;
				FResMediosDPF[i][j].Eficiencia = false;
				FResMediosDPF[i][j].EficienciaSUM = 0.;
				FResMediosDPF[i][j].GastoParedPondSUM = 0.;
				FResMediosDPF[i][j].EficienciaMED = 0.;
				FResMediosDPF[i][j].Porosidad = false;
				FResMediosDPF[i][j].PorosidadSUM = 0.;
				FResMediosDPF[i][j].PorosidadMED = 0.;
				FResMediosDPF[i][j].CoeficienteParticion = false;
				FResMediosDPF[i][j].CoeficienteParticionSUM = 0.;
				FResMediosDPF[i][j].CoeficienteParticionMED = 0.;
				FResMediosDPF[i][j].DiametroUC = false;
				FResMediosDPF[i][j].DiametroUCSUM = 0.;
				FResMediosDPF[i][j].DiametroUCMED = 0.;
				FResMediosDPF[i][j].ShapeFactor = false;
				FResMediosDPF[i][j].ShapeFactorSUM = 0.;
				FResMediosDPF[i][j].ShapeFactorMED = 0.;
				FResMediosDPF[i][j].Kreg1 = false;
				FResMediosDPF[i][j].Kreg1SUM = 0.;
				FResMediosDPF[i][j].Kreg1MED = 0.;
				FResMediosDPF[i][j].Kreg2 = false;
				FResMediosDPF[i][j].Kreg2SUM = 0.;
				FResMediosDPF[i][j].Kreg2MED = 0.;
				FResMediosDPF[i][j].Qreg = false;
				FResMediosDPF[i][j].QregSUM = 0.;
				FResMediosDPF[i][j].QregMED = 0.;
				//FResMediosDPF[i][j].Q1=false;
				//FResMediosDPF[i][j].Q1SUM=0.;
				//FResMediosDPF[i][j].Q1MED=0.;
				//FResMediosDPF[i][j].Q2=false;
				//FResMediosDPF[i][j].Q2SUM=0.;
				//FResMediosDPF[i][j].Q2MED=0.;
				FResMediosDPF[i][j].TasaFraccionMasicaEspecies = false;
				FResMediosDPF[i][j].TasaFraccionSUM = new double[FNumeroEspecies - FIntEGR];
				FResMediosDPF[i][j].TasaFraccionMED = new double[FNumeroEspecies - FIntEGR];
				FResMediosDPF[i][j].FraccionMasicaEspeciesSalida = false;
				FResMediosDPF[i][j].FraccionSalidaSUM = new double[FNumeroEspecies - FIntEGR];
				FResMediosDPF[i][j].FraccionSalidaMED = new double[FNumeroEspecies - FIntEGR];
				for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
					FResMediosDPF[i][j].TasaFraccionSUM[k] = 0.;
					FResMediosDPF[i][j].TasaFraccionMED[k] = 0.;
					FResMediosDPF[i][j].FraccionSalidaSUM[k] = 0.;
					FResMediosDPF[i][j].FraccionSalidaMED[k] = 0.;
				}
				FResMediosDPF[i][j].EficienciaBrown = false;
				FResMediosDPF[i][j].EficienciaBrownSUM = 0.;
				FResMediosDPF[i][j].EficienciaBrownMED = 0.;
				FResMediosDPF[i][j].EficienciaInter = false;
				FResMediosDPF[i][j].EficienciaInterSUM = 0.;
				FResMediosDPF[i][j].EficienciaInterMED = 0.;
				FResMediosDPF[i][j].EficienciaIner = false;
				FResMediosDPF[i][j].EficienciaInerSUM = 0.;
				FResMediosDPF[i][j].EficienciaInerMED = 0.;
				FResMediosDPF[i][j].EficienciaPLSUM = 0.;
				FResMediosDPF[i][j].EficienciaPLMED = 0.;

				fscanf(fich, "%lf %d ", &FResMediosDPF[i][j].Distancia, &NumVars);

				for(int k = 0; k < NumVars; k++) {
					fscanf(fich, "%d ", &TipoVar);
					switch(TipoVar) {
					case 0:
						FResMediosDPF[i][j].VelocidadParedCanalEntrada = true;
						break;
					case 1:
						FResMediosDPF[i][j].VelocidadParedCanalSalida = true;
						break;
					case 2:
						FResMediosDPF[i][j].DPFSootMass = true;
						break;
					case 3:
						FResMediosDPF[i][j].BeamSootMass = true;
						break;
					case 4:
						FResMediosDPF[i][j].CVSootMass = true;
						break;
					case 5:
						FResMediosDPF[i][j].WallSootMass = true;
						break;
					case 6:
						FResMediosDPF[i][j].LayerSootMass = true;
						break;
					case 7:
						FResMediosDPF[i][j].EspesorSoot = true;
						break;
					case 8:
						FResMediosDPF[i][j].Kwall = true;
						break;
					case 9:
						FResMediosDPF[i][j].Ksoot = true;
						break;
					case 10:
						FResMediosDPF[i][j].Eficiencia = true;
						break;
					case 11:
						FResMediosDPF[i][j].Porosidad = true;
						break;
					case 12:
						FResMediosDPF[i][j].CoeficienteParticion = true;
						break;
					case 13:
						FResMediosDPF[i][j].DiametroUC = true;
						break;
					case 14:
						FResMediosDPF[i][j].ShapeFactor = true;
						break;
					case 15:
						FResMediosDPF[i][j].Kreg1 = true;
						break;
					case 16:
						FResMediosDPF[i][j].Kreg2 = true;
						break;
					case 17:
						FResMediosDPF[i][j].Qreg = true;
						break;
					case 18:
						FResMediosDPF[i][j].TemperaturaParedCE = true;
						break;
					case 19:
						FResMediosDPF[i][j].TemperaturaIntermediaPared = true;
						break;
					case 20:
						FResMediosDPF[i][j].TemperaturaParedCS = true;
						break;
					case 21:
						FResMediosDPF[i][j].TemperaturaExternaSuperficie = true;
						break;
					case 22:
						FResMediosDPF[i][j].TemperaturaMediaSuperficie = true;
						break;
					case 23:
						FResMediosDPF[i][j].TemperaturaInternaSuperficie = true;
						break;
					case 24:
						FResMediosDPF[i][j].FraccionMasicaEspeciesSalida = true;
						break;
					case 25:
						FResMediosDPF[i][j].EficienciaBrown = true;
						break;
					case 26:
						FResMediosDPF[i][j].EficienciaInter = true;
						break;
					case 27:
						FResMediosDPF[i][j].EficienciaIner = true;
						break;
					case 28:
						FResMediosDPF[i][j].EficienciaPL = true;
						break;
					case 29:
						FResMediosDPF[i][j].KwallClean = true;
						break;
					case 30:
						FResMediosDPF[i][j].KwallLoaded = true;
						break;

					case 31:
						FResMediosDPF[i][j].KsootIn = true;
						break;
					case 32:
						FResMediosDPF[i][j].KsootDep = true;
						break;
					case 33:
						FResMediosDPF[i][j].EspesorSootIn = true;
						break;
					case 34:
						FResMediosDPF[i][j].LayerSootMassDep = true;
						break;

					default:
						std::cout << "Resultados medios en la DPF " << FNumeroDPF << " no implementados " << std::endl;
					}
				}
			}
		}

		fscanf(fich, "%d %d", &FNumResMediosCE, &FNumResMediosCS);

		fgetpos(fich, &filepos);
		fclose(fich);

		if(FNumResMediosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->LeeResultadosMediosCanalDPF(FNumResMediosCE, FileWAM, filepos);
			}
		}
		if(FNumResMediosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->LeeResultadosMediosCanalDPF(FNumResMediosCS, FileWAM, filepos);
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::LeeResultadosMediosDPF en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CabeceraResultadosMedios(stringstream& medoutput, stEspecies *DatosEspecies) const {
	try {

		std::string Label;
		std::ostringstream TextDist;
		TextDist.precision(8);

		for(int i = 0; i < FNumResMedios; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				TextDist << FResMediosDPF[i][j].Distancia;
				if(FResMediosDPF[i][j].VelocidadParedCanalEntrada) {
					Label = "\t" + PutLabel(802) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(909);
					medoutput << Label.c_str();
					//fprintf(fich,"\tV_PCE_DPF_%d_Haz_%d_a_%5.3f_m(m/s)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].VelocidadParedCanalSalida) {
					Label = "\t" + PutLabel(803) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(909);
					medoutput << Label.c_str();
					//fprintf(fich,"\tV_PCS_DPF_%d_Haz_%d_a_%5.3f_m(m/s)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].DPFSootMass) {
					Label = "\t" + PutLabel(804) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tDPF_soot_mass_DPF_%d(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].BeamSootMass) {
					Label = "\t" + PutLabel(804) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tBeam_soot_mass_DPF_%d_Haz_%d(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].CVSootMass) {
					Label = "\t" + PutLabel(836) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tCV_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].WallSootMass) {
					Label = "\t" + PutLabel(834) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tWall_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].LayerSootMass) {
					Label = "\t" + PutLabel(835) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tLayer_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].LayerSootMassDep) {
					Label = "\t" + PutLabel(846) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					medoutput << Label.c_str();
					//fprintf(fich,"\tLayer_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].EspesorSoot) {
					Label = "\t" + PutLabel(805) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEspesor_Soot_DPF_Haz_%d_a_%5.3f_m(mm)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].EspesorSootIn) {
					Label = "\t" + PutLabel(843) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEspesor_Soot_DPF_Haz_%d_a_%5.3f_m(mm)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].KwallClean) {
					Label = "\t" + PutLabel(841) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].KwallLoaded) {
					Label = "\t" + PutLabel(842) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].Kwall) {
					Label = "\t" + PutLabel(806) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].KsootIn) {
					Label = "\t" + PutLabel(844) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].KsootDep) {
					Label = "\t" + PutLabel(845) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].Ksoot) {
					Label = "\t" + PutLabel(807) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					medoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].Eficiencia) {
					Label = "\t" + PutLabel(808) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].EficienciaBrown) {
					Label = "\t" + PutLabel(837) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEficienciaBrowniana_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);

				}
				if(FResMediosDPF[i][j].EficienciaInter) {
					Label = "\t" + PutLabel(838) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEficienciaInterception_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}

				if(FResMediosDPF[i][j].EficienciaIner) {
					Label = "\t" + PutLabel(839) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEficienciaInercial_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].EficienciaPL) {
					Label = "\t" + PutLabel(840) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tEficienciaInercial_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].Porosidad) {
					Label = "\t" + PutLabel(809) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tPorosidad_DPF_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].CoeficienteParticion) {
					Label = "\t" + PutLabel(810) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
					//fprintf(fich,"\tCoef_Particion_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].DiametroUC) {
					Label = "\t" + PutLabel(811) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					medoutput << Label.c_str();
					//fprintf(fich,"\tDiametroUC_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].ShapeFactor) {
					Label = "\t" + PutLabel(833) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
				}
				if(FResMediosDPF[i][j].Kreg1) {
					Label = "\t" + PutLabel(812) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(922);
					medoutput << Label.c_str();
					//fprintf(fich,"\tKreg_termica_DPF_%d_Haz_%d_a_%5.3f_m(1/s)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].Kreg2) {
					Label = "\t" + PutLabel(813) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(922);
					medoutput << Label.c_str();
					//fprintf(fich,"\tKreg_cata_DPF_Haz_%d_a_%5.3f_m(1/s)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].Qreg) {
					Label = "\t" + PutLabel(814) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(907);
					medoutput << Label.c_str();
					//fprintf(fich,"\tQreg_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
					/*if(FResMediosDPF[i][j].Q1)
					 fprintf(fich,"\tQ_HC_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
					 if(FResMediosDPF[i][j].Q2)
					 fprintf(fich,"\tQ_CO_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);*/
				}
				if(FResMediosDPF[i][j].TemperaturaParedCE) {
					Label = "\t" + PutLabel(815) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					medoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_CE_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].TemperaturaIntermediaPared) {
					Label = "\t" + PutLabel(816) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					medoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_Intermedia_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(FResMediosDPF[i][j].TemperaturaParedCS) {
					Label = "\t" + PutLabel(817) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					medoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_CS_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResMediosDPF[i][j].Distancia);
				}
				if(j + 1 == FNumeroHacesCanales) {
					if(FResMediosDPF[i][j].TemperaturaExternaSuperficie) {
						Label = "\t" + PutLabel(818) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						medoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Ext_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResMediosDPF[i][j].Distancia);
					}
					if(FResMediosDPF[i][j].TemperaturaMediaSuperficie) {
						Label = "\t" + PutLabel(819) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						medoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Med_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResMediosDPF[i][j].Distancia);
					}
					if(FResMediosDPF[i][j].TemperaturaInternaSuperficie) {
						Label = "\t" + PutLabel(820) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						medoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Int_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResMediosDPF[i][j].Distancia);
					}
				}
				/*if(FResMediosDPF[i][j].TasaFraccionMasicaEspecies){
				 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
				 fprintf(fich,"\tTasa_Y_%s_DPF_Haz_%d_a_%5.3f_m(1/s)",DatosEspecies[k].Nombre,FNumeroDPF,k,FResMediosDPF[i][j].Distancia);
				 }
				 }  */
				if(FResMediosDPF[i][j].FraccionMasicaEspeciesSalida) {
					for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
						Label = "\t" + PutLabel(821) + DatosEspecies[k].Nombre + PutLabel(822) + PutLabel(800) + std::to_string(
									FNumeroDPF) + PutLabel(801) + std::to_string(j + 1) + PutLabel(316)
								+ TextDist.str() + PutLabel(317) + PutLabel(901);
						medoutput << Label.c_str();
						//fprintf(fich,"\tY_%s_pared_DPF_Haz_%d_a_%5.3f_m(1/s)",DatosEspecies[k].Nombre,FNumeroDPF,k,FResMediosDPF[i][j].Distancia);
					}
				}

			}
		}

		if(FNumResMediosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->CabeceraResultadosMedios(medoutput, DatosEspecies);
			}
		}
		if(FNumResMediosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->CabeceraResultadosMedios(medoutput, DatosEspecies);
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CabeceraResultadosMedios en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::ImprimeResultadosMedios(stringstream& medoutput) const {
	try {

//FILE *fich=fopen(FileSALIDA,"a");

		for(int i = 0; i < FNumResMedios; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				if(FResMediosDPF[i][j].VelocidadParedCanalEntrada)
					medoutput << "\t" << FResMediosDPF[i][j].VelocidadParedCanalEntradaMED;
				if(FResMediosDPF[i][j].VelocidadParedCanalSalida)
					medoutput << "\t" << FResMediosDPF[i][j].VelocidadParedCanalSalidaMED;
				if(FResMediosDPF[i][j].DPFSootMass)
					medoutput << "\t" << FResMediosDPF[i][j].DPFSootMassMED;
				if(FResMediosDPF[i][j].BeamSootMass)
					medoutput << "\t" << FResMediosDPF[i][j].BeamSootMassMED;
				if(FResMediosDPF[i][j].CVSootMass)
					medoutput << "\t" << FResMediosDPF[i][j].CVSootMassMED;
				if(FResMediosDPF[i][j].WallSootMass)
					medoutput << "\t" << FResMediosDPF[i][j].WallSootMassMED;
				if(FResMediosDPF[i][j].LayerSootMass)
					medoutput << "\t" << FResMediosDPF[i][j].LayerSootMassMED;

				if(FResMediosDPF[i][j].LayerSootMassDep)
					medoutput << "\t" << FResMediosDPF[i][j].LayerSootMassDepMED;

				if(FResMediosDPF[i][j].EspesorSoot)
					medoutput << "\t" << FResMediosDPF[i][j].EspesorSootMED;

				if(FResMediosDPF[i][j].EspesorSootIn)
					medoutput << "\t" << FResMediosDPF[i][j].EspesorSootInMED;

				if(FResMediosDPF[i][j].KwallClean)
					medoutput << "\t" << FResMediosDPF[i][j].KwallCleanMED;
				if(FResMediosDPF[i][j].KwallLoaded)
					medoutput << "\t" << FResMediosDPF[i][j].KwallLoadedMED;
				if(FResMediosDPF[i][j].Kwall)
					medoutput << "\t" << FResMediosDPF[i][j].KwallMED;

				if(FResMediosDPF[i][j].KsootIn)
					medoutput << "\t" << FResMediosDPF[i][j].KsootInMED;
				if(FResMediosDPF[i][j].KsootDep)
					medoutput << "\t" << FResMediosDPF[i][j].KsootDepMED;

				if(FResMediosDPF[i][j].Ksoot)
					medoutput << "\t" << FResMediosDPF[i][j].KsootMED;
				if(FResMediosDPF[i][j].Eficiencia)
					medoutput << "\t" << FResMediosDPF[i][j].EficienciaMED;
				if(FResMediosDPF[i][j].EficienciaBrown)
					medoutput << "\t" << FResMediosDPF[i][j].EficienciaBrownMED;
				if(FResMediosDPF[i][j].EficienciaInter)
					medoutput << "\t" << FResMediosDPF[i][j].EficienciaInterMED;
				if(FResMediosDPF[i][j].EficienciaIner)
					medoutput << "\t" << FResMediosDPF[i][j].EficienciaInerMED;
				if(FResMediosDPF[i][j].EficienciaPL)
					medoutput << "\t" << FResMediosDPF[i][j].EficienciaPLMED;
				if(FResMediosDPF[i][j].Porosidad)
					medoutput << "\t" << FResMediosDPF[i][j].PorosidadMED;
				if(FResMediosDPF[i][j].CoeficienteParticion)
					medoutput << "\t" << FResMediosDPF[i][j].CoeficienteParticionMED;
				if(FResMediosDPF[i][j].DiametroUC)
					medoutput << "\t" << FResMediosDPF[i][j].DiametroUCMED;
				if(FResMediosDPF[i][j].ShapeFactor)
					medoutput << "\t" << FResMediosDPF[i][j].ShapeFactorMED;
				if(FResMediosDPF[i][j].Kreg1)
					medoutput << "\t" << FResMediosDPF[i][j].Kreg1MED;
				if(FResMediosDPF[i][j].Kreg2)
					medoutput << "\t" << FResMediosDPF[i][j].Kreg2MED;
				if(FResMediosDPF[i][j].Qreg)
					medoutput << "\t" << FResMediosDPF[i][j].QregMED;
				/*if(FResMediosDPF[i][j].Q1)
				 fprintf(fich,"\t%lg",FResMediosDPF[i][j].Q1MED);
				 if(FResMediosDPF[i][j].Q2)
				 fprintf(fich,"\t%lg",FResMediosDPF[i][j].Q2MED);*/
				if(FResMediosDPF[i][j].TemperaturaParedCE)
					medoutput << "\t" << FResMediosDPF[i][j].TemperaturaParedCEMED;
				if(FResMediosDPF[i][j].TemperaturaIntermediaPared)
					medoutput << "\t" << FResMediosDPF[i][j].TemperaturaIntermediaParedMED;
				if(FResMediosDPF[i][j].TemperaturaParedCS)
					medoutput << "\t" << FResMediosDPF[i][j].TemperaturaParedCSMED;
				if(j + 1 == FNumeroHacesCanales) {
					if(FResMediosDPF[i][j].TemperaturaExternaSuperficie)
						medoutput << "\t" << FResMediosDPF[i][j].TemperaturaExternaSuperficieMED;
					if(FResMediosDPF[i][j].TemperaturaMediaSuperficie)
						medoutput << "\t" << FResMediosDPF[i][j].TemperaturaMediaSuperficieMED;
					if(FResMediosDPF[i][j].TemperaturaInternaSuperficie)
						medoutput << "\t" << FResMediosDPF[i][j].TemperaturaInternaSuperficieMED;
				}
				/*if(FResMediosDPF[i][j].TasaFraccionMasicaEspecies){
				 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
				 fprintf(fich,"\t%lg",FResMediosDPF[i][j].TasaFraccionMED[k]);
				 }
				 }  */
				if(FResMediosDPF[i][j].FraccionMasicaEspeciesSalida) {
					for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
						medoutput << "\t" << FResMediosDPF[i][j].FraccionSalidaMED[k];
					}
				}

			}
		}

//fclose(fich);

		if(FNumResMediosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->ImprimeResultadosMedios(medoutput);
			}
		}
		if(FNumResMediosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->ImprimeResultadosMedios(medoutput);
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::ImprimeResultadosMedios en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculaResultadosMedios(double theta) {
	int n1 = 0, n2 = 0;
	double dist, d, Vble, Vble2, GastoParedPond1, GastoParedPond2, EficPond1, EficPond2;

	try {

		FTiempoMedSUM += FDeltaTimeDPF;

		for(int i = 0; i < FNumResMedios; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				dist = FResMediosDPF[i][j].Distancia / FCanal[j][0]->getXRef();
				n1 = (int) floor(dist);
				if(n1 >= FCanal[j][0]->getNin() - 1) {
					if(FResMediosDPF[i][j].VelocidadParedCanalEntrada)
						FResMediosDPF[i][j].VelocidadParedCanalEntradaSUM += FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].VelocidadParedCanalSalida)
						FResMediosDPF[i][j].VelocidadParedCanalSalidaSUM += FVelocidadPared[j][1][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].DPFSootMass)
						FResMediosDPF[i][j].DPFSootMassSUM += FDPFSootMass * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].BeamSootMass)
						FResMediosDPF[i][j].BeamSootMassSUM += FBeamSootMass[j] * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].CVSootMass)
						FResMediosDPF[i][j].CVSootMassSUM += FCVSootMass[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].WallSootMass)
						FResMediosDPF[i][j].WallSootMassSUM += FWallSootMass[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].LayerSootMass)
						FResMediosDPF[i][j].LayerSootMassSUM += FLayerSootMass[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;

					if(FResMediosDPF[i][j].LayerSootMassDep)
						FResMediosDPF[i][j].LayerSootMassDepSUM += FLayerSootMassDep[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;

					if(FResMediosDPF[i][j].EspesorSootIn)
						FResMediosDPF[i][j].EspesorSootInSUM += FEspesorSootIn[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;

					if(FResMediosDPF[i][j].EspesorSoot)
						FResMediosDPF[i][j].EspesorSootSUM += FEspesorSoot[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].KwallClean)
						FResMediosDPF[i][j].KwallCleanSUM += FKwallClean[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].KwallLoaded)
						FResMediosDPF[i][j].KwallLoadedSUM += FKwallLoaded[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].Kwall)
						FResMediosDPF[i][j].KwallSUM += FKwall[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;

					if(FResMediosDPF[i][j].KsootIn)
						FResMediosDPF[i][j].KsootInSUM += FKsootIn[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].KsootDep)
						FResMediosDPF[i][j].KsootDepSUM += FKsootDep[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;

					if(FResMediosDPF[i][j].Ksoot)
						FResMediosDPF[i][j].KsootSUM += FKsoot[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].Eficiencia)
						FResMediosDPF[i][j].GastoParedPondSUM += FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1] *
								FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i)
								* FDeltaTimeDPF;
					FResMediosDPF[i][j].EficienciaSUM += FEficiencia[j][FCanal[j][0]->getNin() - 1] *
														 FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1]
														 * FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i) * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].Porosidad)
						FResMediosDPF[i][j].PorosidadSUM += FPorosidad[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].CoeficienteParticion)
						FResMediosDPF[i][j].CoeficienteParticionSUM += FCoeficienteParticion[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].DiametroUC)
						FResMediosDPF[i][j].DiametroUCSUM += FDiametroUnidadColectora[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF * 1000;
					if(FResMediosDPF[i][j].ShapeFactor)
						if(FCalculoFiltrado) {
							FResMediosDPF[i][j].ShapeFactorSUM += FShapeFactor[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
						} else
							FResMediosDPF[i][j].ShapeFactorSUM += FShapeFactorDiscrete;
					if(FResMediosDPF[i][j].Kreg1)
						FResMediosDPF[i][j].Kreg1SUM += FKreg1[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].Kreg2)
						FResMediosDPF[i][j].Kreg2SUM += FKreg2[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].Qreg)
						FResMediosDPF[i][j].QregSUM += FQreg[j][FCanal[j][0]->getNin() - 1] * FDeltaTimeDPF;
					/*if(FResMediosDPF[i][j].Q1)
					 FResMediosDPF[i][j].Q1SUM+=FQ1[j][FCanal[j][0]->getNin()-1]*FDeltaTimeDPF;
					 if(FResMediosDPF[i][j].Q2)
					 FResMediosDPF[i][j].Q2SUM+=FQ2[j][FCanal[j][0]->getNin()-1]*FDeltaTimeDPF;  */
					if(FResMediosDPF[i][j].TemperaturaParedCE)
						FResMediosDPF[i][j].TemperaturaParedCESUM += FTPared[j][FCanal[j][0]->getNin() - 1][2] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].TemperaturaIntermediaPared)
						FResMediosDPF[i][j].TemperaturaIntermediaParedSUM += FTPared[j][FCanal[j][1]->getNin() - 1][1] * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].TemperaturaParedCS)
						FResMediosDPF[i][j].TemperaturaParedCSSUM += FTPared[j][FCanal[j][0]->getNin() - 1][0] * FDeltaTimeDPF;
					if(j + 1 == FNumeroHacesCanales) {
						if(FResMediosDPF[i][j].TemperaturaExternaSuperficie)
							FResMediosDPF[i][j].TemperaturaExternaSuperficieSUM += FTSuperficie[FCanal[j][0]->getNin() - 1][2] * FDeltaTimeDPF;
						if(FResMediosDPF[i][j].TemperaturaMediaSuperficie)
							FResMediosDPF[i][j].TemperaturaMediaSuperficieSUM += FTSuperficie[FCanal[j][0]->getNin() - 1][1] * FDeltaTimeDPF;
						if(FResMediosDPF[i][j].TemperaturaInternaSuperficie)
							FResMediosDPF[i][j].TemperaturaInternaSuperficieSUM += FTSuperficie[FCanal[j][0]->getNin() - 1][0] * FDeltaTimeDPF;
					}
					/* if(FResMediosDPF[i][j].TasaFraccionMasicaEspecies){
					 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
					 FResMediosDPF[i][j].TasaFraccionSUM[k]+=FTasaFraccionMasicaEspecie[j][FCanal[j][0]->getNin()-1][k]*FDeltaTimeDPF;
					 }
					 }  */
					if(FResMediosDPF[i][j].FraccionMasicaEspeciesSalida) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							FResMediosDPF[i][j].FraccionSalidaSUM[k] += FFraccionMasicaEspecieSalida[j][FCanal[j][0]->getNin() - 1][k] *
									FDeltaTimeDPF;
						}
					}
					if(FResMediosDPF[i][j].EficienciaBrown)
						FResMediosDPF[i][j].EficienciaBrownSUM += FEficienciaBrown[j][FCanal[j][0]->getNin() - 1] *
								FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1]
								* FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i) * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].EficienciaInter)
						FResMediosDPF[i][j].EficienciaInterSUM += FEficienciaInter[j][FCanal[j][0]->getNin() - 1] *
								FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1]
								* FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i) * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].EficienciaIner)
						FResMediosDPF[i][j].EficienciaInerSUM += FEficienciaIner[j][FCanal[j][0]->getNin() - 1] *
								FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1]
								* FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i) * FDeltaTimeDPF;
					if(FResMediosDPF[i][j].EficienciaPL)
						FResMediosDPF[i][j].EficienciaPLSUM += FEficienciaPL[j][FCanal[j][0]->getNin() - 1] *
															   FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1]
															   * FAreaVCCanalEntrada[j][FCanal[j][0]->getNin() - 1] * FCanal[j][0]->GetDensidad(i) * FDeltaTimeDPF;

				} else {
					n2 = n1 + 1;
					d = dist - (double) n1;
					if(FResMediosDPF[i][j].VelocidadParedCanalEntrada) {
						Vble = Interpola(FVelocidadPared[j][0][n1], FVelocidadPared[j][0][n2], 1., d);
						FResMediosDPF[i][j].VelocidadParedCanalEntradaSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].VelocidadParedCanalSalida) {
						Vble = Interpola(FVelocidadPared[j][1][n1], FVelocidadPared[j][1][n2], 1., d);
						FResMediosDPF[i][j].VelocidadParedCanalSalidaSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].DPFSootMass) {
						FResMediosDPF[i][j].DPFSootMassSUM += FDPFSootMass * FDeltaTimeDPF * 1000;
					}
					if(FResMediosDPF[i][j].BeamSootMass) {
						FResMediosDPF[i][j].BeamSootMassSUM += FBeamSootMass[j] * FDeltaTimeDPF * 1000;
					}
					if(FResMediosDPF[i][j].CVSootMass) {
						Vble = Interpola(FCVSootMass[j][n1], FCVSootMass[j][n2], 1., d);
						FResMediosDPF[i][j].CVSootMassSUM += Vble * FDeltaTimeDPF * 1000;
					}
					if(FResMediosDPF[i][j].WallSootMass) {
						Vble = Interpola(FWallSootMass[j][n1], FWallSootMass[j][n2], 1., d);
						FResMediosDPF[i][j].WallSootMassSUM += Vble * FDeltaTimeDPF * 1000;
					}
					if(FResMediosDPF[i][j].LayerSootMass) {
						Vble = Interpola(FLayerSootMass[j][n1], FLayerSootMass[j][n2], 1., d);
						FResMediosDPF[i][j].LayerSootMassSUM += Vble * FDeltaTimeDPF * 1000;

					}
					if(FResMediosDPF[i][j].LayerSootMassDep) {
						Vble = Interpola(FLayerSootMassDep[j][n1], FLayerSootMassDep[j][n2], 1., d);
						FResMediosDPF[i][j].LayerSootMassDepSUM += Vble * FDeltaTimeDPF * 1000;

					}
					if(FResMediosDPF[i][j].EspesorSootIn) {
						Vble = Interpola(FEspesorSootIn[j][n1], FEspesorSootIn[j][n2], 1., d);
						FResMediosDPF[i][j].EspesorSootInSUM += Vble * FDeltaTimeDPF;

					}
					if(FResMediosDPF[i][j].EspesorSoot) {
						Vble = Interpola(FEspesorSoot[j][n1], FEspesorSoot[j][n2], 1., d);
						FResMediosDPF[i][j].EspesorSootSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].KwallClean) {
						Vble = Interpola(FKwallClean[j][n1], FKwallClean[j][n2], 1., d);
						FResMediosDPF[i][j].KwallCleanSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].KwallLoaded) {
						Vble = Interpola(FKwallLoaded[j][n1], FKwallLoaded[j][n2], 1., d);
						FResMediosDPF[i][j].KwallLoadedSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].Kwall) {
						Vble = Interpola(FKwall[j][n1], FKwall[j][n2], 1., d);
						FResMediosDPF[i][j].KwallSUM += Vble * FDeltaTimeDPF;

					}
					if(FResMediosDPF[i][j].KsootIn) {
						Vble = Interpola(FKsootIn[j][n1], FKsootIn[j][n2], 1., d);
						FResMediosDPF[i][j].KsootInSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].KsootDep) {
						Vble = Interpola(FKsootDep[j][n1], FKsootDep[j][n2], 1., d);
						FResMediosDPF[i][j].KsootDepSUM += Vble * FDeltaTimeDPF;

					}
					if(FResMediosDPF[i][j].Ksoot) {
						Vble = Interpola(FKsoot[j][n1], FKsoot[j][n2], 1., d);
						FResMediosDPF[i][j].KsootSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].Eficiencia) {
						EficPond1 = FEficiencia[j][n1] * FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] * FCanal[j][0]->GetDensidad(n1);
						EficPond2 = FEficiencia[j][n2] * FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] * FCanal[j][0]->GetDensidad(n2);
						Vble = Interpola(EficPond1, EficPond2, 1., d);
						GastoParedPond1 = FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] * FCanal[j][0]->GetDensidad(n1);
						GastoParedPond2 = FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] * FCanal[j][0]->GetDensidad(n2);
						Vble2 = Interpola(GastoParedPond1, GastoParedPond2, 1., d);
						FResMediosDPF[i][j].EficienciaSUM += Vble * FDeltaTimeDPF;
						FResMediosDPF[i][j].GastoParedPondSUM += Vble2 * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].Porosidad) {
						Vble = Interpola(FPorosidad[j][n1], FPorosidad[j][n2], 1., d);
						FResMediosDPF[i][j].PorosidadSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].CoeficienteParticion) {
						Vble = Interpola(FCoeficienteParticion[j][n1], FCoeficienteParticion[j][n2], 1., d);
						FResMediosDPF[i][j].CoeficienteParticionSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].DiametroUC) {
						Vble = Interpola(FDiametroUnidadColectora[j][n1], FDiametroUnidadColectora[j][n2], 1., d);
						FResMediosDPF[i][j].DiametroUCSUM += Vble * FDeltaTimeDPF * 1000;
					}
					if(FResMediosDPF[i][j].ShapeFactor) {
						if(FCalculoFiltrado) {
							Vble = Interpola(FShapeFactor[j][n1], FShapeFactor[j][n2], 1., d);
							FResMediosDPF[i][j].ShapeFactorSUM += Vble * FDeltaTimeDPF;
						} else
							FResMediosDPF[i][j].ShapeFactorSUM += FShapeFactorDiscrete;
					}
					if(FResMediosDPF[i][j].Kreg1) {
						Vble = Interpola(FKreg1[j][n1], FKreg1[j][n2], 1., d);
						FResMediosDPF[i][j].Kreg1SUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].Kreg2) {
						Vble = Interpola(FKreg2[j][n1], FKreg2[j][n2], 1., d);
						FResMediosDPF[i][j].Kreg2SUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].Qreg) {
						Vble = Interpola(FQreg[j][n1], FQreg[j][n2], 1., d);
						FResMediosDPF[i][j].QregSUM += Vble * FDeltaTimeDPF;
						/* }if(FResMediosDPF[i][j].Q1){
						 Vble=Interpola(FQ1[j][n1],FQ1[j][n2],1.,d);
						 FResMediosDPF[i][j].Q1SUM+=Vble*FDeltaTimeDPF;
						 }if(FResMediosDPF[i][j].Q2){
						 Vble=Interpola(FQ2[j][n1],FQ2[j][n2],1.,d);
						 FResMediosDPF[i][j].Q2SUM+=Vble*FDeltaTimeDPF;     */
					}
					if(FResMediosDPF[i][j].TemperaturaParedCE) {
						Vble = Interpola(FTPared[j][n1][2], FTPared[j][n2][2], 1., d);
						FResMediosDPF[i][j].TemperaturaParedCESUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].TemperaturaIntermediaPared) {
						Vble = Interpola(FTPared[j][n1][1], FTPared[j][n2][1], 1., d);
						FResMediosDPF[i][j].TemperaturaIntermediaParedSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].TemperaturaParedCS) {
						Vble = Interpola(FTPared[j][n1][0], FTPared[j][n2][0], 1., d);
						FResMediosDPF[i][j].TemperaturaParedCSSUM += Vble * FDeltaTimeDPF;
					}
					if(j + 1 == FNumeroHacesCanales) {
						if(FResMediosDPF[i][j].TemperaturaExternaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][2], FTSuperficie[n2][2], 1., d);
							FResMediosDPF[i][j].TemperaturaExternaSuperficieSUM += Vble * FDeltaTimeDPF;
						}
						if(FResMediosDPF[i][j].TemperaturaMediaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][1], FTSuperficie[n2][1], 1., d);
							FResMediosDPF[i][j].TemperaturaMediaSuperficieSUM += Vble * FDeltaTimeDPF;
						}
						if(FResMediosDPF[i][j].TemperaturaInternaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][0], FTSuperficie[n2][0], 1., d);
							FResMediosDPF[i][j].TemperaturaInternaSuperficieSUM += Vble * FDeltaTimeDPF;
						}
					}
					/* if(FResMediosDPF[i][j].TasaFraccionMasicaEspecies){
					 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
					 Vble=Interpola(FTasaFraccionMasicaEspecie[j][n1][k],FTasaFraccionMasicaEspecie[j][n2][k],1.,d);
					 FResMediosDPF[i][j].TasaFraccionSUM[k]+=Vble*FDeltaTimeDPF;
					 }
					 }*/
					if(FResMediosDPF[i][j].FraccionMasicaEspeciesSalida) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							Vble = Interpola(FFraccionMasicaEspecieSalida[j][n1][k], FFraccionMasicaEspecieSalida[j][n2][k], 1., d);
							FResMediosDPF[i][j].FraccionSalidaSUM[k] += Vble * FDeltaTimeDPF;
						}
					}
					if(FResMediosDPF[i][j].EficienciaBrown) {
						EficPond1 = FEficienciaBrown[j][n1] * FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] *
									FCanal[j][0]->GetDensidad(n1);
						EficPond2 = FEficienciaBrown[j][n2] * FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] *
									FCanal[j][0]->GetDensidad(n2);
						Vble = Interpola(EficPond1, EficPond2, 1., d);
						FResMediosDPF[i][j].EficienciaBrownSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].EficienciaInter) {
						EficPond1 = FEficienciaInter[j][n1] * FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] *
									FCanal[j][0]->GetDensidad(n1);
						EficPond2 = FEficienciaInter[j][n2] * FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] *
									FCanal[j][0]->GetDensidad(n2);
						Vble = Interpola(EficPond1, EficPond2, 1., d);
						FResMediosDPF[i][j].EficienciaInterSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].EficienciaIner) {
						EficPond1 = FEficienciaIner[j][n1] * FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] * FCanal[j][0]->GetDensidad(
										n1);
						EficPond2 = FEficienciaIner[j][n2] * FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] * FCanal[j][0]->GetDensidad(
										n2);
						Vble = Interpola(EficPond1, EficPond2, 1., d);
						FResMediosDPF[i][j].EficienciaInerSUM += Vble * FDeltaTimeDPF;
					}
					if(FResMediosDPF[i][j].EficienciaPL) {
						EficPond1 = FEficienciaPL[j][n1] * FVelocidadPared[j][0][n1] * FAreaVCCanalEntrada[j][n1] * FCanal[j][0]->GetDensidad(
										n1);
						EficPond2 = FEficienciaPL[j][n2] * FVelocidadPared[j][0][n2] * FAreaVCCanalEntrada[j][n2] * FCanal[j][0]->GetDensidad(
										n2);
						Vble = Interpola(EficPond1, EficPond2, 1., d);
						FResMediosDPF[i][j].EficienciaPLSUM += Vble * FDeltaTimeDPF;

					}
				}
			}
		}

		if(theta > FControlResMed * FAnguloTotalCiclo) {

			for(int i = 0; i < FNumResMedios; i++) {
				for(int j = 0; j < FNumeroHacesCanales; j++) {
					if(FResMediosDPF[i][j].VelocidadParedCanalEntrada) {
						FResMediosDPF[i][j].VelocidadParedCanalEntradaMED = FResMediosDPF[i][j].VelocidadParedCanalEntradaSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].VelocidadParedCanalEntradaSUM = 0.;
					}
					if(FResMediosDPF[i][j].VelocidadParedCanalSalida) {
						FResMediosDPF[i][j].VelocidadParedCanalSalidaMED = FResMediosDPF[i][j].VelocidadParedCanalSalidaSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].VelocidadParedCanalSalidaSUM = 0.;
					}
					if(FResMediosDPF[i][j].DPFSootMass) {
						FResMediosDPF[i][j].DPFSootMassMED = FResMediosDPF[i][j].DPFSootMassSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].DPFSootMassSUM = 0.;
					}
					if(FResMediosDPF[i][j].BeamSootMass) {
						FResMediosDPF[i][j].BeamSootMassMED = FResMediosDPF[i][j].BeamSootMassSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].BeamSootMassSUM = 0.;
					}
					if(FResMediosDPF[i][j].CVSootMass) {
						FResMediosDPF[i][j].CVSootMassMED = FResMediosDPF[i][j].CVSootMassSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].CVSootMassSUM = 0.;
					}
					if(FResMediosDPF[i][j].WallSootMass) {
						FResMediosDPF[i][j].WallSootMassMED = FResMediosDPF[i][j].WallSootMassSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].WallSootMassSUM = 0.;
					}
					if(FResMediosDPF[i][j].LayerSootMass) {
						FResMediosDPF[i][j].LayerSootMassMED = FResMediosDPF[i][j].LayerSootMassSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].LayerSootMassSUM = 0.;

					}
					if(FResMediosDPF[i][j].LayerSootMassDep) {
						FResMediosDPF[i][j].LayerSootMassDepMED = FResMediosDPF[i][j].LayerSootMassDepSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].LayerSootMassDepSUM = 0.;

					}
					if(FResMediosDPF[i][j].EspesorSootIn) {
						FResMediosDPF[i][j].EspesorSootInMED = FResMediosDPF[i][j].EspesorSootInSUM / FTiempoMedSUM * 1000;
						FResMediosDPF[i][j].EspesorSootInSUM = 0.;

					}
					if(FResMediosDPF[i][j].EspesorSoot) {
						FResMediosDPF[i][j].EspesorSootMED = FResMediosDPF[i][j].EspesorSootSUM / FTiempoMedSUM * 1000;
						FResMediosDPF[i][j].EspesorSootSUM = 0.;
					}
					if(FResMediosDPF[i][j].KwallClean) {
						FResMediosDPF[i][j].KwallCleanMED = FResMediosDPF[i][j].KwallCleanSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KwallCleanSUM = 0.;
					}
					if(FResMediosDPF[i][j].KwallLoaded) {
						FResMediosDPF[i][j].KwallLoadedMED = FResMediosDPF[i][j].KwallLoadedSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KwallLoadedSUM = 0.;
					}
					if(FResMediosDPF[i][j].Kwall) {
						FResMediosDPF[i][j].KwallMED = FResMediosDPF[i][j].KwallSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KwallSUM = 0.;

					}
					if(FResMediosDPF[i][j].KsootIn) {
						FResMediosDPF[i][j].KsootInMED = FResMediosDPF[i][j].KsootInSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KsootSUM = 0.;
					}
					if(FResMediosDPF[i][j].KsootDep) {
						FResMediosDPF[i][j].KsootDepMED = FResMediosDPF[i][j].KsootDepSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KsootDepSUM = 0.;

					}
					if(FResMediosDPF[i][j].Ksoot) {
						FResMediosDPF[i][j].KsootMED = FResMediosDPF[i][j].KsootSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].KsootSUM = 0.;
					}
					if(FResMediosDPF[i][j].Eficiencia) {
						FResMediosDPF[i][j].EficienciaMED = FResMediosDPF[i][j].EficienciaSUM / FResMediosDPF[i][j].GastoParedPondSUM;
						FResMediosDPF[i][j].EficienciaSUM = 0.;
					}
					if(FResMediosDPF[i][j].Porosidad) {
						FResMediosDPF[i][j].PorosidadMED = FResMediosDPF[i][j].PorosidadSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].PorosidadSUM = 0.;
					}
					if(FResMediosDPF[i][j].CoeficienteParticion) {
						FResMediosDPF[i][j].CoeficienteParticionMED = FResMediosDPF[i][j].CoeficienteParticionSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].CoeficienteParticionSUM = 0.;
					}
					if(FResMediosDPF[i][j].DiametroUC) {
						FResMediosDPF[i][j].DiametroUCMED = FResMediosDPF[i][j].DiametroUCSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].DiametroUCSUM = 0.;
					}
					if(FResMediosDPF[i][j].ShapeFactor) {
						if(FCalculoFiltrado) {
							FResMediosDPF[i][j].ShapeFactorMED = FResMediosDPF[i][j].ShapeFactorSUM / FTiempoMedSUM;
							FResMediosDPF[i][j].ShapeFactorSUM = 0.;
						} else
							FResMediosDPF[i][j].ShapeFactorMED = FShapeFactorDiscrete;
					}
					if(FResMediosDPF[i][j].Kreg1) {
						FResMediosDPF[i][j].Kreg1MED = FResMediosDPF[i][j].Kreg1SUM / FTiempoMedSUM;
						FResMediosDPF[i][j].Kreg1SUM = 0.;
					}
					if(FResMediosDPF[i][j].Kreg2) {
						FResMediosDPF[i][j].Kreg2MED = FResMediosDPF[i][j].Kreg2SUM / FTiempoMedSUM;
						FResMediosDPF[i][j].Kreg2SUM = 0.;
					}
					if(FResMediosDPF[i][j].Qreg) {
						FResMediosDPF[i][j].QregMED = FResMediosDPF[i][j].QregSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].QregSUM = 0.;
						/* }if(FResMediosDPF[i][j].Q1){
						 FResMediosDPF[i][j].Q1MED=FResMediosDPF[i][j].Q1SUM/FTiempoMedSUM;
						 FResMediosDPF[i][j].Q1SUM=0.;
						 }if(FResMediosDPF[i][j].Q2){
						 FResMediosDPF[i][j].Q2MED=FResMediosDPF[i][j].Q2SUM/FTiempoMedSUM;
						 FResMediosDPF[i][j].Q2SUM=0.; */
					}
					if(FResMediosDPF[i][j].TemperaturaParedCE) {
						FResMediosDPF[i][j].TemperaturaParedCEMED = FResMediosDPF[i][j].TemperaturaParedCESUM / FTiempoMedSUM;
						FResMediosDPF[i][j].TemperaturaParedCESUM = 0.;
					}
					if(FResMediosDPF[i][j].TemperaturaIntermediaPared) {
						FResMediosDPF[i][j].TemperaturaIntermediaParedMED = FResMediosDPF[i][j].TemperaturaIntermediaParedSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].TemperaturaIntermediaParedSUM = 0.;
					}
					if(FResMediosDPF[i][j].TemperaturaParedCS) {
						FResMediosDPF[i][j].TemperaturaParedCSMED = FResMediosDPF[i][j].TemperaturaParedCSSUM / FTiempoMedSUM;
						FResMediosDPF[i][j].TemperaturaParedCSSUM = 0.;
					}
					if(j + 1 == FNumeroHacesCanales) {
						if(FResMediosDPF[i][j].TemperaturaExternaSuperficie) {
							FResMediosDPF[i][j].TemperaturaExternaSuperficieMED = FResMediosDPF[i][j].TemperaturaExternaSuperficieSUM /
									FTiempoMedSUM;
							FResMediosDPF[i][j].TemperaturaExternaSuperficieSUM = 0.;
						}
						if(FResMediosDPF[i][j].TemperaturaMediaSuperficie) {
							FResMediosDPF[i][j].TemperaturaMediaSuperficieMED = FResMediosDPF[i][j].TemperaturaMediaSuperficieSUM / FTiempoMedSUM;
							FResMediosDPF[i][j].TemperaturaMediaSuperficieSUM = 0.;
						}
						if(FResMediosDPF[i][j].TemperaturaInternaSuperficie) {
							FResMediosDPF[i][j].TemperaturaInternaSuperficieMED = FResMediosDPF[i][j].TemperaturaInternaSuperficieSUM /
									FTiempoMedSUM;
							FResMediosDPF[i][j].TemperaturaInternaSuperficieSUM = 0.;
						}
					}
					/*if(FResMediosDPF[i][j].TasaFraccionMasicaEspecies){
					 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
					 FResMediosDPF[i][j].TasaFraccionMED[k]=FResMediosDPF[i][j].TasaFraccionSUM[k]/FTiempoMedSUM;
					 FResMediosDPF[i][j].TasaFraccionSUM[k]=0.;
					 }
					 }*/
					if(FResMediosDPF[i][j].FraccionMasicaEspeciesSalida) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							FResMediosDPF[i][j].FraccionSalidaMED[k] = FResMediosDPF[i][j].FraccionSalidaSUM[k] / FTiempoMedSUM;
							FResMediosDPF[i][j].FraccionSalidaSUM[k] = 0.;
						}
					}
					if(FResMediosDPF[i][j].EficienciaBrown) {
						FResMediosDPF[i][j].EficienciaBrownMED = FResMediosDPF[i][j].EficienciaBrownSUM / FResMediosDPF[i][j].GastoParedPondSUM;
						FResMediosDPF[i][j].EficienciaBrownSUM = 0.;
					}
					if(FResMediosDPF[i][j].EficienciaInter) {
						FResMediosDPF[i][j].EficienciaInterMED = FResMediosDPF[i][j].EficienciaInterSUM / FResMediosDPF[i][j].GastoParedPondSUM;
						FResMediosDPF[i][j].EficienciaInterSUM = 0.;
					}
					if(FResMediosDPF[i][j].EficienciaIner) {
						FResMediosDPF[i][j].EficienciaInerMED = FResMediosDPF[i][j].EficienciaInerSUM / FResMediosDPF[i][j].GastoParedPondSUM;
						FResMediosDPF[i][j].EficienciaInerSUM = 0.;
					}
					if(FResMediosDPF[i][j].EficienciaPL) {
						FResMediosDPF[i][j].EficienciaPLMED = FResMediosDPF[i][j].EficienciaPLSUM / FResMediosDPF[i][j].GastoParedPondSUM;
						FResMediosDPF[i][j].EficienciaPLSUM = 0.;
						FResMediosDPF[i][j].GastoParedPondSUM = 0.;
					}
				}
			}
			FTiempoMedSUM = 0.;
			FControlResMed = FControlResMed + 1.;
		}

		if(FNumResMediosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->CalculaResultadosMedios(theta);
			}
		}
		if(FNumResMediosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->CalculaResultadosMedios(theta);
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculaResultadosMedios en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::LeeResultadosInstantaneosDPF(const char *FileWAM, fpos_t &filepos) {
	int NumVars = 0, TipoVar = 0;
	try {
		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d ", &FNumResInstantaneos);
		FResInstantaneosDPF = new stResInstantDPF*[FNumResInstantaneos];
		for(int i = 0; i < FNumResInstantaneos; i++) {
			FResInstantaneosDPF[i] = new stResInstantDPF[FNumeroHacesCanales];
		}

		for(int i = 0; i < FNumResInstantaneos; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada = false;
				FResInstantaneosDPF[i][j].VelocidadParedCanalEntradaINS = 0.;
				FResInstantaneosDPF[i][j].VelocidadParedCanalSalida = false;
				FResInstantaneosDPF[i][j].VelocidadParedCanalSalidaINS = 0.;
				FResInstantaneosDPF[i][j].DPFSootMass = false;
				FResInstantaneosDPF[i][j].DPFSootMassINS = 0.;
				FResInstantaneosDPF[i][j].BeamSootMass = false;
				FResInstantaneosDPF[i][j].BeamSootMassINS = 0.;
				FResInstantaneosDPF[i][j].CVSootMass = false;
				FResInstantaneosDPF[i][j].CVSootMassINS = 0.;
				FResInstantaneosDPF[i][j].WallSootMass = false;
				FResInstantaneosDPF[i][j].WallSootMassINS = 0.;
				FResInstantaneosDPF[i][j].LayerSootMass = false;
				FResInstantaneosDPF[i][j].LayerSootMassINS = 0.;

				FResInstantaneosDPF[i][j].LayerSootMassDep = false;
				FResInstantaneosDPF[i][j].LayerSootMassDepINS = 0.;

				FResInstantaneosDPF[i][j].EspesorSootIn = false;
				FResInstantaneosDPF[i][j].EspesorSootInINS = 0.;

				FResInstantaneosDPF[i][j].EspesorSoot = false;
				FResInstantaneosDPF[i][j].EspesorSootINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaParedCE = false;
				FResInstantaneosDPF[i][j].TemperaturaParedCEINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaIntermediaPared = false;
				FResInstantaneosDPF[i][j].TemperaturaIntermediaParedINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaParedCS = false;
				FResInstantaneosDPF[i][j].TemperaturaParedCSINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie = false;
				FResInstantaneosDPF[i][j].TemperaturaExternaSuperficieINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie = false;
				FResInstantaneosDPF[i][j].TemperaturaMediaSuperficieINS = 0.;
				FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie = false;
				FResInstantaneosDPF[i][j].TemperaturaInternaSuperficieINS = 0.;
				FResInstantaneosDPF[i][j].KwallClean = false;
				FResInstantaneosDPF[i][j].KwallCleanINS = 0.;
				FResInstantaneosDPF[i][j].KwallLoaded = false;
				FResInstantaneosDPF[i][j].KwallLoadedINS = 0.;
				FResInstantaneosDPF[i][j].Kwall = false;
				FResInstantaneosDPF[i][j].KwallINS = 0.;

				FResInstantaneosDPF[i][j].KsootIn = false;
				FResInstantaneosDPF[i][j].KsootInINS = 0.;
				FResInstantaneosDPF[i][j].KsootDep = false;
				FResInstantaneosDPF[i][j].KsootDepINS = 0.;

				FResInstantaneosDPF[i][j].Ksoot = false;
				FResInstantaneosDPF[i][j].KsootINS = 0.;
				FResInstantaneosDPF[i][j].Eficiencia = false;
				FResInstantaneosDPF[i][j].EficienciaINS = 0.;
				FResInstantaneosDPF[i][j].Porosidad = false;
				FResInstantaneosDPF[i][j].PorosidadINS = 0.;
				FResInstantaneosDPF[i][j].CoeficienteParticion = false;
				FResInstantaneosDPF[i][j].CoeficienteParticionINS = 0.;
				FResInstantaneosDPF[i][j].DiametroUC = false;
				FResInstantaneosDPF[i][j].DiametroUCINS = 0.;
				FResInstantaneosDPF[i][j].Kreg1 = false;
				FResInstantaneosDPF[i][j].Kreg1INS = 0.;
				FResInstantaneosDPF[i][j].Kreg2 = false;
				FResInstantaneosDPF[i][j].Kreg2INS = 0.;
				FResInstantaneosDPF[i][j].Qreg = false;
				FResInstantaneosDPF[i][j].QregINS = 0.;
				FResInstantaneosDPF[i][j].Q1 = false;
				FResInstantaneosDPF[i][j].Q1INS = 0.;
				FResInstantaneosDPF[i][j].Q2 = false;
				FResInstantaneosDPF[i][j].Q2INS = 0.;
				FResInstantaneosDPF[i][j].TasaFraccionMasicaEspecies = false;
				FResInstantaneosDPF[i][j].TasaFraccionINS = new double[FNumeroEspecies - FIntEGR];
				FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida = false;
				FResInstantaneosDPF[i][j].FraccionSalidaINS = new double[FNumeroEspecies - FIntEGR];
				FResInstantaneosDPF[i][j].ShapeFactor = false;
				FResInstantaneosDPF[i][j].ShapeFactorINS = 0.;
				for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
					FResInstantaneosDPF[i][j].TasaFraccionINS[k] = 0.;
					FResInstantaneosDPF[i][j].FraccionSalidaINS[k] = 0.;
				}
				FResInstantaneosDPF[i][j].EficienciaBrown = false;
				FResInstantaneosDPF[i][j].EficienciaBrownINS = 0.;
				FResInstantaneosDPF[i][j].EficienciaInter = false;
				FResInstantaneosDPF[i][j].EficienciaInterINS = 0.;
				FResInstantaneosDPF[i][j].EficienciaIner = false;
				FResInstantaneosDPF[i][j].EficienciaInerINS = 0.;
				FResInstantaneosDPF[i][j].EficienciaPL = false;
				FResInstantaneosDPF[i][j].EficienciaPLINS = 0.;

				fscanf(fich, "%lf %d ", &FResInstantaneosDPF[i][j].Distancia, &NumVars);

				for(int k = 0; k < NumVars; k++) {
					fscanf(fich, "%d ", &TipoVar);
					switch(TipoVar) {
					case 0:
						FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada = true;
						break;
					case 1:
						FResInstantaneosDPF[i][j].VelocidadParedCanalSalida = true;
						break;
					case 2:
						FResInstantaneosDPF[i][j].DPFSootMass = true;
						break;
					case 3:
						FResInstantaneosDPF[i][j].BeamSootMass = true;
						break;
					case 4:
						FResInstantaneosDPF[i][j].CVSootMass = true;
						break;
					case 5:
						FResInstantaneosDPF[i][j].WallSootMass = true;
						break;
					case 6:
						FResInstantaneosDPF[i][j].LayerSootMass = true;
						break;
					case 7:
						FResInstantaneosDPF[i][j].EspesorSoot = true;
						break;
					case 8:
						FResInstantaneosDPF[i][j].Kwall = true;
						break;
					case 9:
						FResInstantaneosDPF[i][j].Ksoot = true;
						break;
					case 10:
						FResInstantaneosDPF[i][j].Eficiencia = true;
						break;
					case 11:
						FResInstantaneosDPF[i][j].Porosidad = true;
						break;
					case 12:
						FResInstantaneosDPF[i][j].CoeficienteParticion = true;
						break;
					case 13:
						FResInstantaneosDPF[i][j].DiametroUC = true;
						break;
					case 14:
						FResInstantaneosDPF[i][j].ShapeFactor = true;
						break;
					case 15:
						FResInstantaneosDPF[i][j].Kreg1 = true;
						break;
					case 16:
						FResInstantaneosDPF[i][j].Kreg2 = true;
						break;
					case 17:
						FResInstantaneosDPF[i][j].Qreg = true;
						break;
					case 18:
						FResInstantaneosDPF[i][j].TemperaturaParedCE = true;
						break;
					case 19:
						FResInstantaneosDPF[i][j].TemperaturaIntermediaPared = true;
						break;
					case 20:
						FResInstantaneosDPF[i][j].TemperaturaParedCS = true;
						break;
					case 21:
						FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie = true;
						break;
					case 22:
						FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie = true;
						break;
					case 23:
						FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie = true;
						break;
					case 24:
						FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida = true;
						break;
					case 25:
						FResInstantaneosDPF[i][j].EficienciaBrown = true;
						break;
					case 26:
						FResInstantaneosDPF[i][j].EficienciaInter = true;
						break;
					case 27:
						FResInstantaneosDPF[i][j].EficienciaIner = true;
						break;
					case 28:
						FResInstantaneosDPF[i][j].EficienciaPL = true;
						break;
					case 29:
						FResInstantaneosDPF[i][j].KwallClean = true;
						break;
					case 30:
						FResInstantaneosDPF[i][j].KwallLoaded = true;
						break;

					case 31:
						FResInstantaneosDPF[i][j].KsootIn = true;
						break;
					case 32:
						FResInstantaneosDPF[i][j].KsootDep = true;
						break;
					case 33:
						FResInstantaneosDPF[i][j].EspesorSootIn = true;
						break;
					case 34:
						FResInstantaneosDPF[i][j].LayerSootMassDep = true;
						break;

					default:
						std::cout << "Resultados instantaneos en la DPF " << FNumeroDPF << " no implementados " << std::endl;
					}
				}
			}
		}

		fscanf(fich, "%d %d", &FNumResInstantaneosCE, &FNumResInstantaneosCS);

		fgetpos(fich, &filepos);
		fclose(fich);

		if(FNumResInstantaneosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->LeeResultadosInstantaneosCanalDPF(FNumResInstantaneosCE, FileWAM, filepos);
			}
		}
		if(FNumResInstantaneosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->LeeResultadosInstantaneosCanalDPF(FNumResInstantaneosCS, FileWAM, filepos);
			}
		}

//fgetpos(fich, &filepos);
//fclose(fich);
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::LeeResultadosInstantaneosDPF en la DPF: " << FNumeroDPF << std::endl;

		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CabeceraResultadosInstantaneos(stringstream& insoutput, stEspecies *DatosEspecies) const {
	try {

//FILE *fich=fopen(FileSALIDA,"a");
//float Dist;
		std::string Label;
		std::ostringstream TextDist;
		TextDist.precision(8);

		for(int i = 0; i < FNumResInstantaneos; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				TextDist << FResInstantaneosDPF[i][j].Distancia;
				if(FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada) {
					Label = "\t" + PutLabel(802) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(909);
					insoutput << Label.c_str();
					//fprintf(fich,"\tV_PCE_DPF_%d_Haz_%d_a_%5.3f_m(m/s)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].VelocidadParedCanalSalida) {
					Label = "\t" + PutLabel(803) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(909);
					insoutput << Label.c_str();
					//fprintf(fich,"\tV_PCS_DPF_%d_Haz_%d_a_%5.3f_m(m/s)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].DPFSootMass) {
					Label = "\t" + PutLabel(804) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tSoot_mass_DPF_%d_(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].BeamSootMass) {
					Label = "\t" + PutLabel(804) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tSoot_mass_DPF_%d_Haz_%d(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].CVSootMass) {
					Label = "\t" + PutLabel(836) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tCV_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].WallSootMass) {
					Label = "\t" + PutLabel(834) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tWall_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].LayerSootMass) {
					Label = "\t" + PutLabel(835) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tLayer_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].LayerSootMassDep) {
					Label = "\t" + PutLabel(846) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(919);
					insoutput << Label.c_str();
					//fprintf(fich,"\tLayer_soot_mass_DPF_%d_Haz_%d_a_%5.3f_m(kg)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].EspesorSootIn) {
					Label = "\t" + PutLabel(843) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEspesor_Soot_DPF_Haz_%d_a_%5.3f_m(mm)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].EspesorSoot) {
					Label = "\t" + PutLabel(805) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEspesor_Soot_DPF_Haz_%d_a_%5.3f_m(mm)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].KwallClean) {
					Label = "\t" + PutLabel(841) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].KwallLoaded) {
					Label = "\t" + PutLabel(842) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].Kwall) {
					Label = "\t" + PutLabel(806) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_wall_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].KsootIn) {
					Label = "\t" + PutLabel(844) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].KsootDep) {
					Label = "\t" + PutLabel(845) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}

				if(FResInstantaneosDPF[i][j].Ksoot) {
					Label = "\t" + PutLabel(807) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(915);
					insoutput << Label.c_str();
					//fprintf(fich,"\tK_soot_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].Eficiencia) {
					Label = "\t" + PutLabel(808) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].EficienciaBrown) {
					Label = "\t" + PutLabel(837) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].EficienciaInter) {
					Label = "\t" + PutLabel(838) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].EficienciaIner) {
					Label = "\t" + PutLabel(839) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].EficienciaPL) {
					Label = "\t" + PutLabel(840) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tEficiencia_DPF_%d_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].Porosidad) {
					Label = "\t" + PutLabel(809) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tPorosidad_DPF_Haz_%d_a_%5.3f_m(-)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].CoeficienteParticion) {
					Label = "\t" + PutLabel(810) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
					//fprintf(fich,"\tCoef_Particion_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].DiametroUC) {
					Label = "\t" + PutLabel(811) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(920);
					insoutput << Label.c_str();
					//fprintf(fich,"\tDiametroUC_DPF_%d_Haz_%d_a_%5.3f_m(m-2)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].ShapeFactor) {
					Label = "\t" + PutLabel(833) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
				}
				if(FResInstantaneosDPF[i][j].Kreg1) {
					Label = "\t" + PutLabel(812) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(922);
					insoutput << Label.c_str();
					//fprintf(fich,"\tKreg_termica_DPF_%d_Haz_%d_a_%5.3f_m(1/s)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].Kreg2) {
					Label = "\t" + PutLabel(813) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(922);
					insoutput << Label.c_str();
					//fprintf(fich,"\tKreg_cata_DPF_Haz_%d_a_%5.3f_m(1/s)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].Qreg) {
					Label = "\t" + PutLabel(814) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(907);
					insoutput << Label.c_str();
					//fprintf(fich,"\tQreg_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
					/*if(FResInstantaneosDPF[i][j].Q1)
					 fprintf(fich,"\tQ_HC_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
					 if(FResInstantaneosDPF[i][j].Q2)
					 fprintf(fich,"\tQ_CO_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);*/
				}
				if(FResInstantaneosDPF[i][j].TemperaturaParedCE) {
					Label = "\t" + PutLabel(815) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					insoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_CE_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].TemperaturaIntermediaPared) {
					Label = "\t" + PutLabel(816) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					insoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_Intermedia_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(FResInstantaneosDPF[i][j].TemperaturaParedCS) {
					Label = "\t" + PutLabel(817) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
								j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317) + PutLabel(910);
					insoutput << Label.c_str();
					//fprintf(fich,"\tT_Pared_CS_DPF_%d_Haz_%d_a_%5.3f_m(J)",FNumeroDPF,j,FResInstantaneosDPF[i][j].Distancia);
				}
				if(j + 1 == FNumeroHacesCanales) {
					if(FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie) {
						Label = "\t" + PutLabel(818) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						insoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Ext_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResInstantaneosDPF[i][j].Distancia);
					}
					if(FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie) {
						Label = "\t" + PutLabel(819) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						insoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Med_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResInstantaneosDPF[i][j].Distancia);
					}
					if(FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie) {
						Label = "\t" + PutLabel(820) + PutLabel(800) + std::to_string(FNumeroDPF) + PutLabel(801) + std::to_string(
									j + 1) + PutLabel(316) + TextDist.str() + PutLabel(317)
								+ PutLabel(910);
						insoutput << Label.c_str();
						//fprintf(fich,"\tT_Superficial_Int_DPF_%d_a_%5.3f_m(J)",FNumeroDPF,FResInstantaneosDPF[i][j].Distancia);
					}
				}
				/*if(FResInstantaneosDPF[i][j].TasaFraccionMasicaEspecies){
				 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
				 fprintf(fich,"\tTasa_Y_%s_DPF_Haz_%d_a_%5.3f_m(1/s)",DatosEspecies[k].Nombre,FNumeroDPF,k,FResInstantaneosDPF[i][j].Distancia);
				 }
				 } */
				if(FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida) {
					for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
						Label = "\t" + PutLabel(821) + DatosEspecies[k].Nombre + PutLabel(822) + PutLabel(800) + std::to_string(
									FNumeroDPF) + PutLabel(801) + std::to_string(j + 1) + PutLabel(316)
								+ TextDist.str() + PutLabel(317) + PutLabel(901);
						insoutput << Label.c_str();
						//fprintf(fich,"\tY_%s_pared_DPF_Haz_%d_a_%5.3f_m(1/s)",DatosEspecies[k].Nombre,FNumeroDPF,k,FResInstantaneosDPF[i][j].Distancia);
					}
				}

			}
		}

//fclose(fich);

		if(FNumResInstantaneosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->CabeceraResultadosInstantaneos(insoutput, DatosEspecies);
			}
		}
		if(FNumResInstantaneosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->CabeceraResultadosInstantaneos(insoutput, DatosEspecies);
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CabeceraResultadosInstantaneos en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::ImprimeResultadosInstantaneos(stringstream& insoutput) const {
	try {

//FILE *fich=fopen(FileSALIDA,"a");

		for(int i = 0; i < FNumResInstantaneos; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				if(FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada)
					insoutput << "\t" << FResInstantaneosDPF[i][j].VelocidadParedCanalEntradaINS;
				if(FResInstantaneosDPF[i][j].VelocidadParedCanalSalida)
					insoutput << "\t" << FResInstantaneosDPF[i][j].VelocidadParedCanalSalidaINS;
				if(FResInstantaneosDPF[i][j].DPFSootMass)
					insoutput << "\t" << FResInstantaneosDPF[i][j].DPFSootMassINS;
				if(FResInstantaneosDPF[i][j].BeamSootMass)
					insoutput << "\t" << FResInstantaneosDPF[i][j].BeamSootMassINS;
				if(FResInstantaneosDPF[i][j].CVSootMass)
					insoutput << "\t" << FResInstantaneosDPF[i][j].CVSootMassINS;
				if(FResInstantaneosDPF[i][j].WallSootMass)
					insoutput << "\t" << FResInstantaneosDPF[i][j].WallSootMassINS;
				if(FResInstantaneosDPF[i][j].LayerSootMass)
					insoutput << "\t" << FResInstantaneosDPF[i][j].LayerSootMassINS;

				if(FResInstantaneosDPF[i][j].LayerSootMassDep)
					insoutput << "\t" << FResInstantaneosDPF[i][j].LayerSootMassDepINS;

				if(FResInstantaneosDPF[i][j].EspesorSootIn)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EspesorSootInINS;

				if(FResInstantaneosDPF[i][j].EspesorSoot)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EspesorSootINS;
				if(FResInstantaneosDPF[i][j].KwallClean)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KwallCleanINS;
				if(FResInstantaneosDPF[i][j].KwallLoaded)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KwallLoadedINS;
				if(FResInstantaneosDPF[i][j].Kwall)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KwallINS;

				if(FResInstantaneosDPF[i][j].KsootIn)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KsootInINS;
				if(FResInstantaneosDPF[i][j].KsootDep)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KsootDepINS;

				if(FResInstantaneosDPF[i][j].Ksoot)
					insoutput << "\t" << FResInstantaneosDPF[i][j].KsootINS;
				if(FResInstantaneosDPF[i][j].Eficiencia)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EficienciaINS;
				if(FResInstantaneosDPF[i][j].EficienciaBrown)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EficienciaBrownINS;
				if(FResInstantaneosDPF[i][j].EficienciaInter)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EficienciaInterINS;
				if(FResInstantaneosDPF[i][j].EficienciaIner)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EficienciaInerINS;
				if(FResInstantaneosDPF[i][j].EficienciaPL)
					insoutput << "\t" << FResInstantaneosDPF[i][j].EficienciaPLINS;
				if(FResInstantaneosDPF[i][j].Porosidad)
					insoutput << "\t" << FResInstantaneosDPF[i][j].PorosidadINS;
				if(FResInstantaneosDPF[i][j].CoeficienteParticion)
					insoutput << "\t" << FResInstantaneosDPF[i][j].CoeficienteParticionINS;
				if(FResInstantaneosDPF[i][j].DiametroUC)
					insoutput << "\t" << FResInstantaneosDPF[i][j].DiametroUCINS;
				if(FResInstantaneosDPF[i][j].ShapeFactor)
					insoutput << "\t" << FResInstantaneosDPF[i][j].ShapeFactorINS;
				if(FResInstantaneosDPF[i][j].Kreg1)
					insoutput << "\t" << FResInstantaneosDPF[i][j].Kreg1INS;
				if(FResInstantaneosDPF[i][j].Kreg2)
					insoutput << "\t" << FResInstantaneosDPF[i][j].Kreg2INS;
				if(FResInstantaneosDPF[i][j].Qreg)
					insoutput << "\t" << FResInstantaneosDPF[i][j].QregINS;
				/*if(FResInstantaneosDPF[i][j].Q1)
				 fprintf(fich,"\t%lg",FResInstantaneosDPF[i][j].Q1INS);
				 if(FResInstantaneosDPF[i][j].Q2)
				 fprintf(fich,"\t%lg",FResInstantaneosDPF[i][j].Q2INS);*/
				if(FResInstantaneosDPF[i][j].TemperaturaParedCE)
					insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaParedCEINS;
				if(FResInstantaneosDPF[i][j].TemperaturaIntermediaPared)
					insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaIntermediaParedINS;
				if(FResInstantaneosDPF[i][j].TemperaturaParedCS)
					insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaParedCSINS;
				if(j + 1 == FNumeroHacesCanales) {
					if(FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie)
						insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaExternaSuperficieINS;
					if(FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie)
						insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaMediaSuperficieINS;
					if(FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie)
						insoutput << "\t" << FResInstantaneosDPF[i][j].TemperaturaInternaSuperficieINS;
				}
				/*if(FResInstantaneosDPF[i][j].TasaFraccionMasicaEspecies){
				 for(int k=0;k<FNumeroEspecies-FIntEGR;k++){
				 fprintf(fich,"\t%lg",FResInstantaneosDPF[i][j].TasaFraccionINS[k]);
				 }
				 }*/
				if(FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida) {
					for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
						insoutput << "\t" << FResInstantaneosDPF[i][j].FraccionSalidaINS[k];
					}
				}

			}
		}

		if(FNumResInstantaneosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->ImprimeResultadosInstantaneos(insoutput);
			}
		}
		if(FNumResInstantaneosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->ImprimeResultadosInstantaneos(insoutput);
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::ImprimeResultadosInstantaneos en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculaResultadosInstantaneos() {
	int n1 = 0, n2 = 0;
	double dist = 0., d = 0., Vble = 0.;

	try {
		for(int i = 0; i < FNumResInstantaneos; i++) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				dist = FResInstantaneosDPF[i][j].Distancia / FCanal[j][0]->getXRef();
				n1 = (int) floor(dist);
				if(n1 >= FCanal[j][0]->getNin() - 1) {
					if(FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada)
						FResInstantaneosDPF[i][j].VelocidadParedCanalEntradaINS = FVelocidadPared[j][0][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].VelocidadParedCanalSalida)
						FResInstantaneosDPF[i][j].VelocidadParedCanalSalidaINS = FVelocidadPared[j][1][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].DPFSootMass)
						FResInstantaneosDPF[i][j].DPFSootMassINS = FDPFSootMass * 1000;
					if(FResInstantaneosDPF[i][j].BeamSootMass)
						FResInstantaneosDPF[i][j].BeamSootMassINS = FBeamSootMass[j] * 1000;
					if(FResInstantaneosDPF[i][j].CVSootMass)
						FResInstantaneosDPF[i][j].CVSootMassINS = FCVSootMass[j][FCanal[j][0]->getNin() - 1] * 1000;
					if(FResInstantaneosDPF[i][j].WallSootMass)
						FResInstantaneosDPF[i][j].WallSootMassINS = FWallSootMass[j][FCanal[j][0]->getNin() - 1] * 1000;
					if(FResInstantaneosDPF[i][j].LayerSootMass)
						FResInstantaneosDPF[i][j].LayerSootMassINS = FLayerSootMass[j][FCanal[j][0]->getNin() - 1] * 1000;

					if(FResInstantaneosDPF[i][j].LayerSootMassDep)
						FResInstantaneosDPF[i][j].LayerSootMassDepINS = FLayerSootMassDep[j][FCanal[j][0]->getNin() - 1] * 1000;

					if(FResInstantaneosDPF[i][j].EspesorSootIn)
						FResInstantaneosDPF[i][j].EspesorSootInINS = FEspesorSootIn[j][FCanal[j][0]->getNin() - 1] * 1000;

					if(FResInstantaneosDPF[i][j].EspesorSoot)
						FResInstantaneosDPF[i][j].EspesorSootINS = FEspesorSoot[j][FCanal[j][0]->getNin() - 1] * 1000;
					if(FResInstantaneosDPF[i][j].KwallClean)
						FResInstantaneosDPF[i][j].KwallCleanINS = FKwallClean[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].KwallLoaded)
						FResInstantaneosDPF[i][j].KwallLoadedINS = FKwallLoaded[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Kwall)
						FResInstantaneosDPF[i][j].KwallINS = FKwall[j][FCanal[j][0]->getNin() - 1];

					if(FResInstantaneosDPF[i][j].KsootIn)
						FResInstantaneosDPF[i][j].KsootInINS = FKsootIn[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].KsootDep)
						FResInstantaneosDPF[i][j].KsootDepINS = FKsootDep[j][FCanal[j][0]->getNin() - 1];

					if(FResInstantaneosDPF[i][j].Ksoot)
						FResInstantaneosDPF[i][j].KsootINS = FKsoot[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Eficiencia)
						FResInstantaneosDPF[i][j].EficienciaINS = FEficiencia[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Porosidad)
						FResInstantaneosDPF[i][j].PorosidadINS = FPorosidad[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].CoeficienteParticion)
						FResInstantaneosDPF[i][j].CoeficienteParticionINS = FCoeficienteParticion[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].DiametroUC)
						FResInstantaneosDPF[i][j].DiametroUCINS = FDiametroUnidadColectora[j][FCanal[j][0]->getNin() - 1] * 1000;
					if(FResInstantaneosDPF[i][j].ShapeFactor)
						if(FCalculoFiltrado) {
							FResInstantaneosDPF[i][j].ShapeFactorINS = FShapeFactor[j][FCanal[j][0]->getNin() - 1];
						} else
							FResInstantaneosDPF[i][j].ShapeFactorINS = FShapeFactorDiscrete;
					if(FResInstantaneosDPF[i][j].Kreg1)
						FResInstantaneosDPF[i][j].Kreg1INS = FKreg1[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Kreg2)
						FResInstantaneosDPF[i][j].Kreg2INS = FKreg2[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Qreg)
						FResInstantaneosDPF[i][j].QregINS = FQreg[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Q1)
						FResInstantaneosDPF[i][j].Q1INS = FQ1[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].Q2)
						FResInstantaneosDPF[i][j].Q2INS = FQ2[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].TemperaturaParedCE)
						FResInstantaneosDPF[i][j].TemperaturaParedCEINS = FTPared[j][FCanal[j][0]->getNin() - 1][2];
					if(FResInstantaneosDPF[i][j].TemperaturaIntermediaPared)
						FResInstantaneosDPF[i][j].TemperaturaIntermediaParedINS = FTPared[j][FCanal[j][1]->getNin() - 1][1];
					if(FResInstantaneosDPF[i][j].TemperaturaParedCS)
						FResInstantaneosDPF[i][j].TemperaturaParedCSINS = FTPared[j][FCanal[j][0]->getNin() - 1][0];
					if(j + 1 == FNumeroHacesCanales) {
						if(FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie)
							FResInstantaneosDPF[i][j].TemperaturaExternaSuperficieINS = FTSuperficie[FCanal[j][0]->getNin() - 1][2];
						if(FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie)
							FResInstantaneosDPF[i][j].TemperaturaMediaSuperficieINS = FTSuperficie[FCanal[j][0]->getNin() - 1][1];
						if(FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie)
							FResInstantaneosDPF[i][j].TemperaturaInternaSuperficieINS = FTSuperficie[FCanal[j][0]->getNin() - 1][0];
					}
					if(FResInstantaneosDPF[i][j].TasaFraccionMasicaEspecies) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							FResInstantaneosDPF[i][j].TasaFraccionINS[k] = FTasaFraccionMasicaEspecie[j][FCanal[j][0]->getNin() - 1][k];
						}
					}
					if(FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							FResInstantaneosDPF[i][j].FraccionSalidaINS[k] = FFraccionMasicaEspecieSalida[j][FCanal[j][0]->getNin() - 1][k];
						}
					}
					if(FResInstantaneosDPF[i][j].EficienciaBrown)
						FResInstantaneosDPF[i][j].EficienciaBrownINS = FEficienciaBrown[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].EficienciaInter)
						FResInstantaneosDPF[i][j].EficienciaInterINS = FEficienciaInter[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].EficienciaIner)
						FResInstantaneosDPF[i][j].EficienciaInerINS = FEficienciaIner[j][FCanal[j][0]->getNin() - 1];
					if(FResInstantaneosDPF[i][j].EficienciaPL)
						FResInstantaneosDPF[i][j].EficienciaPLINS = FEficienciaPL[j][FCanal[j][0]->getNin() - 1];

				} else {
					n2 = n1 + 1;
					d = dist - (double) n1;
					if(FResInstantaneosDPF[i][j].VelocidadParedCanalEntrada) {
						Vble = Interpola(FVelocidadPared[j][0][n1], FVelocidadPared[j][0][n2], 1., d);
						FResInstantaneosDPF[i][j].VelocidadParedCanalEntradaINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].VelocidadParedCanalSalida) {
						Vble = Interpola(FVelocidadPared[j][1][n1], FVelocidadPared[j][1][n2], 1., d);
						FResInstantaneosDPF[i][j].VelocidadParedCanalSalidaINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].DPFSootMass) {
						FResInstantaneosDPF[i][j].DPFSootMassINS = FDPFSootMass * 1000;
					}
					if(FResInstantaneosDPF[i][j].BeamSootMass) {
						FResInstantaneosDPF[i][j].BeamSootMassINS = FBeamSootMass[j] * 1000;
					}
					if(FResInstantaneosDPF[i][j].CVSootMass) {
						Vble = Interpola(FCVSootMass[j][n1], FCVSootMass[j][n2], 1., d);
						FResInstantaneosDPF[i][j].CVSootMassINS = Vble * 1000;
					}
					if(FResInstantaneosDPF[i][j].WallSootMass) {
						//Vble=Interpola(FMasaSootLayer[j][n1],FMasaSootLayer[j][n2],1.,d);
						//FResInstantaneosDPF[i][j].MasaSootINS=Vble;
						//FResInstantaneosDPF[i][j].MasaSootINS=FMasaSootTotal;
						Vble = Interpola(FWallSootMass[j][n1], FWallSootMass[j][n2], 1., d);
						//FResInstantaneosDPF[i][j].MasaSootINS=Vble;
						FResInstantaneosDPF[i][j].WallSootMassINS = Vble * 1000;
					}
					if(FResInstantaneosDPF[i][j].LayerSootMass) {
						Vble = Interpola(FLayerSootMass[j][n1], FLayerSootMass[j][n2], 1., d);
						FResInstantaneosDPF[i][j].LayerSootMassINS = Vble * 1000;

					}
					if(FResInstantaneosDPF[i][j].LayerSootMassDep) {
						Vble = Interpola(FLayerSootMassDep[j][n1], FLayerSootMassDep[j][n2], 1., d);
						FResInstantaneosDPF[i][j].LayerSootMassDepINS = Vble * 1000;

					}
					if(FResInstantaneosDPF[i][j].EspesorSootIn) {
						Vble = Interpola(FEspesorSootIn[j][n1], FEspesorSootIn[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EspesorSootInINS = Vble * 1000;

					}
					if(FResInstantaneosDPF[i][j].EspesorSoot) {
						Vble = Interpola(FEspesorSoot[j][n1], FEspesorSoot[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EspesorSootINS = Vble * 1000;
					}
					if(FResInstantaneosDPF[i][j].KwallClean) {
						Vble = Interpola(FKwallClean[j][n1], FKwallClean[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KwallCleanINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].KwallLoaded) {
						Vble = Interpola(FKwallLoaded[j][n1], FKwallLoaded[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KwallLoadedINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Kwall) {
						Vble = Interpola(FKwall[j][n1], FKwall[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KwallINS = Vble;

					}
					if(FResInstantaneosDPF[i][j].KsootIn) {
						Vble = Interpola(FKsootIn[j][n1], FKsootIn[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KsootInINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].KsootDep) {
						Vble = Interpola(FKsootDep[j][n1], FKsootDep[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KsootDepINS = Vble;

					}
					if(FResInstantaneosDPF[i][j].Ksoot) {
						Vble = Interpola(FKsoot[j][n1], FKsoot[j][n2], 1., d);
						FResInstantaneosDPF[i][j].KsootINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Eficiencia) {
						Vble = Interpola(FEficiencia[j][n1], FEficiencia[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EficienciaINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Porosidad) {
						Vble = Interpola(FPorosidad[j][n1], FPorosidad[j][n2], 1., d);
						FResInstantaneosDPF[i][j].PorosidadINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].CoeficienteParticion) {
						Vble = Interpola(FCoeficienteParticion[j][n1], FCoeficienteParticion[j][n2], 1., d);
						FResInstantaneosDPF[i][j].CoeficienteParticionINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].DiametroUC) {
						Vble = Interpola(FDiametroUnidadColectora[j][n1], FDiametroUnidadColectora[j][n2], 1., d);
						FResInstantaneosDPF[i][j].DiametroUCINS = Vble * 1000.;
					}
					if(FResInstantaneosDPF[i][j].ShapeFactor) {
						if(FCalculoFiltrado) {
							Vble = Interpola(FShapeFactor[j][n1], FShapeFactor[j][n2], 1., d);
							FResInstantaneosDPF[i][j].ShapeFactorINS = Vble;
						} else
							FResInstantaneosDPF[i][j].ShapeFactorINS = FShapeFactorDiscrete;
					}
					if(FResInstantaneosDPF[i][j].Kreg1) {
						Vble = Interpola(FKreg1[j][n1], FKreg1[j][n2], 1., d);
						FResInstantaneosDPF[i][j].Kreg1INS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Kreg2) {
						Vble = Interpola(FKreg2[j][n1], FKreg2[j][n2], 1., d);
						FResInstantaneosDPF[i][j].Kreg2INS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Qreg) {
						Vble = Interpola(FQreg[j][n1], FQreg[j][n2], 1., d);
						FResInstantaneosDPF[i][j].QregINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Q1) {
						Vble = Interpola(FQ1[j][n1], FQ1[j][n2], 1., d);
						FResInstantaneosDPF[i][j].Q1INS = Vble;
					}
					if(FResInstantaneosDPF[i][j].Q2) {
						Vble = Interpola(FQ2[j][n1], FQ2[j][n2], 1., d);
						FResInstantaneosDPF[i][j].Q2INS = Vble;
					}
					if(FResInstantaneosDPF[i][j].TemperaturaParedCE) {
						Vble = Interpola(FTPared[j][n1][2], FTPared[j][n2][2], 1., d);
						FResInstantaneosDPF[i][j].TemperaturaParedCEINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].TemperaturaIntermediaPared) {
						Vble = Interpola(FTPared[j][n1][1], FTPared[j][n2][1], 1., d);
						FResInstantaneosDPF[i][j].TemperaturaIntermediaParedINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].TemperaturaParedCS) {
						Vble = Interpola(FTPared[j][n1][0], FTPared[j][n2][0], 1., d);
						FResInstantaneosDPF[i][j].TemperaturaParedCSINS = Vble;
					}
					if(j + 1 == FNumeroHacesCanales) {
						if(FResInstantaneosDPF[i][j].TemperaturaExternaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][2], FTSuperficie[n2][2], 1., d);
							FResInstantaneosDPF[i][j].TemperaturaExternaSuperficieINS = Vble;
						}
						if(FResInstantaneosDPF[i][j].TemperaturaMediaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][1], FTSuperficie[n2][1], 1., d);
							FResInstantaneosDPF[i][j].TemperaturaMediaSuperficieINS = Vble;
						}
						if(FResInstantaneosDPF[i][j].TemperaturaInternaSuperficie) {
							Vble = Interpola(FTSuperficie[n1][0], FTSuperficie[n2][0], 1., d);
							FResInstantaneosDPF[i][j].TemperaturaInternaSuperficieINS = Vble;
						}
					}
					if(FResInstantaneosDPF[i][j].TasaFraccionMasicaEspecies) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							Vble = Interpola(FTasaFraccionMasicaEspecie[j][n1][k], FTasaFraccionMasicaEspecie[j][n2][k], 1., d);
							FResInstantaneosDPF[i][j].TasaFraccionINS[k] = Vble;
						}
					}
					if(FResInstantaneosDPF[i][j].FraccionMasicaEspeciesSalida) {
						for(int k = 0; k < FNumeroEspecies - FIntEGR; k++) {
							Vble = Interpola(FFraccionMasicaEspecieSalida[j][n1][k], FFraccionMasicaEspecieSalida[j][n2][k], 1., d);
							FResInstantaneosDPF[i][j].FraccionSalidaINS[k] = Vble;
						}
					}
					if(FResInstantaneosDPF[i][j].EficienciaBrown) {
						Vble = Interpola(FEficienciaBrown[j][n1], FEficienciaBrown[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EficienciaBrownINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].EficienciaInter) {
						Vble = Interpola(FEficienciaInter[j][n1], FEficienciaInter[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EficienciaInterINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].EficienciaIner) {
						Vble = Interpola(FEficienciaIner[j][n1], FEficienciaIner[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EficienciaInerINS = Vble;
					}
					if(FResInstantaneosDPF[i][j].EficienciaPL) {
						Vble = Interpola(FEficienciaPL[j][n1], FEficienciaPL[j][n2], 1., d);
						FResInstantaneosDPF[i][j].EficienciaPLINS = Vble;
					}
				}
			}
		}

		if(FNumResInstantaneosCE != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][0]->CalculaResultadosInstantaneos();
			}
		}
		if(FNumResInstantaneosCS != 0) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				FCanal[j][1]->CalculaResultadosInstantaneos();
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculaResultadosInstantaneos en la DPF: " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TDPF::Interpola(double vizq, double vder, double axid, double xif) {
	double xx = 0., yy = 0.;
	double ret_val = 0.;

	try {

		xx = vder - vizq;
		if(axid != 0.) {
			yy = (xx / axid) * xif;
			ret_val = vizq + yy;
		} else {
			std::cout << "ERROR: valores entrada Interpolacion en la DPF " << FNumeroDPF << std::endl;
		}
		return ret_val;
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::Interpolacion  en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception("");
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoEstabilidadDPF() {
	double TimeMinAnterior = 0.;

	try {

		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FCanal[j][0]->EstabilidadMetodoCalculo();
			FCanal[j][1]->EstabilidadMetodoCalculo();
		}

		FTime0DPF = FTime1DPF;
		FTime1DPF = FCanal[0][0]->getTime1();
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			if(FCanal[j][0]->getTime1() <= FTime1DPF) {
				FTime1DPF = FCanal[j][0]->getTime1();
			}
			if(FCanal[j][1]->getTime1() <= FTime1DPF) {
				FTime1DPF = FCanal[j][1]->getTime1();
			}
		}
		FDeltaTimeDPF = FTime1DPF - FTime0DPF;
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FCanal[j][0]->putTime1(FTime1DPF);
			FCanal[j][0]->putDeltaTime(FDeltaTimeDPF);
			FCanal[j][1]->putTime1(FTime1DPF);
			FCanal[j][1]->putDeltaTime(FDeltaTimeDPF);
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoEstabilidadDPF  en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::AjustaPaso(double TiempoFinPaso) {
	try {
		FTime1DPF = TiempoFinPaso;
		FDeltaTimeDPF = FTime1DPF - FTime0DPF;
		for(int j = 0; j < FNumeroHacesCanales; j++) {
			FCanal[j][0]->putTime1(FTime1DPF);
			FCanal[j][0]->putDeltaTime(FDeltaTimeDPF);
			FCanal[j][1]->putTime1(FTime1DPF);
			FCanal[j][1]->putDeltaTime(FDeltaTimeDPF);
		}

	} catch(exception &N) {
		std::cout << "ERROR: TDPF::AjustaPaso en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::CalculoSubmodelos() {
	try {

		CalculoKsoot();
//CalculoEspesorSoot();
		CalculoVelocidadPared();

		if(FCalculoRegeneracion || FCalculoFiltrado) {
			SubmodeloFiltrado();
			//SubmodeloRegeneracion();
		}

		if(!FCalculoRegeneracion && !FCalculoFiltrado) {
			for(int j = 0; j < FNumeroHacesCanales; j++) {
				for(int i = 0; i < FCanal[j][0]->getNin(); i++) {
					if(FCanal[j][0]->getNodoInicialFuente() <= i) {
						FFraccionMasicaEspecieSalida[j][i][0] = FCanal[j][0]->GetYespecie(i, 0); // O2
						FFraccionMasicaEspecieSalida[j][i][1] = FCanal[j][0]->GetYespecie(i, 1); // CO2
						FFraccionMasicaEspecieSalida[j][i][2] = FCanal[j][0]->GetYespecie(i, 2); // H2O
						FFraccionMasicaEspecieSalida[j][i][3] = FCanal[j][0]->GetYespecie(i, 3);    // HC
						FFraccionMasicaEspecieSalida[j][i][4] = FCanal[j][0]->GetYespecie(i, 4); // Soot
						FFraccionMasicaEspecieSalida[j][i][5] = FCanal[j][0]->GetYespecie(i, 5);   // NOx
						FFraccionMasicaEspecieSalida[j][i][6] = FCanal[j][0]->GetYespecie(i, 6); // CO
						if(FNumeroEspecies == 9) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7);   // N2
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // EGR
						} else if(FNumeroEspecies == 10) {
							FFraccionMasicaEspecieSalida[j][i][7] = FCanal[j][0]->GetYespecie(i, 7); // Combustible
							FFraccionMasicaEspecieSalida[j][i][8] = FCanal[j][0]->GetYespecie(i, 8);   // N2
							FFraccionMasicaEspecieSalida[j][i][9] = FCanal[j][0]->GetYespecie(i, 9);   // EGR
						}
					}
				}
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::CalculoSubmodelos en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TDPF::ComunicacionTubos(TCondicionContorno **CC, int numCC) {
	try {

		for(int i = 0; i < numCC; i++) {
			if(CC[i]->getTipoCC() == nmPipeToPlenumConnection) {
				if(FCanal[0][0]->GetNumeroDeposito() == dynamic_cast<TCCDeposito*>(CC[i])->getNumeroDeposito()
				   && !dynamic_cast<TCCDeposito*>(CC[i])->getUnionDPF()) {
					FTuboEntradaDPF = CC[i]->GetTuboExtremo(0).Pipe;
					if(CC[i]->GetTuboExtremo(0).TipoExtremo == nmLeft) {
						FNodoTuboEntrada = 0;
					} else {
						FNodoTuboEntrada = CC[i]->GetTuboExtremo(0).Pipe->getNin() - 1;
					}

				} else if(FCanal[0][1]->GetNumeroDeposito() == dynamic_cast<TCCDeposito*>(CC[i])->getNumeroDeposito()
						  && !dynamic_cast<TCCDeposito*>(CC[i])->getUnionDPF()) {
					FTuboSalidaDPF = CC[i]->GetTuboExtremo(0).Pipe;
					if(CC[i]->GetTuboExtremo(0).TipoExtremo == nmLeft) {
						FNodoTuboSalida = 0;
					} else {
						FNodoTuboSalida = CC[i]->GetTuboExtremo(0).Pipe->getNin() - 1;
					}
				}
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TDPF::ComunicacionTubos en la DPF " << FNumeroDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#pragma package(smart_init)

