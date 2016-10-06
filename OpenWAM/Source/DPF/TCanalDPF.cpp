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

#include "TCanalDPF.h"
#include "TCondicionContorno.h"
#include "TCCDeposito.h"
#include "TBloqueMotor.h"
#include "TDPF.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TCanalDPF::TCanalDPF(int NumeroEspecies, int j, nmTipoCalculoEspecies CalculoEspecies, nmCalculoGamma CalculoGamma,
					 bool HayEGR, TDPF *DPF, int numerocanal, int NumeroDPF) {

	FDPF = DPF;
	FNumDPF = NumeroDPF;
	FNumeroHaz = j + 1;
	FNumeroCanal = numerocanal;
	if(FNumeroCanal == 0) {
		FTipoCanal = nmCanalEntrada;
	} else if(FNumeroCanal == 1) {
		FTipoCanal = nmCanalSalida;
	}
	FF = 28.454;

	FNodoInicialFuenteCE = 0;

	FHayEGR = HayEGR;
	if(FHayEGR)
		FIntEGR = 0;
	else
		FIntEGR = 1;

	FCicloTubo = 0;

//FTranstermico=false;

	FNumeroEspecies = NumeroEspecies;
	FCalculoEspecies = CalculoEspecies;
	FCalculoGamma = CalculoGamma;
	/*switch(CalculoEspecies){
	 case 0: FCalculoEspecies=nmCalculoSimple; break;
	 case 1: FCalculoEspecies=nmCalculoCompleto; break;
	 }

	 switch(CalculoGamma){
	 case 0: FCalculoGamma=nmGammaConstante; break;
	 case 1: FCalculoGamma=nmComposicion; break;
	 case 2: FCalculoGamma=nmComposicionTemperatura; break;
	 }*/
	FNumEcuaciones = 3 + (FNumeroEspecies - 1 - FIntEGR);

	FComposicionInicial = NULL;
	FFraccionMasicaEspecie = NULL;
	FFraccionMasicaEspecie1 = NULL;
	FTasaFraccionMasicaEspecie = NULL;
	FFraccionMasicaSalida = NULL;
	FFraccionMasicaCC = NULL;
	FVelocidadCC = NULL;
	FDensidadCC = NULL;
	FAreaCC = NULL;
	FGamma = NULL;
	FRMezcla = NULL;
	FGamma1 = NULL;
	FGamma3 = NULL;
	FGamma4 = NULL;
	FGamma5 = NULL;
	FGamma6 = NULL;

	FDExtTramo = NULL;
	FLadoCanalTramo = NULL;
	FLadoMenorCanalTramo = NULL;
	FLadoMayorCanalTramo = NULL;
	FLadoCanal = NULL;
	FLTramo = NULL;
	FDiametroTubo = NULL;
	FDiametroD12 = NULL;
	FDiametroS12 = NULL;
	FLadoCanalD12 = NULL;
	FLadoCanalS12 = NULL;
	FDerLin = NULL;
	FDerLin12 = NULL;
	FDerLinArea = NULL;
	FDerLinArea12 = NULL;
	FArea = NULL;
	FArea12 = NULL;
	FPresion0 = NULL;
	FAsonido0 = NULL;
	FVelocidad0 = NULL;
	FPresion1 = NULL;
	FAsonido1 = NULL;
	FVelocidad1 = NULL;
	FUt = NULL;
	FU0 = NULL;
	FU1 = NULL;
	FU12 = NULL;
	FW = NULL;
	FV1 = NULL;
	FV2 = NULL;
	FUfct0 = NULL;
	FUfct1 = NULL;
	FUfctd = NULL;
	FUfctad = NULL;
	Ffl = NULL;
	FdU = NULL;
	FDeltaFCTd = NULL;
	FflU = NULL;
	FaU = NULL;
	FCoefTurbulencia = NULL;
	ResultadosMedios = NULL;
	ResultInstantaneos = NULL;
	FNumResMedios = 0.;
	FNumResInstant = 0.;
	FNumDistSensores = 0.;
	FIntercooler = false;
	FMod.Modelo = nmLaxWendroff;
	FMod.SubModelo = nmNinguno;
	FMod.OpcionSubModelo = nmNinguna;
	FMod.FormulacionLeyes = nmSinArea;
	/*FTVD=NULL;
	 FTVDdU=NULL;
	 FTVDpp=NULL;
	 FTVDpn=NULL;
	 FTVDrp=NULL;
	 FTVDrn=NULL;
	 FTVDphp=NULL;
	 FTVDphn=NULL;
	 FTVDGp=NULL;
	 FTVDGn=NULL; */
	FCourantLocal = NULL;
	Fhi = NULL;
	Fhe = NULL;
	Frho = NULL;
	FRe = NULL;
	FTVD.Bmas = NULL;
	FTVD.Bvector = NULL;
	FTVD.Bmen = NULL;
	FTVD.Qmatrix = NULL;
	FTVD.Pmatrix = NULL;
	FTVD.gflux = NULL;
	FTVD.Alpha = NULL;
	FTVD.Beta = NULL;
	FTVD.DeltaU = NULL;
	FTVD.DeltaB = NULL;
	FTVD.DeltaW = NULL;
	FTVD.hLandaD = NULL;
	FTVD.LandaD = NULL;
	FTVD.Phi = NULL;
	FTVD.R = NULL;
	FTVD.W = NULL;

	FViscosidadDinamica = NULL;
	Fqreg = NULL;
	Fq_reac1 = NULL;
	Fq_reac2 = NULL;
	FH0Pared = NULL;
	FEspesorSoot = NULL;
	FEficiencia = NULL;
	FRreg1 = NULL;
	FRreg2 = NULL;
	FSupEspecifica = NULL;
	FLongitudVC = NULL;
	FTPared = NULL;
	FContador = 0;
	FVelPro = NULL;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TCanalDPF::~TCanalDPF() {
	if(FComposicionInicial)
		delete[] FComposicionInicial;
	if(FFraccionMasicaEspecie != NULL) {
		for(int i = 0; i < FNin; i++)
			delete[] FFraccionMasicaEspecie[i];
		delete[] FFraccionMasicaEspecie;
	}
	if(FFraccionMasicaEspecie1 != NULL) {
		for(int i = 0; i < FNin; i++)
			delete[] FFraccionMasicaEspecie1[i];
		delete[] FFraccionMasicaEspecie1;
	}
	if(FTasaFraccionMasicaEspecie != NULL) {
		for(int i = 0; i < FNin; i++)
			delete[] FTasaFraccionMasicaEspecie[i];
		delete[] FTasaFraccionMasicaEspecie;
	}
	if(FFraccionMasicaSalida != NULL) {
		for(int i = 0; i < FNin; i++)
			delete[] FFraccionMasicaSalida[i];
		delete[] FFraccionMasicaSalida;
	}
	if(FFraccionMasicaCC != NULL) {
		for(int i = 0; i < 1; i++)
			delete[] FFraccionMasicaCC[i];
		delete[] FFraccionMasicaCC;
	}

	if(FVelocidadCC != NULL)
		delete[] FVelocidadCC;
	if(FDensidadCC != NULL)
		delete[] FDensidadCC;
	if(FAreaCC != NULL)
		delete[] FAreaCC;
	if(FGamma != NULL)
		delete[] FGamma;
	if(FRMezcla != NULL)
		delete[] FRMezcla;
	if(FGamma1 != NULL)
		delete[] FGamma1;
	if(FGamma3 != NULL)
		delete[] FGamma3;
	if(FGamma4 != NULL)
		delete[] FGamma4;
	if(FGamma5 != NULL)
		delete[] FGamma5;
	if(FGamma6 != NULL)
		delete[] FGamma6;

	if(FDExtTramo != NULL)
		delete[] FDExtTramo;
	if(FLadoCanalTramo != NULL)
		delete[] FLadoCanalTramo;
	if(FLadoMenorCanalTramo != NULL)
		delete[] FLadoMenorCanalTramo;
	if(FLadoMayorCanalTramo != NULL)
		delete[] FLadoMayorCanalTramo;
	if(FLadoCanal != NULL)
		delete[] FLadoCanal;
	if(FLTramo != NULL)
		delete[] FLTramo;
	if(FDiametroTubo != NULL)
		delete[] FDiametroTubo;
	if(FDiametroD12 != NULL)
		delete[] FDiametroD12;
	if(FDiametroS12 != NULL)
		delete[] FDiametroS12;
	if(FLadoCanalD12 != NULL)
		delete[] FLadoCanalD12;
	if(FLadoCanalS12 != NULL)
		delete[] FLadoCanalS12;
	if(FDerLin != NULL)
		delete[] FDerLin;
	if(FDerLin12 != NULL)
		delete[] FDerLin12;
	if(FDerLinArea != NULL)
		delete[] FDerLinArea;
	if(FDerLinArea12 != NULL)
		delete[] FDerLinArea12;
	if(FArea != NULL)
		delete[] FArea;
	if(FArea12 != NULL)
		delete[] FArea12;
	if(FPresion0 != NULL)
		delete[] FPresion0;
	if(FAsonido0 != NULL)
		delete[] FAsonido0;
	if(FVelocidad0 != NULL)
		delete[] FVelocidad0;
	if(FPresion1 != NULL)
		delete[] FPresion1;
	if(FAsonido1 != NULL)
		delete[] FAsonido1;
	if(FVelocidad1 != NULL)
		delete[] FVelocidad1;
	if(Fhi != NULL)
		delete[] Fhi;
	if(Fhe != NULL)
		delete[] Fhe;
	if(Frho != NULL)
		delete[] Frho;
	if(FRe != NULL)
		delete[] FRe;
	if(FTPared != NULL)
		delete[] FTPared;
	if(FVelPro != NULL)
		delete[] FVelPro;

	if(FUt != NULL) {
		for(int i = 0; i < 3; i++)
			delete[] FUt[i];
		delete[] FUt;
	}
	if(FU0 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FU0[i];
		delete[] FU0;
	}
	if(FU1 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FU1[i];
		delete[] FU1;
	}
	if(FU12 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FU12[i];
		delete[] FU12;
	}
	if(FW != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FW[i];
		delete[] FW;
	}
	if(FV1 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FV1[i];
		delete[] FV1;
	}
	if(FV2 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FV2[i];
		delete[] FV2;
	}
	if(FUfct0 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FUfct0[i];
		delete[] FUfct0;
	}
	if(FUfct1 != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FUfct1[i];
		delete[] FUfct1;
	}
	if(FUfctd != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FUfctd[i];
		delete[] FUfctd;
	}
	if(FUfctad != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FUfctad[i];
		delete[] FUfctad;
	}
	if(Ffl != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] Ffl[i];
		delete[] Ffl;
	}
	if(FdU != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FdU[i];
		delete[] FdU;
	}
	if(FDeltaFCTd != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FDeltaFCTd[i];
		delete[] FDeltaFCTd;
	}
	if(FflU != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FflU[i];
		delete[] FflU;
	}
	if(FaU != NULL) {
		for(int i = 0; i < FNumEcuaciones; i++)
			delete[] FaU[i];
		delete[] FaU;
	}
	if(FCoefTurbulencia != NULL)
		delete[] FCoefTurbulencia;

	for(int i = 0; i < FNumResMedios; i++) {
		if(ResultadosMedios[i].FraccionSUM != NULL)
			delete[] ResultadosMedios[i].FraccionSUM;
		if(ResultadosMedios[i].FraccionMED != NULL)
			delete[] ResultadosMedios[i].FraccionMED;
	}
	if(ResultadosMedios != NULL)
		delete[] ResultadosMedios;
	for(int i = 0; i < FNumResInstant; i++) {
		if(ResultInstantaneos[i].FraccionINS != NULL)
			delete[] ResultInstantaneos[i].FraccionINS;
	}
	if(ResultInstantaneos != NULL)
		delete[] ResultInstantaneos;

	/*if(FTVD!=NULL){
	 for (int i=0;i<3;i++)delete[] FTVD[i];
	 delete[] FTVD;
	 }
	 if(FTVDdU!=NULL){
	 for (int i=0;i<3;i++)delete[] FTVDdU[i];
	 delete[] FTVDdU;
	 }
	 if(FTVDpp!=NULL)delete[] FTVDpp;
	 if(FTVDpn!=NULL)delete[] FTVDpn;
	 if(FTVDphp!=NULL)delete[] FTVDphp;
	 if(FTVDphn!=NULL)delete[] FTVDphn;
	 if(FTVDrp!=NULL)delete[] FTVDrp;
	 if(FTVDrn!=NULL)delete[] FTVDrn;
	 if(FTVDGp!=NULL)delete[] FTVDGp;
	 if(FTVDGn!=NULL)delete[] FTVDGn; */
	if(FCourantLocal != NULL)
		delete[] FCourantLocal;

	for(int i = 0; i < FNumEcuaciones; i++) {
		for(int k = 0; k < FNumEcuaciones; k++) {
			delete[] FTVD.Pmatrix[i][k];
			delete[] FTVD.Qmatrix[i][k];
		}
		delete[] FTVD.Pmatrix[i];
		delete[] FTVD.Qmatrix[i];
	}
	for(int i = 0; i < FNumEcuaciones; i++) {
		delete[] FTVD.Bmas[i];
		delete[] FTVD.Bmen[i];
		delete[] FTVD.Bvector[i];
		delete[] FTVD.gflux[i];
		delete[] FTVD.Alpha[i];
		delete[] FTVD.Beta[i];
		delete[] FTVD.DeltaU[i];
		delete[] FTVD.DeltaB[i];
		delete[] FTVD.DeltaW[i];
		delete[] FTVD.hLandaD[i];
		delete[] FTVD.LandaD[i];
		delete[] FTVD.Phi[i];
		delete[] FTVD.W[i];
		delete[] FTVD.R[i];
	}
	delete[] FTVD.Pmatrix;
	delete[] FTVD.Qmatrix;
	delete[] FTVD.Bmas;
	delete[] FTVD.Bmen;
	delete[] FTVD.Bvector;
	delete[] FTVD.gflux;
	delete[] FTVD.Alpha;
	delete[] FTVD.Beta;
	delete[] FTVD.DeltaU;
	delete[] FTVD.DeltaB;
	delete[] FTVD.DeltaW;
	delete[] FTVD.hLandaD;
	delete[] FTVD.LandaD;
	delete[] FTVD.Phi;
	delete[] FTVD.W;
	delete[] FTVD.R;

	if(FViscosidadDinamica != NULL)
		delete[] FViscosidadDinamica;
	if(Fqreg != NULL)
		delete[] Fqreg;
	if(Fq_reac1 != NULL)
		delete[] Fq_reac1;
	if(Fq_reac2 != NULL)
		delete[] Fq_reac2;
	if(FH0Pared != NULL)
		delete[] FH0Pared;
	if(FEspesorSoot != NULL)
		delete[] FEspesorSoot;
	if(FEficiencia != NULL)
		delete[] FEficiencia;
	if(FRreg1 != NULL)
		delete[] FRreg1;
	if(FRreg2 != NULL)
		delete[] FRreg2;
	if(FSupEspecifica != NULL)
		delete[] FSupEspecifica;
	if(FLongitudVC != NULL)
		delete[] FLongitudVC;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::LeeDatosGeneralesCanal(const char *FileWAM, fpos_t &filepos, int NodoIzquierdo, int NodoDerecho) {
	try {
		int TipTC = 0, Nodo = 0;
		int metodo[4];
		double c_cese = 0., fracciontotal = 0.;

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d %d %d ", &Nodo, &FNTramos, &FNumeroParesdeCanales);

		if(Nodo == NodoIzquierdo) {
			FNodoIzq = NodoIzquierdo;
			FNodoDer = 0;
		} else if(Nodo == NodoDerecho) {
			FNodoDer = NodoDerecho;
			FNodoIzq = 0;
		}
		fscanf(fich, "%lf %lf %lf ", &FTini, &FPini, &FVelMedia);
		fscanf(fich, " %lf ", &FCoefAjusFric);

		if(FDPF->getTctpt() == 0) {
			fscanf(fich, " %lf %lf %lf ", &Ftmax, &Ftmin, &FlowTcoef);
		}

		TipTC = 2;
		switch(TipTC) {
		case 1:
			FTipoTransCal = nmTuboAdmision;
			break;
		case 2:
			FTipoTransCal = nmTuboEscape;
			break;
		case 3:
			FTipoTransCal = nmPipaEscape;
			break;
		case 4:
			FTipoTransCal = nmPipaAdmision;
			break;
		}

		FComposicionInicial = new double[FNumeroEspecies - FIntEGR];
		for(int i = 0; i < FNumeroEspecies - 1; i++) {
			fscanf(fich, "%lf ", &FComposicionInicial[i]);
			fracciontotal += FComposicionInicial[i];

		}
		if(FHayEGR) {
			if(FCalculoEspecies == nmCalculoCompleto) {
				if(FComposicionInicial[0] > 0.2)
					FComposicionInicial[FNumeroEspecies - 1] = 0.;
				else
					FComposicionInicial[FNumeroEspecies - 1] = 1.;
			} else {
				if(FComposicionInicial[0] > 0.5)
					FComposicionInicial[FNumeroEspecies - 1] = 1.;
				else
					FComposicionInicial[FNumeroEspecies - 1] = 0.;
			}
		}

		if(fracciontotal != 1.) {
			std::cout << "ERROR: La fracci�n m�sica total no puede ser distinta de 1. Repasa la lectura en el haz  " <<
					  FNumeroHaz << std::endl;
			throw Exception(" ");
		}

		fscanf(fich, "%lf ", &FMallado);

		fscanf(fich, "%d ", &metodo[0]);

		if(metodo[0] == 0) {
			FMod.Modelo = nmLaxWendroff;
			fscanf(fich, "%d ", &metodo[1]);
			if(metodo[1] == 0) {
				// Lax&Wendroff
				metodo[3] = 1;
				switch(metodo[3]) {
				case 0:
					FMod.FormulacionLeyes = nmSinArea;
					break;
				case 1:
					FMod.FormulacionLeyes = nmConArea;
					break;
				}
			}

			if(metodo[1] == 1) {
				// Lax&Wendroff + FCT
				FMod.FormulacionLeyes = nmConArea;
				FMod.SubModelo = nmFCT;
				fscanf(fich, "%d ", &metodo[2]);
				switch(metodo[2]) {
				case 0:
					FMod.OpcionSubModelo = nmDDNAD;
					FMod.Difusion = nmDamping;
					FMod.Antidifusion = nmNaive;
					break;
				case 1:
					FMod.OpcionSubModelo = nmDDPAD;
					FMod.Difusion = nmDamping;
					FMod.Antidifusion = nmPhoenical;
					break;
				case 2:
					FMod.OpcionSubModelo = nmDDEAD;
					FMod.Difusion = nmDamping;
					FMod.Antidifusion = nmExplicit;
					break;
				case 3:
					FMod.OpcionSubModelo = nmDSNAD;
					FMod.Difusion = nmSmoothing;
					FMod.Antidifusion = nmNaive;
					break;
				case 4:
					FMod.OpcionSubModelo = nmDSPAD;
					FMod.Difusion = nmSmoothing;
					FMod.Antidifusion = nmPhoenical;
					break;
				case 5:
					FMod.OpcionSubModelo = nmDSEAD;
					FMod.Difusion = nmSmoothing;
					FMod.Antidifusion = nmPhoenical;
					break;
				}
			}
		} else if(metodo[0] == 2) {
			FMod.Modelo = nmTVD;
			FMod.FormulacionLeyes = nmConArea;
		}

		/*if(metodo[0] == 1){
		 fscanf(fich,"%lf ",c_cese);
		 printf(" CE-SE con c = %lf\n",c_cese);
		 printf("WARNING: No calcular conductos c�nicos con el CE-SE");

		 FMod.Modelo=nmCESE;
		 FCcese=c_cese;
		 } */

		fscanf(fich, "%lf ", &FCourant);
		printf(" con N. Courant %lf\n\n", FCourant);

		fgetpos(fich, &filepos);
		fclose(fich);

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::LeeDatosGeneralesCanal en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::LeeDatosGeometricosCanal(const char *FileWAM, fpos_t &filepos) {

	FDExtTramo = new double[FNTramos + 1];
	FLTramo = new double[FNTramos + 1];
	FLadoCanalTramo = new double[FNTramos + 1];
	FLadoMenorCanalTramo = new double[FNTramos + 1];
	FLadoMayorCanalTramo = new double[FNTramos + 1];

	try {

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		FTipoMallado = nmDistancia;

		fscanf(fich, "%d ", &FTipoSeccionCanal);

		switch(FTipoSeccionCanal) {
		case 0:
			FSeccionCanal = nmCuadrada;
			break;
		case 1:
			FSeccionCanal = nmCircular;
			break;
		case 2:
			FSeccionCanal = nmRectangular;
			break;
		case 3:
			FSeccionCanal = nmTriangular;
			break;
		}

		if(FSeccionCanal == nmCuadrada) {
			fscanf(fich, "%lf ", &FLadoCanalTramo[0]);
			FDExtTramo[0] = __cons::SQR_4_Pi * FLadoCanalTramo[0];
			for(int i = 1; i <= FNTramos; i++) {
				fscanf(fich, "%lf %lf ", &FLTramo[i], &FLadoCanalTramo[i]);
				FDExtTramo[i] = __cons::SQR_4_Pi * FLadoCanalTramo[i];
			}
		} else if(FSeccionCanal == nmCircular) {
			fscanf(fich, "%lf ", &FDExtTramo[0]);
			for(int i = 1; i <= FNTramos; i++) {
				fscanf(fich, "%lf %lf ", &FLTramo[i], &FDExtTramo[i]);
			}
		} else if(FSeccionCanal == nmRectangular) {
			fscanf(fich, "%lf %lf ", &FLadoMenorCanalTramo[0], &FLadoMayorCanalTramo[0]);
			FDExtTramo[0] = sqrt(FLadoMenorCanalTramo[0] * FLadoMayorCanalTramo[0] * __cons::_4_Pi);
			for(int i = 1; i <= FNTramos; i++) {
				fscanf(fich, "%lf %lf ", &FLTramo[i], &FLadoMenorCanalTramo[i], &FLadoMayorCanalTramo[i]);
				FDExtTramo[i] = pow(__cons::_4_Pi * FLadoMenorCanalTramo[i] * FLadoMayorCanalTramo[i], 0.5);
			}
		} else if(FSeccionCanal == nmTriangular) {
			fscanf(fich, "%lf ", &FLadoCanalTramo[0]);
			FDExtTramo[0] = __cons::SQR_4_Pi * FLadoCanalTramo[0] * 1.316074012952492; // sqrt(tan(60))=1.316074012952492
			for(int i = 1; i <= FNTramos; i++) {
				fscanf(fich, "%lf %lf ", &FLTramo[i], &FLadoCanalTramo[i]);
				FDExtTramo[i] = __cons::SQR_4_Pi * FLadoCanalTramo[i] * 1.316074012952492; // sqrt(tan(60))=1.316074012952492
			}
		}

		CalculoPuntosMalla();

		fgetpos(fich, &filepos);
		fclose(fich);

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::LeeDatosGeometricosCanal en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculoPuntosMalla() {

	try {
		double *FLTotalTramo;
		double xx = 0.;

		FLTotalTramo = new double[FNTramos + 1];

		FLTramo[0] = 0.;
		FLTotalTramo[0] = FLTramo[0];

		for(int i = 1; i <= FNTramos; i++) {
			FLTotalTramo[i] = FLTotalTramo[i - 1] + FLTramo[i];
		}
		FLongitudTotal = FLTotalTramo[FNTramos];

		if(FTipoMallado == nmDistancia)
			xx = FLongitudTotal / FMallado + 0.5;
		if(FTipoMallado != nmDistancia) {
			std::cout << "WARNING: No se ha definido una forma correcta de mallado en el Canal: " << FNumeroHaz << " de la DPF " <<
					  FNumDPF << std::endl;
		}

		FNin = (int) xx + 1;
		if(FNin < 3)
			FNin = 3;

		FXref = FLongitudTotal / (double)(FNin - 1);

		int ii = 0;
		FDiametroTubo = new double[FNin];
		FArea = new double[FNin];
		FLadoCanal = new double[FNin];
		for(int i = 1; i < FNin - 1; i++) {
			do {
				ii = ii + 1;
			} while(FLTotalTramo[ii] < i * FXref);
			ii = ii - 1;
			double r1 = i * FXref - FLTotalTramo[ii];
			FLadoCanal[i] = InterpolaTubo(FLadoCanalTramo[ii], FLadoCanalTramo[ii + 1], FLTramo[ii + 1], r1);
			FDiametroTubo[i] = InterpolaTubo(FDExtTramo[ii], FDExtTramo[ii + 1], FLTramo[ii + 1], r1);
			//FArea[i]=Pi*FDiametroTubo[i]*FDiametroTubo[i]/4.;
			FArea[i] = pow(FLadoCanal[i], 2.);
		}
		FDiametroTubo[0] = FDExtTramo[0];
		FLadoCanal[0] = FLadoCanalTramo[0];
//FArea[0]=Pi*FDiametroTubo[0]*FDiametroTubo[0]/4.;
		FArea[0] = pow(FLadoCanal[0], 2.);

		FDiametroTubo[FNin - 1] = FDExtTramo[FNTramos];
		FLadoCanal[FNin - 1] = FLadoCanalTramo[FNTramos];
		FArea[FNin - 1] = pow(FLadoCanal[FNin - 1], 2.);
//FArea[FNin-1]=Pi*FDiametroTubo[FNin-1]*FDiametroTubo[FNin-1]/4.;

		FVelPro = new double[FNin];
		FDiametroD12 = new double[FNin];
		FDiametroS12 = new double[FNin];
		FArea12 = new double[FNin];
		FLadoCanalD12 = new double[FNin];
		FLadoCanalS12 = new double[FNin];
		for(int i = 0; i < FNin - 1; i++) {
			FDiametroD12[i] = (FDiametroTubo[i + 1] + FDiametroTubo[i]) / 2.;
			FDiametroS12[i] = sqrt((pow(FDiametroTubo[i + 1], 2.) + pow(FDiametroTubo[i], 2.)) / 2.);
			FArea12[i] = (FArea[i + 1] + FArea[i]) / 2.;
			FLadoCanalD12[i] = (FLadoCanal[i + 1] + FLadoCanal[i]) / 2.;
			FLadoCanalS12[i] = sqrt((pow(FLadoCanal[i + 1], 2.) + pow(FLadoCanal[i], 2.)) / 2.);
		}

		delete[] FLTotalTramo;

		NodoTerminoFuente();

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculoPuntosMalla en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::NodoTerminoFuente() {
	try {

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin; i++) {
				if(FDPF->getLongitudTapon() > i * FXref) {
					FNodoInicialFuente = i + 1;
					FNodoFinalFuente = FNin - 1;
				}
			}
			FDistanciaInterpolacion = FNodoInicialFuente * FXref - FDPF->getLongitudTapon();
		} else if(FTipoCanal == nmCanalSalida) {
			int i = 0;
			while((FDPF->getLongitudTapon() + i * FXref) < FDPF->getLongitudEfec() && i < FNin) {
				FNodoInicialFuente = 0;
				FNodoFinalFuente = i;
				i++;
			}
			FDistanciaInterpolacion = FDPF->getLongitudEfec() - FDPF->getLongitudTapon() - FNodoFinalFuente * FXref;
			FNodoInicialFuenteCE = FDPF->GetCanal(FNumeroHaz - 1, 0)->getNodoInicialFuente();
		}
		if(FDistanciaInterpolacion < 1e-5) {
			FDistanciaInterpolacion = 0.;
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::NodoTerminoFuente en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ComunicacionCanal_CC(TCondicionContorno **CC) {
	try {

// Habra que revisar las CC susceptibles de estar en la trampa (Uni�n a dep�sito,UET,
// extremo cerrado).

		if(FNodoDer == 0) {
			for(int i = 0; i < CC[FNodoIzq - 1]->getNumeroTubosCC(); i++) {
				if(FNumeroHaz - 1 == CC[FNodoIzq - 1]->GetTuboExtremo(i).NumeroHaz) {
					FTuboCCNodoIzq = i;
				}
			}
			FNumeroDeposito = dynamic_cast<TCCDeposito*>(CC[FNodoIzq - 1])->getNumeroDeposito();
		}

		if(FNodoIzq == 0) {
			for(int i = 0; i < CC[FNodoDer - 1]->getNumeroTubosCC(); i++) {
				if(FNumeroHaz - 1 == CC[FNodoDer - 1]->GetTuboExtremo(i).NumeroHaz) {
					FTuboCCNodoDer = i;
				}
			}
			FNumeroDeposito = dynamic_cast<TCCDeposito*>(CC[FNodoDer - 1])->getNumeroDeposito();
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ComunicacionCanal_CC en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::InterpolaTubo(double vizq, double vder, double axid, double xif) {
	double xx = 0., yy = 0.;
	double ret_val = 0.;

	try {

		xx = vder - vizq;
		if(axid != 0.) {
			yy = (xx / axid) * xif;
			ret_val = vizq + yy;
		} else {
			std::cout << "ERROR: valores entrada InterpolaTubo en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		}
		return ret_val;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::InterpolaTubo en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::IniciaVariablesFundamentalesCanalDPF() {
	try {
		double RMezclaIni = 0., CpMezclaIni = 0., CvMezclaIni = 0., GammaIni = 0.;

		FViscosidadDinamica = new double[FNin];
		Fqreg = new double[FNin];
		Fq_reac1 = new double[FNin];
		Fq_reac2 = new double[FNin];
		FH0Pared = new double[FNin];
		FEspesorSoot = new double[FNin];
		FEficiencia = new double[FNin];
		FRreg1 = new double[FNin];
		FRreg2 = new double[FNin];
		FSupEspecifica = new double[FNin];
		FLongitudVC = new double[FNin];

		for(int i = 0; i < FNin; i++) {
			FViscosidadDinamica[i] = 1.4615e-6 * pow(__units::degCToK(FTini),
									 1.5) / (__units::degCToK(FTini) + 110.4); //0.0000172*pow((FTini+273.)/293,0.74);   // Konstandopoulos 2003-01-0846
			FVelPro[i] = 5.;
			FEspesorSoot[i] = 0.;
			FEficiencia[i] = 0.;
		}

		FPresion0 = new double[FNin];
		FAsonido0 = new double[FNin];
		FVelocidad0 = new double[FNin];
		FPresion1 = new double[FNin];
		FAsonido1 = new double[FNin];
		FVelocidad1 = new double[FNin];
		FFraccionMasicaEspecie = new double*[FNin];
		FFraccionMasicaEspecie1 = new double*[FNin];
		FTasaFraccionMasicaEspecie = new double*[FNin];
		FFraccionMasicaSalida = new double*[FNin];
		Frho = new double[FNin];
		for(int i = 0; i < FNin; i++)
			FFraccionMasicaEspecie[i] = new double[FNumeroEspecies - FIntEGR];
		for(int i = 0; i < FNin; i++)
			FFraccionMasicaEspecie1[i] = new double[FNumeroEspecies - FIntEGR];
		for(int i = 0; i < FNin; i++)
			FTasaFraccionMasicaEspecie[i] = new double[FNumeroEspecies - FIntEGR];
		for(int i = 0; i < FNin; i++)
			FFraccionMasicaSalida[i] = new double[FNumeroEspecies - FIntEGR];
		FFraccionMasicaCC = new double*[2];
		for(int i = 0; i < 2; i++)
			FFraccionMasicaCC[i] = new double[FNumeroEspecies - FIntEGR];
		FVelocidadCC = new double[2];
		FDensidadCC = new double[2];
		FAreaCC = new double[2];
		FGamma = new double[FNin];
		FGamma1 = new double[FNin];
		FGamma3 = new double[FNin];
		FGamma4 = new double[FNin];
		FGamma5 = new double[FNin];
		FGamma6 = new double[FNin];
		FRMezcla = new double[FNin];
		FTPared = new double[FNin];
		Fhi = new double[FNin];
		Fhe = new double[FNin];
		FRe = new double[FNin];
		FCoefTurbulencia = new double[FNin];

		FUt = new double*[3];
		for(int i = 0; i < 3; i++)
			FUt[i] = new double[FNin];
		FU0 = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FU0[i] = new double[FNin];
		FU1 = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FU1[i] = new double[FNin];
		FU12 = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FU12[i] = new double[FNin];
		FW = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FW[i] = new double[FNin];
		FV1 = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FV1[i] = new double[FNin];
		FV2 = new double*[FNumEcuaciones];
		for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
			FV2[i] = new double[FNin];
		if(FMod.SubModelo == nmFCT) {
			FUfct0 = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FUfct0[i] = new double[FNin + 1];
			FUfct1 = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FUfct1[i] = new double[FNin + 1];
			FUfctd = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FUfctd[i] = new double[FNin + 1];
			FUfctad = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FUfctad[i] = new double[FNin + 1];
			Ffl = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				Ffl[i] = new double[FNin + 1];
			FdU = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FdU[i] = new double[FNin + 1];
			FDeltaFCTd = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FDeltaFCTd[i] = new double[FNin + 1];
			FflU = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FflU[i] = new double[FNin + 1];
			FaU = new double*[FNumEcuaciones];
			for(int i = 0; i < 3 + (FNumeroEspecies - 1); i++)
				FaU[i] = new double[FNin + 1];
		}

		if(FMod.Modelo == nmTVD)
			DimensionaTVD();

		FDerLin = new double[FNin - 1];
		FDerLin12 = new double[FNin - 1];

		FDerLinArea = new double[FNin - 1];
		FDerLinArea12 = new double[FNin - 1];

		FTime0 = 0.;
		FTime1 = 0.;
		FDeltaTime = 0.;

		for(int i = 0; i < FNin - 1; i++) {
			FDerLin[i] = DerLinF(FDiametroTubo[i], FDiametroTubo[i + 1], FXref);
			FDerLinArea[i] = DerLinFArea(FArea[i], FArea[i + 1], FXref);
		}

		for(int i = 1; i < FNin - 1; i++) {
			FDerLin12[i] = DerLinF(FDiametroD12[i], FDiametroD12[i - 1], FXref);
			FDerLinArea12[i] = DerLinFArea(FArea12[i - 1], FArea12[i], FXref);
		}

// Calculo de Gamma y R para la composici�n inicial.
		if(FCalculoEspecies == nmCalculoCompleto) {

			RMezclaIni = CalculoCompletoRMezcla(FComposicionInicial[0], FComposicionInicial[1], FComposicionInicial[2], 0,
												FCalculoGamma, nmMEP);
			CpMezclaIni = CalculoCompletoCpMezcla(FComposicionInicial[0], FComposicionInicial[1], FComposicionInicial[2], 0,
												  __units::degCToK(FTini), FCalculoGamma, nmMEP);
			GammaIni = CalculoCompletoGamma(RMezclaIni, CpMezclaIni, FCalculoGamma);

		} else if(FCalculoEspecies == nmCalculoSimple) {

			RMezclaIni = CalculoSimpleRMezcla(FComposicionInicial[0], 0, FCalculoGamma, nmMEP);
			CvMezclaIni = CalculoSimpleCvMezcla(__units::degCToK(FTini), FComposicionInicial[0], 0, FCalculoGamma, nmMEP);
			GammaIni = CalculoSimpleGamma(RMezclaIni, CvMezclaIni, FCalculoGamma);

		}

		for(int i = 0; i < FNin; i++) {

			for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
				FFraccionMasicaEspecie[i][j] = FComposicionInicial[j];
				FTasaFraccionMasicaEspecie[i][j] = 0.;
				FFraccionMasicaSalida[i][j] = 0.;
			}
			FRMezcla[i] = RMezclaIni;
			FGamma[i] = GammaIni;
			FGamma1[i] = __Gamma::G1(FGamma[i]);
			FGamma3[i] = __Gamma::G3(FGamma[i]);
			FGamma4[i] = __Gamma::G4(FGamma[i]);
			FGamma5[i] = __Gamma::G5(FGamma[i]);
			FGamma6[i] = __Gamma::G6(FGamma[i]);
			FCoefTurbulencia[i] = 1.;
			Fhi[i] = 0.;
			Fhe[i] = 0.;
			FRe[i] = 0.;

			FPresion0[i] = FPini;
			FAsonido0[i] = sqrt(__units::degCToK(FTini) * FGamma[i] * FRMezcla[i]) / __cons::ARef;
			FVelocidad0[i] = FVelMedia / __cons::ARef;
			FPresion1[i] = FPini;
			FAsonido1[i] = sqrt(__units::degCToK(FTini) * FGamma[i] * FRMezcla[i]) / __cons::ARef;
			FVelocidad1[i] = FVelMedia / __cons::ARef;
			Frho[i] = __units::BarToPa(FPresion0[i]) / FRMezcla[i] / __units::degCToK(FTini);

			if(FMod.FormulacionLeyes == nmConArea) {
				Transforma1Area(FVelocidad0[i], FAsonido0[i], FPresion0[i], FU0, FArea[i], FGamma[i], FGamma1[i],
								FFraccionMasicaEspecie[i], i);

				Transforma1Area(FVelocidad1[i], FAsonido1[i], FPresion1[i], FU1, FArea[i], FGamma[i], FGamma1[i],
								FFraccionMasicaEspecie[i], i);
			} else {
				std::cout << "ERROR: El tipo de formulacion de las leyes de conservacion no esta bien definido" << std::endl;
				throw Exception("");
			}
		}

		for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
			for(int i = 0; i < 2; i++) {
				FFraccionMasicaCC[i][j] = FComposicionInicial[j];
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::IniciaVariablesFundamentalesCanalDPF en el haz " << FNumeroHaz << "de la DPF " <<
				  FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ActualizaPropiedadesGas() {
	try {
		double CvMezcla = 0., CpMezcla = 0., temperatura = 0.;

		for(int i = 0; i < FNin; i++) {
			temperatura = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);

			// Calculo de Gamma y R a partir de la composici�n en cada nodo.
			if(FCalculoEspecies == nmCalculoCompleto) {

				FRMezcla[i] = CalculoCompletoRMezcla(FFraccionMasicaEspecie[i][0], FFraccionMasicaEspecie[i][1],
													 FFraccionMasicaEspecie[i][2], 0, FCalculoGamma, nmMEP);
				CpMezcla = CalculoCompletoCpMezcla(FFraccionMasicaEspecie[i][0], FFraccionMasicaEspecie[i][1],
												   FFraccionMasicaEspecie[i][2], 0, temperatura, FCalculoGamma, nmMEP);
				FGamma[i] = CalculoCompletoGamma(FRMezcla[i], CpMezcla, FCalculoGamma);

			} else if(FCalculoEspecies == nmCalculoSimple) {

				FRMezcla[i] = CalculoSimpleRMezcla(FFraccionMasicaEspecie[i][0], 0, FCalculoGamma, nmMEP);
				CvMezcla = CalculoSimpleCvMezcla(temperatura, FFraccionMasicaEspecie[i][0], 0, FCalculoGamma, nmMEP);
				FGamma[i] = CalculoSimpleGamma(FRMezcla[i], CvMezcla, FCalculoGamma);
			}

			FGamma1[i] = __Gamma::G1(FGamma[i]);
			FGamma3[i] = __Gamma::G3(FGamma[i]);
			FGamma4[i] = __Gamma::G4(FGamma[i]);
			FGamma5[i] = __Gamma::G5(FGamma[i]);
			FGamma6[i] = __Gamma::G6(FGamma[i]);

			Frho[i] = __units::BarToPa(FPresion0[i]) / FRMezcla[i] / temperatura;
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ActualizaPropiedadesGas en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Transforma1Area(double v, double a, double p, double **U, double Area, double Gamma, double Gamma1,
								double *Yespecie, int i) {
	try {

		U[0][i] = Gamma * __units::BarToPa(p) / pow(a * __cons::ARef, 2) * Area;
		U[1][i] = U[0][i] * v * __cons::ARef;
		U[2][i] = Area * __units::BarToPa(p) / Gamma1 + U[1][i] * v * __cons::ARef / 2.0;
		for(int j = 3; j < 3 + (FNumeroEspecies - 2); j++) {
			U[j][i] = U[0][i] * Yespecie[j - 3];
		}
		if(FHayEGR)
			U[3 + (FNumeroEspecies - 2)][i] = U[0][i] * Yespecie[FNumeroEspecies - 1];

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Transforma1Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Transforma2Area(double *v, double *a, double *p, double **U, double Area, double Gamma, double Gamma1,
								double *Yespecie, int i) {
	try {
		double fraccionmasicaacum = 0.;

		*v = U[1][i] / U[0][i] / __cons::ARef;
		*p = (U[2][i] - U[1][i] * *v * __cons::ARef / 2.0) * Gamma1 / Area;
		*a = sqrt(Gamma * *p * Area / U[0][i] / __cons::ARef2);
		*p = __units::PaToBar(*p);

// Soluci�n del Transporte de Especies Qu�micas.
		fraccionmasicaacum = 0.;
		for(int j = 0; j < FNumeroEspecies - 2; j++) {
			Yespecie[j] = U[j + 3][i] / U[0][i];
			fraccionmasicaacum += Yespecie[j];
		}
		Yespecie[FNumeroEspecies - 2] = 1. - fraccionmasicaacum;
		if(FHayEGR)
			Yespecie[FNumeroEspecies - 1] = U[FNumeroEspecies - 2 + 3][i] / U[0][i];

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Transforma2Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Transforma3Area(double **Ufct, double **U, double Area, double Gamma, double Gamma1, double Gamma6,
								int i) {

	try {

		Ufct[0][i] = U[1][i];    // Gasto

		Ufct[1][i] = Gamma * U[2][i] / U[0][i] - pow(U[1][i], 2.) * Gamma1 / 2. / U[0][i] / U[0][i];
		Ufct[2][i] = ((U[2][i] - pow(U[1][i], 2.) / U[0][i] / 2.) * Gamma1)
					 * pow((1. + (Gamma1 / 2.) * pow(U[1][i], 2.) / U[0][i] / U[0][i] / (Gamma * ((U[2][i] - pow(U[1][i],
							 2.) / U[0][i] / 2.) * Gamma1) / U[0][i])), Gamma * Gamma6) / Area;

		for(int j = 3; j < FNumEcuaciones; j++) {
			Ufct[j][i] = U[j][i] * U[1][i] / U[0][i];  // Gasto de cada especie.

		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Transforma3Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Transforma4Area(double **U1, double **Ufctd, double Area, double Gamma, double Gamma1, double Gamma3,
								double Gamma4, double Gamma6, int i) {
	double error = 0., fu = 0., dfu = 0., vel = 0., vel1 = 0.;
	double v, a, p, *Y;
	bool peta = false;
	double pruebadefuego1 = 0., pruebadefuego2 = 0.;

	try {

		Y = new double[FNumeroEspecies - 1];
//secante = false;
		vel = Ufctd[0][i] / (U1[0][i]);
		error = 1.;

//Newton Raphson
		while(error > 0.00000001) {

			fu = vel - Ufctd[0][i] / (Area * Ufctd[2][i] * Gamma4) * pow(2. * Ufctd[1][i],
					(Gamma / Gamma1)) * pow((2. * Ufctd[1][i] - pow(vel, 2.)), -Gamma6);

			dfu = 1. - Ufctd[0][i] * vel / (Area * Ufctd[2][i] * Gamma) * pow(2. * Ufctd[1][i],
					(Gamma / Gamma1)) * pow((2. * Ufctd[1][i] - pow(vel, 2.)), -Gamma / Gamma1);

			vel1 = vel - fu / dfu;
			error = fabs(vel1 - vel);
			vel = vel1;
		}

		if(!peta) {
			v = vel / __cons::ARef;
			a = pow(Gamma1 * (Ufctd[1][i] - (vel * vel) / 2.), 0.5) / __cons::ARef;
			p = Ufctd[2][i] / pow((1 + Gamma3 * pow(v / a, 2.)), Gamma / Gamma1) / 1.e5;
			for(int j = 0; j < FNumeroEspecies - 1 - FIntEGR; j++) {
				if(FTipoCanal == nmCanalSalida && i * FXref < 0.015) {
					Y[j] = FFraccionMasicaSalida[i][j];
				} else if(Ufctd[0][i] > 0.00000001) {
					Y[j] = Ufctd[j + 3][i] / Ufctd[0][i]; // Fracci�n m�sica de cada especie.
				} else {
					Y[j] = FFraccionMasicaEspecie[i][j];
				}
			}

			Transforma1Area(v, a, p, U1, Area, Gamma, Gamma1, Y, i);
		}

		delete[] Y;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Transforma4Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::EstabilidadMetodoCalculo() {
	double VTotalMax = 0., VTotalNodo = 0.;

	try {

		VTotalMax = (FAsonido0[0] + fabs(FVelocidad0[0])) * __cons::ARef;
		for(int i = 0; i < FNin; i++) {
			VTotalNodo = (FAsonido0[i] + fabs(FVelocidad0[i])) * __cons::ARef;
			if(VTotalNodo > VTotalMax)
				VTotalMax = VTotalNodo;
		}
		FDeltaTime = FCourant * FXref / VTotalMax;

		if(FMod.Modelo == nmTVD) {
			TVD_Estabilidad();
		}

		FTime0 = FTime1;
		FTime1 = FTime0 + FDeltaTime;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::EstabilidadMetodoCalculo en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaVariablesFundamentales() {
	try {
		if(FMod.Modelo == nmLaxWendroff && FMod.FormulacionLeyes == nmConArea) {
			LaxWendroffArea();
			if(FMod.SubModelo == nmFCT) {
				FluxCorrectedTransport();
			}
		} else if(FMod.Modelo == nmTVD) {
			TVD_Limitador();
		} else {
			std::cout << "ERROR: Metodo de calculo no implementado" << std::endl;
			throw Exception("");
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaVariablesFundamentales en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::FluxCorrectedTransport() {
	double c1 = 0., c2 = 0., c3 = 0., c4 = 0., sign = 0.;

	try {
		if(FNin > 3) {
			for(int k = 0; k < FNumEcuaciones; k++) {
				FU1[k][0] = FU0[k][0];
				FU1[k][FNin - 1] = FU0[k][FNin - 1];
			}
			// Transformaci�n de variables
			for(int i = 0; i < FNin; i++) {
				Transforma3Area(FUfct0, FU0, FArea[i], FGamma[i], FGamma1[i], FGamma6[i], i);

				Transforma3Area(FUfct1, FU1, FArea[i], FGamma[i], FGamma1[i], FGamma6[i], i);
			}

			// Etapa de Difusion
			if(FMod.Difusion == nmDamping) {
				for(int i = 0; i < FNin - 1; i++) {
					for(int k = 0; k < FNumEcuaciones; k++) {
						Ffl[k][i] = (FUfct0[k][i + 1] - FUfct0[k][i]) / 8.;
					}
				}
			} else if(FMod.Difusion == nmSmoothing) {
				for(int i = 0; i < FNin - 1; i++) {
					for(int k = 0; k < FNumEcuaciones; k++) {
						Ffl[k][i] = (FUfct1[k][i + 1] - FUfct1[k][i]) / 8.;
					}
				}
			} else {
				std::cout << "ERROR: Metodo de Difusion mal definido para el FCT en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
						  std::endl;
				throw Exception("");
			}

			for(int i = 1; i < FNin - 1; i++) {
				for(int k = 0; k < FNumEcuaciones; k++) {
					FdU[k][i] = Ffl[k][i] - Ffl[k][i - 1];
					FUfctd[k][i] = FUfct1[k][i] + FdU[k][i];
				}
			}

			for(int k = 0; k < FNumEcuaciones; k++) {
				FUfctd[k][0] = FUfct1[k][0];
				FUfctd[k][FNin - 1] = FUfct1[k][FNin - 1];
			}

			// Etapa de Antidifusion
			if(FMod.Antidifusion == nmExplicit) {
				for(int i = 1; i < FNin; i++) {
					for(int k = 0; k < FNumEcuaciones; k++) {
						FDeltaFCTd[k][i] = FUfctd[k][i] - FUfctd[k][i - 1];
					}
				}
			} else if(FMod.Antidifusion == nmNaive) {
				for(int i = 1; i < FNin; i++) {
					for(int k = 0; k < FNumEcuaciones; k++) {
						FDeltaFCTd[k][i] = FUfct0[k][i] - FUfct0[k][i - 1];
					}
				}
			} else if(FMod.Antidifusion == nmPhoenical) {
				for(int i = 1; i < FNin; i++) {
					for(int k = 0; k < FNumEcuaciones; k++) {
						FDeltaFCTd[k][i] = FUfct1[k][i] - FUfct1[k][i - 1];
					}
				}
			} else {
				std::cout << "ERROR: Metodo de Antidifusion mal definido para el FCT en el haz " << FNumeroHaz << "de la DPF " <<
						  FNumDPF << std::endl;
				throw Exception("");
			}

			for(int k = 0; k < FNumEcuaciones; k++) {
				FDeltaFCTd[k][0] = FDeltaFCTd[k][1];
				FDeltaFCTd[k][FNin] = FDeltaFCTd[k][FNin - 1];
			}

			for(int i = 1; i < FNin - 2; i++) {
				for(int k = 0; k < FNumEcuaciones; k++) {
					if(FDeltaFCTd[k][i + 1] >= 0.) {
						sign = 1.;
					} else {
						sign = -1.;
					}
					c1 = 5. * sign * FDeltaFCTd[k][i] / 8.;
					c2 = fabs(FDeltaFCTd[k][i + 1]) / 8.;
					c3 = 5. * sign * FDeltaFCTd[k][i + 2] / 8.;
					if(c1 < c2)
						c4 = c1;
					else
						c4 = c2;
					if(c3 < c4)
						c4 = c3;
					if(c4 > 0.)
						FflU[k][i] = sign * c4;
					else
						FflU[k][i] = 0.;
				}
			}

			//Extremo Izquierdo
			for(int k = 0; k < FNumEcuaciones; k++) {
				if(FDeltaFCTd[k][1] >= 0.) {
					sign = 1.;
				} else {
					sign = -1.;
				}
				c2 = fabs(FDeltaFCTd[k][1]) / 8.;
				c3 = 5. * sign * FDeltaFCTd[k][2] / 8.;
				if(c2 < c3)
					c4 = c2;
				else
					c4 = c3;
				if(c4 > 0.)
					FflU[k][0] = sign * c4;
				else
					FflU[k][0] = 0.;
			}
			//Extremo Derecho
			for(int k = 0; k < FNumEcuaciones; k++) {
				if(FDeltaFCTd[k][FNin - 2] >= 0.) {
					sign = 1.;
				} else {
					sign = -1.;
				}
				c1 = 5. * sign * FDeltaFCTd[k][FNin - 2] / 8.;
				c2 = fabs(FDeltaFCTd[k][FNin - 1]) / 8.;
				if(c1 < c2)
					c4 = c1;
				else
					c4 = c2;
				if(c4 > 0.)
					FflU[k][FNin - 2] = sign * c4;
				else
					FflU[k][FNin - 2] = 0.;
			}

			for(int i = 1; i < FNin - 1; i++) {
				for(int k = 0; k < FNumEcuaciones; k++) {
					FaU[k][i] = -FflU[k][i] + FflU[k][i - 1];
					FUfctad[k][i] = FUfctd[k][i] + FaU[k][i];
				}
				/*for(int k=3;k<FNumEcuaciones;k++){
				 FUfctad[k][i]=FUfct1[k][i];
				 } */
				Transforma4Area(FU1, FUfctad, FArea[i], FGamma[i], FGamma1[i], FGamma3[i], FGamma4[i], FGamma6[i], i);
			}

		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::FluxCorrectedTransport en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::LaxWendroffArea() {
	try {
		int Nodos = 0;
		double *VelocidadPared;
		double x1, x2, x3, x4, *hi12, *rho12, *Re12, *TPared12, *CoefTurbulencia12, *Gamma12, *Rmezcla12, *Gamma1_12,
			   *VelocidadPared12, *ViscosidadDinamica12, *q_reac1_12, *q_reac2_12, *H0Pared12,/***TasaFraccionMasicaEspecie12*/
			   **FraccionMasicaSalida12, *EspesorSoot12, *Eficiencia12, *Rreg1_12, *Rreg2_12, *SupEspecifica12, *LongitudVC12,
			   **FraccionMasica12;

		hi12 = new double[FNin - 1];
		rho12 = new double[FNin - 1];
		Re12 = new double[FNin - 1];
		TPared12 = new double[FNin - 1];
		CoefTurbulencia12 = new double[FNin - 1];
		Gamma12 = new double[FNin - 1];
		Rmezcla12 = new double[FNin - 1];
		Gamma1_12 = new double[FNin - 1];
		VelocidadPared = new double[FNin];
		VelocidadPared12 = new double[FNin - 1];
		ViscosidadDinamica12 = new double[FNin - 1];
		q_reac1_12 = new double[FNin - 1];
		q_reac2_12 = new double[FNin - 1];
		H0Pared12 = new double[FNin - 1];
		EspesorSoot12 = new double[FNin - 1];
		Eficiencia12 = new double[FNin - 1];
		Rreg1_12 = new double[FNin - 1];
		Rreg2_12 = new double[FNin - 1];
		SupEspecifica12 = new double[FNin - 1];
		LongitudVC12 = new double[FNin - 1];
		FraccionMasica12 = new double*[FNin - 1];
		for(int i = 0; i < FNin - 1; i++)
			FraccionMasica12[i] = new double[FNumeroEspecies - FIntEGR];
		FraccionMasicaSalida12 = new double*[FNin - 1];
		for(int i = 0; i < FNin - 1; i++)
			FraccionMasicaSalida12[i] = new double[FNumeroEspecies - FIntEGR];

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin; i++) {
				VelocidadPared[i] = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i);
				FEspesorSoot[i] = FDPF->GetEspesorSoot(FNumeroHaz - 1, i);
				FArea[i] = pow(FLadoCanal[i] - 2 * FEspesorSoot[i], 2.);
				FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 2);
				double temp = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
				FH0Pared[i] = FRMezcla[i] * FGamma[i] / FGamma1[i] * temp + pow(VelocidadPared[i], 2.) / 2;

			}
		} else {
			for(int i = 0; i < FNin; i++) {
				if(FNodoInicialFuenteCE == 0) {
					VelocidadPared[i] = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i);
					double temp = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
					FH0Pared[i] = FRMezcla[i] * FGamma[i] / (FGamma[i] - 1) * temp + pow(VelocidadPared[i], 2.) / 2;

					FEficiencia[i] = FDPF->GetEficiencia(FNumeroHaz - 1, i);
					FRreg1[i] = FDPF->GetRreg1(FNumeroHaz - 1, i);
					FRreg2[i] = FDPF->GetRreg2(FNumeroHaz - 1, i);
					Fq_reac1[i] = 0.;
					Fq_reac2[i] = 0.;
					FSupEspecifica[i] = FDPF->GetSupEspecifica(FNumeroHaz - 1, i);
					FLongitudVC[i] = FDPF->GetLongitudVC(FNumeroHaz - 1, i);
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 0);
					for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
						FFraccionMasicaSalida[i][j] = FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, i, j);
					}
				} else if(FNodoInicialFuenteCE + i <= FNin - 1) {
					if(i == 0) {
						VelocidadPared[i] = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i);
						FH0Pared[i] = FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i)
									  - (FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i + 1) - FDPF->GetCanal(FNumeroHaz - 1,
											  0)->GetH0Pared(FNodoInicialFuenteCE + i)) / FXref
									  * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FEficiencia[i] = FDPF->GetEficiencia(FNumeroHaz - 1, FNodoInicialFuenteCE + i)
										 - (FDPF->GetEficiencia(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1) - FDPF->GetEficiencia(FNumeroHaz - 1,
												 FNodoInicialFuenteCE + i)) / FXref
										 * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FRreg1[i] = FDPF->GetRreg1(FNumeroHaz - 1, FNodoInicialFuenteCE + i)
									- (FDPF->GetRreg1(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1) - FDPF->GetRreg1(FNumeroHaz - 1,
											FNodoInicialFuenteCE + i)) / FXref
									* FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FRreg2[i] = FDPF->GetRreg2(FNumeroHaz - 1, FNodoInicialFuenteCE + i)
									- (FDPF->GetRreg2(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1) - FDPF->GetRreg2(FNumeroHaz - 1,
											FNodoInicialFuenteCE + i)) / FXref
									* FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						Fq_reac1[i] = 0.;
						Fq_reac2[i] = 0.;
						FSupEspecifica[i] = FDPF->GetSupEspecifica(FNumeroHaz - 1, FNodoInicialFuenteCE + i)
											- (FDPF->GetSupEspecifica(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1) - FDPF->GetSupEspecifica(FNumeroHaz - 1,
													FNodoInicialFuenteCE + i)) / FXref
											* FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FLongitudVC[i] = FDPF->GetLongitudVC(FNumeroHaz - 1, FNodoInicialFuenteCE + i)
										 - (FDPF->GetLongitudVC(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1) - FDPF->GetLongitudVC(FNumeroHaz - 1,
												 FNodoInicialFuenteCE + i)) / FXref
										 * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0)
									 - (FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1, 0) - FDPF->GetTPared(FNumeroHaz - 1,
											 FNodoInicialFuenteCE + i, 0)) / FXref
									 * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();

						for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
							FFraccionMasicaSalida[i][j] = FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i, j)
														  - (FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i,
																  j + 1) - FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i, j)) / FXref
														  * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						}
						/*FH0Pared[i]=FDPF->GetCanal(FNumeroHaz-1,0)->GetH0Pared(FNodoInicialFuenteCE+i);
						 FEficiencia[i]=FDPF->GetEficiencia(FNumeroHaz-1,FNodoInicialFuenteCE+i);
						 FRreg1[i]=FDPF->GetRreg1(FNumeroHaz-1,FNodoInicialFuenteCE+i);
						 FRreg2[i]=FDPF->GetRreg2(FNumeroHaz-1,FNodoInicialFuenteCE+i);
						 Fq_reac1[i]=0.;
						 Fq_reac2[i]=0.;
						 FSupEspecifica[i]=FDPF->GetSupEspecifica(FNumeroHaz-1,FNodoInicialFuenteCE+i);
						 FLongitudVC[i]=FDPF->GetLongitudVC(FNumeroHaz-1,FNodoInicialFuenteCE+i);
						 FTPared[i]=FDPF->GetTPared(FNumeroHaz-1,FNodoInicialFuenteCE+i,0);

						 for(int j=0;j<FNumeroEspecies-FIntEGR;j++){
						 FFraccionMasicaSalida[i][j]=FDPF->GetFraccionMasicaSalida(FNumeroHaz-1,FNodoInicialFuenteCE+i,j);
						 }  */
					} else {
						VelocidadPared[i] = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i);
						FH0Pared[i] = InterpolaTubo(FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE - 1 + i),
													FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i), 1., (FXref - FDPF->GetCanal(FNumeroHaz - 1,
															0)->getDistanciaInterpolacion()) / FXref);
						FEficiencia[i] = InterpolaTubo(FDPF->GetEficiencia(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i),
													   FDPF->GetEficiencia(FNumeroHaz - 1, FNodoInicialFuenteCE + i), 1.,
													   (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						FRreg1[i] = InterpolaTubo(FDPF->GetRreg1(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i), FDPF->GetRreg1(FNumeroHaz - 1,
												  FNodoInicialFuenteCE + i), 1.,
												  (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						FRreg2[i] = InterpolaTubo(FDPF->GetRreg2(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i), FDPF->GetRreg2(FNumeroHaz - 1,
												  FNodoInicialFuenteCE + i), 1.,
												  (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						Fq_reac1[i] = 0.;
						Fq_reac2[i] = 0.;
						FSupEspecifica[i] = InterpolaTubo(FDPF->GetSupEspecifica(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i),
														  FDPF->GetSupEspecifica(FNumeroHaz - 1, FNodoInicialFuenteCE + i), 1.,
														  (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						FLongitudVC[i] = InterpolaTubo(FDPF->GetLongitudVC(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i),
													   FDPF->GetLongitudVC(FNumeroHaz - 1, FNodoInicialFuenteCE + i), 1.,
													   (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						FTPared[i] = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i, 0),
												   FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0), 1.,
												   (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);

						for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
							FFraccionMasicaSalida[i][j] = InterpolaTubo(FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i,
														  j),
														  FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i, j), 1., (FXref - FDPF->GetCanal(FNumeroHaz - 1,
																  0)->getDistanciaInterpolacion()) / FXref);
						}
					}
				} else {
					FH0Pared[i] = 0.;
					VelocidadPared[i] = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i);
					FEficiencia[i] = 0.;
					FRreg1[i] = 0.;
					FRreg2[i] = 0.;
					FSupEspecifica[i] = FDPF->GetSupEspecifica(FNumeroHaz - 1, FNin - 1);
					Fq_reac1[i] = 0.;
					Fq_reac2[i] = 0.;
					if(i == FNin - 1) {
						FLongitudVC[i] = FDPF->GetLongitudVC(FNumeroHaz - 1, FNin - 1);
					} else {
						FLongitudVC[i] = FDPF->GetLongitudVC(FNumeroHaz - 1, FNin - 2);
					}
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNin - 1, 0);
					for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
						FFraccionMasicaSalida[i][j] = 0.;
					}
				}
			}
		}

// Paso Primero
		Nodos = FNin;

		CalculaFlujo(FU0, FW, FGamma, FGamma1, Nodos);

		CalculaFuente1Area(FU0, FV1, FArea, FGamma1,
						   Nodos); // Hay ke modificar el valor del area en funci�n del espesor de Soot

		CalculaFuente2Area(FU0, FV2, FDiametroTubo, FCoefTurbulencia, Fhi, Frho, FTPared, FGamma, FRMezcla, FGamma1, FLadoCanal,
						   VelocidadPared, FViscosidadDinamica, Fq_reac1, Fq_reac2,
						   FFraccionMasicaSalida, FH0Pared, Nodos, FArea, FEspesorSoot, FEficiencia, FRreg1, FRreg2, FSupEspecifica, FLongitudVC,
						   FFraccionMasicaEspecie);

		for(int i = 0; i < FNin - 1; i++) {
			FDerLinArea[i] = DerLinFArea(FArea[i], FArea[i + 1], FXref);
			for(int j = 0; j < FNumEcuaciones; j++) {
				x1 = (FU0[j][i] + FU0[j][i + 1]) / 2.;
				x2 = -FDeltaTime / 2. / FXref * (FW[j][i + 1] - FW[j][i]);
				x3 = -FDeltaTime / 4. * (FV1[j][i + 1] + FV1[j][i]) * FDerLinArea[i];
				x4 = -FDeltaTime / 4. * (FV2[j][i + 1] + FV2[j][i]);

				if(FTipoCanal == nmCanalSalida && i == FNodoFinalFuente) {
					if(j != 1) {
						x4 = -FDeltaTime / 8. * (FV2[j][i]);
					}
				} else if(FTipoCanal == nmCanalEntrada && i == FNodoInicialFuente - 1) {
					if(j != 1) {
						x4 = -FDeltaTime / 8. * (FV2[j][i + 1]);
					}
				}
				FU12[j][i] = x1 + x2 + x3 + x4;
			}
		}

// Paso segundo
		for(int i = 0; i < FNin - 1; i++) {
			hi12[i] = (Fhi[i] + Fhi[i + 1]) / 2.;
			rho12[i] = (Frho[i] + Frho[i + 1]) / 2.;
			Re12[i] = (FRe[i] + FRe[i + 1]) / 2.;
			TPared12[i] = (FTPared[i] + FTPared[i + 1]) / 2.;
			CoefTurbulencia12[i] = (FCoefTurbulencia[i] + FCoefTurbulencia[i + 1]) / 2.;
			Gamma12[i] = (FGamma[i] + FGamma[i + 1]) / 2.;
			Rmezcla12[i] = (FRMezcla[i] + FRMezcla[i + 1]) / 2.;
			Gamma1_12[i] = (FGamma1[i] + FGamma1[i + 1]) / 2.;
			VelocidadPared12[i] = (VelocidadPared[i] + VelocidadPared[i + 1]) / 2.;
			ViscosidadDinamica12[i] = (FViscosidadDinamica[i] + FViscosidadDinamica[i + 1]) / 2.;
			q_reac1_12[i] = (Fq_reac1[i] + Fq_reac1[i + 1]) / 2.;
			q_reac2_12[i] = (Fq_reac2[i] + Fq_reac2[i + 1]) / 2.;
			H0Pared12[i] = (FH0Pared[i] + FH0Pared[i + 1]) / 2.;
			EspesorSoot12[i] = (FEspesorSoot[i] + FEspesorSoot[i + 1]) / 2.;
			Eficiencia12[i] = (FEficiencia[i] + FEficiencia[i + 1]) / 2.;
			Rreg1_12[i] = (FRreg1[i] + FRreg1[i + 1]) / 2.;
			Rreg2_12[i] = (FRreg2[i] + FRreg2[i + 1]) / 2.;
			SupEspecifica12[i] = (FSupEspecifica[i] + FSupEspecifica[i + 1]) / 2.;
			LongitudVC12[i] = (FLongitudVC[i] + FLongitudVC[i + 1]) / 2.;
			FArea12[i] = (FArea[i] + FArea[i + 1]) / 2.;
			for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
				FraccionMasicaSalida12[i][j] = (FFraccionMasicaSalida[i][j] + FFraccionMasicaSalida[i + 1][j]) / 2.;
				FraccionMasica12[i][j] = (FFraccionMasicaEspecie[i][j] + FFraccionMasicaEspecie[i + 1][j]) / 2.;
			}
			if(FTipoCanal == nmCanalSalida && i == FNodoFinalFuente) {
				VelocidadPared12[i] = VelocidadPared[i];
				H0Pared12[i] = FH0Pared[i];
				Eficiencia12[i] = FEficiencia[i];
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					FraccionMasicaSalida12[i][j] = FFraccionMasicaSalida[i][j];
					FraccionMasica12[i][j] = FFraccionMasicaEspecie[i][j];
				}
			} else if(FTipoCanal == nmCanalEntrada && i == FNodoInicialFuente - 1) {
				VelocidadPared12[i] = VelocidadPared[i + 1];
				H0Pared12[i] = FH0Pared[i + 1];
				Eficiencia12[i] = FEficiencia[i + 1];
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					FraccionMasicaSalida12[i][j] = FFraccionMasicaSalida[i + 1][j];
					FraccionMasica12[i][j] = FFraccionMasicaEspecie[i + 1][j];
				}
			}

		}
		Nodos = FNin - 1;
		CalculaFlujo(FU12, FW, Gamma12, Gamma1_12, Nodos);

		CalculaFuente1Area(FU12, FV1, FArea12, Gamma1_12, Nodos);

		CalculaFuente2Area(FU12, FV2, FDiametroD12, CoefTurbulencia12, hi12, rho12, TPared12, Gamma12, Rmezcla12, Gamma1_12,
						   FLadoCanalD12, VelocidadPared12, ViscosidadDinamica12, q_reac1_12,
						   q_reac2_12, FraccionMasicaSalida12, H0Pared12, Nodos, FArea12, EspesorSoot12, Eficiencia12, Rreg1_12, Rreg2_12,
						   SupEspecifica12, LongitudVC12, FraccionMasica12);

		for(int i = 1; i < FNin - 1; i++) {
			FDerLinArea12[i] = DerLinFArea(FArea12[i - 1], FArea12[i], FXref);
			for(int j = 0; j < FNumEcuaciones; j++) {
				x1 = FU0[j][i];
				x2 = -FDeltaTime / FXref * (FW[j][i] - FW[j][i - 1]);
				x3 = -FDeltaTime / 2. * (FV1[j][i] + FV1[j][i - 1]) * FDerLinArea12[i];
				x4 = -FDeltaTime / 2. * (FV2[j][i] + FV2[j][i - 1]);
				FU1[j][i] = x1 + x2 + x3 + x4;
			}

		}

// Liberaci�n de memoria en vectores locales.
		delete[] hi12;
		delete[] rho12;
		delete[] Re12;
		delete[] TPared12;
		delete[] CoefTurbulencia12;
		delete[] Gamma12;
		delete[] Rmezcla12;
		delete[] Gamma1_12;
		delete[] q_reac1_12;
		delete[] q_reac2_12;
		delete[] ViscosidadDinamica12;
		delete[] VelocidadPared;
		delete[] VelocidadPared12;
		delete[] H0Pared12;
		for(int i = 0; i < FNin - 1; i++)
			delete[] FraccionMasicaSalida12[i];
		delete[] FraccionMasicaSalida12;
		delete[] EspesorSoot12;
		delete[] Eficiencia12;
		delete[] Rreg1_12;
		delete[] Rreg2_12;
		delete[] SupEspecifica12;
		delete[] LongitudVC12;
		for(int i = 0; i < FNin - 1; i++) {
			delete[] FraccionMasica12[i];
		}
		delete[] FraccionMasica12;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::LaxWendroffArea en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaFlujo(double **U, double **W, double *Gamma, double *Gamma1, int Nodos) {
	try {
		for(int i = 0; i < Nodos; i++) {
			if(fabs(U[1][i]) > 1e-100) {
				W[0][i] = U[1][i];
				W[1][i] = U[2][i] * Gamma1[i] - (Gamma[i] - 3.0) * pow(U[1][i], 2) / U[0][i] / 2.;
				W[2][i] = Gamma[i] * U[2][i] * U[1][i] / U[0][i] - Gamma1[i] * pow(U[1][i], 3) / (pow(U[0][i], 2) * 2.);
				for(int j = 3; j < FNumEcuaciones; j++) {
					W[j][i] = U[j][i] * U[1][i] / U[0][i];
				}
			} else {
				W[0][i] = 0.;
				W[1][i] = U[2][i] * Gamma1[i];
				W[2][i] = 0.;
				for(int j = 3; j < FNumEcuaciones; j++) {
					W[j][i] = 0.;
				}
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaFlujo en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaFuente1(double **U, double **V1, double *Gamma, double *Gamma1, int Nodos) {
	try {
		for(int i = 0; i < Nodos; i++) {
			V1[0][i] = U[1][i];
			V1[1][i] = pow(U[1][i], 2.) / U[0][i];
			V1[2][i] = Gamma[i] * U[2][i] * U[1][i] / U[0][i] - Gamma1[i] * pow(U[1][i], 3.) / (pow(U[0][i], 2.) * 2.);
			for(int j = 3; j < FNumEcuaciones; j++) {
				V1[j][i] = U[j][i] * U[1][i] / U[0][i];
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaFuente1 en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaFuente1Area(double **U, double **V1, double *Area, double *Gamma1, int Nodos) {
	double p = 0.;
	try {
		for(int i = 0; i < Nodos; i++) {
			p = (U[2][i] - pow(U[1][i], 2.) / U[0][i] / 2.0) * Gamma1[i] / 1.e5 / Area[i];

			V1[0][i] = 0.;
			V1[1][i] = -1.e5 * p;
			V1[2][i] = 0.;
			for(int j = 3; j < FNumEcuaciones; j++) {
				V1[j][i] = 0.;
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaFuente1Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaFuente2(double **U, double **V2, double *diame, double *CoefTurbulencia, double *hi, double *rho,
							   double *TempParedTubo, double *Gamma, double *Rmezcla, double *Gamma1,
							   double *LadoCanal, double *VelocidadPared, double *ViscosidadDinamica, double *q_reac1, double *q_reac2,
							   double **FraccionMasicaSalida, double *H0Pared, int Nodos, double *EspesorSoot,
							   double *Eficiencia, double *Rreg1, double *Rreg2, double *SupEspecifica, double *LongitudVC, double **FraccionMasica) {
	double v = 0., a = 0., p = 0., tgas = 0., q = 0.;

	try {
		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < Nodos; i++) {
				// paso de las variables en funcion de la velocidad,asonido y presion
				v = U[1][i] / U[0][i] / __cons::ARef;
				p = (U[2][i] - U[1][i] * v * __cons::ARef / 2.0) * Gamma1[i];
				a = sqrt(Gamma[i] * p / U[0][i] / __cons::ARef2);
				p = __units::PaToBar(p);
				if(v > 1e200 || v < -1e200) {
					std::cout << "ERROR: Valor de velocidad no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " <<
							  i << std::endl;
					throw Exception("Error Velociad");
				}
				if(p > 1e200 || p < 0) {
					std::cout << "ERROR: Valor de presion no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " << i
							  << std::endl;
					throw Exception("Error presion");
				}
				if(a > 1e200 || a < 0) {
					std::cout << "ERROR: Valor de velocidad del sonido no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
							  " nodo " << i << std::endl;
					throw Exception("Error velocidad del sonido");
				}

				// determinacion factor calor
				if(FDPF->getCoefAjusTC() == 0) {
					q = 0.;
				} else {
					tgas = a * a * __cons::ARef2 / (Gamma[i] * Rmezcla[i]);
					TransmisionCalor(tgas, LadoCanal[i] - 2 * EspesorSoot[i], &q, CoefTurbulencia[i], hi[i], rho[i], TempParedTubo[i]);
					q = q * FDPF->getCoefAjusTC();
				}

				V2[0][i] = 4 / (LadoCanal[i] - 2 * EspesorSoot[i]) * rho[i] * VelocidadPared[i];
				V2[1][i] = FCoefAjusFric * FF * ViscosidadDinamica[i] * v / pow(LadoCanal[i], 2);
				V2[2][i] = -U[0][i] * (q - 4 * H0Pared[i] * VelocidadPared[i] / (LadoCanal[i] - 2 * EspesorSoot[i]));
				for(int j = 3; j < FNumEcuaciones; j++) {
					V2[j][i] = 4 / (LadoCanal[i] - 2 * EspesorSoot[i]) * rho[i] * VelocidadPared[i] * U[j][i] / U[0][i];
				}
			}
		} else if(FTipoCanal == nmCanalSalida) {
			for(int i = 0; i < Nodos; i++) {
				// paso de las variables en funcion de la velocidad,asonido y presion
				v = U[1][i] / U[0][i] / __cons::ARef;
				p = (U[2][i] - U[1][i] * v * __cons::ARef / 2.0) * Gamma1[i];
				a = sqrt(Gamma[i] * p / U[0][i] / __cons::ARef2);
				p = __units::PaToBar(p);
				if(v > 1e200 || v < -1e200) {
					std::cout << "ERROR: Valor de velocidad no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " <<
							  i << std::endl;
					throw Exception("Error Velociad");
				}
				if(p > 1e200 || p < 0) {
					std::cout << "ERROR: Valor de presion no v�lido ven el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " << i
							  << std::endl;
					throw Exception("Error presion");
				}
				if(a > 1e200 || a < 0) {
					std::cout << "ERROR: Valor de velocidad del sonido no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
							  " nodo " << i << std::endl;
					throw Exception("Error velocidad del sonido");
				}

				// determinacion factor calor
				if(FDPF->getCoefAjusTC() == 0) {
					q = 0.;
				} else {
					tgas = a * a * __cons::ARef2 / (Gamma[i] * Rmezcla[i]);
					TransmisionCalor(tgas, LadoCanal[i], &q, CoefTurbulencia[i], hi[i], rho[i], TempParedTubo[i]);
					q = q * FDPF->getCoefAjusTC();
				}

				V2[0][i] = -4 * LadoCanal[i] * rho[i] * VelocidadPared[i];
				V2[1][i] = FCoefAjusFric * FF * ViscosidadDinamica[i] * v / pow(LadoCanal[i], 2);
				V2[2][i] = -U[0][i] * (q + q_reac1[i] + q_reac2[i] + 4 * H0Pared[i] * VelocidadPared[i] / LadoCanal[i]);
				for(int j = 3; j < FNumEcuaciones; j++) {
					//V2[j][i]=-4/LadoCanal[i]*rho[i]*VelocidadPared[i]*U[j][i]/U[0][i]+U[0][i]*TasaFraccionMasicaEspecie[j-3][i];
					V2[j][i] = -4 * LadoCanal[i] * rho[i] * VelocidadPared[i] * FraccionMasicaSalida[j - 3][i];
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaFuente2 en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaFuente2Area(double **U, double **V2, double *diame, double *CoefTurbulencia, double *hi,
								   double *rho, double *TempParedTubo, double *Gamma, double *Rmezcla, double *Gamma1,
								   double *LadoCanal, double *VelocidadPared, double *ViscosidadDinamica, double *q_reac1, double *q_reac2,
								   double **FraccionMasicaSalida, double *H0Pared, int Nodos, double *Area,
								   double *EspesorSoot, double *Eficiencia, double *Rreg1, double *Rreg2, double *SupEspecifica, double *LongitudVC,
								   double **FraccionMasica) {
	double v = 0., a = 0., p = 0., tgas = 0., q = 0.;

	try {
		/*for(int i=0;i<Nodos;i++){
		 // paso de las variables en funcion de la velocidad,asonido y presion
		 v=U[1][i]/U[0][i]/__cons::ARef;
		 p=(U[2][i]-U[1][i]*v*__cons::ARef/2.0)*Gamma1[i]/1e5/Area[i];
		 a=sqrt(Gamma[i]*p*1e5*Area[i]/U[0][i]/__cons::ARef/__cons::ARef);

		 diame=sqrt(Area[i]*4/Pi);
		 // calculo factor friccion
		 if(v==0. || FCoefAjusFric==0){
		 g=0.;
		 }else{
		 Colebrook(FFriccion,diame,&f,Re[i]);
		 g=f*v*v*v*__cons::ARef*__cons::ARef/fabs(v)*2/diame*FCoefAjusFric;
		 }

		 // determinacion factor calor
		 if(FCoefAjusTC==0){   //  q=0
		 q=0.;
		 }
		 else{
		 tgas=a*a*__cons::ARef*__cons::ARef/(Gamma[i]*Rmezcla[i]);
		 TransmisionCalor(tgas,diame,&q,CoefTurbulencia[i],hi[i],rho[i],TempParedTubo[i]);
		 q=q*FCoefAjusTC;
		 }

		 // asignacion a V2
		 V2[0][i]=0.0;
		 V2[1][i]= U[0][i]*g;
		 V2[2][i]= -U[0][i]*q;
		 for(int j=3;j<FNumEcuaciones;j++){
		 V2[j][i]=0.;
		 }
		 } */

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < Nodos; i++) {
				// paso de las variables en funcion de la velocidad,asonido y presion
				v = U[1][i] / U[0][i] / __cons::ARef;
				p = (U[2][i] - U[1][i] * v * __cons::ARef / 2.0) * Gamma1[i] / Area[i];
				a = sqrt(Gamma[i] * p * Area[i] / U[0][i] / __cons::ARef2);
				p = __units::PaToBar(p);
				if(v > 1e200 || v < -1e200) {
					std::cout << "ERROR: Valor de velocidad no v�lido en el haz " << FNumeroHaz << " de la DPF " << FNumDPF << " nodo " <<
							  i << std::endl;
					throw Exception("Error Velociad");
				}
				if(p > 1e200 || p < 0) {
					std::cout << "ERROR: Valor de presion no v�lido en el haz " << FNumeroHaz << " de la DPF " << FNumDPF << " nodo " << i
							  << std::endl;
					throw Exception("Error presion");
				}
				if(a > 1e200 || a < 0) {
					std::cout << "ERROR: Valor de velocidad del sonido no v�lido en el haz " << FNumeroHaz << " de la DPF " << FNumDPF <<
							  " nodo " << i << std::endl;
					throw Exception("Error velocidad del sonido");
				}

				// determinacion factor calor
				if(FDPF->getCoefAjusTC() == 0) {
					q = 0.;
				} else {
					tgas = a * a * __cons::ARef2 / (Gamma[i] * Rmezcla[i]);
					TransmisionCalor(tgas, LadoCanal[i] - 2 * EspesorSoot[i], &q, CoefTurbulencia[i], hi[i], rho[i], TempParedTubo[i]);

					if(FDPF->getTctpt() == 0) {
						if(tgas < __units::degCToK(Ftmax) && tgas > __units::degCToK(Ftmin)) {
							q = q * ((tgas - __units::degCToK(Ftmin)) * (FDPF->getCoefAjusTC() - FlowTcoef) / (Ftmax - Ftmin) + FlowTcoef);
						} else if(tgas > __units::degCToK(Ftmax)) {
							q = q * FDPF->getCoefAjusTC();
						} else {
							q = q * FlowTcoef;
						}
					} else {
						q = q * FDPF->getCoefAjusTC();
						//q = q*0.02; apaño momentaneo hecho con Pedro antes de enterarnos que habia que inicializar EspesorSoot a 0
					}
				}

				V2[0][i] = 4 * (LadoCanal[i] - 2 * EspesorSoot[i]) * rho[i] * VelocidadPared[i];
				V2[1][i] = FCoefAjusFric * FF * ViscosidadDinamica[i] * v * __cons::ARef;
				V2[2][i] = -U[0][i] * (q - 4 * H0Pared[i] * VelocidadPared[i] * (LadoCanal[i] - 2 * EspesorSoot[i]) / Area[i]);
				for(int j = 3; j < FNumEcuaciones; j++) {
					V2[j][i] = 4 * (LadoCanal[i] - 2 * EspesorSoot[i]) * rho[i] * VelocidadPared[i] * U[j][i] / U[0][i];
				}
			}
		} else if(FTipoCanal == nmCanalSalida) {
			for(int i = 0; i < Nodos; i++) {
				// paso de las variables en funcion de la velocidad,asonido y presion
				v = U[1][i] / U[0][i] / __cons::ARef;
				p = (U[2][i] - U[1][i] * v * __cons::ARef / 2.0) * Gamma1[i] / Area[i];
				a = sqrt(Gamma[i] * p * Area[i] / U[0][i] / __cons::ARef2);
				p = __units::PaToBar(p);
				if(v > 1e200 || v < -1e200) {
					std::cout << "ERROR: Valor de velocidad no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " <<
							  i << std::endl;
					throw Exception("Error Velociad");
				}
				if(p > 1e200 || p < 0) {
					std::cout << "ERROR: Valor de presion no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << " nodo " << i
							  << std::endl;
					throw Exception("Error presion");
				}
				if(a > 1e200 || a < 0) {
					std::cout << "ERROR: Valor de velocidad del sonido no v�lido en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
							  " nodo " << i << std::endl;
					throw Exception("Error velocidad del sonido");
				}

				// determinacion factor calor
				if(FDPF->getCoefAjusTC() == 0) {
					q = 0.;
				} else {
					tgas = a * a * __cons::ARef2 / (Gamma[i] * Rmezcla[i]);
					TransmisionCalor(tgas, LadoCanal[i], &q, CoefTurbulencia[i], hi[i], rho[i], TempParedTubo[i]);

					if(FDPF->getTctpt() == 0) {
						if(tgas < __units::degCToK(Ftmax) && tgas > __units::degCToK(Ftmin)) {
							q = q * ((tgas - __units::degCToK(Ftmin)) * (FDPF->getCoefAjusTC() - FlowTcoef) / (Ftmax - Ftmin) + FlowTcoef);
						} else if(tgas > __units::degCToK(Ftmax)) {
							q = q * FDPF->getCoefAjusTC();
						} else {
							q = q * FlowTcoef;
						}
					} else {
						q = q * FDPF->getCoefAjusTC();
					}
				}

				V2[0][i] = -4 * LadoCanal[i] * rho[i] * VelocidadPared[i];
				V2[1][i] = FCoefAjusFric * FF * ViscosidadDinamica[i] * v * __cons::ARef;
				V2[2][i] = -U[0][i] * (q + q_reac1[i] + q_reac2[i] + 4 * H0Pared[i] * VelocidadPared[i] * LadoCanal[i] / Area[i]);
				for(int j = 3; j < FNumEcuaciones; j++) {
					V2[j][i] = -4 * LadoCanal[i] * rho[i] * VelocidadPared[i] * FraccionMasicaSalida[i][j - 3];
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaFuente2Area en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Colebrook(double rug, double dia, double *f, double Re) {
	double u = 0., p = 0., t = 0., rho = 0., mu = 0., re = 0., temp = 0.;

	try {
		if(Re < 2000.0 && Re > 1.0) {
			*f = 32 / Re;
		} else {
			if(Re < 1.0) {
				Re = 1.0;
				*f = 32. / Re;
			} else {
				if(Re < 4000.0) {
					Re = 4000;
				}
				temp = (rug / (3700.0 * dia) + 5.74 / pow(Re, 0.9));
				temp = pow(log10(temp), 2.0);
				*f = 0.0625 / temp;
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Colebrook en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::TransmisionCalor(double tgas, double lado, double *q, double CoefTurbulencia, double hi, double rho,
								 double Tw) {

	try {

		if(hi != 0) {
			switch(FTipoTransCal) {
			case nmPipaEscape:
				*q = 4. * hi * (__units::degCToK(Tw) - tgas) / rho / lado;
				break;

			case nmTuboEscape:
				*q = 4. * hi * (__units::degCToK(Tw) - tgas) / rho / lado * CoefTurbulencia;
				break;

			case nmTuboAdmision:
				*q = 4. * hi * (__units::degCToK(Tw) - tgas) / rho / lado;
				break;

			default:
				std::cout << "WARNING: Tipo de transmision de calor no valida en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
						  std::endl;
				*q = 0.;
			}
		} else
			*q = 0.;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::TransmisionCalor en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::DerLinF(double d1, double d2, double xref) {
	double dm = 0., ret_val = 0.;

	try {

		dm = (d1 + d2) / 2.0;
		ret_val = 2. * (d2 - d1) / dm / xref;
		return ret_val;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::DerLinF en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::DerLinFArea(double area1, double area2, double xref) {
	double ret_val = 0.;

	try {
		ret_val = (area2 - area1) / xref;
		return ret_val;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::DerLinFArea en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ActualizaValoresNuevos(TCondicionContorno **CC) {
	try {

		double a = 0., v = 0., p = 0.;
		double LandaIzq = 0., BetaIzq = 0., EntropiaIzq = 0.;
		double LandaDer = 0., BetaDer = 0., EntropiaDer = 0.;
		double *YIzq, *YDer;

		YIzq = new double[FNumeroEspecies - FIntEGR];
		YDer = new double[FNumeroEspecies - FIntEGR];

		if(FNodoDer == 0) {
			LandaIzq = CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).Landa;
			BetaIzq = CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).Beta;
			EntropiaIzq = CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).Entropia;
			for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
				YIzq[i] = CC[FNodoIzq - 1]->GetFraccionMasicaEspecie(i);
			}
		} else {
			LandaIzq = FLandaExtremoCerrado;
			BetaIzq = FBetaExtremoCerrado;
			EntropiaIzq = FEntropiaExtremoCerrado;
			for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
				YIzq[i] = FFraccionMasicaCC[0][i];
			}
		}

		TransformaContorno(&LandaIzq, &BetaIzq, &EntropiaIzq, &a, &v, &p, 1, FGamma1[0], FGamma3[0], FGamma4[0], FGamma5[0]);

		if(FMod.FormulacionLeyes == nmConArea)
			Transforma1Area(v, a, p, FU1, FArea[0], FGamma[0], FGamma1[0], YIzq, 0);

		if(FNodoIzq == 0) {
			LandaDer = CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).Landa;
			BetaDer = CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).Beta;
			EntropiaDer = CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).Entropia;
			for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
				YDer[i] = CC[FNodoDer - 1]->GetFraccionMasicaEspecie(i);
			}
		} else {
			LandaDer = FLandaExtremoCerrado;
			BetaDer = FBetaExtremoCerrado;
			EntropiaDer = FEntropiaExtremoCerrado;
			for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
				YDer[i] = FFraccionMasicaCC[1][i];
			}
		}

		TransformaContorno(&LandaDer, &BetaDer, &EntropiaDer, &a, &v, &p, 1, FGamma1[FNin - 1], FGamma3[FNin - 1],
						   FGamma4[FNin - 1], FGamma5[FNin - 1]);

		if(FMod.FormulacionLeyes == nmConArea)
			Transforma1Area(v, a, p, FU1, FArea[FNin - 1], FGamma[FNin - 1], FGamma1[FNin - 1], YDer, FNin - 1);

		if(FMod.FormulacionLeyes == nmConArea) {
			for(int i = 0; i < FNin; i++) {
				FVelocidad1[i] = FVelocidad0[i];
				FAsonido1[i] = FAsonido0[i];
				FPresion1[i] = FPresion0[i];
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					FFraccionMasicaEspecie1[i][j] = FFraccionMasicaEspecie[i][j];
				}
				Transforma2Area(&FVelocidad0[i], &FAsonido0[i], &FPresion0[i], FU1, FArea[i], FGamma[i], FGamma1[i],
								FFraccionMasicaEspecie[i], i);

				for(int k = 0; k < FNumEcuaciones; k++) {
					FU0[k][i] = FU1[k][i];
				}
			}
		}

		delete[] YIzq;
		delete[] YDer;

	}

	catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ActualizaValoresNuevos en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::TransformaContorno(double *L, double *B, double *E, double *a, double *v, double *p, int modo,
								   double Gamma1, double Gamma3, double Gamma4, double Gamma5) {
	try {
		if(modo == 0) {
			*L = (*a + Gamma3 * *v);
			*B = (*a - Gamma3 * *v);
			*E = *a / pow(*p, Gamma5);
		} else {
			*a = (*L + *B) / 2.;
			*v = (*L - *B) / Gamma1;
			*p = pow(*a / *E, Gamma4);
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::TransformaContorno en el haz " << FNumeroHaz << "de la DPF " << FNumDPF << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ReduccionFlujoSubsonico() {
	double Machx = 0., Machy = 0., Velocidady = 0., Sonidoy = 0.;

	try {

		for(int i = 0; i < FNin; i++) {
			Machx = FVelocidad0[i] / FAsonido0[i];
			if(-1. >= Machx || Machx > 1.) {
				Machy = Machx / fabs(Machx) * sqrt((pow(Machx, 2) + 2. / FGamma1[i]) / (FGamma4[i] * pow(Machx, 2) - 1.));
				Sonidoy = FAsonido0[i] * sqrt((FGamma1[i] / 2. * pow(Machx, 2) + 1.) / (FGamma1[i] / 2. * pow(Machy, 2) + 1.));

				Velocidady = Sonidoy * Machy;
				FAsonido0[i] = Sonidoy;
				FVelocidad0[i] = Velocidady;
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ReduccionFlujoSubsonico en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ReduccionFlujoSubsonicoFCT() {
	double Machx = 0., Machy = 0., Velocidady = 0., Sonidoy = 0.;
	double velocidad = 0., asonido = 0., presion = 0.;

	try {

		for(int i = 1; i < FNin - 1; i++) {
			velocidad = FU1[1][i] / FU1[0][i] / __cons::ARef;
			presion = (FU1[2][i] - FU1[1][i] * velocidad * __cons::ARef / 2.0) * FGamma1[i] / FArea[i];
			asonido = sqrt(FGamma[i] * presion * FArea[i] / FU1[0][i] / __cons::ARef / __cons::ARef);
			Machx = velocidad / asonido;
			if(-1. >= Machx || Machx > 1.) {
				Machy = Machx / fabs(Machx) * sqrt((pow(Machx, 2) + 2. / FGamma1[i]) / (FGamma4[i] * pow(Machx, 2) - 1.));
				Sonidoy = asonido * sqrt((FGamma1[i] / 2. * pow(Machx, 2) + 1.) / (FGamma1[i] / 2. * pow(Machy, 2) + 1.));

				Velocidady = Sonidoy * Machy;
				asonido = Sonidoy;
				velocidad = Velocidady;
				FU1[0][i] = FGamma[i] * presion / pow(asonido * __cons::ARef, 2) * FArea[i];
				FU1[1][i] = FU1[0][i] * velocidad * __cons::ARef;
				FU1[2][i] = FArea[i] * presion / FGamma1[i] + FU1[1][i] * velocidad * __cons::ARef / 2.0;

			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ReduccionFlujoSubsonicoFCT en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::LeeResultadosMediosCanalDPF(int NumResMedios, const char *FileWAM, fpos_t &filepos) {
	int NumVars = 0, TipoVar = 0;

	try {

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		FNumResMedios = NumResMedios;
		ResultadosMedios = new stResMediosTubo[FNumResMedios];
		FTiempoMedSUM = 0.;
		FControlResMed = 1;

		for(int i = 0; i < FNumResMedios; i++) {
			ResultadosMedios[i].TemperaturaGas = false;
			ResultadosMedios[i].TemperaturaGasSUM = 0.;
			ResultadosMedios[i].TemperaturaGasMED = 0;
			ResultadosMedios[i].Pressure = false;
			ResultadosMedios[i].PresionSUM = 0.;
			ResultadosMedios[i].PresionMED = 0.;
			ResultadosMedios[i].Velocity = false;
			ResultadosMedios[i].VelocidadSUM = 0.;
			ResultadosMedios[i].VelocidadMED = 0.;
			ResultadosMedios[i].Massflow = false;
			ResultadosMedios[i].GastoSUM = 0.;
			ResultadosMedios[i].GastoMED = 0.;
			ResultadosMedios[i].TemperaturaInternaPared = false;
			ResultadosMedios[i].TemperaturaInternaParedSUM = 0.;
			ResultadosMedios[i].TemperaturaInternaParedMED = 0.;
			ResultadosMedios[i].TemperaturaIntermediaPared = false;
			ResultadosMedios[i].TemperaturaIntermediaParedSUM = 0.;
			ResultadosMedios[i].TemperaturaIntermediaParedMED = 0.;
			ResultadosMedios[i].TemperaturaExternaPared = false;
			ResultadosMedios[i].TemperaturaExternaParedSUM = 0.;
			ResultadosMedios[i].TemperaturaExternaParedMED = 0.;
			ResultadosMedios[i].NITmedio = false;
			ResultadosMedios[i].NITmedioSUM = 0;
			ResultadosMedios[i].NITmedioMED = 0;
			ResultadosMedios[i].CoefPelInterior = false;
			ResultadosMedios[i].CoefPelInteriorSUM = 0.;
			ResultadosMedios[i].CoefPelInteriorMED = 0.;
			ResultadosMedios[i].FraccionMasicaEspecies = false;
			ResultadosMedios[i].FraccionSUM = new double[FNumeroEspecies - FIntEGR];
			ResultadosMedios[i].FraccionMED = new double[FNumeroEspecies - FIntEGR];
			for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
				ResultadosMedios[i].FraccionSUM[j] = 0.;
				ResultadosMedios[i].FraccionMED[j] = 0.;
			}

			fscanf(fich, "%lf %d ", &ResultadosMedios[i].Distancia, &NumVars);

			for(int j = 0; j < NumVars; j++) {
				fscanf(fich, "%d ", &TipoVar);
				switch(TipoVar) {
				case 0:
					ResultadosMedios[i].TemperaturaGas = true;
					break;
				case 1:
					ResultadosMedios[i].Pressure = true;
					break;
				case 2:
					ResultadosMedios[i].Velocity = true;
					break;
				case 3:
					ResultadosMedios[i].Massflow = true;
					break;
				case 4:
					ResultadosMedios[i].CoefPelInterior = true;
					break;
				case 5:
					ResultadosMedios[i].FraccionMasicaEspecies = true;
					break;
				default:
					std::cout << "WARNING: El tipo de variable seleccionada para la salida de resultados medios no es v�lida" <<
							  std::endl;
				}

			}

		}

		fgetpos(fich, &filepos);
		fclose(fich);

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::LeeResultadosMedios en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CabeceraResultadosMedios(stringstream& medoutput, stEspecies *DatosEspecies) const {
	try {
		std::string TipoCanal;
		std::string Label;
		std::ostringstream TextDist;
		TextDist.precision(8);

		if(FTipoCanal == nmCanalEntrada) {
			TipoCanal = new char[3];
			TipoCanal = "_Ent";
		} else {
			TipoCanal = new char[3];
			TipoCanal = "_Sal";
		}
		for(int i = 0; i < FNumResMedios; i++) {
			TextDist << ResultadosMedios[i].Distancia;
			if(ResultadosMedios[i].TemperaturaGas) {
				Label = "\t" + PutLabel(823) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(910);
				medoutput << Label.c_str();
			}
			if(ResultadosMedios[i].Pressure) {
				Label = "\t" + PutLabel(824) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(908);
				medoutput << Label.c_str();
			}
			if(ResultadosMedios[i].Velocity) {
				Label = "\t" + PutLabel(825) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(909);
				medoutput << Label.c_str();
			}
			if(ResultadosMedios[i].Massflow) {
				Label = "\t" + PutLabel(826) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(904);
				medoutput << Label.c_str();
			}
			if(ResultadosMedios[i].CoefPelInterior) {
				Label = "\t" + PutLabel(827) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(911);
				medoutput << Label.c_str();
			}
			if(ResultadosMedios[i].FraccionMasicaEspecies) {
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					Label = "\t" + PutLabel(821) + DatosEspecies[j].Nombre + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(
								801) + std::to_string(FNumeroHaz) + TipoCanal + PutLabel(316)
							+ TextDist.str() + PutLabel(317) + PutLabel(901);
					medoutput << Label.c_str();
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CabeceraResultadosMedios en el haz " << FNumeroHaz << "de la DPF " << FNumDPF <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ImprimeResultadosMedios(stringstream& medoutput) const {
	try {

		for(int i = 0; i < FNumResMedios; i++) {
			if(ResultadosMedios[i].TemperaturaGas)
				medoutput << "\t" << ResultadosMedios[i].TemperaturaGasMED;
			if(ResultadosMedios[i].Pressure)
				medoutput << "\t" << ResultadosMedios[i].PresionMED;
			if(ResultadosMedios[i].Velocity)
				medoutput << "\t" << ResultadosMedios[i].VelocidadMED;
			if(ResultadosMedios[i].Massflow)
				medoutput << "\t" << ResultadosMedios[i].GastoMED;
			if(ResultadosMedios[i].CoefPelInterior)
				medoutput << "\t" << ResultadosMedios[i].CoefPelInteriorMED;
			if(ResultadosMedios[i].FraccionMasicaEspecies) {
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					medoutput << "\t" << ResultadosMedios[i].FraccionMED[j];
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ImprimeResultadosMedios en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::LeeResultadosInstantaneosCanalDPF(int NumResInstant, const char *FileWAM, fpos_t &filepos) {
	int NumVars = 0, TipoVar = 0;

	try {

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

//fscanf(fich,"%d ",&FNumResInstant);  // Puntos del tubo en los que se piden resultados instant�neos
		FNumResInstant = NumResInstant;
		ResultInstantaneos = new stResInstantTubo[FNumResInstant];

		for(int i = 0; i < FNumResInstant; i++) {
			ResultInstantaneos[i].Pressure = false;
			ResultInstantaneos[i].PresionINS = 0.;
			ResultInstantaneos[i].Velocity = false;
			ResultInstantaneos[i].VelocidadINS = 0.;
			ResultInstantaneos[i].TemperaturaGas = false;
			ResultInstantaneos[i].TemperaturaGasINS = 0.;
			ResultInstantaneos[i].FlujoMasico = false;
			ResultInstantaneos[i].FlujoMasicoINS = 0.;
			ResultInstantaneos[i].VelocidadDerecha = false;
			ResultInstantaneos[i].VelocidadDerechaINS = 0.;
			ResultInstantaneos[i].VelocidadIzquierda = false;
			ResultInstantaneos[i].VelocidadIzquierdaINS = 0.;
			ResultInstantaneos[i].PresionDerecha = false;
			ResultInstantaneos[i].PresionDerechaINS = 0.;
			ResultInstantaneos[i].PresionIzquierda = false;
			ResultInstantaneos[i].PresionIzquierdaINS = 0.;
			ResultInstantaneos[i].NIT = false;
			ResultInstantaneos[i].NITINS = 0.;
			ResultInstantaneos[i].TemperaturaInternaPared = false;
			ResultInstantaneos[i].TemperaturaInternaParedINS = 0.;
			ResultInstantaneos[i].TemperaturaIntermediaPared = false;
			ResultInstantaneos[i].TemperaturaIntermediaParedINS = 0.;
			ResultInstantaneos[i].TemperaturaExternaPared = false;
			ResultInstantaneos[i].TemperaturaExternaParedINS = 0.;
			ResultInstantaneos[i].CoefPelInterior = false;
			ResultInstantaneos[i].CoefPelInteriorINS = 0.;
			ResultInstantaneos[i].FraccionMasicaEspecies = false;
			ResultInstantaneos[i].FraccionINS = new double[FNumeroEspecies - FIntEGR];
			ResultInstantaneos[i].Gamma = false;
			for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
				ResultInstantaneos[i].FraccionINS[j] = 0.;
			}

			fscanf(fich, "%lf %d ", &ResultInstantaneos[i].Distancia, &NumVars);

			for(int j = 0; j < NumVars; j++) {
				fscanf(fich, "%d ", &TipoVar);
				switch(TipoVar) {
				case 0:
					ResultInstantaneos[i].Pressure = true;
					break;
				case 1:
					ResultInstantaneos[i].Velocity = true;
					break;
				case 2:
					ResultInstantaneos[i].TemperaturaGas = true;
					break;
				case 3:
					ResultInstantaneos[i].FlujoMasico = true;
					break;
				case 4:
					ResultInstantaneos[i].VelocidadDerecha = true;
					break;
				case 5:
					ResultInstantaneos[i].VelocidadIzquierda = true;
					break;
				case 6:
					ResultInstantaneos[i].PresionDerecha = true;
					break;
				case 7:
					ResultInstantaneos[i].PresionIzquierda = true;
					break;
				case 8:
					ResultInstantaneos[i].CoefPelInterior = true;
					break;
				case 9:
					ResultInstantaneos[i].FraccionMasicaEspecies = true;
					break;
				case 10:
					ResultInstantaneos[i].Gamma = true;
					break;
				default:
					std::cout << "WARNING: El tipo de variable seleccionada para la salida de resultados instantaneos no es v�lida" <<
							  std::endl;
				}
			}
		}

		fgetpos(fich, &filepos);
		fclose(fich);

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::LeeResultadosInstantaneosCanalDPF en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CabeceraResultadosInstantaneos(stringstream& insoutput, stEspecies *DatosEspecies) const {
	try {
		std::string TipoCanal;
		std::string Label;
		std::ostringstream TextDist;
		TextDist.precision(8);

		if(FTipoCanal == nmCanalEntrada) {
			TipoCanal = new char[3];
			TipoCanal = "_Ent";
		} else {
			TipoCanal = new char[3];
			TipoCanal = "_Sal";
		}

		for(int i = 0; i < FNumResInstant; i++) {
			TextDist << ResultInstantaneos[i].Distancia;
			if(ResultInstantaneos[i].Pressure) {
				Label = "\t" + PutLabel(824) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(908);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].Velocity) {
				Label = "\t" + PutLabel(825) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(909);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].TemperaturaGas) {
				Label = "\t" + PutLabel(823) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(910);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].FlujoMasico) {
				Label = "\t" + PutLabel(826) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(904);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].VelocidadDerecha) {
				Label = "\t" + PutLabel(829) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(909);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].VelocidadIzquierda) {
				Label = "\t" + PutLabel(830) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(909);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].PresionDerecha) {
				Label = "\t" + PutLabel(831) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(908);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].PresionIzquierda) {
				Label = "\t" + PutLabel(832) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(908);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].CoefPelInterior) {
				Label = "\t" + PutLabel(827) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(911);
				insoutput << Label.c_str();
			}
			if(ResultInstantaneos[i].FraccionMasicaEspecies) {
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					Label = "\t" + PutLabel(821) + DatosEspecies[j].Nombre + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(
								801) + std::to_string(FNumeroHaz) + TipoCanal + PutLabel(316)
							+ TextDist.str() + PutLabel(317) + PutLabel(901);
					insoutput << Label.c_str();
				}
			}
			if(ResultInstantaneos[i].Gamma) {
				Label = "\t" + PutLabel(828) + PutLabel(800) + std::to_string(FNumDPF) + PutLabel(801) + std::to_string(
							FNumeroHaz) + TipoCanal + PutLabel(316) + TextDist.str() + PutLabel(317)
						+ PutLabel(901);
				insoutput << Label.c_str();
			}

		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CabeceraResultadosInstantaneos en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::ImprimeResultadosInstantaneos(stringstream& insoutput) const {
	try {

		for(int i = 0; i < FNumResInstant; i++) {
			if(ResultInstantaneos[i].Pressure)
				insoutput << "\t" << ResultInstantaneos[i].PresionINS;
			if(ResultInstantaneos[i].Velocity)
				insoutput << "\t" << ResultInstantaneos[i].VelocidadINS;
			if(ResultInstantaneos[i].TemperaturaGas)
				insoutput << "\t" << ResultInstantaneos[i].TemperaturaGasINS;
			if(ResultInstantaneos[i].FlujoMasico)
				insoutput << "\t" << ResultInstantaneos[i].FlujoMasicoINS;
			if(ResultInstantaneos[i].VelocidadDerecha)
				insoutput << "\t" << ResultInstantaneos[i].VelocidadDerechaINS;
			if(ResultInstantaneos[i].VelocidadIzquierda)
				insoutput << "\t" << ResultInstantaneos[i].VelocidadIzquierdaINS;
			if(ResultInstantaneos[i].PresionDerecha)
				insoutput << "\t" << ResultInstantaneos[i].PresionDerechaINS;
			if(ResultInstantaneos[i].PresionIzquierda)
				insoutput << "\t" << ResultInstantaneos[i].PresionIzquierdaINS;
			if(ResultInstantaneos[i].CoefPelInterior)
				insoutput << "\t" << ResultInstantaneos[i].CoefPelInteriorINS;
			if(ResultInstantaneos[i].FraccionMasicaEspecies) {
				for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
					insoutput << "\t" << ResultInstantaneos[i].FraccionINS[j];
				}
			}
			if(ResultInstantaneos[i].Gamma)
				insoutput << "\t" << ResultInstantaneos[i].GammaINS;

		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::ResultadosInstantaneos en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaResultadosMedios(double theta) {
	double dist, Temp, Dens, Area, Gto1, Gto2, Vble, d, Rmezcla, Gamma, GastoPonderacion, a;
	int n1 = 0, n2 = 0;

	try {

		FTiempoMedSUM += FDeltaTime;

		for(int i = 0; i < FNumResMedios; i++) {
			dist = ResultadosMedios[i].Distancia / FXref;
			n1 = (int) floor(dist);
			if(n1 >= FNin - 1) {
				Temp = __units::KTodegC(pow(FAsonido0[FNin - 1] * __cons::ARef, 2.) / (FGamma[FNin - 1] * FRMezcla[FNin - 1]));

				if(ResultadosMedios[i].TemperaturaGas || ResultadosMedios[i].FraccionMasicaEspecies) {
					Dens = __units::BarToPa(FPresion0[FNin - 1]) / FRMezcla[FNin - 1] / __units::degCToK(Temp);
					//Area=pow(FDiametroTubo[FNin-1],2.)*Pi/4.;
					GastoPonderacion = FVelocidad0[FNin - 1] * __cons::ARef * FArea[FNin - 1] * Dens;
					ResultadosMedios[i].PonderacionSUM += GastoPonderacion * FDeltaTime;
					ResultadosMedios[i].GastoPonderacionSUM += GastoPonderacion;
				}
				if(ResultadosMedios[i].TemperaturaGas)
					ResultadosMedios[i].TemperaturaGasSUM += Temp * GastoPonderacion;
				if(ResultadosMedios[i].Pressure)
					ResultadosMedios[i].PresionSUM += FPresion0[FNin - 1] * FDeltaTime;
				if(ResultadosMedios[i].Velocity)
					ResultadosMedios[i].VelocidadSUM += FVelocidad0[FNin - 1] * __cons::ARef * FDeltaTime;
				if(ResultadosMedios[i].Massflow) {
					Dens = __units::BarToPa(FPresion0[FNin - 1]) / FRMezcla[FNin - 1] / __units::degCToK(Temp);
					//Area=pow(FDiametroTubo[FNin-1],2.)*Pi/4.;
					ResultadosMedios[i].GastoSUM += FVelocidad0[FNin - 1] * __cons::ARef * FArea[FNin - 1] * Dens * FDeltaTime;
				}
				if(ResultadosMedios[i].TemperaturaInternaPared)
					ResultadosMedios[i].TemperaturaInternaParedSUM += FTPTubo[0][FNin - 1] * FDeltaTime;
				if(ResultadosMedios[i].TemperaturaIntermediaPared)
					ResultadosMedios[i].TemperaturaIntermediaParedSUM += FTPTubo[1][FNin - 1] * FDeltaTime;
				if(ResultadosMedios[i].TemperaturaExternaPared)
					ResultadosMedios[i].TemperaturaExternaParedSUM += FTPTubo[2][FNin - 1] * FDeltaTime;
				if(ResultadosMedios[i].NITmedio) {
					double nit = CalculaNIT(FAsonido0[FNin - 1], FVelocidad0[FNin - 1], FPresion0[FNin - 1], FDiametroTubo[FNin - 1],
											FGamma[FNin - 1], FRMezcla[FNin - 1]);
					ResultadosMedios[i].NITmedioSUM += nit * FDeltaTime;
				}
				if(ResultadosMedios[i].CoefPelInterior) {
					ResultadosMedios[i].CoefPelInteriorSUM += Fhi[FNin - 1] * FDeltaTime;
				}
				if(ResultadosMedios[i].FraccionMasicaEspecies) {
					for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
						ResultadosMedios[i].FraccionSUM[j] += FFraccionMasicaEspecie[FNin - 1][j] * GastoPonderacion * FDeltaTime;
					}
				}
			} else {
				n2 = n1 + 1;
				d = dist - (double) n1;

				a = InterpolaTubo(FAsonido0[n1], FAsonido0[n2], 1., d);
				Gamma = InterpolaTubo(FGamma[n1], FGamma[n2], 1., d);
				Rmezcla = InterpolaTubo(FRMezcla[n1], FRMezcla[n2], 1., d);
				Temp = __units::KTodegC(pow(a * __cons::ARef, 2.) / (Gamma * Rmezcla));

				if(ResultadosMedios[i].TemperaturaGas || ResultadosMedios[i].FraccionMasicaEspecies) {
					Dens = FGamma[n1] * __units::BarToPa(FPresion0[n1]) / pow(FAsonido0[n1] * __cons::ARef, 2.);
					//Area=pow(FDiametroTubo[n1],2.)*Pi/4.;
					Gto1 = FVelocidad0[n1] * __cons::ARef * Dens * FArea[n1];
					Dens = FGamma[n2] * __units::BarToPa(FPresion0[n2]) / pow(FAsonido0[n2] * __cons::ARef, 2.);
					//Area=pow(FDiametroTubo[n2],2.)*Pi/4.;
					Gto2 = FVelocidad0[n2] * __cons::ARef * Dens * FArea[n2];
					GastoPonderacion = InterpolaTubo(Gto1, Gto2, 1., d);
					ResultadosMedios[i].GastoPonderacionSUM += GastoPonderacion;
					ResultadosMedios[i].PonderacionSUM += GastoPonderacion * FDeltaTime;
				}

				if(ResultadosMedios[i].TemperaturaGas)
					ResultadosMedios[i].TemperaturaGasSUM += Temp * GastoPonderacion;
				if(ResultadosMedios[i].Pressure) {
					Vble = InterpolaTubo(FPresion0[n1], FPresion0[n2], 1., d);
					ResultadosMedios[i].PresionSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].Velocity) {
					Vble = InterpolaTubo(FVelocidad0[n1], FVelocidad0[n2], 1., d);
					ResultadosMedios[i].VelocidadSUM += Vble * __cons::ARef * FDeltaTime;
				}
				if(ResultadosMedios[i].Massflow) {
					Dens = FGamma[n1] * __units::BarToPa(FPresion0[n1]) / pow(FAsonido0[n1] * __cons::ARef, 2.);
					//Area=pow(FDiametroTubo[n1],2.)*Pi/4.;
					Gto1 = FVelocidad0[n1] * __cons::ARef * Dens * FArea[n1];
					Dens = FGamma[n2] * __units::BarToPa(FPresion0[n2]) / pow(FAsonido0[n2] * __cons::ARef, 2.);
					//Area=pow(FDiametroTubo[n2],2.)*Pi/4.;
					Gto2 = FVelocidad0[n2] * __cons::ARef * Dens * FArea[n2];
					Vble = InterpolaTubo(Gto1, Gto2, 1., d);
					ResultadosMedios[i].GastoSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].TemperaturaInternaPared) {
					Vble = InterpolaTubo(FTPTubo[0][n1], FTPTubo[0][n2], 2., d);
					ResultadosMedios[i].TemperaturaInternaParedSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].TemperaturaIntermediaPared) {
					Vble = InterpolaTubo(FTPTubo[1][n1], FTPTubo[1][n2], 2., d);
					ResultadosMedios[i].TemperaturaIntermediaParedSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].TemperaturaExternaPared) {
					Vble = InterpolaTubo(FTPTubo[2][n1], FTPTubo[2][n2], 2., d);
					ResultadosMedios[i].TemperaturaExternaParedSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].NITmedio) {
					//double a=InterpolaTubo(FAsonido0[n1],FAsonido0[n2],1.,d);
					double v = InterpolaTubo(FVelocidad0[n1], FVelocidad0[n2], 1., d);
					double p = InterpolaTubo(FPresion0[n1], FPresion0[n2], 1., d);
					double diam = InterpolaTubo(FDiametroTubo[n1], FDiametroTubo[n2], 1., d);

					double nit = CalculaNIT(a, v, p, diam, Gamma, Rmezcla);
					ResultadosMedios[i].NITmedioSUM += nit * FDeltaTime;
				}
				if(ResultadosMedios[i].CoefPelInterior) {
					Vble = InterpolaTubo(Fhi[n1], Fhi[n2], 1., d);
					ResultadosMedios[i].CoefPelInteriorSUM += Vble * FDeltaTime;
				}
				if(ResultadosMedios[i].FraccionMasicaEspecies) {
					for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
						Vble = InterpolaTubo(FFraccionMasicaEspecie[n1][j], FFraccionMasicaEspecie[n2][j], 1., d);
						ResultadosMedios[i].FraccionSUM[j] += Vble * GastoPonderacion * FDeltaTime;
					}
				}

			}
		}

		if(theta > FControlResMed * FDPF->getAnguloTotalCiclo()) {

			for(int i = 0; i < FNumResMedios; i++) {
				if(ResultadosMedios[i].Pressure) {
					ResultadosMedios[i].PresionMED = ResultadosMedios[i].PresionSUM / FTiempoMedSUM;
					ResultadosMedios[i].PresionSUM = 0.;
				}
				if(ResultadosMedios[i].TemperaturaGas) {
					ResultadosMedios[i].TemperaturaGasMED = ResultadosMedios[i].TemperaturaGasSUM / ResultadosMedios[i].GastoPonderacionSUM;
					ResultadosMedios[i].TemperaturaGasSUM = 0.;
					ResultadosMedios[i].GastoPonderacionSUM = 0.;
				}
				if(ResultadosMedios[i].Velocity) {
					ResultadosMedios[i].VelocidadMED = ResultadosMedios[i].VelocidadSUM / FTiempoMedSUM;
					ResultadosMedios[i].VelocidadSUM = 0.;
				}
				if(ResultadosMedios[i].Massflow) {
					ResultadosMedios[i].GastoMED = ResultadosMedios[i].GastoSUM / FTiempoMedSUM;
					ResultadosMedios[i].GastoSUM = 0.;
				}
				if(ResultadosMedios[i].TemperaturaInternaPared) {
					ResultadosMedios[i].TemperaturaInternaParedMED = ResultadosMedios[i].TemperaturaInternaParedSUM / FTiempoMedSUM;
					ResultadosMedios[i].TemperaturaInternaParedSUM = 0.;
				}
				if(ResultadosMedios[i].TemperaturaIntermediaPared) {
					ResultadosMedios[i].TemperaturaIntermediaParedMED = ResultadosMedios[i].TemperaturaIntermediaParedSUM / FTiempoMedSUM;
					ResultadosMedios[i].TemperaturaIntermediaParedSUM = 0.;
				}
				if(ResultadosMedios[i].TemperaturaExternaPared) {
					ResultadosMedios[i].TemperaturaExternaParedMED = ResultadosMedios[i].TemperaturaExternaParedSUM / FTiempoMedSUM;
					ResultadosMedios[i].TemperaturaExternaParedSUM = 0.;
				}
				if(ResultadosMedios[i].NITmedio) {
					ResultadosMedios[i].NITmedioMED = ResultadosMedios[i].NITmedioSUM / FTiempoMedSUM;
					ResultadosMedios[i].NITmedioSUM = 0.;
				}
				if(ResultadosMedios[i].CoefPelInterior) {
					ResultadosMedios[i].CoefPelInteriorMED = ResultadosMedios[i].CoefPelInteriorSUM / FTiempoMedSUM;
					ResultadosMedios[i].CoefPelInteriorSUM = 0.;
				}
				if(ResultadosMedios[i].FraccionMasicaEspecies) {
					for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
						if(ResultadosMedios[i].PonderacionSUM == 0.) {
							ResultadosMedios[i].FraccionMED[j] = 0.;
						} else {
							ResultadosMedios[i].FraccionMED[j] = ResultadosMedios[i].FraccionSUM[j] / ResultadosMedios[i].PonderacionSUM;
						}
						ResultadosMedios[i].FraccionSUM[j] = 0.;
					}
					ResultadosMedios[i].PonderacionSUM = 0.;
				}
			}
			FTiempoMedSUM = 0.;
			FControlResMed = FControlResMed + 1.;
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaResultadosMedios en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaResultadosInstantaneos() {
	double dist = 0., d = 0.;
	int n1 = 0, n2 = 0;

	try {

		if(FNumResInstant != 0) {
			for(int i = 0; i < FNumResInstant; i++) {
				dist = ResultInstantaneos[i].Distancia / FXref;
				n1 = (int) floor(dist);
				if(n1 >= FNin - 1) {
					if(ResultInstantaneos[i].Pressure)
						ResultInstantaneos[i].PresionINS = FPresion0[FNin - 1];
					if(ResultInstantaneos[i].Velocity)
						ResultInstantaneos[i].VelocidadINS = FVelocidad0[FNin - 1] * __cons::ARef;
					if(ResultInstantaneos[i].TemperaturaGas) {
						double temp = __units::KTodegC(pow(FAsonido0[FNin - 1] * __cons::ARef, 2.) / (FGamma[FNin - 1] * FRMezcla[FNin - 1]));
						ResultInstantaneos[i].TemperaturaGasINS = temp;
					}
					if(ResultInstantaneos[i].FlujoMasico) {
						//double area=pow(FDiametroTubo[FNin-1],2)*Pi/4.;
						double dens = FGamma[FNin - 1] * __units::BarToPa(FPresion0[FNin - 1]) / pow(FAsonido0[FNin - 1] * __cons::ARef, 2.);
						double gto = FArea[FNin - 1] * dens * FVelocidad0[FNin - 1] * __cons::ARef;
						ResultInstantaneos[i].FlujoMasicoINS = gto;
					}
					if(ResultInstantaneos[i].VelocidadDerecha) {
						double ason = FAsonido0[FNin - 1] * __cons::ARef;
						double vel = FGamma1[FNin - 1] / 2 * FVelocidad0[FNin - 1] * __cons::ARef;
						double Aa = ason / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]) / __cons::ARef;
						double VelDer = ((ason + vel) / __cons::ARef - Aa) / FGamma1[FNin - 1] * __cons::ARef;
						ResultInstantaneos[i].VelocidadDerechaINS = VelDer;
					}
					if(ResultInstantaneos[i].VelocidadIzquierda) {
						double ason = FAsonido0[FNin - 1] * __cons::ARef;
						double vel = FGamma1[FNin - 1] / 2 * FVelocidad0[FNin - 1] * __cons::ARef;
						double Aa = ason / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]) / __cons::ARef;
						double VelIzq = -((ason - vel) / __cons::ARef - Aa) / FGamma1[FNin - 1] * __cons::ARef;
						ResultInstantaneos[i].VelocidadIzquierdaINS = VelIzq;
					}
					if(ResultInstantaneos[i].PresionDerecha) {
						double ason = FAsonido0[FNin - 1] * __cons::ARef;
						double vel = FGamma1[FNin - 1] / 2 * FVelocidad0[FNin - 1] * __cons::ARef;
						double Aa = ason / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]) / __cons::ARef;
						double PreDer = pow(((ason + vel) / __cons::ARef / Aa + 1) / 2., FGamma4[FNin - 1]);
						ResultInstantaneos[i].PresionDerechaINS = PreDer;
					}
					if(ResultInstantaneos[i].PresionIzquierda) {
						double ason = FAsonido0[FNin - 1] * __cons::ARef;
						double vel = FGamma1[FNin - 1] / 2 * FVelocidad0[FNin - 1] * __cons::ARef;
						double Aa = ason / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]) / __cons::ARef;
						double PreIzq = pow(((ason - vel) / __cons::ARef / Aa + 1) / 2., FGamma4[FNin - 1]);
						ResultInstantaneos[i].PresionIzquierdaINS = PreIzq;
					}
					if(ResultInstantaneos[i].NIT) {
						double nit = CalculaNIT(FAsonido0[FNin - 1], FVelocidad0[FNin - 1], FPresion0[FNin - 1], FDiametroTubo[FNin - 1],
												FGamma[FNin - 1], FRMezcla[FNin - 1]);
						ResultInstantaneos[i].NITINS = nit;
					}
					if(ResultInstantaneos[i].TemperaturaInternaPared)
						ResultInstantaneos[i].TemperaturaInternaParedINS = FTPTubo[0][FNin - 1];
					if(ResultInstantaneos[i].TemperaturaIntermediaPared)
						ResultInstantaneos[i].TemperaturaIntermediaParedINS = FTPTubo[1][FNin - 1];
					if(ResultInstantaneos[i].TemperaturaExternaPared)
						ResultInstantaneos[i].TemperaturaExternaParedINS = FTPTubo[2][FNin - 1];
					if(ResultInstantaneos[i].CoefPelInterior)
						ResultInstantaneos[i].CoefPelInteriorINS = Fhi[FNin - 1];
					if(ResultInstantaneos[i].FraccionMasicaEspecies) {
						for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
							ResultInstantaneos[i].FraccionINS[j] = FFraccionMasicaEspecie[FNin - 1][j];
						}
					}
					if(ResultInstantaneos[i].Gamma)
						ResultInstantaneos[i].GammaINS = FGamma[FNin - 1];
				} else {
					n2 = n1 + 1;
					d = dist - (double) n1;
					if(ResultInstantaneos[i].Pressure) {
						double pres = InterpolaTubo(FPresion0[n1], FPresion0[n2], 1., d);
						ResultInstantaneos[i].PresionINS = pres;
					}
					if(ResultInstantaneos[i].Velocity) {
						double vel = InterpolaTubo(FVelocidad0[n1], FVelocidad0[n2], 1., d);
						ResultInstantaneos[i].VelocidadINS = vel * __cons::ARef;
					}
					if(ResultInstantaneos[i].TemperaturaGas) {
						double temp1 = pow(FAsonido0[n1] * __cons::ARef, 2.) / (FGamma[n1] * FRMezcla[n1]);
						double temp2 = pow(FAsonido0[n2] * __cons::ARef, 2.) / (FGamma[n2] * FRMezcla[n2]);
						double temp = InterpolaTubo(temp1, temp2, 1., d);
						ResultInstantaneos[i].TemperaturaGasINS = __units::KTodegC(temp);
					}
					if(ResultInstantaneos[i].FlujoMasico) {
						//double area1=pow(FDiametroTubo[n1],2)*Pi/4.;
						double dens1 = FGamma[n1] * __units::BarToPa(FPresion0[n1]) / pow(FAsonido0[n1] * __cons::ARef, 2.);
						double gto1 = FArea[n1] * dens1 * FVelocidad0[n1] * __cons::ARef;
						//double area2=pow(FDiametroTubo[n2],2)*Pi/4.;
						double dens2 = FGamma[n2] * __units::BarToPa(FPresion0[n2]) / pow(FAsonido0[n2] * __cons::ARef, 2.);
						double gto2 = FArea[n2] * dens2 * FVelocidad0[n2] * __cons::ARef;
						double gto = InterpolaTubo(gto1, gto2, 1., d);
						ResultInstantaneos[i].FlujoMasicoINS = gto;
					}
					if(ResultInstantaneos[i].VelocidadDerecha) {
						double ason1 = FAsonido0[n1] * __cons::ARef;
						double vel1 = FGamma1[n1] / 2 * FVelocidad0[n1] * __cons::ARef;
						double Aa1 = ason1 / pow(FPresion0[n1], FGamma5[n1]) / __cons::ARef;
						double VelDer1 = ((ason1 + vel1) / __cons::ARef - Aa1) / FGamma1[n1] * __cons::ARef;
						double ason2 = FAsonido0[n2] * __cons::ARef;
						double vel2 = FGamma1[n2] / 2 * FVelocidad0[n2] * __cons::ARef;
						double Aa2 = ason2 / pow(FPresion0[n2], FGamma5[n2]) / __cons::ARef;
						double VelDer2 = ((ason2 + vel2) / __cons::ARef - Aa2) / FGamma1[n2] * __cons::ARef;
						double VelDer = InterpolaTubo(VelDer1, VelDer2, 1., d);
						ResultInstantaneos[i].VelocidadDerechaINS = VelDer;
					}
					if(ResultInstantaneos[i].VelocidadIzquierda) {
						double ason1 = FAsonido0[n1] * __cons::ARef;
						double vel1 = FGamma1[n1] / 2 * FVelocidad0[n1] * __cons::ARef;
						double Aa1 = ason1 / pow(FPresion0[n1], FGamma5[n1]) / __cons::ARef;
						double VelIzq1 = -((ason1 - vel1) / __cons::ARef - Aa1) / FGamma1[n1] * __cons::ARef;
						double ason2 = FAsonido0[n2] * __cons::ARef;
						double vel2 = FGamma1[n2] / 2 * FVelocidad0[n2] * __cons::ARef;
						double Aa2 = ason2 / pow(FPresion0[n2], FGamma5[n2]) / __cons::ARef;
						double VelIzq2 = -((ason2 - vel2) / __cons::ARef - Aa2) / FGamma1[n2] * __cons::ARef;
						double VelIzq = InterpolaTubo(VelIzq1, VelIzq2, 1., d);
						ResultInstantaneos[i].VelocidadIzquierdaINS = VelIzq;
					}
					if(ResultInstantaneos[i].PresionDerecha) {
						double ason1 = FAsonido0[n1] * __cons::ARef;
						double vel1 = FGamma1[n1] / 2 * FVelocidad0[n1] * __cons::ARef;
						double Aa1 = ason1 / pow(FPresion0[n1], FGamma5[n1]) / __cons::ARef;
						double PreDer1 = pow(((ason1 + vel1) / __cons::ARef / Aa1 + 1) / 2., FGamma4[n1]);
						double ason2 = FAsonido0[n2] * __cons::ARef;
						double vel2 = FGamma1[n2] / 2 * FVelocidad0[n2] * __cons::ARef;
						double Aa2 = ason2 / pow(FPresion0[n2], FGamma5[n2]) / __cons::ARef;
						double PreDer2 = pow(((ason2 + vel2) / __cons::ARef / Aa2 + 1) / 2., FGamma4[n2]);
						double PreDer = InterpolaTubo(PreDer1, PreDer2, 1., d);
						ResultInstantaneos[i].PresionDerechaINS = PreDer;
					}
					if(ResultInstantaneos[i].PresionIzquierda) {
						double ason1 = FAsonido0[n1] * __cons::ARef;
						double vel1 = FGamma1[n1] / 2 * FVelocidad0[n1] * __cons::ARef;
						double Aa1 = ason1 / pow(FPresion0[n1], FGamma5[n1]) / __cons::ARef;
						double PreIzq1 = pow(((ason1 - vel1) / __cons::ARef / Aa1 + 1) / 2., FGamma4[n1]);
						double ason2 = FAsonido0[n2] * __cons::ARef;
						double vel2 = FGamma1[n2] / 2 * FVelocidad0[n2] * __cons::ARef;
						double Aa2 = ason2 / pow(FPresion0[n2], FGamma5[n2]) / __cons::ARef;
						double PreIzq2 = pow(((ason2 - vel2) / __cons::ARef / Aa2 + 1) / 2., FGamma4[n2]);
						double PreIzq = InterpolaTubo(PreIzq1, PreIzq2, 1., d);
						ResultInstantaneos[i].PresionIzquierdaINS = PreIzq;
					}
					if(ResultInstantaneos[i].NIT) {
						double nit1 = CalculaNIT(FAsonido0[n1], FVelocidad0[n1], FPresion0[n1], FDiametroTubo[n1], FGamma[n1], FRMezcla[n1]);
						double nit2 = CalculaNIT(FAsonido0[n2], FVelocidad0[n2], FPresion0[n2], FDiametroTubo[n2], FGamma[n2], FRMezcla[n2]);
						double nit = InterpolaTubo(nit1, nit2, 1., d);
						ResultInstantaneos[i].NITINS = nit;
					}
					if(ResultInstantaneos[i].TemperaturaInternaPared) {
						double TP = InterpolaTubo(FTPTubo[0][n1], FTPTubo[0][n2], 1., d);
						ResultInstantaneos[i].TemperaturaInternaParedINS = TP;
					}
					if(ResultInstantaneos[i].TemperaturaIntermediaPared) {
						double TP = InterpolaTubo(FTPTubo[1][n1], FTPTubo[1][n2], 1., d);
						ResultInstantaneos[i].TemperaturaIntermediaParedINS = TP;
					}
					if(ResultInstantaneos[i].TemperaturaExternaPared) {
						double TP = InterpolaTubo(FTPTubo[2][n1], FTPTubo[2][n2], 1., d);
						ResultInstantaneos[i].TemperaturaExternaParedINS = TP;
					}
					if(ResultInstantaneos[i].CoefPelInterior) {
						double hi = InterpolaTubo(Fhi[n1], Fhi[n2], 1., d);
						ResultInstantaneos[i].CoefPelInteriorINS = hi;
					}
					if(ResultInstantaneos[i].FraccionMasicaEspecies) {
						for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
							double Fraccion = InterpolaTubo(FFraccionMasicaEspecie[n1][j], FFraccionMasicaEspecie[n2][j], 1., d);
							ResultInstantaneos[i].FraccionINS[j] = Fraccion;
						}
					}
					if(ResultInstantaneos[i].Gamma) {
						double gamma = InterpolaTubo(FGamma[n1], FGamma[n2], 1., d);
						ResultInstantaneos[i].GammaINS = gamma;
					}
				}
			}

		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaResultadosInstantaneos en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::CalculaNIT(double a, double v, double p, double d, double Gamma, double R) {
	try {
		double kp = 0.;
		double tem0 = 0., pre0 = 0., area = 0., gto = 0., nit = 0., V = 0., A = 0., V2 = 0., A2 = 0.;

		kp = Gamma * R / (Gamma - 1);
		V = v * __cons::ARef;
		V2 = V * V;
		A = a * __cons::ARef;
		A2 = A * A;
		tem0 = A2 / (Gamma * R) + V2 / 2. / kp;
		pre0 = __units::BarToPa(p) * pow((1 + V2 / 2. / kp / A2), (kp / 287.));
		area = __geom::Circle_area(d);
		gto = Gamma * __units::BarToPa(p) * area * V / A2;
		nit = gto * kp * tem0 * (1 - pow(pre0 / 100000., (-287. / kp)));
		return nit;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaNIT en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaCoeficientePeliculaExterior(double pamb) {
	try {

		double dtem = 0., temed = 0., rhog = 0., viscext = 0., Re = 0., Pr = 0., cond = 0., vel = 0., n1 = 0., n2 = 0., L = 0.;

		if(FDPF->getTipoCalcTempPared() != nmTempConstante && FDPF->getCoefAjusTC() != 0) {
			for(int i = 0; i < FNin; i++) {
				dtem = fabs(__units::degCToK(FDPF->GetTSuperficie(i, 1)) - FDPF->getTExt());
				temed = (__units::degCToK(FDPF->GetTSuperficie(i, 1)) + FDPF->getTExt()) / 2.;

				// Calculo de las caracter�sticas del refrigerante
				switch(FDPF->getTipoRefrig()) {
				case nmAire:
					// Condiciones del aire a temperatura media (temed) y 0.1m/s
					rhog = __units::BarToPa(pamb) / __R::Air / temed; // Densidad del aire atmosf�rico (se considera R=cte=287)
					viscext = 1.4615e-6 * pow(temed, 1.5) / (temed + 110.4);
					Re = rhog * 0.7 * FDPF->getDiametroExt() / viscext;
					Pr = 0.7;
					cond = -8.39061e-09 * pow(temed, 2) + 7.05256e-05 * temed + 6.51528e-03;
					break;
				case nmAgua:
					// Condiciones de agua saturada a (temed) y 1m/s
					rhog = 980;
					viscext = -2.632351E-09 * pow(temed, 3.) + 2.737629E-06 * pow(temed, 2.) - 9.530709E-04 * temed + 1.114642E-01;
					Re = rhog * 1. * FDPF->getDiametroExt() / viscext;
					Pr = -2.022269E-05 * pow(temed, 3.) + 2.106518E-02 * pow(temed, 2.) - 7.340298E+00 * temed + 8.581110E+02;
					cond = 9.496332E-09 * pow(temed, 3.) - 1.707697E-05 * pow(temed, 3.) + 9.183462E-03 * temed - 8.626578E-01;
					break;
				default:
					std::cout << "WARNING: Tipo de refrigeraci�n mal definida en el tubo: " << FNumeroHaz << std::endl;
				}

				// Calculo del coeficiente de pel�cula exterior
				// Termino de conveccion de Churchill Bernstein
				if((2e4 < Re) && (Re < 4e5)) {
					n1 = 1 / 2;
					n2 = 1;
				} else {
					n1 = 5 / 8;
					n2 = 4 / 5;
				}
				Fhe[i] = 0.3
						 + 0.62 * pow(Re, 0.5) * pow(Pr, 0.333333) / pow(1 + pow(0.4 / Pr, 0.666666), 0.25) * pow(1 + pow(Re / 282000, n1),
								 n2) * cond
						 / (FDPF->GetDiametroExtHaz(FNumeroHaz - 1) + 2 * FDPF->getEspesorAislante() + 2 * FDPF->getEspesorMetal() + 2 *
							FDPF->getEspesorAire());
				// Termino de radiaci�n
				if(dtem != 0.) {
					Fhe[i] = Fhe[i] + __cons::Sigma * FDPF->getEmisividad() * (pow(__units::degCToK(FDPF->GetTSuperficie(i, 1)),
							 4) - pow(FDPF->getTExt(), 4.)) / dtem;
				}
				Fhe[i] = Fhe[i] * FDPF->getCoefExt();

			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TTubo::CalculaCoeficientePeliculaExterior en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaCoeficientePeliculaInterior() {
	try {
		double Tg = 0., cesp = 0., viscgas = 0., cond = 0., zzz = 0., czz = 0., uq1 = 0., Nu = 0.;

		if(FDPF->getCoefAjusTC() != 0) {
			for(int i = 0; i < FNin; i++) {
				Tg = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
				cesp = FGamma[i] * FRMezcla[i] / FGamma1[i];
				FViscosidadDinamica[i] = 1.4615e-6 * pow(Tg, 1.5) / (Tg + 110.4);
				cond = FViscosidadDinamica[i] * cesp / 0.709;
				if(FTime1 != 0) {
					zzz = 0.013 / (FTime1 - FTime0);
					czz = 2 / (zzz + 1);
					uq1 = fabs(FVelocidad0[i] * __cons::ARef);
					FVelPro[i] = czz * uq1 + (1 - czz) * FVelPro[i];
					//FVelPro[i]=5.;
				}
				if(FTipoCanal == nmCanalEntrada) {
					FRe[i] = Frho[i] * FVelPro[i] * (FLadoCanal[i] - 2 * FEspesorSoot[i]) / FViscosidadDinamica[i];
					//Fhi[i]=3.08*cond/FLadoCanal[i]*FDPF->getCoefAjusTC();
					//Nu=0.571*pow(FRe[i]*(FLadoCanal[i]-2*FEspesorSoot[i])/FDPF->getLongitudEfec(),2./3.)*0.2;
					Nu = 0.571 * pow(FRe[i] * (FLadoCanal[i] - 2 * FEspesorSoot[i]) / FDPF->getLongitudEfec(), 2. / 3.);
					//Fhi[i]=Nu*cond/(FLadoCanal[i]-2*FEspesorSoot[i])*FDPF->getCoefAjusTC();
					Fhi[i] = Nu * cond / (FLadoCanal[i] - 2 * FEspesorSoot[i]);
				} else {
					FRe[i] = Frho[i] * FVelPro[i] * FLadoCanal[i] / FViscosidadDinamica[i];
					//Fhi[i]=3.08*cond/FLadoCanal[i]*FDPF->getCoefAjusTC();
					//Nu=0.571*pow(FRe[i]*FLadoCanal[i]/FDPF->getLongitudEfec(),2./3.)*0.2;
					Nu = 0.571 * pow(FRe[i] * FLadoCanal[i] / FDPF->getLongitudEfec(), 2. / 3.);
					//Fhi[i]=Nu*cond/FLadoCanal[i]*FDPF->getCoefAjusTC();
					Fhi[i] = Nu * cond / FLadoCanal[i];
				}
			}
		} else if(FCoefAjusFric != 0) {
			for(int i = 0; i < FNin; i++) {
				Tg = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
				FViscosidadDinamica[i] = 1.4615e-6 * pow(Tg, 1.5) / (Tg + 110.4);
				if(FTime1 != 0) {
					zzz = 0.013 / (FTime1 - FTime0);
					czz = 2 / (zzz + 1);
					uq1 = fabs(FVelocidad0[i] * __cons::ARef);
					FVelPro[i] = czz * uq1 + (1 - czz) * FVelPro[i];
				}
				FRe[i] = Frho[i] * FVelPro[i] * FLadoCanal[i] / FViscosidadDinamica[i];
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaCoeficientePeliculaInterior en el haz " << FNumeroHaz << "de la DPF" <<
				  std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::SalidaGeneralTubos(stEspecies *DatosEspecies) const {
	try {
		printf("__________________________________\n");
		printf("RESULTADOS MEDIOS DE LOS CONDUCTOS\n");
		printf("__________________________________\n\n\n");
		for(int i = 0; i < FNumResMedios; i++) {
			std::cout << "Valores medios del punto " << i << " situado en el conducto num " << FNumeroHaz;
			std::cout << " a " << ResultadosMedios[i].Distancia << " metros del extremos izquierdo" << std::endl;
			std::cout << "Temperatura media = " << ResultadosMedios[i].TemperaturaGasMED << " C" << std::endl;
			std::cout << "Presion media     = " << ResultadosMedios[i].PresionMED << " bar" << std::endl;
			std::cout << "Velocidad media   = " << ResultadosMedios[i].VelocidadMED << " m/s" << std::endl;
			std::cout << "Gasto Medio       = " << ResultadosMedios[i].GastoMED << " kg/s" << std::endl;
			std::cout << "NIT Medio         = " << ResultadosMedios[i].NITmedioMED << " Watios" << std::endl << std::endl;
			std::cout << "Fraccion M�sica Media de " << DatosEspecies[i].Nombre << ": " << ResultadosMedios[i].FraccionMED[i] <<
					  std::endl;
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::SalidaGeneralTubos en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::Maximo(double x, double y) {
	if(x < y)
		return y;
	else
		return x;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::Minimo(double x, double y) {
	if(x < y)
		return x;
	else
		return y;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::AjustaPaso(double TiempoFinPaso) {
	try {
		FTime1 = TiempoFinPaso;
		FDeltaTime = FTime1 - FTime0;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::AjustaPaso en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaCaracteristicasExtremos(TCondicionContorno **CC, double DeltaTiempo) {
	try {
		double temp1 = 0., temp2 = 0.;
		nmPipeEnd TipoExtremo;
		if(FNodoDer == 0) {
			if(FVelocidad0[0] <= 0) {
				CC[FNodoIzq - 1]->PutEntropia(FTuboCCNodoIzq,
											  Interpola_Entropia(CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).TipoExtremo, DeltaTiempo));
			}
			CC[FNodoIzq - 1]->PutBeta(FTuboCCNodoIzq,
									  Interpola_Caracteristica(CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).Entropia, 1, 0, DeltaTiempo,
											  CC[FNodoIzq - 1]->GetTuboExtremo(FTuboCCNodoIzq).TipoExtremo));

			// Nodo derecho: extremo cerrado
			if(FVelocidad0[FNin - 1] >= 0) {
				TipoExtremo = nmRight;
				//FEntropiaExtremoCerrado=Interpola_Entropia(TipoExtremo,DeltaTiempo);
				temp1 = pow(FAsonido0[FNin - 2] * __cons::ARef, 2.) / (FGamma[FNin - 2] * FRMezcla[FNin - 2]);
				temp2 = pow(FAsonido0[FNin - 3] * __cons::ARef, 2.) / (FGamma[FNin - 3] * FRMezcla[FNin - 3]);
				FEntropiaExtremoCerrado = sqrt((temp1 - (temp2 - temp1) / 2) * FGamma[FNin - 1] * FRMezcla[FNin - 1]) / __cons::ARef /
										  pow(FPresion0[FNin - 1], FGamma5[FNin - 1]);
			}
			FLandaExtremoCerrado = Interpola_Caracteristica(FEntropiaExtremoCerrado, -1, FNin - 1, DeltaTiempo, TipoExtremo);
			FBetaExtremoCerrado = FLandaExtremoCerrado; //Cerrado
			//FBetaExtremoCerrado=FEntropiaExtremoCerrado; // Anecoico
		}

		if(FNodoIzq == 0) {
			if(FVelocidad0[FNin - 1] >= 0) {
				CC[FNodoDer - 1]->PutEntropia(FTuboCCNodoDer,
											  Interpola_Entropia(CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).TipoExtremo, DeltaTiempo));
			}
			CC[FNodoDer - 1]->PutLanda(FTuboCCNodoDer,
									   Interpola_Caracteristica(CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).Entropia, -1, FNin - 1, DeltaTiempo,
											   CC[FNodoDer - 1]->GetTuboExtremo(FTuboCCNodoDer).TipoExtremo));

			// Nodo izquierdo: extremo cerrado
			if(FVelocidad0[0] <= 0) {
				TipoExtremo = nmLeft;
				//FEntropiaExtremoCerrado=Interpola_Entropia(TipoExtremo,DeltaTiempo);
				temp1 = pow(FAsonido0[1] * __cons::ARef, 2.) / (FGamma[1] * FRMezcla[1]);
				temp2 = pow(FAsonido0[2] * __cons::ARef, 2.) / (FGamma[2] * FRMezcla[2]);
				FEntropiaExtremoCerrado = sqrt((temp1 - (temp2 - temp1) / 2) * FGamma[0] * FRMezcla[0]) / __cons::ARef / pow(
											  FPresion0[0], FGamma5[0]);
			}
			FBetaExtremoCerrado = Interpola_Caracteristica(FEntropiaExtremoCerrado, 1, 0, DeltaTiempo, TipoExtremo);
			FLandaExtremoCerrado = FBetaExtremoCerrado; // Cerrado
			//FLandaExtremoCerrado=FEntropiaExtremoCerrado; // Anecoico
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaCaracteristicasExtremos en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::Interpola_Entropia(nmPipeEnd TipoExtremoTubo, double DeltaTiempo) {
	try {

		int signo = 1;
		int extremo = 0;
		int indiceCC = 0;

		if(TipoExtremoTubo == nmRight) {        // Righ pipe end
			signo = -1;
			extremo = FNin - 1;
			indiceCC = 1;
		}
		double dtdx = DeltaTiempo / FXref;
		int ind = extremo;
		double entropia = 0.;
		double velocidadp = 0.;
		double asonidop = 0.;

		if(DeltaTiempo < 1e-15 || FVelocidad0[extremo] == 0.0) {

			Calculo_Entropia(&entropia, &velocidadp, extremo, 0., signo, DeltaTiempo, indiceCC, TipoExtremoTubo);

		} else {
			int ind1 = ind + signo;

			stPathOrigin PathOrigin(FU0[0][ind], FU0[1][ind], FU0[0][ind1], FU0[1][ind1], dtdx, signo);

			double dist = zbrent(PathOrigin, 0., 1., 1e-5);

			Calculo_Entropia(&entropia, &velocidadp, ind, dist, signo, DeltaTiempo, indiceCC, TipoExtremoTubo);
		}
		return entropia;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Interpola_Entropia en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Calculo_Entropia(double *entropia, double *velocidadp, int ind, double dist, int signo,
								 double DeltaTiempo, int indiceCC, nmPipeEnd TipoExtremoTubo) {
	try {
		int ind1 = ind + signo;

		double w0 = FU0[0][ind];
		double w1 = FU0[1][ind];
		double w2 = FU0[2][ind];
		double gammap = FGamma[ind];
		double Rmezclap = FRMezcla[ind];
		double espesorsootp = 0.;
		double tptubop;
		if(FTipoCanal == nmCanalEntrada) {
			tptubop = FDPF->GetTPared(FNumeroHaz - 1, ind, 2);
			espesorsootp = FEspesorSoot[ind];
		} else {
			tptubop = FDPF->GetTPared(FNumeroHaz - 1, ind, 0);
		}
		double hip = Fhi[ind];
		double coffp = FCoefTurbulencia[ind];
		double ladocanalp = FLadoCanal[ind];
		double VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind);
		double VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind);
		double viscdinamicap = FViscosidadDinamica[ind];
		double areap = pow2(ladocanalp - 2 * espesorsootp);
		for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
			FFraccionMasicaCC[indiceCC][j] = FFraccionMasicaEspecie[ind][j];
		}

		if(dist > 0. || dist < 1.) {
			w0 = Interpola(FU0[0][ind], FU0[0][ind1], 1., dist);
			w1 = Interpola(FU0[1][ind], FU0[1][ind1], 1., dist);
			w2 = Interpola(FU0[2][ind], FU0[2][ind1], 1., dist);
			gammap = Interpola(FGamma[ind], FGamma[ind1], 1., dist);
			Rmezclap = Interpola(FRMezcla[ind], FRMezcla[ind1], 1., dist);
			ladocanalp = InterpolaTubo(FLadoCanal[ind], FLadoCanal[ind1], 1., dist);
			if(FTipoCanal == nmCanalEntrada) {
				espesorsootp = InterpolaTubo(FEspesorSoot[ind], FEspesorSoot[ind1], 1., dist);
			}
			for(int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
				FFraccionMasicaCC[indiceCC][j] = InterpolaTubo(FFraccionMasicaEspecie[ind][j], FFraccionMasicaEspecie[ind1][j], 1.,
												 dist);
			}
			if(FTipoCanal == nmCanalEntrada && TipoExtremoTubo == nmRight) {
				VelocidadParedEntradap = InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind),
													   FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind1), 1., dist);
				VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
										1)->getNodoFinalFuente());
			} else if(FTipoCanal == nmCanalSalida && TipoExtremoTubo == nmRight) {
				if(FXref <= FDPF->getLongitudTapon()) {
					VelocidadParedEntradap = 0.;
					VelocidadParedSalidap = 0.;
				} else {
					VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
											 0)->getNodoFinalFuente());
					VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
											1)->getNodoFinalFuente());
				}
			} else if(FTipoCanal == nmCanalEntrada && TipoExtremoTubo == nmLeft) {
				if(FXref <= FDPF->getLongitudTapon()) {
					VelocidadParedEntradap = 0.;
					VelocidadParedSalidap = 0.;
				} else {
					VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
											 0)->getNodoInicialFuente());
					VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
											1)->getNodoInicialFuente());
				}
			} else if(FTipoCanal == nmCanalSalida && TipoExtremoTubo == nmLeft) {
				VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
										 0)->getNodoInicialFuente());
				VelocidadParedSalidap = InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind),
													  FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind1), 1., dist);
			}
			//VelocidadParedEntradap=InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz-1,0,ind),FDPF->GetVelocidadPared(FNumeroHaz-1,0,ind1),1.,dist);
			//VelocidadParedSalidap=InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz-1,1,ind),FDPF->GetVelocidadPared(FNumeroHaz-1,1,ind1),1.,dist);
			viscdinamicap = InterpolaTubo(FViscosidadDinamica[ind], FViscosidadDinamica[ind1], 1., dist);
			if(FDPF->getCoefAjusTC() != 0 || FCoefAjusFric != 0) {
				if(FTipoCanal == nmCanalEntrada) {
					tptubop = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, ind, 2), FDPF->GetTPared(FNumeroHaz - 1, ind1, 2), 1., dist);
				} else {
					tptubop = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, ind, 0), FDPF->GetTPared(FNumeroHaz - 1, ind1, 0), 1., dist);
				}
				hip = InterpolaTubo(Fhi[ind], Fhi[ind1], 1., dist);
				coffp = InterpolaTubo(FCoefTurbulencia[ind], FCoefTurbulencia[ind1], 1., dist);
			}
			areap = pow2(ladocanalp - 2 * espesorsootp);
		}
		double gamma1p = __Gamma::G1(gammap);
		double gamma3p = __Gamma::G3(gammap);
		double gamma5p = __Gamma::G5(gammap);
//*velocidadp=InterpolaTubo(FVelocidad0[ind],FVelocidad0[ind1],1.,dist);
		*velocidadp = w1 / w0 / __cons::ARef;
//double presionp=InterpolaTubo(FPresion0[ind],FPresion0[ind1],1.,dist);
		double rhop = w0 / areap;
//double rhop=InterpolaTubo(Frho[ind],Frho[ind1],1.,dist);
		double asonidop = sqrt(gammap * gamma1p * (w2 / w0 - pow2(*velocidadp * __cons::ARef) / 2)) / __cons::ARef;
		double presionp = __units::PaToBar((w2 - pow2(w1) / 2. / w0) * gamma1p / areap);
//double asonidop=InterpolaTubo(FAsonido0[ind],FAsonido0[ind1],1.,dist);
		double entropiap = asonidop / pow(presionp, gamma5p);

		*entropia = entropiap;

		/*variacion de la entropia debida a la transmision del calor*/
		/*------------------------------------------*/
		if(FDPF->getCoefAjusTC() != 0) {
			double q = 0;
			double tgasp = pow2(asonidop * __cons::ARef) / (gammap * Rmezclap);

			TransmisionCalor(tgasp, ladocanalp - 2 * espesorsootp, &q, coffp, hip, rhop, tptubop);

			//Las siguientes expresiones est�n en la Tesis de Corber�n. P�gina 23
			double dacal = gamma3p * entropiap * q * FDPF->getCoefAjusTC() * DeltaTiempo / pow2(asonidop * __cons::ARef);

			*entropia += dacal;
		}

		/*variacion de la entropia debida al t�rmino de friccion*/
		/*---------------------------------------*/
		double velabs = fabs(*velocidadp);
		if(velabs != 0. & FCoefAjusFric != 0) {
			double dafric = gamma1p * FCoefAjusFric * entropiap * FF * viscdinamicap * DeltaTiempo * pow2(
								velabs / asonidop) / areap;
			*entropia += dafric;
		}

		double deltat = 0.;
		if((TipoExtremoTubo == nmLeft && FNodoDer == 0) || (TipoExtremoTubo == nmRight && FNodoIzq == 0)) {
			if(velabs != 0 && dist * FXref > FDPF->getLongitudTapon()) {
				deltat = (dist * FXref - FDPF->getLongitudTapon()) / ((*velocidadp) * __cons::ARef);
			}
		} else {
			deltat = DeltaTiempo;
		}

		/*variaci�n debida a la variaci�n de momento por el flujo que atraviesa la pared porosa */
		/*---------------------------------------*/
		double var_momento = 0.;
		if(FTipoCanal == nmCanalEntrada) {
			var_momento = -gamma1p * entropiap * pow2(*velocidadp) / pow2(asonidop) * 2. * VelocidadParedEntradap /
						  (ladocanalp - 2. * espesorsootp) * deltat;
		} else {
			var_momento = gamma1p * entropiap * pow2(*velocidadp) / pow2(asonidop) * 2. * VelocidadParedSalidap / ladocanalp *
						  deltat;
		}
		*entropia += var_momento;

		/*variaci�n debida a la variaci�n de energ�a por el flujo que atraviesa la pared porosa */
		/*---------------------------------------*/
		double var_energia = 0.;
		if(FTipoCanal == nmCanalEntrada) {
			var_energia = gamma1p * deltat * entropiap * VelocidadParedEntradap / (ladocanalp - 2. * espesorsootp) / pow2(
							  asonidop) * (pow2(*velocidadp) - pow2(VelocidadParedEntradap / __cons::ARef));
		} else {
			var_energia = gamma1p * deltat * entropiap * VelocidadParedSalidap / ladocanalp / pow2(asonidop) * (-pow2(
							  *velocidadp) + pow2(VelocidadParedSalidap / __cons::ARef));
		}
		*entropia += var_energia;
	}

	catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Calculo_Entropia en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::Interpola_Caracteristica(double entropia, int signo, int extremo, double DeltaTiempo,
		nmPipeEnd TipoExtremoTubo) {
	try {

		double dtdx = DeltaTiempo / FXref;
		int ind = extremo;
		double caracteristica = 0.;
		double velocidadp = 0.;
		double asonidop = 0.;

		if(DeltaTiempo < 1e-15) {
			Calculo_Caracteristica(&caracteristica, &velocidadp, &asonidop, ind, 0., signo, entropia, DeltaTiempo, TipoExtremoTubo);
		} else {
			int ind1 = ind + signo;

			stCharOrigin CharOrigin(FU0[0][ind], FU0[1][ind], FU0[2][ind], FU0[0][ind1], FU0[1][ind1], FU0[2][ind1], FGamma[ind],
									FGamma[ind1], dtdx, signo);

			double dist = zbrent(CharOrigin, 0., 1., 1e-5);

			Calculo_Caracteristica(&caracteristica, &velocidadp, &asonidop, ind, dist, signo, entropia, DeltaTiempo,
								   TipoExtremoTubo);
		}
		return caracteristica;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Interpola_Caracteristica en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::Calculo_Caracteristica(double *caracteristica, double *velocidadp, double *asonidop, int ind,
									   double dist, int signo, double entropia, double DeltaTiempo, nmPipeEnd TipoExtremoTubo) {
	try {

		double w0 = FU0[0][ind];
		double w1 = FU0[1][ind];
		double w2 = FU0[2][ind];

		int ind1 = ind + signo;

		double gammap = FGamma[ind];
		double Rmezclap = FRMezcla[ind];
		double espesorsootp = 0.;
		double tptubop;
		if(FTipoCanal == nmCanalEntrada) {
			tptubop = FDPF->GetTPared(FNumeroHaz - 1, ind, 2);
			espesorsootp = FEspesorSoot[ind];
		} else {
			tptubop = FDPF->GetTPared(FNumeroHaz - 1, ind, 0);
		}
		double hip = Fhi[ind];
		double coffp = FCoefTurbulencia[ind];
		double ladocanalp = FLadoCanal[ind];
		double VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind);
		double VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind);
		double viscdinamicap = FViscosidadDinamica[ind];
		double areap = pow2(ladocanalp - 2 * espesorsootp);

		if(dist < 1. && dist > 0) {
			w0 = Interpola(FU0[0][ind], FU0[0][ind1], 1., dist);
			w1 = Interpola(FU0[1][ind], FU0[1][ind1], 1., dist);
			w2 = Interpola(FU0[2][ind], FU0[2][ind1], 1., dist);
			gammap = Interpola(FGamma[ind], FGamma[ind1], 1., dist);
			Rmezclap = Interpola(FRMezcla[ind], FRMezcla[ind1], 1., dist);
			ladocanalp = InterpolaTubo(FLadoCanal[ind], FLadoCanal[ind1], 1., dist);
			if(FTipoCanal == nmCanalEntrada && TipoExtremoTubo == nmRight) {
				VelocidadParedEntradap = InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind),
													   FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, ind1), 1., dist);
				VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
										1)->getNodoFinalFuente());
			} else if(FTipoCanal == nmCanalSalida && TipoExtremoTubo == nmRight) {
				if(FXref <= FDPF->getLongitudTapon()) {
					VelocidadParedEntradap = 0.;
					VelocidadParedSalidap = 0.;
				} else {
					VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
											 0)->getNodoFinalFuente());
					VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
											1)->getNodoFinalFuente());
				}
			} else if(FTipoCanal == nmCanalEntrada && TipoExtremoTubo == nmLeft) {
				if(FXref <= FDPF->getLongitudTapon()) {
					VelocidadParedEntradap = 0.;
					VelocidadParedSalidap = 0.;
				} else {
					VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
											 0)->getNodoInicialFuente());
					VelocidadParedSalidap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, FDPF->GetCanal(FNumeroHaz - 1,
											1)->getNodoInicialFuente());
				}
			} else if(FTipoCanal == nmCanalSalida && TipoExtremoTubo == nmLeft) {
				VelocidadParedEntradap = FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, FDPF->GetCanal(FNumeroHaz - 1,
										 0)->getNodoInicialFuente());
				VelocidadParedSalidap = InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind),
													  FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, ind1), 1., dist);
			}
			//VelocidadParedEntradap=InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz-1,0,ind),FDPF->GetVelocidadPared(FNumeroHaz-1,0,ind1),1.,dist);
			//VelocidadParedSalidap=InterpolaTubo(FDPF->GetVelocidadPared(FNumeroHaz-1,1,ind),FDPF->GetVelocidadPared(FNumeroHaz-1,1,ind1),1.,dist);
			viscdinamicap = InterpolaTubo(FViscosidadDinamica[ind], FViscosidadDinamica[ind1], 1., dist);
			if(FDPF->getCoefAjusTC() != 0 || FCoefAjusFric != 0) {
				if(FTipoCanal == nmCanalEntrada) {
					tptubop = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, ind, 2), FDPF->GetTPared(FNumeroHaz - 1, ind, 2), 1., dist);
					espesorsootp = InterpolaTubo(FEspesorSoot[ind], FEspesorSoot[ind1], 1., dist);
				} else {
					tptubop = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, ind, 0), FDPF->GetTPared(FNumeroHaz - 1, ind1, 0), 1., dist);
				}
				hip = InterpolaTubo(Fhi[ind], Fhi[ind1], 1., dist);
				coffp = InterpolaTubo(FCoefTurbulencia[ind], FCoefTurbulencia[ind1], 1., dist);
			}
			areap = pow2(ladocanalp - 2 * espesorsootp);
		}
		double gamma1p = __Gamma::G1(gammap);
		double gamma3p = __Gamma::G3(gammap);
		double gamma5p = __Gamma::G5(gammap);
		double rhop = w0 / areap;
		*velocidadp = w1 / w0 / __cons::ARef;
		*asonidop = sqrt(gammap * gamma1p * (w2 / w0 - pow2(*velocidadp * __cons::ARef) / 2)) / __cons::ARef;
//*velocidadp=InterpolaTubo(FVelocidad0[ind],FVelocidad0[ind1],1.,dist);
//double rhop=InterpolaTubo(Frho[ind],Frho[ind1],1.,dist);
//*asonidop=InterpolaTubo(FAsonido0[ind],FAsonido0[ind1],1.,dist);
		*caracteristica = (*asonidop - signo * gamma3p * *velocidadp);

//Las siguientes expresiones se pueden encontrar en la Tesis de Corber�n
//P�gina 22
		/*variacion debida a la transmision del calor*/
		/*------------------------------------------*/
		if(FDPF->getCoefAjusTC() != 0) {
			double q = 0;

			double tgasp = pow2(*asonidop * __cons::ARef) / (gammap * Rmezclap);

			TransmisionCalor(tgasp, ladocanalp - 2 * espesorsootp, &q, coffp, hip, rhop, tptubop);

			double dacal = gamma3p * gamma1p * q * FDPF->getCoefAjusTC() * DeltaTiempo / (pow2(__cons::ARef) * *asonidop);

			*caracteristica += dacal;
		}

		/*variacion debida a la variacion entropia*/
		/*----------------------------------------*/
		double presionp = __units::PaToBar((w2 - pow2(w1) / 2. / w0) * gamma1p / areap);
		double entropiap = *asonidop / pow(presionp, gamma5p);
		double increentropia = entropia - entropiap;
		double daen = *asonidop * increentropia / entropiap;
		*caracteristica += daen;

		/*variacion debida al cambio de seccion*/
		/*-------------------------------------*/
		double daar = 0.;
		if(signo == 1) {
			daar = -gamma3p * *asonidop * *velocidadp * pow2(FLadoCanal[ind1] - FLadoCanal[ind] - 2 *
					(FEspesorSoot[ind1] - FEspesorSoot[ind])) * DeltaTiempo * __cons::ARef
				   / (pow2(ladocanalp - 2 * espesorsootp) * FXref);
		} else if(signo == -1) {
			daar = -gamma3p * *asonidop * *velocidadp * pow(FLadoCanal[ind] - FLadoCanal[ind1] - 2 *
					(FEspesorSoot[ind] - FEspesorSoot[ind1]), 2.) * DeltaTiempo * __cons::ARef
				   / (pow2(ladocanalp - 2 * espesorsootp) * FXref);
		}
		*caracteristica += daar;

		/*variacion debida al termino de friccion*/
		/*---------------------------------------*/
		double velabs = fabs(*velocidadp);
		if(velabs != 0. & FCoefAjusFric != 0) {
			double dafric = gamma3p * FCoefAjusFric * (signo + gamma1p * *velocidadp / *asonidop) * FF * viscdinamicap * DeltaTiempo
							* *velocidadp / pow2(ladocanalp - 2. * espesorsootp);

			*caracteristica += dafric;
		}

		double deltat = 0.;
		if((TipoExtremoTubo == nmLeft && FNodoDer == 0) || (TipoExtremoTubo == nmRight && FNodoIzq == 0)) {
			if(velabs != 0 && dist * FXref > FDPF->getLongitudTapon()) {
				deltat = (dist * FXref - FDPF->getLongitudTapon()) / ((*asonidop - signo * *velocidadp) * __cons::ARef);
			}
		} else {
			deltat = DeltaTiempo;
		}

		/*variaci�n debida a la variaci�n de momento por el flujo que atraviesa la pared porosa */
		/*---------------------------------------*/
		double var_momento = 0.;
		if(FTipoCanal == nmCanalEntrada) {
			var_momento = gamma1p * 2. * VelocidadParedEntradap * *velocidadp / (ladocanalp - 2. * espesorsootp) *
						  (-gamma1p * *velocidadp / (*asonidop) - signo) * deltat;
		} else {
			var_momento = gamma1p * 2. * VelocidadParedSalidap * *velocidadp / ladocanalp * (gamma1p * *velocidadp /
						  (*asonidop) + signo) * deltat;
		}
		*caracteristica += var_momento;

		/*variaci�n debida a la variaci�n de energ�a por el flujo que atraviesa la pared porosa */
		/*---------------------------------------*/
		double var_energia = 0.;
		if(FTipoCanal == nmCanalEntrada) {
			var_energia = pow(gamma1p, 2.) * deltat * VelocidadParedEntradap / (ladocanalp - 2. * espesorsootp) / (*asonidop) *
						  (pow(*velocidadp, 2.) - pow(VelocidadParedEntradap / __cons::ARef, 2.));
		} else {
			var_energia = pow(gamma1p, 2.) * deltat * VelocidadParedSalidap / ladocanalp / (*asonidop) * (-pow(*velocidadp,
						  2.) + pow(VelocidadParedSalidap / __cons::ARef, 2.));
		}
		*caracteristica += var_energia;

		/*variaci�n debida al t�rmino de flujo a trav�s de la pared*/ // Ojo, que la velocidad en la pared est� ya con dimensiones (m/s)
		double var_flujo = 0.;
		if(FTipoCanal == nmCanalEntrada) {
			var_flujo = -2 * gamma1p * deltat * VelocidadParedEntradap * *asonidop / (ladocanalp - 2 * espesorsootp);
		} else {
			var_flujo = 2 * gamma1p * deltat * VelocidadParedSalidap * *asonidop / ladocanalp;
		}

		*caracteristica += var_flujo;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::Calculo_Caracteristica en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::InicializaCaracteristicas(TCondicionContorno **CC) {
	try {
		if(FNodoDer == 0) {
			CC[FNodoIzq - 1]->PutLanda(FTuboCCNodoIzq, FAsonido0[0] + FGamma3[0] * FVelocidad0[0]);
			CC[FNodoIzq - 1]->PutBeta(FTuboCCNodoIzq, FAsonido0[0] - FGamma3[0] * FVelocidad0[0]);
			CC[FNodoIzq - 1]->PutEntropia(FTuboCCNodoIzq, FAsonido0[0] / pow(FPresion0[0], FGamma5[0]));

			FLandaExtremoCerrado = FAsonido0[FNin - 1] + FGamma3[FNin - 1] * FVelocidad0[FNin - 1];
			FBetaExtremoCerrado = FAsonido0[FNin - 1] - FGamma3[FNin - 1] * FVelocidad0[FNin - 1];
			FEntropiaExtremoCerrado = FAsonido0[FNin - 1] / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]);
		}
		if(FNodoIzq == 0) {

			FLandaExtremoCerrado = FAsonido0[0] + FGamma3[0] * FVelocidad0[0];
			FBetaExtremoCerrado = FAsonido0[0] - FGamma3[0] * FVelocidad0[0];
			FEntropiaExtremoCerrado = FAsonido0[0] / pow(FPresion0[0], FGamma5[0]);

			CC[FNodoDer - 1]->PutLanda(FTuboCCNodoDer, FAsonido0[FNin - 1] + FGamma3[FNin - 1] * FVelocidad0[FNin - 1]);
			CC[FNodoDer - 1]->PutBeta(FTuboCCNodoDer, FAsonido0[FNin - 1] - FGamma3[FNin - 1] * FVelocidad0[FNin - 1]);
			CC[FNodoDer - 1]->PutEntropia(FTuboCCNodoDer, FAsonido0[FNin - 1] / pow(FPresion0[FNin - 1], FGamma5[FNin - 1]));

		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::InicializaCaracteristicas en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaB() {
	try {
		double p = 0., f = 0., tgas = 0., g = 0., q = 0., diamemed = 0.;
		double Rm = 0., Rm1 = 0., gamma = 0., gamma1 = 0., Vmed = 0., H1 = 0., H2 = 0., Amed = 0., Hmed = 0., rhoAmed = 0.;
		double qrA1, qrA2, qrAmed, Remed, Rmed, himed, rhomed, cturbmed, twallmed;
		double LadoCanalmed, EspesorSootmed, VelocidadParedmed, H0Paredmed, ViscosidadDinamicamed, friccionmed;

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin - 1; i++) {
				FTVD.Bvector[0][i] = 0.;
				FTVD.Bvector[1][i] = 0.;
				FTVD.Bvector[2][i] = 0.;

				Remed = 0.5 * FRe[i] + 0.5 * FRe[i + 1];
				Rm = sqrtRhoA[i + 1] / sqrtRhoA[i + 1];
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i + 1] + FGamma[i]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = (Frho[i] + Rm * Frho[i + 1]) / Rm1;
				LadoCanalmed = (FLadoCanal[i] + Rm * FLadoCanal[i + 1]) / Rm1;
				EspesorSootmed = (FEspesorSoot[i] + Rm * FEspesorSoot[i + 1]) / Rm1;
				VelocidadParedmed = (FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i) + Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 0,
									 i + 1)) / Rm1;
				if(VelocidadParedmed != 0) {
					printf("ya\n");
				}
				ViscosidadDinamicamed = (FViscosidadDinamica[i] + Rm * FViscosidadDinamica[i + 1]) / Rm1;

				FTVD.Bvector[0][i] = 4. * (LadoCanalmed - 2 * EspesorSootmed) * rhomed * VelocidadParedmed * FXref;
				H0Paredmed = (FH0Pared[i] + Rm * FH0Pared[i + 1]) / Rm1;

				if(FArea[i] != FArea[i + 1] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i + 1] + FVelocidad0[i]) / Rm1;
					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef2 * FAsonido0[i] * FAsonido0[i] / gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i + 1] * FVelocidad0[i + 1] + __cons::ARef2 * FAsonido0[i + 1] * FAsonido0[i +
							1] / gamma1;
					Hmed = (Rm * H2 + H1) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					rhoAmed = sqrtRhoA[i + 1] * sqrtRhoA[i + 1];
					p = rhomed * Amed * Amed / gamma;
				}

				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(FArea[i] != FArea[i + 1]) {
						FTVD.Bvector[1][i] = p * (FArea[i] - FArea[i + 1]);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bvector[1][i] = p * (FArea[i] - FArea[i + 1]) + FXref * friccionmed;
				}
				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrA1=(FArea[i] * Frho[i])*(-4*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i]);
					//qrA2=(FArea[i+1] * Frho[i+1])*(-4*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i+1]);
					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = -4 * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed;

					FTVD.Bvector[2][i] = -FXref * qrAmed;
				} else {
					diamemed = 0.5 * FDiametroTubo[i] + 0.5 * FDiametroTubo[i + 1];
					Rmed = 0.5 * FRMezcla[i] + 0.5 * FRMezcla[i + 1];
					himed = 0.5 * Fhi[i] + 0.5 * Fhi[i + 1];
					cturbmed = 0.5 * FCoefTurbulencia[i] + 0.5 * FCoefTurbulencia[i + 1];
					twallmed = 0.5 * FTPared[i] + 0.5 * FTPared[i + 1];

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();

					//qrA1=(FArea[i] * Frho[i])*(q-4*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i]);
					//qrA2=(FArea[i+1] * Frho[i+1])*(q-4*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i+1]);

					//qrAmed = 0.5*(qrA1 + qrA2);
					//qrAmed=((0.5*FArea[i]+ 0.5*FArea[i+1]) * rhomed)*(q*rhoAmed-4*H0Paredmed*VelocidadParedmed*LadoCanalmed/(0.5*FArea[i]+ 0.5*FArea[i+1]));
					qrAmed = q * rhoAmed - 4 * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed;
					FTVD.Bvector[2][i] = -FXref * qrAmed;
				}
			}
		} else if(FTipoCanal == nmCanalSalida) {
			for(int i = 0; i < FNin - 1; i++) {
				FTVD.Bvector[0][i] = 0.;
				FTVD.Bvector[1][i] = 0.;
				FTVD.Bvector[2][i] = 0.;

				Remed = 0.5 * FRe[i] + 0.5 * FRe[i + 1];
				Rm = sqrtRhoA[i + 1] / sqrtRhoA[i + 1];
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i + 1] + FGamma[i]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = 0.5 * Frho[i] + 0.5 * Frho[i + 1];
				LadoCanalmed = 0.5 * FLadoCanal[i] + 0.5 * FLadoCanal[i + 1];
				VelocidadParedmed = (FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i) + Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 1,
									 i + 1)) / Rm1;
				ViscosidadDinamicamed = 0.5 * FViscosidadDinamica[i] + 0.5 * FViscosidadDinamica[i + 1];

				FTVD.Bvector[0][i] = -4. * LadoCanalmed * rhomed * VelocidadParedmed * FXref;

				H0Paredmed = (FH0Pared[i] + Rm * FH0Pared[i + 1]) / Rm1;

				if(FArea[i] != FArea[i + 1] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i + 1] + FVelocidad0[i]) / Rm1;

					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef * FAsonido0[i] * __cons::ARef * FAsonido0[i] /
						 gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i + 1] * FVelocidad0[i + 1] + __cons::ARef2 * FAsonido0[i + 1] * FAsonido0[i +
							1] / gamma1;
					Hmed = (Rm * H2 + H1) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					rhoAmed = sqrtRhoA[i + 1] * sqrtRhoA[i + 1];
					p = rhomed * Amed * Amed / gamma;
				}

				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(FArea[i] != FArea[i + 1]) {
						FTVD.Bvector[1][i] = p * (FArea[i] - FArea[i + 1]);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bvector[1][i] = p * (FArea[i] - FArea[i + 1]) + FXref * friccionmed;
				}
				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrA1=(FArea[i] * Frho[i])*(4*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i]);
					//qrA2=(FArea[i+1] * Frho[i+1])*(4*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i+1]);
					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = 4 * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed;

					FTVD.Bvector[2][i] = -FXref * qrAmed;
				} else {
					diamemed = 0.5 * FDiametroTubo[i] + 0.5 * FDiametroTubo[i + 1];
					Rmed = 0.5 * FRMezcla[i] + 0.5 * FRMezcla[i + 1];
					himed = 0.5 * Fhi[i] + 0.5 * Fhi[i + 1];
					cturbmed = 0.5 * FCoefTurbulencia[i] + 0.5 * FCoefTurbulencia[i + 1];
					twallmed = 0.5 * FTPared[i] + 0.5 * FTPared[i + 1];

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();
					//qrA1=(FArea[i] * Frho[i])*(q+4*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i]);
					//qrA2=(FArea[i+1] * Frho[i+1])*(q+4*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i+1]);
					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = q * rhoAmed + 4 * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed;

					FTVD.Bvector[2][i] = -FXref * qrAmed;
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaB en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaBmen() {
	try {
		double p, f, tgas, g, q, diamemed;
		double rhoAmed, Rm, Rm1, H1, H2, gamma, gamma1, Hmed, Vmed, Amed, B;
		double qrA1, qrA2, qrAmed, Remed, Rmed, himed, rhomed, cturbmed, twallmed;
		double LadoCanalmed, EspesorSootmed, VelocidadParedmed, H0Paredmed, ViscosidadDinamicamed, friccionmed, Area12;

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 1; i < FNin; i++) {
				FTVD.Bmen[0][i] = 0.;
				FTVD.Bmen[1][i] = 0.;
				FTVD.Bmen[2][i] = 0.;
				rhoAmed = sqrt(sqrtRhoA[i - 1] * pow3(sqrtRhoA[i]));
				B = FU0[0][i] + rhoAmed + sqrtRhoA[i] * sqrtRhoA[i - 1];
				Rm = B / sqrt(pow3(sqrtRhoA[i - 1]) * sqrtRhoA[i]);
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i] + FGamma[i - 1]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = sqrt(Frho[i] * sqrt(Frho[i] * Frho[i - 1]));
				LadoCanalmed = (Rm * FLadoCanal[i] + FLadoCanal[i - 1]) / Rm1;
				EspesorSootmed = (Rm * FEspesorSoot[i] + FEspesorSoot[i - 1]) / Rm1;
				VelocidadParedmed = (Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i) + FDPF->GetVelocidadPared(FNumeroHaz - 1, 0,
									 i - 1)) / Rm1;
				ViscosidadDinamicamed = (Rm * FViscosidadDinamica[i] + FViscosidadDinamica[i - 1]) / Rm1;

				FTVD.Bmen[0][i] = 4. * (LadoCanalmed - 2. * EspesorSootmed) * rhomed * VelocidadParedmed * FXref;

				Area12 = (FArea[i - 1] + FArea[i]) / 2.;
				H0Paredmed = (Rm * FH0Pared[i] + FH0Pared[i - 1]) / Rm1;

				if(Area12 != FArea[i] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i] + FVelocidad0[i - 1]) / Rm1;

					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef2 * FAsonido0[i] * FAsonido0[i] / gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i - 1] * FVelocidad0[i - 1] + __cons::ARef2 * FAsonido0[i - 1] * FAsonido0[i -
							1] / gamma1;
					Hmed = (Rm * H1 + H2) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					rhoAmed = sqrt(sqrtRhoA[i - 1] * pow3(sqrtRhoA[i]));
					p = rhomed * Amed * Amed / gamma;
				}

				Remed = (Rm * FRe[i] + FRe[i - 1]) / Rm1;
				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(Area12 != FArea[i]) {
						FTVD.Bmen[1][i] = p * (Area12 - FArea[i]);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bmen[1][i] = p * (Area12 - FArea[i]) + FXref * friccionmed;
				}
				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrA1=(FArea[i] * Frho[i])*(-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i]);
					//qrA2=(Area12 * rhomed)*(-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/Area12);

					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = -4. * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed;
					FTVD.Bmen[2][i] = -FXref * qrAmed;
				} else {
					diamemed = (Rm * FDiametroTubo[i] + FDiametroTubo[i - 1]) / Rm1;
					Rmed = (Rm * FRMezcla[i] + FRMezcla[i - 1]) / Rm1;
					himed = (Rm * Fhi[i] + Fhi[i - 1]) / Rm1;
					cturbmed = (Rm * FCoefTurbulencia[i] + FCoefTurbulencia[i - 1]) / Rm1;
					twallmed = (Rm * FTPared[i] + FTPared[i - 1]) / Rm1;

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();
					//qrA1=(FArea[i] * Frho[i])*(q-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/FArea[i]);
					//qrA2=(Area12 * rhomed)*(q-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/Area12);
					//qrAmed = 0.5*(qrA1 + qrA2);
					//qrAmed =(Area12 * rhomed)*(q-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/Area12);
					qrAmed = q * rhoAmed - 4. * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed;
					FTVD.Bmen[2][i] = -FXref * qrAmed;
				}

			}
		} else if(FTipoCanal == nmCanalSalida) {
			for(int i = 1; i < FNin; i++) {
				FTVD.Bmen[0][i] = 0.;
				FTVD.Bmen[1][i] = 0.;
				FTVD.Bmen[2][i] = 0.;

				rhoAmed = sqrt(sqrtRhoA[i - 1] * pow3(sqrtRhoA[i]));
				B = FU0[0][i] + rhoAmed + sqrtRhoA[i] * sqrtRhoA[i - 1];
				Rm = B / sqrt(pow3(sqrtRhoA[i - 1]) * sqrtRhoA[i]);
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i] + FGamma[i - 1]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = sqrt(Frho[i] * sqrt(Frho[i] * Frho[i - 1]));
				LadoCanalmed = (Rm * FLadoCanal[i] + FLadoCanal[i - 1]) / Rm1;
				VelocidadParedmed = (Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i) + FDPF->GetVelocidadPared(FNumeroHaz - 1, 1,
									 i - 1)) / Rm1;
				ViscosidadDinamicamed = (Rm * FViscosidadDinamica[i] + FViscosidadDinamica[i - 1]) / Rm1;

				FTVD.Bmen[0][i] = -4. * LadoCanalmed * rhomed * VelocidadParedmed * FXref;

				Area12 = (FArea[i - 1] + FArea[i]) / 2.;
				H0Paredmed = (Rm * FH0Pared[i] + FH0Pared[i - 1]) / Rm1;

				if(Area12 != FArea[i] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i] + FVelocidad0[i - 1]) / Rm1;

					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef * __cons::ARef * FAsonido0[i] * FAsonido0[i] /
						 gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i - 1] * FVelocidad0[i - 1] + __cons::ARef2 * FAsonido0[i - 1] * FAsonido0[i -
							1] / gamma1;
					Hmed = (Rm * H1 + H2) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					rhoAmed = sqrt(sqrtRhoA[i - 1] * pow3(sqrtRhoA[i]));
					p = rhomed * Amed * Amed / gamma;
				}

				Remed = (Rm * FRe[i] + FRe[i - 1]) / Rm1;

				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(Area12 != FArea[i]) {
						FTVD.Bmen[1][i] = p * (Area12 - FArea[i]);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bmen[1][i] = p * (Area12 - FArea[i]) + FXref * friccionmed;
				}

				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrA1=(FArea[i] * Frho[i])*(4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i]);
					//qrA2=(Area12 * rhomed)*(4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/Area12);
					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = 4. * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed;
					FTVD.Bmen[2][i] = -FXref * qrAmed;

				} else {
					diamemed = (Rm * FDiametroTubo[i] + FDiametroTubo[i - 1]) / Rm1;
					Rmed = (Rm * FRMezcla[i] + FRMezcla[i - 1]) / Rm1;
					himed = (Rm * Fhi[i] + Fhi[i - 1]) / Rm1;
					cturbmed = (Rm * FCoefTurbulencia[i] + FCoefTurbulencia[i - +1]) / Rm1;
					twallmed = (Rm * FTPared[i] + FTPared[i - 1]) / Rm1;

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();
					//qrA1=(FArea[i] * Frho[i])*(q+4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/FArea[i]);
					//qrA2=(Area12 * rhomed)*(q+4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/Area12);
					//qrAmed = 0.5*(qrA1 + qrA2);
					qrAmed = q * rhoAmed + 4. * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed;

					FTVD.Bmen[2][i] = -FXref * qrAmed;
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaBmas en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaBmas() {
	try {
		double v = 0., a = 0., p = 0., f = 0., tgas = 0., g = 0., q = 0., diamemed = 0.;
		double rhoAmed = 0., Amed = 0., Vmed = 0., Hmed = 0., H1 = 0., H2 = 0., Rm = 0., Rm1 = 0., gamma = 0., gamma1 = 0.,
			   B = 0.;
		double qrA1, qrA2, qrAmed, Remed, Rmed, himed, rhomed, cturbmed, twallmed;
		double LadoCanalmed, EspesorSootmed, VelocidadParedmed, H0Paredmed, ViscosidadDinamicamed, friccionmed, Area12;

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin - 1; i++) {
				FTVD.Bmas[0][i] = 0.;
				FTVD.Bmas[1][i] = 0.;
				FTVD.Bmas[2][i] = 0.;

				rhoAmed = sqrt(sqrtRhoA[i + 1] * pow3(sqrtRhoA[i]));
				B = FU0[0][i] + rhoAmed + sqrtRhoA[i] * sqrtRhoA[i + 1];
				Rm = B / sqrt(pow3(sqrtRhoA[i + 1]) * sqrtRhoA[i]);
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i] + FGamma[i + 1]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = sqrt(Frho[i] * sqrt(Frho[i] * Frho[i + 1]));
				LadoCanalmed = (Rm * FLadoCanal[i] + FLadoCanal[i + 1]) / Rm1;
				EspesorSootmed = (Rm * FEspesorSoot[i] + FEspesorSoot[i + 1]) / Rm1;
				VelocidadParedmed = (Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i) + FDPF->GetVelocidadPared(FNumeroHaz - 1, 0,
									 i + 1)) / Rm1;
				ViscosidadDinamicamed = (Rm * FViscosidadDinamica[i] + FViscosidadDinamica[i + 1]) / Rm1;

				FTVD.Bmas[0][i] = 4. * (LadoCanalmed - 2. * EspesorSootmed) * rhomed * VelocidadParedmed * FXref;

				Area12 = (FArea[i] + FArea[i + 1]) / 2.;
				H0Paredmed = (Rm * FH0Pared[i] + FH0Pared[i + 1]) / Rm1;

				if(Area12 != FArea[i] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i] + FVelocidad0[i + 1]) / Rm1;

					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef2 * FAsonido0[i] * FAsonido0[i] / gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i + 1] * FVelocidad0[i + 1] + __cons::ARef * __cons::ARef * FAsonido0[i + 1] *
						 FAsonido0[i + 1] / gamma1;
					Hmed = (Rm * H1 + H2) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					rhoAmed = sqrt(sqrtRhoA[i + 1] * pow3(sqrtRhoA[i]));
					p = rhomed * Amed * Amed / gamma;
				}

				Remed = (Rm * FRe[i] + FRe[i + 1]) / Rm1;
				p = rhomed * Amed * Amed / gamma;
				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(Area12 != FArea[i]) {
						FTVD.Bmas[1][i] = p * (FArea[i] - Area12);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bmas[1][i] = p * (FArea[i] - Area12) + FXref * friccionmed;
				}
				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrAmed =(Area12 * rhomed)*(-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/Area12);
					qrAmed = -4. * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed / FXref;
					FTVD.Bmen[2][i] = -FXref * qrAmed;

				} else {
					diamemed = (Rm * FDiametroTubo[i] + FDiametroTubo[i + 1]) / Rm1;
					Rmed = (Rm * FRMezcla[i] + FRMezcla[i + 1]) / Rm1;
					himed = (Rm * Fhi[i] + Fhi[i + 1]) / Rm1;
					cturbmed = (Rm * FCoefTurbulencia[i] + FCoefTurbulencia[i + 1]) / Rm1;
					twallmed = (Rm * FTPared[i] + FTPared[i + 1]) / Rm1;

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();
					//qrAmed =(Area12 * rhomed)*(q-4.*H0Paredmed*VelocidadParedmed*(LadoCanalmed-2*EspesorSootmed)/Area12);
					qrAmed = q * rhoAmed - 4. * H0Paredmed * VelocidadParedmed * (LadoCanalmed - 2 * EspesorSootmed) * rhomed / FXref;
					FTVD.Bmen[2][i] = -FXref * qrAmed;
				}

			}
		}
		if(FTipoCanal == nmCanalSalida) {
			for(int i = 0; i < FNin - 1; i++) {
				FTVD.Bmas[0][i] = 0.;
				FTVD.Bmas[1][i] = 0.;
				FTVD.Bmas[2][i] = 0.;

				rhoAmed = sqrt(sqrtRhoA[i + 1] * pow3(sqrtRhoA[i]));
				B = FU0[0][i] + rhoAmed + sqrtRhoA[i] * sqrtRhoA[i + 1];
				Rm = B / sqrt(pow3(sqrtRhoA[i + 1]) * sqrtRhoA[i]);
				Rm1 = Rm + 1;
				gamma = (Rm * FGamma[i] + FGamma[i + 1]) / Rm1;
				gamma1 = gamma - 1;

				rhomed = sqrt(Frho[i] * sqrt(Frho[i] * Frho[i + 1]));
				LadoCanalmed = (Rm * FLadoCanal[i] + FLadoCanal[i + 1]) / Rm1;
				VelocidadParedmed = (Rm * FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i) + FDPF->GetVelocidadPared(FNumeroHaz - 1, 1,
									 i + 1)) / Rm1;
				ViscosidadDinamicamed = (Rm * FViscosidadDinamica[i] + FViscosidadDinamica[i + 1]) / Rm1;

				FTVD.Bmas[0][i] = -4. * LadoCanalmed * rhomed * VelocidadParedmed * FXref;

				Area12 = (FArea[i + 1] + FArea[i]) / 2.;
				H0Paredmed = (Rm * FH0Pared[i] + FH0Pared[i + 1]) / Rm1;

				if(Area12 != FArea[i] || FCoefAjusFric != 0 || FDPF->getCoefAjusTC() != 0) {
					Vmed = __cons::ARef * (Rm * FVelocidad0[i] + FVelocidad0[i + 1]) / Rm1;

					H1 = 0.5 * __cons::ARef2 * FVelocidad0[i] * FVelocidad0[i] + __cons::ARef2 * FAsonido0[i] * FAsonido0[i] / gamma1;
					H2 = 0.5 * __cons::ARef2 * FVelocidad0[i + 1] * FVelocidad0[i + 1] + __cons::ARef2 * FAsonido0[i + 1] * FAsonido0[i +
							1] / gamma1;
					Hmed = (Rm * H1 + H2) / Rm1;
					Amed = sqrt(gamma1 * (Hmed - 0.5 * Vmed * Vmed));
					p = rhomed * Amed * Amed / gamma;
				}
				Remed = (Rm * FRe[i] + FRe[i + 1]) / Rm1;

				if(Remed <= 1e-6 || FCoefAjusFric == 0) {
					if(Area12 != FArea[i]) {
						FTVD.Bmas[1][i] = p * (FArea[i] - Area12);
					}
				} else {
					friccionmed = FCoefAjusFric * FF * ViscosidadDinamicamed * Vmed;
					FTVD.Bmas[1][i] = p * (FArea[i] - Area12) + FXref * friccionmed;
				}
				if(Remed <= 1e-6 || FDPF->getCoefAjusTC() == 0) {
					//qrAmed =(Area12 * rhomed)*(4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/Area12);
					qrAmed = 4. * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed / FXref;
					FTVD.Bmas[2][i] = -FXref * qrAmed;
				} else {
					diamemed = (Rm * FDiametroTubo[i] + FDiametroTubo[i + 1]) / Rm1;
					Rmed = (Rm * FRMezcla[i] + FRMezcla[i + 1]) / Rm1;
					himed = (Rm * Fhi[i] + Fhi[i + 1]) / Rm1;
					cturbmed = (Rm * FCoefTurbulencia[i] + FCoefTurbulencia[i + 1]) / Rm1;
					twallmed = (Rm * FTPared[i] + FTPared[i + 1]) / Rm1;

					tgas = Amed * Amed / gamma / Rmed;

					TransmisionCalor(tgas, diamemed, &q, cturbmed, himed, rhomed, twallmed);
					q = q * FDPF->getCoefAjusTC();
					//qrAmed =(Area12 * rhomed)*(q+4.*H0Paredmed*VelocidadParedmed*LadoCanalmed/Area12);
					qrAmed = q * rhoAmed + 4. * H0Paredmed * VelocidadParedmed * LadoCanalmed * rhomed / FXref;
					FTVD.Bmas[2][i] = -FXref * qrAmed;
				}
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaBmas en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::CalculaMatrizJacobiana() {
	try {
		double Rmed = 0., Rmed1 = 0., Vmed = 0., Hmed = 0., Amed = 0., gamma = 0., gamma1 = 0., gamma2 = 0., H1 = 0., H2 = 0.;
		double *Ymed;

		Ymed = new double[FNumeroEspecies - 1 - FIntEGR];

		for(int i = 0; i < FNin - 1; i++) {
			RoeConstants();

			// Calculo de las medias de Roe
			Rmed = sqrtRhoA[i + 1] / sqrtRhoA[i + 1];
			Rmed1 = Rmed + 1;
			gamma = (Rmed * FGamma[i + 1] + FGamma[i]) / Rmed1;
			gamma1 = gamma - 1;
			gamma2 = gamma1 / 2;
			H1 = 0.5 * FVelocidad0[i] * FVelocidad0[i] * pow(__cons::ARef, 2.) + pow(FAsonido0[i] * __cons::ARef, 2.) / gamma1;
			H2 = 0.5 * FVelocidad0[i + 1] * FVelocidad0[i + 1] * pow(__cons::ARef, 2.) + pow(FAsonido0[i + 1] * __cons::ARef,
					2.) / gamma1;
			Vmed = (Rmed * FVelocidad0[i + 1] + FVelocidad0[i]) * __cons::ARef / Rmed1;
			Hmed = (Rmed * H2 + H1) / Rmed1;
			Amed = pow(gamma1 * (Hmed - 0.5 * Vmed * Vmed), 0.5);
			for(int j = 0; j < FNumeroEspecies - 1 - FIntEGR; j++) {
				Ymed[j] = (Rmed * FFraccionMasicaEspecie[i + 1][j] + FFraccionMasicaEspecie[i][j]) / Rmed1;
			}

			FTVD.Pmatrix[1][0][i] = Vmed - Amed;
			FTVD.Pmatrix[1][1][i] = Vmed;
			FTVD.Pmatrix[1][2][i] = Vmed + Amed;

			FTVD.Pmatrix[2][0][i] = Hmed - Vmed * Amed;
			FTVD.Pmatrix[2][1][i] = Vmed * Vmed / 2.;
			FTVD.Pmatrix[2][2][i] = Hmed + Vmed * Amed;

			// Calculo de la matriz Q (autovectores a la derecha)
			FTVD.Qmatrix[0][0][i] = Vmed / Amed / 2. + gamma2 * Vmed * Vmed / Amed / Amed / 2.;
			FTVD.Qmatrix[0][1][i] = -0.5 / Amed - gamma2 * Vmed / Amed / Amed;
			FTVD.Qmatrix[0][2][i] = gamma2 / Amed / Amed;

			FTVD.Qmatrix[1][0][i] = 1. - gamma2 * Vmed * Vmed / Amed / Amed;
			FTVD.Qmatrix[1][1][i] = gamma1 * Vmed / Amed / Amed;
			FTVD.Qmatrix[1][2][i] = -gamma1 / Amed / Amed;

			FTVD.Qmatrix[2][0][i] = -Vmed / Amed / 2. + gamma2 / 2. * Vmed * Vmed / Amed / Amed;
			FTVD.Qmatrix[2][1][i] = 0.5 / Amed - gamma2 * Vmed / Amed / Amed;
			FTVD.Qmatrix[2][2][i] = gamma2 / Amed / Amed;

			// Valores propios de la matriz jacobiana
			FTVD.Alpha[0][i] = FTVD.Pmatrix[1][0][i];
			FTVD.Alpha[1][i] = Vmed;
			FTVD.Alpha[2][i] = FTVD.Pmatrix[1][2][i];

			for(int j = 3; j < FNumEcuaciones; j++) {
				FTVD.Alpha[j][i] = Vmed;

			}
		}
		delete[] Ymed;
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::CalculaMatrizJacobiana en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::TVD_Estabilidad() {
	try {
		double VTotalMax = 0.;
		double VLocal = 0., DeltaT_tvd = 0.;
		double B1 = 0., B2 = 0., FraccionMasicaSalida1 = 0., FraccionMasicaSalida2 = 0., temp = 0.;

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin; i++) {
				FEspesorSoot[i] = FDPF->GetEspesorSoot(FNumeroHaz - 1, i);
				FArea[i] = pow(FLadoCanal[i] - 2 * FEspesorSoot[i], 2.);
				FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 2);
				temp = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
				FH0Pared[i] = FRMezcla[i] * FGamma[i] / FGamma1[i] * temp + pow(FVelocidad0[i] * __cons::ARef, 2.) / 2;
			}
		} else {
			for(int i = 0; i < FNin; i++) {
				if(FNodoInicialFuenteCE == 0) {
					FH0Pared[i] = FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(i);
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 0);
				} else if(FNodoInicialFuenteCE + i <= FNin - 1) {
					if(i == 0) {
						FH0Pared[i] = FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i)
									  - (FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i + 1) - FDPF->GetCanal(FNumeroHaz - 1,
											  0)->GetH0Pared(FNodoInicialFuenteCE + i)) / FXref
									  * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0)
									 - (FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1, 0) - FDPF->GetTPared(FNumeroHaz - 1,
											 FNodoInicialFuenteCE + i, 0)) / FXref
									 * FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion();
						//FH0Pared[i]=FDPF->GetCanal(FNumeroHaz-1,0)->GetH0Pared(FNodoInicialFuenteCE+i);
						//FTPared[i]=FDPF->GetTPared(FNumeroHaz-1,FNodoInicialFuenteCE+i,0);
					} else {
						FH0Pared[i] = InterpolaTubo(FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE - 1 + i),
													FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i), 1., (FXref - FDPF->GetCanal(FNumeroHaz - 1,
															0)->getDistanciaInterpolacion()) / FXref);
						FTPared[i] = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i, 0),
												   FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0), 1.,
												   (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);

					}
				} else {
					FH0Pared[i] = 0.;
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNin - 1, 0);
				}
			}
		}

		CalculaFlujo(FU0, FTVD.W, FGamma, FGamma1, FNin);

		CalculaMatrizJacobiana();

		CalculaB();

		for(int i = 0; i < FNin - 1; i++) {
			for(int k = 0; k < 3; ++k) {
				FTVD.DeltaU[k][i] = FTVD.Qmatrix[k][0][i] * (FU0[0][i + 1] - FU0[0][i]) + FTVD.Qmatrix[k][1][i] *
									(FU0[1][i + 1] - FU0[1][i]) + FTVD.Qmatrix[k][2][i] * (FU0[2][i + 1] - FU0[2][i]);
				/*FTVD.DeltaB[k][i]=FTVD.Qmatrix[k][0][i]*FTVD.Bvector[0][i]+
				 FTVD.Qmatrix[k][1][i]*FTVD.Bvector[1][i]+
				 FTVD.Qmatrix[k][2][i]*FTVD.Bvector[2][i];*/
				FTVD.DeltaB[k][i] = FTVD.Qmatrix[k][0][i] * (FTVD.Bvector[0][i + 1] - FTVD.Bvector[0][i]) + FTVD.Qmatrix[k][1][i] *
									(FTVD.Bvector[1][i + 1] - FTVD.Bvector[1][i])
									+ FTVD.Qmatrix[k][2][i] * (FTVD.Bvector[2][i + 1] - FTVD.Bvector[2][i]);
				FTVD.DeltaW[k][i] = FTVD.Qmatrix[k][0][i] * (FTVD.W[0][i + 1] - FTVD.W[0][i] + FTVD.Bvector[0][i]) +
									FTVD.Qmatrix[k][1][i] * (FTVD.W[1][i + 1] - FTVD.W[1][i] + FTVD.Bvector[1][i])
									+ FTVD.Qmatrix[k][2][i] * (FTVD.W[2][i + 1] - FTVD.W[2][i] + FTVD.Bvector[2][i]);
				if(fabs(FTVD.DeltaU[k][i]) < 1e-3) {
					FTVD.Beta[k][i] = 0.;
				} else {
					FTVD.Beta[k][i] = FTVD.DeltaB[k][i] / FTVD.DeltaU[k][i];
				}
				if(FTVD.Alpha[k][i] + FTVD.Beta[k][i] != 0) {
					if((VLocal = fabs(FTVD.Alpha[k][i] + FTVD.Beta[k][i])) > VTotalMax) {
						VTotalMax = VLocal;
					}
				}
			}
			for(int k = 3; k < FNumEcuaciones; k++) {
				if(FTipoCanal == nmCanalEntrada) {
					B1 = 4. * (FLadoCanal[i] - 2 * FDPF->GetEspesorSoot(FNumeroHaz - 1,
							   i)) * Frho[i] * FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i) * FFraccionMasicaEspecie[i][k - 3] * FXref;
					B2 = 4. * (FLadoCanal[i + 1] - 2 * FDPF->GetEspesorSoot(FNumeroHaz - 1,
							   i + 1)) * Frho[i + 1] * FDPF->GetVelocidadPared(FNumeroHaz - 1, 0, i + 1)
						 * FFraccionMasicaEspecie[i + 1][k - 3] * FXref;
				} else if(FTipoCanal == nmCanalSalida) {
					if(FNodoInicialFuenteCE == 0) {
						FraccionMasicaSalida1 = FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, i, k - 3);
						FraccionMasicaSalida2 = FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, i, k - 3);
					} else if(FNodoInicialFuenteCE + i <= FNin - 1) {
						FraccionMasicaSalida1 = InterpolaTubo(FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i,
															  k - 3),
															  FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i, k - 3), 1.,
															  (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						if(FNodoInicialFuenteCE + i == FNin - 1) {
							FraccionMasicaSalida2 = 0.;
						} else {
							FraccionMasicaSalida2 = InterpolaTubo(FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i, k - 3),
																  FDPF->GetFraccionMasicaSalida(FNumeroHaz - 1, FNodoInicialFuenteCE + i + 1, k - 3), 1.,
																  (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
						}
					} else {
						FraccionMasicaSalida1 = 0.;
						FraccionMasicaSalida2 = 0.;
					}
					B1 = 4. * FLadoCanal[i] * Frho[i] * FDPF->GetVelocidadPared(FNumeroHaz - 1, 1, i) * FraccionMasicaSalida1 * FXref;
					B2 = 4. * FLadoCanal[i + 1] * Frho[i + 1] * FDPF->GetVelocidadPared(FNumeroHaz - 1, 1,
							i + 1) * FraccionMasicaSalida2 * FXref;
				}
				FTVD.DeltaB[k][i] = B2 - B1;
				FTVD.DeltaU[k][i] = FU0[0][i + 1] - FU0[0][i];
				if(fabs(FTVD.DeltaU[k][i]) < 1e-3) {
					FTVD.Beta[k][i] = 0.;
				} else {
					FTVD.Beta[k][i] = FTVD.DeltaB[k][i] / FTVD.DeltaU[k][i];
				}
				if(FTVD.Alpha[k][i] + FTVD.Beta[k][i] != 0) {
					if((VLocal = fabs(FTVD.Alpha[k][i] + FTVD.Beta[k][i])) > VTotalMax) {
						VTotalMax = VLocal;
					}
				}
			}
		}
		DeltaT_tvd = FCourant * FXref / VTotalMax;
		if(DeltaT_tvd < FDeltaTime)
			FDeltaTime = DeltaT_tvd;

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::TVD_Estabilidad en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::TVD_Limitador() {
	try {
		double temp = 0.;
		double dtdx = FDeltaTime / FXref;

		if(FTipoCanal == nmCanalEntrada) {
			for(int i = 0; i < FNin; i++) {
				FEspesorSoot[i] = FDPF->GetEspesorSoot(FNumeroHaz - 1, i);
				FArea[i] = pow(FLadoCanal[i] - 2 * FEspesorSoot[i], 2.);
				FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 2);
				temp = pow(FAsonido0[i] * __cons::ARef, 2.) / (FGamma[i] * FRMezcla[i]);
				FH0Pared[i] = FRMezcla[i] * FGamma[i] / FGamma1[i] * temp + pow(FVelocidad0[i] * __cons::ARef, 2.) / 2;
			}
		} else {
			for(int i = 0; i < FNin; i++) {
				if(FNodoInicialFuenteCE == 0) {
					FH0Pared[i] = FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(i);
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, i, 0);
				} else if(FNodoInicialFuenteCE + i <= FNin - 1) {
					if(i == 0) {
						FH0Pared[i] = FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i);
						FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0);
					} else {
						FH0Pared[i] = InterpolaTubo(FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE - 1 + i),
													FDPF->GetCanal(FNumeroHaz - 1, 0)->GetH0Pared(FNodoInicialFuenteCE + i), 1., (FXref - FDPF->GetCanal(FNumeroHaz - 1,
															0)->getDistanciaInterpolacion()) / FXref);
						FTPared[i] = InterpolaTubo(FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE - 1 + i, 0),
												   FDPF->GetTPared(FNumeroHaz - 1, FNodoInicialFuenteCE + i, 0), 1.,
												   (FXref - FDPF->GetCanal(FNumeroHaz - 1, 0)->getDistanciaInterpolacion()) / FXref);
					}
				} else {
					FH0Pared[i] = 0.;
					FTPared[i] = FDPF->GetTPared(FNumeroHaz - 1, FNin - 1, 0);
					;
				}
			}
		}

		CalculaBmas();

		CalculaBmen();

		for(int i = 0; i < FNin - 1; ++i) {
			// Calculo de las variables en el sistema de coordenadas Q
			for(int k = 0; k < FNumEcuaciones; k++) {
				FTVD.LandaD[k][i] = dtdx * (FTVD.Alpha[k][i] + FTVD.Beta[k][i]);
				if(FTVD.LandaD[k][i] >= 0.) {
					FTVD.hLandaD[k][i] = 1;
				} else {
					FTVD.hLandaD[k][i] = -1;
				}
			}
		}
		for(int i = 1; i < FNin - 2; ++i) {
			for(int k = 0; k < 3; k++) {
				if(fabs(FTVD.DeltaW[k][i]) < 0.0001 || fabs(FTVD.LandaD[k][i]) == 1) {
					FTVD.R[k][i] = 1.;
				} else {
					FTVD.R[k][i] = ((double) FTVD.hLandaD[k][i - FTVD.hLandaD[k][i]] - FTVD.LandaD[k][i - FTVD.hLandaD[k][i]]) *
								   (FTVD.DeltaW[k][i - FTVD.hLandaD[k][i]])
								   / (((double) FTVD.hLandaD[k][i] - FTVD.LandaD[k][i]) * FTVD.DeltaW[k][i]);
				}
			}
			for(int k = 3; k < FNumEcuaciones; k++) {
				if(fabs(FTVD.W[k][i] - FTVD.W[k][i + 1]) < 0.0001 || fabs(FTVD.LandaD[k][i]) == 1) {
					FTVD.R[k][i] = 1.;
				} else {
					FTVD.R[k][i] = ((double) FTVD.hLandaD[k][i - FTVD.hLandaD[k][i]] - FTVD.LandaD[k][i - FTVD.hLandaD[k][i]])
								   * (FTVD.W[k][i + 1 - FTVD.hLandaD[k][i]] - FTVD.W[k][i - FTVD.hLandaD[k][i]]) / (((double) FTVD.hLandaD[k][i] -
										   FTVD.LandaD[k][i]) * (FTVD.W[k][i + 1] - FTVD.W[k][i]));
				}
			}
		}
		for(int k = 0; k < FNumEcuaciones; k++) {
			FTVD.R[k][0] = 0.;
			FTVD.R[k][FNin - 2] = 0.;

		}
		for(int i = 0; i < FNin - 1; ++i) {
			for(int k = 0; k < 3; k++) {
				FTVD.Phi[k][i] = (double) FTVD.hLandaD[k][i] - Limita(FTVD.R[k][i]) * ((double) FTVD.hLandaD[k][i] - FTVD.LandaD[k][i]);
			}
		}

		for(int i = 0; i < FNin - 1; ++i) {
			for(int k = 0; k < 3; ++k) {
				FTVD.gflux[k][i] = 0.5
								   * (FTVD.W[k][i] + FTVD.W[k][i + 1] - FTVD.Bmas[k][i] + FTVD.Bmen[k][i + 1]
									  - (FTVD.Pmatrix[k][0][i] * FTVD.Phi[0][i] * FTVD.DeltaW[0][i] + FTVD.Pmatrix[k][1][i] * FTVD.Phi[1][i] *
										 FTVD.DeltaW[1][i]
										 + FTVD.Pmatrix[k][2][i] * FTVD.Phi[2][i] * FTVD.DeltaW[2][i]));
			}
			for(int k = 3; k < FNumEcuaciones; k++) {
				FTVD.gflux[k][i] = 0.5 * (FTVD.W[k][i] + FTVD.W[k][i + 1] - (double) FTVD.hLandaD[k][i] *
										  (FTVD.W[k][i + 1] - FTVD.W[k][i]))
								   + 0.5 * Limita(FTVD.R[k][i]) * ((double) FTVD.hLandaD[k][i] - FTVD.LandaD[k][i]) * (FTVD.W[k][i + 1] - FTVD.W[k][i]);
			}
		}
		for(int i = 1; i < FNin - 1; ++i) {
			for(int k = 0; k < 3; k++) {
				FU1[k][i] = FU0[k][i] - dtdx * ((FTVD.gflux[k][i] - FTVD.gflux[k][i - 1]) + (FTVD.Bmen[k][i] + FTVD.Bmas[k][i]));
			}
			for(int k = 3; k < FNumEcuaciones; k++) {
				FU1[k][i] = FU0[k][i] - dtdx * (FTVD.gflux[k][i] - FTVD.gflux[k][i - 1]);
			}
		}
	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::TVD_Limitador en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TCanalDPF::Limita(double r) {
//------Van Albada
//double ret_val=(r*r+r)/(1+r*r);

//------Van Leer
//double ret_val=(fabs(r)+r)/(1.+fabs(r));

//------Minmod
	double ret_val = Minimo(2 * fabs(r), 1.);
//double ret_val1=Minimo(2*r,1.);
//double ret_val=Maximo(0,ret_val1);

//------Limitador de Roe (minmod)
//double ret_val1=Minimo(r,1);
//double ret_val=Maximo(0,ret_val1);

//------Limitador Roe Superbee
//double ret_val1=Minimo(2*r,1);
//double ret_val2=Minimo(r,2);
//double ret_val3=Maximo(0,ret_val1);
//double ret_val=Maximo(ret_val2,ret_val3);

//------Limitador de Osher (phi=1.5)
//double ret_val1=Minimo(r,2);
//double ret_val=Maximo(0,ret_val1);

//------Limitador de Sweby (phi=1.5)
//double ret_val1=Minimo(2*r,1);
//double ret_val2=Minimo(r,2);
//double ret_val3=Maximo(0,ret_val1);
//double ret_val=Maximo(ret_val2,ret_val3);

//------Limitador de Ospre (WAterson & Deconinck)
//double ret_val=1.5*(r*r+r)/(r*r+r+1.);

//------Limitador de MC de van Leer
//double ret_val1=Minimo(0.5+0.5*r,2);
//double ret_val2=Minimo(ret_val1,2*r);
//double ret_val=Maximo(0,ret_val2);

//------Limitador de UMIST
//double ret_val1=Minimo(0.75+0.25*r,2);
//double ret_val2=Minimo(ret_val1,0.25+0.75*r);
//double ret_val3=Minimo(ret_val2,2*r);
//double ret_val=Maximo(0,ret_val3);

//------Limitador de Koren
//double ret_val1=Minimo((1+2*r)/3,2);
//double ret_val2=Minimo(ret_val1,2*r);
//double ret_val=Maximo(0,ret_val2);

	return ret_val;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::DimensionaTVD() {
	try {

		FTVD.Bmas = new double*[FNumEcuaciones];
		FTVD.Bvector = new double*[FNumEcuaciones];
		FTVD.Bmen = new double*[FNumEcuaciones];
		FTVD.gflux = new double*[FNumEcuaciones];
		FTVD.Alpha = new double*[FNumEcuaciones];
		FTVD.Beta = new double*[FNumEcuaciones];
		FTVD.DeltaU = new double*[FNumEcuaciones];
		FTVD.DeltaB = new double*[FNumEcuaciones];
		FTVD.DeltaW = new double*[FNumEcuaciones];
		FTVD.hLandaD = new int*[FNumEcuaciones];
		FTVD.LandaD = new double*[FNumEcuaciones];
		FTVD.Phi = new double*[FNumEcuaciones];
		FTVD.R = new double*[FNumEcuaciones];
		FTVD.W = new double*[FNumEcuaciones];
		FTVD.Qmatrix = new double**[FNumEcuaciones];
		FTVD.Pmatrix = new double**[FNumEcuaciones];

		sqrtRhoA = new double[FNin];

		for(int k = 0; k < FNumEcuaciones; ++k) {
			FTVD.Bmas[k] = new double[FNin];
			FTVD.Bvector[k] = new double[FNin];
			FTVD.Bmen[k] = new double[FNin];
			FTVD.gflux[k] = new double[FNin];
			FTVD.Alpha[k] = new double[FNin];
			FTVD.Beta[k] = new double[FNin];
			FTVD.DeltaU[k] = new double[FNin];
			FTVD.DeltaB[k] = new double[FNin];
			FTVD.DeltaW[k] = new double[FNin];
			FTVD.hLandaD[k] = new int[FNin];
			FTVD.LandaD[k] = new double[FNin];
			FTVD.Phi[k] = new double[FNin];
			FTVD.R[k] = new double[FNin];
			FTVD.W[k] = new double[FNin];
		}

		for(int k = 0; k < FNumEcuaciones; ++k) {
			FTVD.Qmatrix[k] = new double*[FNumEcuaciones];
			FTVD.Pmatrix[k] = new double*[FNumEcuaciones];

			for(int i = 0; i < FNumEcuaciones; i++) {
				FTVD.Qmatrix[k][i] = new double[FNin];
				FTVD.Pmatrix[k][i] = new double[FNin];
			}
		}
		for(int i = 0; i < FNin; i++) {

			for(int k = 3; k < FNumEcuaciones; k++) {
				for(int j = 3; j < FNumEcuaciones; j++) {
					FTVD.Pmatrix[k][j][i] = 0.;
				}
			}
			FTVD.Pmatrix[0][0][i] = 1.;
			FTVD.Pmatrix[0][1][i] = 1.;
			FTVD.Pmatrix[0][2][i] = 1.;
			for(int j = 3; j < FNumEcuaciones; j++) {
				FTVD.Pmatrix[j][j][i] = 1.;
			}

			for(int k = 3; k < FNumEcuaciones; k++) {
				for(int j = 3; j < FNumEcuaciones; j++) {
					FTVD.Qmatrix[k][j][i] = 0.;
				}
			}
			for(int j = 3; j < FNumEcuaciones; j++) {
				FTVD.Qmatrix[j][j][i] = 1.;
			}
		}

	} catch(exception &N) {
		std::cout << "ERROR: TCanalDPF::DimensionaTVD en el haz " << FNumeroHaz << "de la DPF" << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TCanalDPF::RoeConstants() {

	sqrtRhoA[0] = sqrt(FU1[0][0]);
	for(int i = 0; i < FNin - 1; i++) {
		sqrtRhoA[i + 1] = sqrt(FU1[0][i + 1]);
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#pragma package(smart_init)

