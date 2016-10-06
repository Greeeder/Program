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

#include "TConcentrico.h"
#include "TTubo.h"
#include "TDPF.h"

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentrico::TConcentrico() {

	FTubo = NULL;
	FDPF = NULL;

	FSUMTime = 0.;
	FCicloTubo = 0;
	FCicloActual = 0;
	FNumeroNodosTemp = 3;

	FTg = NULL;
	FTpantant = NULL;
	FTpantpos = NULL;
	FTpant0 = NULL;
	FTpant1 = NULL;
	FTpant2 = NULL;
	FTPared = NULL;
	FTParedAnt = NULL;

	FRg_int = NULL;
	FRg_int_ext = NULL;
	FRg_ext_int = NULL;
	FR_ext = NULL;
	FR_int_radiacion = NULL;
	FR_int_RadExt = NULL;
	FR_int_AxiAnt = NULL;
	FR_int_AxiPos = NULL;
	FR_int_RadInt = NULL;
	FR_ext_RadExt = NULL;
	FR_ext_AxiAnt = NULL;
	FR_ext_AxiPos = NULL;
	FR_ext_RadInt = NULL;

	FCapIntExt = NULL;
	FCapIntMed = NULL;
	FCapIntInt = NULL;
	FCapExtExt = NULL;
	FCapExtMed = NULL;
	FCapExtInt = NULL;

	FSUMTPTuboIntPro = NULL;
	FSUMTPTuboExtPro = NULL;

	FHayDPF = false;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

TConcentrico::~TConcentrico() {

	for(int j = 0; j < FNumeroTubos; j++) {
		for(int i = 0; i < FNumeroNodosTemp; i++) {
			delete[] FTParedAnt[j][i];
			delete[] FTPared[j][i];
		}
		delete[] FTParedAnt[j];
		delete[] FTPared[j];
	}

	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 2; j++) {
			delete[] FSUMTPTuboExtPro[i][j];
			delete[] FSUMTPTuboExtPro[i][j];
		}
		delete[] FSUMTPTuboIntPro[i];
		delete[] FSUMTPTuboExtPro[i];
	}

	delete[] FTParedAnt;
	delete[] FTPared;
	delete[] FTg;
	delete[] FTpantant;
	delete[] FTpantpos;
	delete[] FTpant0;
	delete[] FTpant1;
	delete[] FTpant2;
	delete[] FRg_int;
	delete[] FRg_int_ext;
	delete[] FRg_ext_int;
	delete[] FR_ext;
	delete[] FR_int_radiacion;
	delete[] FR_int_RadExt;
	delete[] FR_int_AxiAnt;
	delete[] FR_int_AxiPos;
	delete[] FR_int_RadInt;
	delete[] FR_ext_RadExt;
	delete[] FR_ext_AxiAnt;
	delete[] FR_ext_AxiPos;
	delete[] FR_ext_RadInt;

	delete[] FCapIntExt;
	delete[] FCapIntMed;
	delete[] FCapIntInt;
	delete[] FCapExtExt;
	delete[] FCapExtMed;
	delete[] FCapExtInt;

	delete[] FSUMTPTuboIntPro;
	delete[] FSUMTPTuboExtPro;

	delete[] FTubo;
	delete[] FDPF;

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentrico::PutTiempoDPF(double Tiempo) {
	try {
		FDPF[0]->putTime1DPF(Tiempo);
		FDPF[0]->putDeltaTimeDPF(FDPF[0]->getTime1DPF() - FDPF[0]->getTime0DPF());

		for(int k = 0; k < FDPF[0]->getNumeroHacesCanales(); k++) {
			(FDPF[0]->GetCanal(k, 0))->putTime1(FDPF[0]->getTime1DPF());
			(FDPF[0]->GetCanal(k, 0))->putDeltaTime(FDPF[0]->getDeltaTimeDPF());
			(FDPF[0]->GetCanal(k, 1))->putTime1(FDPF[0]->getTime1DPF());
			(FDPF[0]->GetCanal(k, 1))->putDeltaTime(FDPF[0]->getDeltaTimeDPF());
		}

	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Escritura Time1 en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TConcentrico::GetTiempo(int i) {
	try {
		return FTubo[i]->getTime1();
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Time1 en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TConcentrico::GetTiempoDPF() {
	try {
		return FDPF[0]->getTime1DPF();
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Time1DPF en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TConcentrico::PutTiempo(int i, double Tiempo) {
	try {
		FTubo[i]->PutTime1(Tiempo);
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Escritura Time1 en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int TConcentrico::GetNumTubo(int i) {
	try {
		return FTubo[i]->getNumeroTubo();
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Numero de tubo en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int TConcentrico::GetNumDPF() {
	try {
		return FNumDPFInterna;
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Numero de DPF en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

int TConcentrico::GetNumTuboExterno() {
	try {
		return FNumTuboExterno;
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Numero de tubo externo en el tubo concentrico: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TConcentrico::GetEspesorFluido() {
	try {
		return FEspesor_fluido;
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Espesor de la capa de fluido: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TConcentrico::GetTPared(int j, int k, int i) {
	try {
		return FTPared[j][k][i];
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Temperatura de pared: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

double TConcentrico::GetRg_ext_int(int i) {
	try {
		return FRg_ext_int[i];
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion Resistencia Conveccion Capa de aire con tubo exterior: " << FNumeroConcentrico <<
			 endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

bool TConcentrico::GetHayDPF() {
	try {
		return FHayDPF;
	} catch(exception &N) {
		cout << "ERROR: TConcentrico::Peticion saber si hay DPF: " << FNumeroConcentrico << endl;
		cout << "Tipo de error: " << N.what() << endl;
		throw Exception(N.what());
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#pragma package(smart_init)
