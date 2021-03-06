/*--------------------------------------------------------------------------------*\
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

#include "TTimeControl.h"

//---------------------------------------------------------------------------

__fastcall TTimeControl::TTimeControl() {
	TiempoTotal = 0.;
	TiempoCiclo0 = 0;
	TiempoCiclo1 = 0;
	Procesador = 2.2e9; //Hz
	Pasos = 0;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

__fastcall TTimeControl::~TTimeControl() {
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TTimeControl::Ini() {
	int k = 0;
	asm rdtsc
	asm mov k, eax

	Inicio=k;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TTimeControl::Fin() {
	int k = 0;
	asm rdtsc
	asm mov k, eax

	Final=k;
	TiempoTotal += (double)(Final - Inicio) / Procesador;
	Pasos++;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

void TTimeControl::PromediaCiclo() {
	TiempoCiclo0 = TiempoCiclo1;
	TiempoCiclo1 = TiempoTotal;
	printf("Tiempo total: %g seg.\n", TiempoTotal);
	printf("Tiempo ultimo ciclo: %g seg.\n\n", TiempoCiclo1 - TiempoCiclo0);
	printf("Pasos: %d\n", Pasos);
	printf("Tiempo por Steps: %g\n", TiempoTotal / (double) Pasos);
}
#pragma package(smart_init)
