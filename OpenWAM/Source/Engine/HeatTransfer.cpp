/**
* @file HeatTransfer.h
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @date 30 de mar. de 2016
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
* This file solve heat transfer in the cylinder.
*/


#include <iostream>
#include <fstream>
#include <sstream>
#include "Constantes.h"
#include "Math.h"
#include "HeatTransfer.h"
#include "TCrankMechanism.h"

using namespace std;

void HeatTransfer::Inicializar(TCrankMechanism *CrankMechanism, stInicializar Inicializar)
{
	_inicializar = Inicializar;
}

void HeatTransfer::InicializarAlRCA(stInicializarAlRCA InicializarAlRCA)
{
	_inicializarAlRCA = InicializarAlRCA;
}

void HeatTransfer::InicializarPorAngulo(stInicializarPorAngulo InicializarPorAngulo)
{
	_inicializarPorAngulo = InicializarPorAngulo;
}


double HeatTransfer::Xp(double x, double ratio_ctms)
{
	//curva propuesta por Antonio Gil
	return ratio_ctms + (1 / (pow((cosh(x / 100)), 40) + (ratio_ctms / (1 - ratio_ctms))));

}


void HeatTransfer::Calcula_Calor_Woschni(double T, double p)
{

}

double HeatTransfer::Heat(double T, double p) {

	Calcula_Calor_Woschni(T, p);
	return _QRAD + _QW;
}

void HeatTransfer::ReadData(xml_node node)
{

}
