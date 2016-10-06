/**
* @file HeatTransferCombMECSwirl.h
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
* This file solve heat transfer in CIE with swirl.
*/

#include <string>
#include <vector>
//#include "Comun.h"
#include "HeatTransfer.h"


#ifndef __HEATTRANSFERCOMBMECSWIRL_H
#define __HEATTRANSFERCOMBMECSWIRL_H



class HeatTransferCombMECSwirl : public HeatTransfer
{
public:

	void Calcula_Calor_Woschni(double T, double p);

};

#endif