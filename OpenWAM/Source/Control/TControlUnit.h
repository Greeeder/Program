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
#ifndef TControlUnitH
#define TControlUnitH

#include "Globales.h"
#include "TSensor.h"
//---------------------------------------------------------------------------

/*! This oject is used to control parameter in the engine (rack positions, fuel injected ...) */
class TControlUnit {
  protected:

	std::vector<double> FInputs;
	std::vector<double *> FOutputs;

  public:
	/*! Constructor of the controller class*/
	  TControlUnit();

	/*! Destructor of the controller class*/
	~TControlUnit();

	void AppendOutput_ptr(int i, double *output);

	void setOutput(int i, double value);

};

typedef shared_ptr<TControlUnit> ControlUnit_ptr;


#endif
