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

/**
* @file TTurbocharger.h
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
* This file defines a class that for a group of elements that form the turbocharger.
*/

#ifndef TTurbocharger_hpp
#define TTurbocharger_hpp

#include "TIntegrable.hpp"

class TTurbocharger: public TIntegrableGroup
{
public:

	TTurbocharger();

	virtual ~TTurbocharger();

	/*!
	* \brief Integrates one time-step without updating the state vector.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*/
	virtual void IntegrateWithoutUpdating();

	/*!
	* \brief Integrate the flow, updating the state vector.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*
	* Integrates the flow evolution inside the duct, updating the state vector
	* at the end. Uses an already-set time-step.
	*/
	virtual void Solve();

	/*!
	* \brief Integrate the flow, updating the state vector.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*
	* Integrates the flow evolution inside the duct, updating the state
	* vector at the end. Computes until the time passed is reached.
	*
	* \param t Time at the end of the integration. [s]
	*/
	virtual void Solve(double t);

	/*!
	* \brief Updates the state vector with the already-computed results.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*
	* After solving the system and computing the state vector update, this
	* function is called to update it.
	*/
	virtual void UpdateStateVector();
};

#endif
