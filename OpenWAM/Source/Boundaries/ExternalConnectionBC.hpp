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
 * \file ExternalConnectionBC.hpp
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 *
 * \section LICENSE
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
 * \section DESCRIPTION
 * This file declares a boundary condition for an external connection.
 */

#ifndef ExternalConnectionBC_hpp
#define ExternalConnectionBC_hpp

#include "ConstantConditionsBC.hpp"
#include "TCCExternalConnectionVol.h"


/*!
 * \brief An external connection BC.
 */
class TExternalConnectionBC: public TConstantConditionsBC, public TCCExternalConnectionVol {

	friend void attach_to_external_connection(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid, double u, int ID);

  protected:

	/*!
	 * \brief Constructor, only used from friend functions.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param pipe Pipe attached to this BC.
	 * \param pipe_end Pipe end connected to this BC.
	 * \param pipe_cc Constant conditions pipe attached to this BC.
	 * \param virtual_pipe Virtual pipe used by this BC.
	 */
	TExternalConnectionBC(const Pipe_ptr & pipe, nmPipeEnd pipe_end,
		const ConstantConditionsPipe_ptr & pipe_cc,
		const VirtualPipe_ptr & virtual_pipe);
	
  public:

	/*!
	 * \brief Gets the boundary pressure, temperature and speed.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/21
	 * 
	 * u is negative it the speed goes out of the OpenWAM domain.
	 * If this BC is attached to the left end of a pipe, the speed is negative if it
	 * goes to the left. If this BC is attached to the right end of a pipe, the speed
	 * is negative if it goes to the right.
	 * 
	 * \param p Pressure at the BC. [bar]
	 * \param T Temperature at the BC. [K]
	 * \param u Speed at the BC. [m / s]
	 */
	virtual void LoadNewData(double* p, double* T, double* u);

	/*!
	 * \brief Gets the boundary mass, enthalpy and momentum flows.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * m is negative if the flow goes out of the OpenWAM domain.
	 * 
	 * \param m Mass flow rate. [kg / s]
	 * \param mht Total enthalpy flux. [W]
	 * \param mom Momentum flux. [N / s]
	 */
	virtual void LoadFluxes(double* m, double *mh0, double *mom);

	/*!
	 * \brief Sets the pressure, temperature, speed and time of the BC virtual cell.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/21
	 *
	 * u is negative it the speed goes out of the OpenWAM domain.
	 * If this BC is attached to the left end of a pipe, the speed is negative if it
	 * goes to the left. If this BC is attached to the right end of a pipe, the speed
	 * is negative if it goes to the right.
	 *
	 * \param u BC speed. [m / s]
	 * \param T BC temperature. [K]
	 * \param p BC pressure. [Pa]
	 * \param t BC time. [s]
	 */
	virtual void UpdateCurrentExternalProperties(double u, double T,
		double p, double t);
};

/*!
 * \brief Unique pointer to a TExternalConnectionBC object.
 */
typedef unique_ptr<TExternalConnectionBC> ExternalConnectionBC_ptr;

/*!
 * \brief Joins a pipe to an external connection BC.
 *
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 * \date 2016/03/15
 * 
 * \param pipe Pipe connected to this boundary condition.
 * \param pipe_end Pipe end connected to the boundary condition.
 * \param p Constant pressure at the BC. [Pa]
 * \param T Constant temperature at the BC. [K]
 * \param fluid Pointer to the fluid.
 * \param u Constant speed at the BC. [m / s]
 * \param ID BC ID, in case it is an external code BC.
 */
void attach_to_external_connection(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid, double u = 0, int ID = 1);

#endif
