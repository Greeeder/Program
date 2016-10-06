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

/*!
 * \file Connect0Dto1D.hpp
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
 * This file declares a boundary condition between a 0D elements to a
 * 1D element.
 */

#ifndef Connect0Dto1D_hpp
#define Connect0Dto1D_hpp

#include "ConstantConditionsBC.hpp"

class TBasicPlenum;

/*!
 * \brief A BC with constant total pressure and total temperature.
 */
class TConnect0Dto1D: public TBoundaryCondition {

	friend void attach_to_0Dto1Dconnection(const int id, const TFlowObject_ptr& pipe, nmPipeEnd pipe_end,
	const TFlowObject_ptr& zeroD, TTipoValvula_ptr valve);

  protected:

	  TBasicPlenum* ZeroD;
	  TPipe* OneD;

	  double T_buffer;
	  bool InletFlow;
	/*!
	 * \brief Default constructor.
	 */
	TConnect0Dto1D();

	/*!
	 * \brief Constructor, only used from friend functions.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param pipe Flow object attached to this BC.
	 * \param pipe_end Pipe end connected to this BC.
	 * \param virtual_pipe Virtual pipe used by this BC.
	 */

	TConnect0Dto1D(const TFlowObject_ptr & pipe, nmPipeEnd pipe_end,
		const VirtualPipe_ptr & virtual_pipe);

	/*!
	* \brief Constructor, only used from friend functions.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	* \date 2016/03/15
	*
	* \param pipe Flow object attached to this BC.
	* \param con_index Index of the entry connected to this BC.
	* \param virtual_pipe Virtual pipe used by this BC.
	*/

	TConnect0Dto1D(const TFlowObject_ptr & pipe, int con_index,
		const VirtualPipe_ptr & virtual_pipe);
	
  public:
	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/09
	 * 
	 * \param t Time for which the flux is going to be computed. [s]
	 * \param dt Time-step. [s]
	 * \return The BC flux vector.
	 */
	virtual ColVector Flux(double t, double dt);

	/*!
	* \brief	Gets the fluid at the boundary.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	18/05/2016
	*
	* \return	The fluid pointer.
	*/

	virtual TFluid_ptr FluidBC();

	/*!
	 * \brief Sets the total pressure and temperature of the BC virtual cell.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/09
	 *
	 * \param p BC pressure. [Pa]
	 * \param T BC temperature. [K]
	 */
	void setPT(double p, double T);

	/*!
	* \brief Updates the BC characteristics.
	*
	* \param t Current time. [s]
	* \param dt Time-step. [s]
	*/
	virtual void updateBCCharacteristics(double t, double dt){}
};

/*!
 * \brief Unique pointer to a TPressureBC object.
 */
typedef unique_ptr<TConnect0Dto1D> Connect0Dto1D_ptr;

/*!
 * \brief Joins a pipe to a pressure BC.
 *
 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * \date 2016/05/09
 * 
 * \param pipe Pipe connected to this boundary condition.
 * \param pipe_end Pipe end connected to the boundary condition.
 * \param zeroD 0-dimensional element (cylinder, plenum ...) connecte to this boundary conditions
 * \param entry Entry number of the 0-dimensional element
 * \param fluid Pointer to the fluid
 */
void attach_to_0Dto1Dconnection(const int id, const TFlowObject_ptr& pipe, nmPipeEnd pipe_end,
	const TFlowObject_ptr& zeroD, TTipoValvula_ptr valve);

#endif
