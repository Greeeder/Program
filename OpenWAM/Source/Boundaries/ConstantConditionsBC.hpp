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
 * \file ConstantConditionsBC.hpp
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
 * This file declares a boundary condition with constant pressure,
 * temperature and speed, including some helper functions.
 */

#ifndef ConstantConditionsBC_hpp
#define ConstantConditionsBC_hpp

#include "PipeConnection.hpp"
#include "TConstantConditionsPipe.hpp"

/*!
 * \brief A BC with constant pressure, temperature and speed.
 */
class TConstantConditionsBC: public TPipeConnection {

	friend void attach_to_constant_BC(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, double u, TFluid_ptr fluid);

  protected:
	ConstantConditionsPipe_ptr FCCPipe; //!< Pipe representing the CC BC.

	double FM;
	double FFLUXU;
	double FH;
	double FDT;

	/*!
	 * \brief Default constructor.
	 */
	TConstantConditionsBC();

	/*!
	 * \brief Constructor, only used from friend functions.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param pipe Pipe attached to this BC.
	 * \param pipe_end Pipe end connected to this BC.
	 * \param pipe_cc Constant conditions pipe attached to this BC.
	 * \param virtual_pipe Virtual pipe used by this BC.
	 */
	TConstantConditionsBC(const Pipe_ptr & pipe, nmPipeEnd pipe_end,
		const ConstantConditionsPipe_ptr & pipe_cc,
		const VirtualPipe_ptr & virtual_pipe);
	
  public:
	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param t Time for which the flux is going to be computed. [s]
	 * \param dt Time-step. [s]
	 * \return The BC flux vector.
	 */
	virtual ColVector Flux(double t, double dt);

	/*!
	 * \brief Sets the pressure, temperature and speed of the BC virtual cell.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param p BC pressure. [Pa]
	 * \param T BC temperature. [K]
	 * \param u BC speed. [m / s]
	 */
	void setPTU(double p, double T, double u = 0);
};

/*!
 * \brief Unique pointer to a TConstantConditionsBC object.
 */
typedef unique_ptr<TConstantConditionsBC> ConstantConditionsBC_ptr;

/*!
 * \brief Joins a pipe to a constant conditions BC.
 *
 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * \date 2016/03/15
 * 
 * \param pipe Pipe connected to this boundary condition.
 * \param pipe_end Pipe end connected to the boundary condition.
 * \param p Constant pressure at the BC. [Pa]
 * \param T Constant temperature at the BC. [K]
 * \param u Constant speed at the BC. [m / s]
 * \param fluid Pointer to the fluid.
 */
void attach_to_constant_BC(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, double u, TFluid_ptr fluid);

#endif
