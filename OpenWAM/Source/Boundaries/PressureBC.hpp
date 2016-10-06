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
 * \file PressureBC.hpp
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
 * This file declares a boundary condition with constant pressure
 * and temperature.
 */

#ifndef PressureBC_hpp
#define PressureBC_hpp

#include "ConstantConditionsBC.hpp"

/*!
 * \brief A BC with constant total pressure and total temperature.
 */
class TPressureBC: public TConstantConditionsBC {

	friend void attach_to_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid);

  protected:
	/*!
	 * \brief Default constructor.
	 */
	TPressureBC();

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
	TPressureBC(const Pipe_ptr & pipe, nmPipeEnd pipe_end,
		const ConstantConditionsPipe_ptr & pipe_cc,
		const VirtualPipe_ptr & virtual_pipe);

	double T_buffer; //!< Buffer temperature, to avoid weird oscillations. [K]
	double FCD; //!< Discharge coefficient.
	
  public:
	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/09
	 * 
	 * \param t Time for which the flux is going to be computed. [s]
	 * \param dt Time-step. [s]
	 * \return The BC flux vector.
	 */
	virtual ColVector Flux(double t, double dt);

	/*!
	 * \brief Gets the discharge coefficient.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/06/20
	 * 
	 * \return The discharge coefficient.
	 */
	double getCD() const;

	/*!
	 * \brief Sets the discharge coefficient.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/06/20
	 * 
	 * \param CD Discharge coefficient.
	 */
	void setCD(double CD);

	/*!
	 * \brief Sets the total pressure and temperature of the BC virtual cell.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/09
	 *
	 * \param p BC pressure. [Pa]
	 * \param T BC temperature. [K]
	 */
	void setPT(double p, double T);
};

/*!
 * \brief Unique pointer to a TPressureBC object.
 */
typedef unique_ptr<TPressureBC> PressureBC_ptr;

/*!
 * \brief Joins a pipe to a pressure BC.
 *
 * \author Luis Miguel Garcia-Cuevas González <luiga12@mot.upv.es>
 * \date 2016/05/09
 * 
 * \param pipe Pipe connected to this boundary condition.
 * \param pipe_end Pipe end connected to the boundary condition.
 * \param p Total pressure at the BC. [Pa]
 * \param T Total temperature at the BC. [K]
 * \param fluid Pointer to the fluid
 */
void attach_to_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid);

#endif
