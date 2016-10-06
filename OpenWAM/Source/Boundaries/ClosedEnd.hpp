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
 * @file ClosedEnd.hpp
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
 * This file declares a closed end boundary condition and helper functions.
 */

#ifndef ClosedEnd_hpp
#define ClosedEnd_hpp

#include "BoundaryCondition.hpp"
#include "TPipe.hpp"

/**
 * @brief A closed end.
 */
class TClosedEnd: public TBoundaryCondition {
  protected:

	/**
	 * @brief Updates the BC characteristics.
	 *
	 * @param t Current time. [s]
	 * @param dt Time-step. [s]
	 */
	virtual void updateBCCharacteristics(double t, double dt);

  public:

	/**
	 * @brief Generic constructor
	 */
	TClosedEnd();

	/**
	 * @brief Initialises the boundary condition.
	 *
	 * @param pipe_0 Pipe connected to this boundary condition.
	 * @param pipe_end Pipe end connected to the boundary condition.
	 */
	TClosedEnd(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end);

	/**
	 * @brief Computes the boundary condition flux vector.
	 *
	 * The mass flow and the enthalpy flow for this boundary condition are
	 * exactly equal to 0. The momentum flux only takes into account the
	 * pressure term.
	 *
	 * @param t Time for which the flux is going to be computed. [s]
	 * @param dt Time-step. [s]
	 * @return The BC flux vector.
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
};

/**
 * @brief Makes a closed end in a pipe.
 *
 * @param pipe Pipe connected to this boundary condition.
 * @param pipe_end Pipe end connected to the boundary condition.
 */
void close_pipe_end(const Pipe_ptr & pipe, nmPipeEnd pipe_end);

#endif
