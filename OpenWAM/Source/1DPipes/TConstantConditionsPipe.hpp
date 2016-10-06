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
 * @file TConstantConditionsPipe.hpp
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
 * This file declares a one-dimensional pipe with a constant state vector.
 */

#ifndef TConstantConditionsPipe_hpp
#define TConstantConditionsPipe_hpp

#include "TPipe.hpp"
#include "AllPipeMethods.hpp"

/*!
 * \brief A pipe with a constant state vector.
 */
class TConstantConditionsPipe: public TPipe {

	friend shared_ptr<TConstantConditionsPipe> create_constant_conditions_pipe(
		double p, double T, double u, double D, double l, unsigned int size, TFluid_ptr fluid);
public:
	/**
	 * @brief Returns the maximum allowable time-step.
	 *
	 * Returns the maximum allowable time-step due to stability criteria.
	 *
	 * @return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep();

	/**
	 * @brief Integrate the flow. In this case, only the time advances..
	 *
	 * Integrates the flow evolution inside the duct. As the pipe has
	 * a constant state vector, it only advances in time.
	 */
	virtual void Solve();

	/**
	 * @brief Integrate the flow. In this case, only the time advances.
	 *
	 * Integrates the flow evolution inside the duct, updating the state
	 * vector at the end. Computes until the time passed is reached. As
	 * the pipe has a constant state vector, it only advances in time.
	 *
	 * @param t Time at the end of the integration. [s]
	 */
	virtual void Solve(double t);
};

/**
 * \brief Shared pointer to a TConstantConditionsPipe object.
 */
typedef shared_ptr<TConstantConditionsPipe> ConstantConditionsPipe_ptr;

/*!
 * \brief Creates a pipe with constant and uniform values of
 * pressure, temperature, speed and diameter.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * \date 15/03/2016
 * 
 * \param p Pipe pressure. [Pa]
 * \param T Pipe temperature. [K]
 * \param u Pipe speed. [m / s]
 * \param D Pipe diameter. [m]
 * \param l Pipe length. [m]
 * \param size Number of cells.
 * \param fluid Pointer to the fluid.
 */
ConstantConditionsPipe_ptr create_constant_conditions_pipe(
	double p, double T, double u, double D, double l, unsigned int size, TFluid_ptr fluid);

#endif
