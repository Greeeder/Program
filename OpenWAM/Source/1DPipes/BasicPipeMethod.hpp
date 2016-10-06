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
 * \file BasicPipeMethod.hpp
 * \author Francisco Jose Arnau <farnau@mot.upv.es>
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
 * This file declares a pure virtual pipe computation method.
 */

#ifndef BasicPipeMethod_hpp
#define BasicPipeMethod_hpp

#include "Math_wam.h"
#include <memory>

class TPipe;

/*!
 * \brief A pure virtual pipe computation method prototype.
 *
 * It is used for time-integrating the flow inside a pipe, and for computing
 * some of the flow properties.
 */
class TBasicPipeMethod {
  protected:
	double FMaxCourant; //!< Maximum Courant number.
	std::string FName; //!< Method name.

	TPipe * FPipe; //!< Pipe that uses the method.

	/*!
	 * \brief Computes the maximum allowable time-step.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void ComputeMaxTimeStep() = 0;

  public:
	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TBasicPipeMethod() {};

	/*!
	 * \brief Connects the method to a pipe.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param pipe Pipe to connect to.
	 */
	virtual void Connect(const shared_ptr<TPipe> & pipe) = 0;

	/*!
	 * \brief Gets the Courant number.
	 *
	 * Gets the Courant number. The maximum allowable time-step is limite
	 * due to the Courant-Friedrichs-Lewy condition, so the eigenvector of the
	 * problem travelling at maximum speed (having the maximum eigenvalue)
	 * doesn't travel more than a given fraction of the mesh size, being this
	 * fraction the Courant number.
	 */
	virtual double getCourant() const = 0;

	/*!
	 * \brief Returns the maximum allowable time-step.
	 *
	 * Returns the maximum allowable time-step due to stability criteria,
	 * probably due to the Courant-Friedrichs-Lewy condition.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep() = 0;

	/*!
	 * \brief Returns the already computed fluxes between cells.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \return Fluxes at the cell boundaries.
	 */
	virtual RowArray getFluxes() const = 0;

	/*!
	 * \brief Returns the name of the method.
	 *
	 * \author Luis Miguel Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \return The name of the method.
	 */
	virtual std::string getName() const = 0;

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void IntegrateWithoutUpdating() = 0;

	/*!
	 * \brief Sets the Courant number.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016-05-04
	 * 
	 * Sets the Courant number. The maximum allowable time-step is limited
	 * due to the Courant-Friedrichs-Lewy condition, so the eigenvector of the
	 * problem travelling at maximum speed (having the maximum eigenvalue)
	 * doesn't travel more than a given fraction of the mesh size, being this
	 * fraction the Courant number.
	 * 
	 * \param C Courant number.
	 */
	void setCourant(double C);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given total pressure, total temperature
	 * and flow speed.
	 *
	 * \param p_t Total pressure. [Pa]
	 * \param T_t Total temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	virtual void setPtTtU(double p_t, double T_t, double u) = 0;

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given total pressure, total temperature
	 * and flow speed, one set of values for each node/cell.
	 *
	 * \param p_t Total pressure. [Pa]
	 * \param T_t Total temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	virtual void setPtTtU(const RowVector& p_t, const RowVector& T_t,
		const RowVector& u) = 0;

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given pressure, temperature and flow speed.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	virtual void setPTU(double p, double T, double u) = 0;

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given pressure, temperature and flow speed,
	 * one set of values for each node/cell.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	virtual void setPTU(const RowVector& p, const RowVector& T, const RowVector& u) = 0;

	/*!
	 * \brief Integrate the flow.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct.
	 */
	virtual void Solve() = 0;

	/*!
	 * \brief	Updates the composition at each cell.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/05/2016
	 */

	virtual void UpdateComposition() = 0;

	/*!
	 * \brief Updates the flow variables with the current state vector values.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void UpdateFlowVariables() = 0;

	/*!
	 * \brief Updates R, gamma and company.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void UpdateGasProperties() = 0;
};

#endif
