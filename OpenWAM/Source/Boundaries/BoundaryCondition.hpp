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
 * \file BoundaryCondition.hpp
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
 * This file declares a basic boundary condition.
 */

#ifndef BoundaryCondition_hpp
#define BoundaryCondition_hpp

#include "TPipe.hpp"
#include "TVirtualPipe.hpp"
#include <memory>
#include <utility>

/*!
 * \brief A generic boundary condition.
 */
class TBoundaryCondition {
  protected:
	TPipe * FPipe_0; //!< Pipe connected to this boundary condition.
	TFlowObject *FlowObj; //!< Pipe connected to this boundary condition.
	nmPipeEnd FPipeEnd_0; //!< Pipe end connected to this boundary condition.
	int FCon_Index; //!< Index of the entry connected to this boundary condition.
	unsigned int FCell_0; //!< Cell connected to this boundary condition.
	double FCurrentTime; //!< Current time for the BC. [s]
	ColVector FFlux_0; //!< BC flux, for memoisation purposes.
	double FLambda; //!< Right-travelling characteristic.
	double FBeta; //!< Left-travelling characteristic.
	double FEntropy; //!< Entropy level.
	double FPressure; //!< Pressure at the BC. [Pa]
	double FSpeed; //!< Speed at the BC. [m / s]
	double FSpeedOfSound; //!< Speed of sound at the BC. [m / s]
	double FTemperature; //!< Temperature at the BC. [K]
	double FArea; //!< BC area. [m ** 2]

	TFluid_ptr WorkFluid; //!< Pointer to the working fluid at Boundary Conditions 

	VirtualPipe_ptr FVirtualPipe; //!< Pipe representing the connection.
	unsigned int FPipeNumber; //!< Pipe number for the virtual pipe.

	/*!
	 * \brief Updates the BC characteristics.
	 *
	 * \param t Current time. [s]
	 * \param dt Time-step. [s]
	 */
	virtual void updateBCCharacteristics(double t, double dt) = 0;

  public:
	/*!
	 * \brief Generic constructor
	 */
	TBoundaryCondition();

	/*!
	 * \brief Initialises the boundary condition.
	 *
	 * \param pipe_0 Pipe connected to this boundary condition.
	 * \param pipe_end Pipe end connected to the boundary condition.
	 */
	TBoundaryCondition(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end);

	/*!
	* \brief Initialises the boundary condition.
	*
	* \param pipe_0 Pipe connected to this boundary condition.
	* \param pipe_end Pipe end connected to the boundary condition.
	*/
	TBoundaryCondition(const TFlowObject_ptr & pipe_0, nmPipeEnd pipe_end);

	/*!
	* \brief Initialises the boundary condition.
	*
	* \param zerod 0-Dimensonal object connected to this boundary condition.
	* \param con_id Index of the element entyr connected to the boundary condition.
	*/
	TBoundaryCondition(const TFlowObject_ptr & zerod, int con_id);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TBoundaryCondition();

	/*!
	 * \brief Computes the pressure, temperature and speed at the BC.
	 *
	 * Computes the pressure, temperature and speed at the BC using the
	 * method of characteristics.
	 *
	 * \param t Current time. [s]
	 * \param dt Time-step. [s]
	 */
	void computePTUFromCharacteristics(double t, double dt);

	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * \param t Time for which the flux is going to be computed. [s]
	 * \param dt Time-step. [s]
	 * \return The BC flux vector.
	 */
	virtual ColVector Flux(double t, double dt) = 0;

	/*!
	 * \brief	Gets the fluid at the boundary.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/05/2016
	 *
	 * \return	The fluid pointer.
	 */

	virtual TFluid_ptr FluidBC() = 0;

	/*!
	 * \brief Gets the BC area.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016-05-02
	 * 
	 * \return The BC area. [m ** 2]
	 */
	double getArea() const;

	/*!
	 * \brief Gets the left-travelling characteristic.
	 *
	 * \return Left-travelling characteristic.
	 */
	double getBeta() const;

	/*!
	 * \brief Returns the total enthalpy flow through the BC.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 * 
	 * \return BC total enthalpy flow. [W]
	 */
	double getEnthalpyFlow() const;

	/*!
	 * \brief Gets the entropy level.
	 *
	 * \return The entropy level.
	 */
	double getEntropy() const;

	/*!
	 * \brief	Gets the fluid at the boundary.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/04/2016
	 *
	 * \return	The fluid pointer.
	 */

	TFluid_ptr getFluid() const{ return WorkFluid; }

	/*!
	 * \brief Gets the right-travelling characteristic.
	 *
	 * \return Right-travelling characteristic.
	 */
	double getLambda() const;

	/*!
	 * \brief Returns the mass flow rate through the BC.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 * 
	 * \return BC mass flow rate. [kg / s]
	 */
	double getMassFlow() const;

	/*!
	 * \brief Returns the BC pressure.
	 *
	 * \return BC pressure. [Pa]
	 */
	double getPressure() const;

	/*!
	 * \brief Returns the BC speed.
	 *
	 * \return BC speed. [m / s]
	 */
	double getSpeed() const;

	/*!
	 * \brief Returns the BC temperature.
	 *
	 * \return BC temperature. [K]
	 */
	double getTemperature() const;

	VirtualPipe_ptr getVirtualPipe() const {
		return FVirtualPipe;
	}

	/*!
	 * \brief Sets the left-travelling characteristic.
	 *
	 * \param Left-travelling characteristic.
	 */
	void setBeta(double beta);

	/*!
	 * \brief Sets the entropy level.
	 *
	 * \param The entropy level.
	 */
	void setEntropy(double entropy);

	/*!
	 * \brief Sets the right-travelling characteristic.
	 *
	 * \param Right-travelling characteristic.
	 */
	void setLambda(double lambda);
};

typedef unique_ptr<TBoundaryCondition> BoundaryCondition_ptr;

#endif
