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
 * \file PipeMethod.hpp
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
 * This file declares a pipe computation method.
 */

#ifndef PipeMethod_hpp
#define PipeMethod_hpp

#include "BasicPipeMethod.hpp"
#include "TPipe.hpp"

/*!
 * \brief A pipe computation method prototype.
 *
 * It is used for time-integrating the flow inside a pipe, and for computing
 * some of the flow properties.
 */
class TPipeMethod: public TBasicPipeMethod {
  protected:
	double FMaxTimeStep; ///< Maximum time-step. [s]

	/*!
	 * \brief Computes a characteristic in a point.
	 */
	double ComputeCharacteristic(double& velocidadp, double& asonidop, int ind, double dist, int signo, double entropia,
								 double DeltaTiempo);

	/*!
	 * \brief Computes the entropy in a given point.
	 */
	double ComputeEntropy(double& velocidadp, int ind, double dist, int signo, double DeltaTiempo, int indiceCC);

	/*!
	 * \brief Computes the characteristics at the pipe ends.
	 *
	 * \param dt Time-step. [s]
	 */
	void ComputePipeEndCharacteristics(double dt);

  public:
	/*!
	 * \brief Connects the method to a pipe.
	 *
	 * \param pipe Pipe to connect to.
	 */
	virtual void Connect(const Pipe_ptr & pipe);

	/*!
	 * \brief Gets the Courant number.
	 *
	 * Gets the Courant number. The maximum allowable time-step is limite
	 * due to the Courant-Friedrichs-Lewy condition, so the eigenvector of the
	 * problem travelling at maximum speed (having the maximum eigenvalue)
	 * doesn't travel more than a given fraction of the mesh size, being this
	 * fraction the Courant number.
	 */
	virtual double getCourant() const;

	/*!
	 * \brief Returns the maximum allowable time-step.
	 *
	 * Returns the maximum allowable time-step due to stability criteria,
	 * probably due to the Courant-Friedrichs-Lewy condition.
	 *
	 * \return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep();

	/*!
	 * \brief Returns the name of the method.
	 *
	 * \author Luis Miguel Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \return The name of the method.
	 */
	virtual std::string getName() const;

	/*!
	 * \brief Interpolates the left-travelling characteristic.
	 *
	 * \param pipe_end Pipe end for which beta will be computed.
	 * \param dt Time step. [s]
	 */
	void InterpolateBeta(nmPipeEnd pipe_end, double dt);

	/*!
	 * \brief Interpolates a characteristic.
	 *
	 * \param pipe_end Pipe end for which the characteristic will be computed.
	 * \param dt Time step. [s]
	 */
	double InterpolateCharacteristic(double entropia, int signo, int extremo, double DeltaTiempo);

	/*!
	 * \brief Interpolates the entropy level.
	 *
	 * \param pipe_end Pipe end for which the entropy will be computed.
	 * \param dt Time step. [s]
	 *
	 * \return The entropy level.
	 */
	double InterpolateEntropy(nmPipeEnd pipe_end, double dt);

	/*!
	 * \brief Interpolates the right-travelling characteristic.
	 *
	 * \param pipe_end Pipe end for which lambda will be computed.
	 * \param dt Time step. [s]
	 */
	void InterpolateLambda(nmPipeEnd pipe_end, double dt);

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
	virtual void setPtTtU(double p_t, double T_t, double u);

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
		const RowVector& u);
};

#endif
