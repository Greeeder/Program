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
 * \file PipeConnection.hpp
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
 * This file declares a pipe connection and helper functions.
 */

#ifndef PipeConnection_hpp
#define PipeConnection_hpp

#include "BoundaryCondition.hpp"
#include "TVirtualPipe.hpp"
#include "Godunov.hpp"
#include <utility>

/*!
 * \brief A connection between two pipes.
 */
class TPipeConnection: public TBoundaryCondition {
	friend void attach_pipes(const Pipe_ptr pipe_0, nmPipeEnd pipe_end_0,
		const Pipe_ptr pipe_1, nmPipeEnd pipe_end_1);
  protected:

	/*!
	 * \brief Updates the BC characteristics.
	 *
	 * \param t Current time. [s]
	 * \param dt Time-step. [s]
	 */
	virtual void updateBCCharacteristics(double t, double dt);

  public:

	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 */
	TPipeConnection();

	/*!
	 * \brief Creates a pipe connection object.
	 *
	 * \param pipe_0 Pipe attached to this connection.
	 * \param pipe_end_0 Pipe end attached to this connection.
	 * \param virtual_pipe Virtual pipe used by this connection.
	 * \param pipe_number 0 if it is the first pipe, 1 if it is the second.
	 */
	TPipeConnection(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end_0,
		const VirtualPipe_ptr & virtual_pipe, unsigned int pipe_number);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TPipeConnection();

	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * The fluxes are computed as the fluxes in a virtual pipe with two cells,
	 * one corresponding to one of the attached pipes and the other one to
	 * one cell of the other pipe.
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
	 * \brief Gets the effective area of the connection.
	 * 
	 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return Effective area of the connection. [m ** 2]
	 */
	double getAeff() const;

	/*!
	 * \brief Sets the effective area of the connection.
	 * 
	 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param Aeff Effective area of the connection. [m ** 2]
	 */
	void setAeff(double Aeff);
};

/*!
 * \brief Joins two pipes.
 *
 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * 
 * \param pipe_0 First pipe connected to this boundary condition.
 * \param pipe_end_0 First pipe end connected to the boundary condition.
 * \param pipe_1 Second pipe connected to this boundary condition.
 * \param pipe_end_1 Second pipe end connected to the boundary condition.
 */
void attach_pipes(const Pipe_ptr pipe_0, nmPipeEnd pipe_end_0,
	const Pipe_ptr pipe_1, nmPipeEnd pipe_end_1);

#endif
