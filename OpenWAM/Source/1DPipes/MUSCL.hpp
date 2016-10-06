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
 * \file MUSCL.hpp
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
 * This file declares a MUSCL finite-volume integrator for
 * one-dimensional pipes.
 */

#ifndef MUSCL_hpp
#define MUSCL_hpp

#include "TPipe.hpp"
#include "PipeMethod.hpp"
#include "Math_wam.h"
#include "RiemannSolvers.hpp"
#include "Godunov.hpp"
#include "Limiters.hpp"

/*!
 * \brief A MUSCL integrator.
 *
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * 
 * A MUSCL integrator.  It is second order-accurate in both
 * time and space.
 */
class TMUSCL: public TGodunov {
  protected:
	Limiter_pt Flim; //!< Slope limiter function.
	RowVector Fphi; //!< Slope limiter values.
	RowArray Fratio; //!< Ratio of successive gradients for the state vector.
	RowArray FU12; //!< State vector, prediction step.
	RowArray FF1; //!< Flux vector, prediction step.
	RowArray FF2; //!< Flux vector, correction step.
	RowArray Fphi_l; //!< Slope limiter values, left side.
	RowArray Fphi_r; //!< Slope limiter values, right side.
	RowArray FV3; //!< Extra source terms.

	/*!
	 * \brief Computes the left and right extrapolated state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * The left and right extrapolated state vector is computed using
	 * a slope-limited, linear extrapolation.
	 *
	 * \param U State vector.
	 */
	void ComputeExtrapolatedValues(const RowArray & U);

	void ComputeSlopeLimiter(const RowArray & U);

  public:
	/*!
	 * \brief Default constructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Initialises the integrator with default values.
	 */
	TMUSCL();

	/*!
	 * \brief Constructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Initialises the integrator and attaches it to a pipe.
	 *
	 * \param pipe Pipe to attach to.
	 */
	TMUSCL(const Pipe_ptr & pipe);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TMUSCL();

	/*!
	 * \brief Connects the MUSCL integrator to a pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param pipe Pipe to connect to.
	 */
	virtual void Connect(const Pipe_ptr & pipe);

	/*!
	 * \brief Solves one Euler step.
	 * 
	 * \param U State vector.
	 * \param V Extra sources vector.
	 * \param t Current time. [s]
	 * \param dt Time-step. [s]
	 */
	void EulerStep(const RowArray& U, const RowArray& V, double t, double dt);

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * The updated state vector is computed using a MUSCL scheme, using Heun's
	 * method for the time integration.
	 */
	virtual void IntegrateWithoutUpdating();

	/*!
	 * \brief Sets the slope limiter function.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * The slope limiter function is needed to limit the slope of the linear
	 * extrapolation for computing the values of the state vector at the cell
	 * interfaces, so the method is 2nd order TVD.
	 *
	 * \param limiter Limiter function.
	 */
	void setLimiter(Limiter_pt limiter);

	/*!
	 * \brief Integrate the flow.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct.
	 *
	 * Solves the fluxes between cells, the boundary condition fluxes and
	 * the source terms. Then, it solves the finite-volume problem and updates
	 * the state vector of #FPipe. The attribute TBasicPipe::FCurrentTime of
	 * #FPipe is also updated with the current simulation time.
	 */
	virtual void Solve();
};


/*!
 * \brief Sets TMUSCL as the computational method for a pipe.
 *
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * 
 * After using this function, the pipe will have a TMUSCL method as its
 * computational method and the method will know which pipe will have its
 * ownership. An approximate Riemann solver function is also needed for
 * the TMUSCL method. Finally, a limiter function is used in order to compute
 * the slope-limited extrapolated values at the cell interfaces.
 *
 * \param pipe Pipe that will be computed using a TMUSCL object.
 * \param rs Riemann solver for the TMUSCL object.
 * \param limiter Limiter function for the TMUSCL object.
 */
void use_muscl(const Pipe_ptr & pipe, const RiemannSolver_pt & rs,
	const Limiter_pt & limiter);

#endif
