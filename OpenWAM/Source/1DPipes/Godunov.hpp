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
 * \file Godunov.hpp
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
 * This file declares a Godunov finite-volume integrator for
 * one-dimensional pipes.
 */

#ifndef Godunov_hpp
#define Godunov_hpp

#include "TPipe.hpp"
#include "PipeMethod.hpp"
#include "Math_wam.h"
#include "RiemannSolvers.hpp"

class TVirtualPipe;

/*!
 * \brief A Godunov integrator.
 *
 * A Godunov integrator.  It is first order-accurate in both
 * time and space. It is a density-based, one-dimensional,
 * finite-volume integrator for computing the Euler equations of
 * fluid dynamics.
 */
class TGodunov: public TPipeMethod {
	friend TVirtualPipe;
  protected:
	RowArray FW; ///< Flux vector;
	RowArray FF; ///< Flux vector at cell interfaces.
	RowArray FFY; ///< Flux vector at cell interfaces for each chemical component [kg / s / m ** 2]
	RowArray FV1; ///< Source terms due to area.
	RowArray FV2; ///< Source terms due to friction and heat.
	RiemannSolver_pt Frs; ///< Approximate Riemann solver.
	RowVector Fa_l; ///< Speed of sound, left extrapolation. [m / s]
	RowVector Fa_r; ///< Speed of sound, right extrapolation. [m / s]
	RowVector Fe_l; ///< Internal energy, left extrapolation.
	RowVector Fe_r; ///< Internal energy, left extrapolation.
	RowVector Fp_l; ///< Pressure, left extrapolation. [Pa]
	RowVector Fp_r; ///< Pressure, right extrapolation. [Pa]
	RowVector Frho_l; ///< Density, left extrapolation [kg / m ** 3]
	RowVector Frho_r; ///< Density, right extrapolation. [kg / m ** 3]
	RowArray FU_l; ///< State vector, left extrapolation.
	RowArray FU_r; ///< State vector, right extrapolation.
	RowVector Fu_l; ///< Speed, left extrapolation. [m / s]
	RowVector Fu_r; ///< Speed, right extrapolation. [m / s]
	RowArray FW_l; ///< Flux vector, left extrapolation.
	RowArray FW_r; ///< Flux vector, right extrapolation.

	/*!
	 * \brief Computes the maximum allowable time-step.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void ComputeMaxTimeStep();

  public:
	/*!
	 * \brief Default constructor.
	 *
	 * Initialises the propagator with default values.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	TGodunov();

	/*!
	 * \brief Constructor.
	 *
	 * Initialises the propagator and attaches it to a pipe.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 *
	 * \param pipe Pipe to attach to.
	 */
	TGodunov(const Pipe_ptr & pipe);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TGodunov();

	/*!
	 * \brief Computes the flux vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param U State vector.
	 * \param W Flux vector.
	 * \param Gamma Specific heat capacities ratio.
	 * \param Gamma1 Specific heat capacities ratio minus 1.
	 */
	void ComputeFlux(const RowArray & U, RowArray & W, const RowVector & Gamma, const RowVector & Gamma1);

	/*!
	 * \brief Computes the source terms due to area change.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param U State vector.
	 * \param V Source terms.
	 * \param A Area. [m ** 2]
	 * \param Gamma1 Gamma - 1.
	 */
	void ComputeSource1(const RowArray & U, RowArray & V, const RowVector & A, const RowVector & Gamma1);

	/*!
	 * \brief	Calculate the source terms related with friction and heat transfer.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	16/05/2016
	 *
	 * \param	U		 	State vector.
	 * \param [in,out]	V	Source terms.
	 * \param	D		 	Cell diameter [m].
	 * \param	f		 	Cell friction coefficient [-].
	 * \param	q		 	Heat transfer [W].
	 * \param	dx		 	Meash size [m].
	 */

	void ComputeSource2(const RowArray & U, RowArray & V, const RowVector & D, const RowVector & f, const RowVector & q, const double dx);

	/*!
	 * \brief Computes the source terms due to friction and heat transfer.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param U State vector.
	 * \param V Source terms.
	 * \param A Area. [m ** 2]
	 * \param h Internal heat transfer coefficient. [W / (m ** 2 * K)]
	 * \param rho Density. [kg / (m ** 3)]
	 * \param Re Reynolds number.
	 * \param Twall Wall temperature. [K]
	 * \param Gamma Specific heat capacities ratio.
	 * \param R Gas constant. [J / (kg * K)]
	 * \param Gamma1 Gamma - 1.
	 */
	void ComputeSource2(const RowArray & U, RowArray & V, const RowVector & A, const RowVector & h, const RowVector & rho,
						const RowVector & Re, const RowVector & TWall, const RowVector & Gamma1);

	/*!
	 * \brief Connects the Godunov integrator to a pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param pipe Pipe to connect to.
	 */
	virtual void Connect(const Pipe_ptr & pipe);

	/*!
	 * \brief Returns the already computed fluxes between cells.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Fluxes at the cell boundaries.
	 */
	virtual RowArray getFluxes() const;

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 */
	virtual void IntegrateWithoutUpdating();

	/*!
	 * \brief Integrate the flow one time-step.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct.
	 *
	 * Solves the fluxes between cells, the boundary condition fluxes and
	 * the source terms. Then, it solves the finite-volume problem and updates
	 * the state vector of #FPipe. The attribute TPipe::FCurrentTime of
	 * #FPipe is also updated with the current simulation time.
	 */
	virtual void Solve();

	/*!
	 * \brief Integrates the central cells of the pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * During the integration of the flow, the method computes the flow
	 * between cells and at the pipe ends. This method solves the flow
	 * in the interfaces between cells inside #FPipe.
	 */
	virtual void SolveCentralCells();

	/*!
	 * \brief Sets the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Sets the state vector with a given pressure, temperature and flow speed.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	virtual void setPTU(double p, double T, double u);

	/*!
	 * \brief Sets the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Sets the state vector with a given pressure, temperature and flow speed,
	 * one set of values for each node/cell.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	void setPTU(const RowVector& p, const RowVector& T, const RowVector& u);

	/*!
	 * \brief Sets the approximate Riemann solver.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * The approximate Riemann solver is needed for computing the flux between
	 * cells.
	 *
	 * \param rs Approximate Riemann solver to be set.
	 */
	void setRiemannSolver(RiemannSolver_pt rs);

	/*!
	 * \brief	Updates the composition at each cell.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/05/2016
	 */

	virtual void UpdateComposition();

	/*!
	 * \brief Updates the flow variables with the current state vector values.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void UpdateFlowVariables();

	/*!
	 * \brief Updates R, gamma and company.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void UpdateGasProperties();

	/*!
	 * \brief Updates the state of the left and right-side variables.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Godunov's method uses the state variables at the left and right side
	 * of each interface between cells. These variables are computed using
	 * this method, and are saved in several attributes of the object for
	 * later use.
	 *
	 * Updates the state of #Frho_l, #Frho_r, #FW_l, #FW_r...
	 */
	void UpdateLeftRightVars();
};

/*!
 * \brief Sets TGodunov as the computational method for a pipe.
 *
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * 
 * After using this function, the pipe will have a TGodunov method as its
 * computational method and the method will know which pipe will have its
 * ownership. An approximate Riemann solver function is also needed for
 * the TGodunov method.
 *
 * \param pipe Pipe that will be computed using a TGodunov object.
 * \param rs Riemann solver for the TGodunov object.
 */
void use_godunov(const Pipe_ptr & pipe, const RiemannSolver_pt & rs);

#endif
