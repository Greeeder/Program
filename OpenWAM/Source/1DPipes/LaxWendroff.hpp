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
 * \file LaxWendroff.hpp
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
 * This file declares a Lax Wendroff finite-differences integrator for
 * one-dimensional pipes.
 */

#ifndef LaxWendroff_hpp
#define LaxWendroff_hpp

#include "TPipe.hpp"
#include "PipeMethod.hpp"
#include "Math_wam.h"

/*!
 * \brief A Lax Wendroff integrator.
 *
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * 
 * A two-steps Lax Wendroff integrator.  It is second order-accurate in both
 * time and space.
 */
class TLaxWendroff: public TPipeMethod {
  protected:
	RowVector Fhi12; ///< Internal heat transfer coefficient at half time-step.
	RowVector Frho12; ///< Density at half time-step. [kg / m ** 3]
	RowVector FRe12; ///< Reynolds number at half time-step.
	RowVector FTWPipe12; ///< Wall temperature at half time-step. [K]
	RowVector FGamma12; ///< Heat capacities ratio at half time-step.
	RowVector FR12; ///< Gas constant at half time-step. [J / (kg * K)]
	RowVector FGamma1_12; ///< Gamma - 1.
	RowVector FDerLinArea_12; ///< Area derivative at half time-step. [m]
	RowVector FArea_12; ///< Area at half cell. [m ** 2]
	RowVector FQint_12; ///< Heat transfer at half cell [W]
	RowVector FFric_12; ///< Friction coefficient at half cell [-]

	RowArray Fx1;
	RowArray Fx2;
	RowArray Fx3;
	RowArray Fx4;
	RowArray Fx1_12;
	RowArray Fx2_12;
	RowArray Fx3_12;
	RowArray Fx4_12;
	RowArray FW; ///< Flux vector;
	RowArray FUY;
	RowArray FY;
	RowArray FUY_12;
	RowArray FWY;
	RowArray FVY;
	RowArray FV1; ///< Source terms due to area.
	RowArray FV2; ///< Source terms due to friction and heat.
	RowArray FU_12; ///< State vector at half time-step.
	RowArray FW_12; ///< Flux vector at half time-step.
	RowArray FV1_12; ///< Source terms due to area at half time-step.
	RowArray FV2_12; ///< Source terms due to friction and heat at half time-step.

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
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Initialises the propagator with default values.
	 */
	TLaxWendroff();

	/*!
	 * \brief Constructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Initialises the propagator and attaches it to a pipe.
	 *
	 * \param pipe Pipe to attach to.
	 */
	TLaxWendroff(const Pipe_ptr & pipe);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TLaxWendroff();

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
	 * \brief Computes the source terms due to friction and heat transfer.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param U State vector.
	 * \param V Source terms.
	 * \param A Area. [m ** 2]
	 * \param q Heat transfer [W]
	 * \param f Friction coefficient [-]		
	 */
	void ComputeSource2(const RowArray & U, RowArray & V, const RowVector & A, const RowVector &q,
		const RowVector &f);

	/*!
	 * \brief Connects the Lax Wendroff integrator to a pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param pipe Pipe to connect to.
	 */
	virtual void Connect(const Pipe_ptr & pipe);

	/*!
	 * \brief Returns the already computed fluxes between nodes.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Fluxes between nodes.
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
	 * \brief Integrate the flow.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct.
	 */
	virtual void Solve();

	/*!
	 * \brief Integrates the central nodes of the pipe.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	void SolveCentralNodes();

	void SolveCentralNodesSpecies();

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
	 * \brief Sets the state vector in a given cell.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Sets the state vector with a given pressure, temperature and flow speed,
	 * but only for a given cell.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 * \param i Cell number.
	 */
	virtual void setPTU(double p, double T, double u, unsigned int i);

	/*!
	 * \brief Sets the updated state vector in a given cell.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Sets the updated state vector with a given pressure, temperature and flow speed,
	 * but only for a given cell.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 * \param i Cell number.
	 */
	virtual void setPTU1(double p, double T, double u, unsigned int i);

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
	 * \brief	Updates the composition (to do).
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/05/2016
	 */

	virtual void UpdateComposition(){};

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
};

/*!
 * \brief Sets TLaxWendroff as the computational method for a pipe.
 *
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * 
 * \param pipe Pipe that will be computed using a TGodunov object.
 */
void use_laxwendroff(const Pipe_ptr & pipe);

#endif
