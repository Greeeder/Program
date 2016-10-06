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
 * \file RiemannSolvers.hpp
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
 * This file declares approximate Riemann solvers that can be used with a
 * Godunov finite-volume integrator.
 */

#ifndef RiemannSolvers_hpp
#define RiemannSolvers_hpp

#include "Math_wam.h"

/*!
 * \brief HLL Riemann solver.
 *
 * An approximate solver by Harten, Lax and van Leer.
 *
 * \param U_l Left side extrapolated value for the state vector.
 * \param U_r Right side extrapolated value for the state vector.
 * \param f_l Left side extrapolated value for the flux vector.
 * \param f_r Right side extrapolated value for the flux vector.
 * \param rho_l Left side extrapolated value for the density. [kg / m ** 3]
 * \param rho_r Right side extrapolated value for the density. [kg / m ** 3]
 * \param p_l Left side extrapolated value for the pressure. [Pa]
 * \param p_r Right side extrapolated value for the pressure. [Pa]
 * \param e_l Left side extrapolated value for the specific internal energy.
 * \param e_r Right side extrapolated value for the specific internal energy.
 * \param c_l Left side extrapolated value for the flow speed. [m / s]
 * \param c_r Right side extrapolated value for the flow speed. [m / s]
 * \param a_l Left side extrapolated value for the speed of sound. [m / s]
 * \param a_r Right side extrapolated value for the speed of sound. [m / s]
 * \return An approximated solution to the Riemann problem.
 */
RowArray HLL(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r);

/*!
 * \brief HLLC Riemann solver.
 *
 * An approximate solver by Harten, Lax and van Leer.
 *
 * \param U_l Left side extrapolated value for the state vector.
 * \param U_r Right side extrapolated value for the state vector.
 * \param f_l Left side extrapolated value for the flux vector.
 * \param f_r Right side extrapolated value for the flux vector.
 * \param rho_l Left side extrapolated value for the density. [kg / m ** 3]
 * \param rho_r Right side extrapolated value for the density. [kg / m ** 3]
 * \param p_l Left side extrapolated value for the pressure. [Pa]
 * \param p_r Right side extrapolated value for the pressure. [Pa]
 * \param e_l Left side extrapolated value for the specific internal energy.
 * \param e_r Right side extrapolated value for the specific internal energy.
 * \param c_l Left side extrapolated value for the flow speed. [m / s]
 * \param c_r Right side extrapolated value for the flow speed. [m / s]
 * \param a_l Left side extrapolated value for the speed of sound. [m / s]
 * \param a_r Right side extrapolated value for the speed of sound. [m / s]
 * \return An approximated solution to the Riemann problem.
 */
RowArray HLLC(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r);

/*!
 * \brief KT Riemann solver.
 *
 * The Kurganov and Tadmor central scheme by Kurganov and Tadmor.  Not exactly
 * an approximate Riemann solver, but it can be used instead.
 *
 * Approximates the intercell fluxes by the Rusanov flux.
 *
 * \param U_l Left side extrapolated value for the state vector.
 * \param U_r Right side extrapolated value for the state vector.
 * \param f_l Left side extrapolated value for the flux vector.
 * \param f_r Right side extrapolated value for the flux vector.
 * \param rho_l Left side extrapolated value for the density. [kg / m ** 3]
 * \param rho_r Right side extrapolated value for the density. [kg / m ** 3]
 * \param p_l Left side extrapolated value for the pressure. [Pa]
 * \param p_r Right side extrapolated value for the pressure. [Pa]
 * \param e_l Left side extrapolated value for the specific internal energy.
 * \param e_r Right side extrapolated value for the specific internal energy.
 * \param c_l Left side extrapolated value for the flow speed. [m / s]
 * \param c_r Right side extrapolated value for the flow speed. [m / s]
 * \param a_l Left side extrapolated value for the speed of sound. [m / s]
 * \param a_r Right side extrapolated value for the speed of sound. [m / s]
 * \return An approximated solution to the Riemann problem.
 */
RowArray KT(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r);

/*!
 * \brief Pointer to an approximate Riemann solver.
 *
 * Pointer to an approximate Riemann solver, useful for the TGodunov.
 *
 * \param U_l Left side extrapolated value for the state vector.
 * \param U_r Right side extrapolated value for the state vector.
 * \param f_l Left side extrapolated value for the flux vector.
 * \param f_r Right side extrapolated value for the flux vector.
 * \param rho_l Left side extrapolated value for the density. [kg / m ** 3]
 * \param rho_r Right side extrapolated value for the density. [kg / m ** 3]
 * \param p_l Left side extrapolated value for the pressure. [Pa]
 * \param p_r Right side extrapolated value for the pressure. [Pa]
 * \param e_l Left side extrapolated value for the specific internal energy.
 * \param e_r Right side extrapolated value for the specific internal energy.
 * \param c_l Left side extrapolated value for the flow speed. [m / s]
 * \param c_r Right side extrapolated value for the flow speed. [m / s]
 * \param a_l Left side extrapolated value for the speed of sound. [m / s]
 * \param a_r Right side extrapolated value for the speed of sound. [m / s]
 * \return An approximated solution to the Riemann problem.
 */
typedef RowArray(*RiemannSolver_pt)(const RowArray &, const RowArray &,
	const RowArray &, const RowArray &,
	const RowVector &, const RowVector &, const RowVector &, const RowVector &,
	const RowVector &, const RowVector &, const RowVector &, const RowVector &,
	const RowVector &, const RowVector &);

#endif
