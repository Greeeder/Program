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
 * \file RiemannSolvers.cpp
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
 * This file defines approximate Riemann solvers that can be used with a
 * Godunov finite-volume integrator.
 */

#include "RiemannSolvers.hpp"

RowArray HLL(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r) {
	auto m = U_l.rows();
	auto n = U_l.cols();
	RowVector c_l = (u_l - a_l).min(u_r - a_r);
	RowVector c_r = (u_r + a_r).max(u_l + a_l);
	RowArray F = (f_l.rowwise() * c_r - f_r.rowwise() * c_l
		+ (U_r - U_l).rowwise() * (c_l * c_r)).rowwise()
		/ (c_r - c_l);
	for (auto i = 0; i < n; i++)
	{
		if (c_l(i) > 0.)
		{
			F.col(i) = f_l.col(i);
		}
		else if (c_r(i) < 0.)
		{
			F.col(i) = f_r.col(i);
		}
	}
	return F;
}

RowArray HLLC(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r) {
	auto m = U_l.rows();
	auto n = U_l.cols();
	RowArray F(m, n);
	RowVector c_l = u_l - a_l;
	RowVector c_r = u_r + a_r;
	RowVector c_star = (p_r - p_l + rho_l * u_l * (c_l - u_l)
		- rho_r * u_r * (c_r - u_r))
		/ (rho_l * (c_l - u_l) - rho_r * (c_r - u_r));
	for (auto i = 0; i < n; i++)
	{
		if (c_l(i) > 0.)
		{
			F.col(i) = f_l.col(i);
		}
		else if (c_r(i) < 0.)
		{
			F.col(i) = f_r.col(i);
		}
		else if (c_star(i) > 0.)
		{
			double u_s = c_l(i) * rho_l(i) * (c_l(i) - u_l(i)) / (c_l(i) - c_star(i));
			F.col(i) = f_l.col(i) - c_l(i) * U_l.col(i);
			F(0, i) += u_s;
			F(1, i) += u_s * c_star(i);
			F(2, i) += u_s * (e_l(i) + (c_star(i) - u_l(i))
				* (c_star(i) + p_l(i) / rho_l(i) / (c_l(i) - u_l(i))));
		}
		else
		{
			double u_s = c_r(i) * rho_r(i) * (c_r(i) - u_r(i)) / (c_r(i) - c_star(i));
			F.col(i) = f_r.col(i) - c_r(i) * U_r.col(i);
			F(0, i) += u_s;
			F(1, i) += u_s * c_star(i);
			F(2, i) += u_s * (e_r(i) + (c_star(i) - u_r(i))
				* (c_star(i) + p_r(i) / rho_r(i) / (c_r(i) - u_r(i))));
		}
	}
	return F;
}

RowArray KT(const RowArray & U_l, const RowArray & U_r, const RowArray & f_l, const RowArray & f_r,
	const RowVector & rho_l, const RowVector & rho_r, const RowVector & p_l, const RowVector & p_r,
	const RowVector & e_l, const RowVector & e_r, const RowVector & u_l, const RowVector & u_r,
	const RowVector & a_l, const RowVector & a_r) {
	int m = U_l.rows();
	int n = U_l.cols();
	RowVector c_l = (u_l.abs() + a_l);
	RowVector c_r = (u_r.abs() + a_r);
	RowVector a_max = c_l;
	for(auto i = 0; i < c_l.cols(); i++) {
		if(c_r(i) > a_max(i)) {
			a_max(i) = c_r(i);
		}
	}
	return 0.5 * ((f_l + f_r) - (U_r - U_l).rowwise() * a_max);
}
