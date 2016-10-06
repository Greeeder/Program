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
 * @file TConstantConditionsPipe.cpp
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
 * This file defines a one-dimensional pipe with a constant state vector.
 */

#include "TConstantConditionsPipe.hpp"

double TConstantConditionsPipe::getMaxTimeStep() {
    return 1E10;
}

void TConstantConditionsPipe::Solve(double t) {
    FCurrentTime = t;
}

void TConstantConditionsPipe::Solve() {
    FCurrentTime += FTimeStep;
}

ConstantConditionsPipe_ptr create_constant_conditions_pipe(
	double p, double T, double u, double D, double l, unsigned int size, TFluid_ptr fluid) {
	RowVector x(2);
	RowVector D_vector(2);
	ConstantConditionsPipe_ptr pipe = std::make_shared<TConstantConditionsPipe> ();
	double dx = 0.;
	x(0) = 0.;
	x(1) = l;
	D_vector.setConstant(D);
	pipe->WorkingFluid.resize(1);
	pipe->WorkingFluid[0] = make_shared<TFluid>(fluid.get());
	dx = (x(1) - x(0)) / size;
	pipe->setGeometry(x, dx, D_vector);
	use_godunov(pipe, &KT);
	pipe->setPTU(p, T, u);
	return pipe;
}
