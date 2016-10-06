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
 * \file ConstantConditionsBC.cpp
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
 * This file defines a boundary condition with constant pressure and
 * temperature, including some helper functions.
 */

#include "ConstantConditionsBC.hpp"

TConstantConditionsBC::TConstantConditionsBC() {}

TConstantConditionsBC::TConstantConditionsBC(const Pipe_ptr & pipe,
	nmPipeEnd pipe_end, const ConstantConditionsPipe_ptr & pipe_cc,
	const VirtualPipe_ptr & virtual_pipe):
	TPipeConnection(pipe, pipe_end, virtual_pipe, 0) {

	FCCPipe = pipe_cc;
	FM = 0.;
	FFLUXU = 0.;
	FH = 0.;
	FDT = 0.;
}

ColVector TConstantConditionsBC::Flux(double t, double dt) {
	FCCPipe->setTimeStep(dt);
	FCCPipe->Solve(t);
	TPipeConnection::Flux(t, dt);
	FM += getMassFlow() * dt;
	FFLUXU += FFlux_0(1) * getArea() * dt;
	FH += getEnthalpyFlow() * dt;
	FDT += dt;
	return FFlux_0;
}

void TConstantConditionsBC::setPTU(double p, double T, double u) {
	FCCPipe->setPTU(p, T, u);
}

void attach_to_constant_BC(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, double u, TFluid_ptr fluid) {
	double D = 0.;
	if (pipe_end == nmLeft) {
		D = pipe->getD()(0);
	} else {
		D = pipe->getD().tail(1)(0);
	}
	ConstantConditionsPipe_ptr constant_pipe =
		create_constant_conditions_pipe(p, T, u, D, 1E10, 1, fluid);

	nmPipeEnd vp_end = nmRight;
	VirtualPipe_ptr virtual_pipe = make_shared<TVirtualPipe>(pipe, pipe_end,
		constant_pipe, vp_end);
	unique_ptr<TGodunov> method(new TGodunov());
	method->Connect(virtual_pipe);
	method->setRiemannSolver(&KT);
	virtual_pipe->FMethod = std::move(method);
	virtual_pipe->FMethod->setPTU(1E5, 300, 0);
	ConstantConditionsBC_ptr first(new TConstantConditionsBC(pipe, pipe_end,
		constant_pipe, virtual_pipe));
	if(pipe_end == nmLeft) {
		dynamic_cast<TConstantConditionsBC*>(first.get())->FArea =
			pipe->getArea()(0);
		pipe->setLeftBC(move(first));
	} else {
		dynamic_cast<TConstantConditionsBC*>(first.get())->FArea =
			pipe->getArea().tail(1)(0);
		pipe->setRightBC(move(first));
	}
}
