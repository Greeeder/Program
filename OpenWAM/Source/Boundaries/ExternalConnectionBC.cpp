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
 * \file ExternalConnectionBC.cpp
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

#include "ExternalConnectionBC.hpp"

TExternalConnectionBC::TExternalConnectionBC(const Pipe_ptr & pipe,
	nmPipeEnd pipe_end, const ConstantConditionsPipe_ptr & pipe_cc,
	const VirtualPipe_ptr & virtual_pipe):
	TConstantConditionsBC(pipe, pipe_end, pipe_cc, virtual_pipe) {}

void TExternalConnectionBC::LoadNewData(double *p, double *T, double *u) {
	double m = FM / FDT;
	if (std::abs(m) > 1E-8) {
		*p = FCCPipe->getPressure()(0);
		double h0 = FH / m / FDT;
		double cp = FCCPipe->getcp()(0);
		double R = FCCPipe->getR()(0);
		double A = getArea();
		double a = pow2(m * R / A);
		double b = 2 * pow2(*p) * cp;
		double c = -2 * pow2(*p) * h0;
		double d = b * b - 4. * a * c;
		*T = (-b + Sqrt(d)) / (2 * a);
		*u = m / A / (*p / R / *T);
		if (FPipeEnd_0 == nmRight) {
			*u = -*u;
		}
		*p = __units::PaToBar(*p);
	} else {
		*p = __units::PaToBar(FCCPipe->getPressure()(0));
		*T = FCCPipe->getTemperature()(0);
		*u = 0.;
	}
	if (*T < 0.) {
		*p = __units::PaToBar(FCCPipe->getPressure()(0));
		*T = FCCPipe->getTemperature()(0);
		*u = 0.;
	}
	FM = 0.;
	FFLUXU = 0.;
	FH = 0.;
	FDT = 0.;
}

void TExternalConnectionBC::LoadFluxes(double *m, double *mh0, double *mom)
{
	if (FPipeEnd_0 == nmLeft) {
		*m = FM / FDT;
		*mh0 = FH / FDT;
		*mom = FFLUXU / FDT;
	} else {
		*m = -FM / FDT;
		*mh0 = -FH / FDT;
		*mom = FFLUXU / FDT;
	}
	FM = 0.;
	FFLUXU = 0.;
	FH = 0.;
	FDT = 0.;
}

void TExternalConnectionBC::UpdateCurrentExternalProperties(double u,
	double T, double p, double t) {
	setPTU(p, T, u);
	FPExt = p;
	FTExt = T;
	FUExt = u;
	TBoundaryCondition::FCurrentTime = t;
}

void attach_to_external_connection(const Pipe_ptr pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid, double u, int ID) {
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
	ExternalConnectionBC_ptr first(new TExternalConnectionBC(pipe, pipe_end,
		constant_pipe, virtual_pipe));
	first->FID = ID;
	if(pipe_end == nmLeft) {
		dynamic_cast<TExternalConnectionBC*>(first.get())->FArea =
			pipe->getArea()(0);
		pipe->setLeftBC(move(first));
	} else {
		dynamic_cast<TExternalConnectionBC*>(first.get())->FArea =
			pipe->getArea().tail(1)(0);
		pipe->setRightBC(move(first));
	}
}
