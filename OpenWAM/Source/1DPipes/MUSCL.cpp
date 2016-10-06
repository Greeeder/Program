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
 * @file MUSCL.cpp
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
 * This file defines a MUSCL finite-volume integrator for
 * one-dimensional pipes.
 */

#include "MUSCL.hpp"
#include "BoundaryCondition.hpp"


TMUSCL::TMUSCL() {
	FMaxCourant = 0.9;
	FName = "MUSCL";
}

TMUSCL::TMUSCL(const Pipe_ptr & pipe) {
	FMaxCourant = 0.9;
	FName = "MUSCL";
	Connect(pipe);
	setPTU(1E5, 300, 0.);
}

TMUSCL::~TMUSCL() {}

void TMUSCL::Connect(const Pipe_ptr & pipe) {
	TGodunov::Connect(pipe);
	int m = pipe->FU0.rows();
	int n = pipe->FMx.cols();
	FU12.setZero(m, n);
	FF1.setZero(m, n + 1);
	FF2.setZero(m, n + 1);
	Fphi_l.setZero(m, n - 2);
	Fphi_r.setZero(m, n - 2);
}


void TMUSCL::setLimiter(Limiter_pt limiter) {
	Flim = limiter;
}


void TMUSCL::ComputeExtrapolatedValues(const RowArray & U) {
	int n = U.cols();
	FU_r = U.leftCols(n - 1);
	FU_l = U.rightCols(n - 1);
	FU_r.rightCols(n - 2) += Fphi_r / 2. * (U.middleCols(1, n - 2)
											- U.leftCols(n - 2));
	FU_l.leftCols(n - 2) -= Fphi_l / 2. * (U.rightCols(n - 2)
										   - U.middleCols(1, n - 2));
	Frho_r = FU_r.row(0);
	Frho_l = FU_l.row(0);
	Fu_r = FU_r.row(1) / Frho_r;
	Fu_l = FU_l.row(1) / Frho_l;
	Fp_r = (FU_r.row(2) - FU_r.row(1) * Fu_r / 2.) * (FPipe->FGamma.leftCols(n - 1)
			- 1.);
	Fp_l = (FU_l.row(2) - FU_l.row(1) * Fu_l / 2.) * (FPipe->FGamma.rightCols(n - 1)
			- 1.);
	Fe_r = FU_r.row(2) / FU_r.row(0);
	Fe_l = FU_l.row(2) / FU_l.row(0);
	Fa_r = (FPipe->FGamma.leftCols(n - 1) * Fp_r / Frho_r).sqrt();
	Fa_l = (FPipe->FGamma.rightCols(n - 1) * Fp_l / Frho_l).sqrt();
}

void TMUSCL::ComputeSlopeLimiter(const RowArray & U) {
	int n = U.cols() - 2;
	Fratio = (U.rightCols(n)
			  - U.middleCols(1, n) + 1E-100)
			 / (U.middleCols(1, n)
				- U.leftCols(n) + 1E-100);
	Flim(Fratio, &Fphi_r);
	Flim(1. / Fratio, &Fphi_l);
}

void TMUSCL::EulerStep(const RowArray& U, const RowArray& V, double t, double dt) {
	FPipe->FIsIntegrated = false;
	auto n = U.cols();
	// Compute the extrapolated values of the state vector at the cell boundaries.
	ComputeSlopeLimiter(U);
	ComputeExtrapolatedValues(U);
	// Compute the fluxes of the left and right extrapolated states.
	ComputeFlux(FU_r, FW_r, FPipe->FGamma.leftCols(n - 1),
				FPipe->FGamma1.leftCols(n - 1));
	ComputeFlux(FU_l, FW_l, FPipe->FGamma.rightCols(n - 1),
				FPipe->FGamma1.rightCols(n - 1));
	/*
	 * Solves the Riemann problem at each cell interface, generating FF.
	 * FF contains the fluxes at the cell interfaces.
	 */
	SolveCentralCells();
	// Computes some extra source terms.
	ComputeSource1(U, FV1, FPipe->FArea, FPipe->FGamma1);
	ComputeSource2(U, FV2, FPipe->FDcell, FPipe->FFric, FPipe->FQint, FPipe->FXref);
	// Fill the value of the fluxes at the left and right ends of the duct.
	FF.col(0) = FPipe->FLeftBC->Flux(t, dt);
	FF.col(n) = FPipe->FRightBC->Flux(t, dt);
	/*
	 * With all the flux and source terms, the updated state vector can be computed.
	 * It is computed using an explicit Euler step.
	 */
	for (auto i = 0; i < 3; i++)
	{
		FU12.row(i) = (FF.row(i).head(n) * FPipe->FArea.head(n)
			- FF.row(i).tail(n) * FPipe->FArea.tail(n)
			+ FV1.row(i) + FV2.row(i) + V.row(i)) / FPipe->FVolume * dt
			+ U.row(i);
	}
}

void TMUSCL::IntegrateWithoutUpdating() {
	FPipe->FIsIntegrated = false;
	auto dt = FPipe->getTimeStep();
	// Using Heun's method, so will have to do two EulerStep calls.
	// FV3 is computed only once and kept constant for both steps.
	FV3 = FPipe->FSource->ComputeSource(FPipe->FCurrentTime, dt);
	// First EulerStep call, to get a prediction of the state vector update.
	EulerStep(FPipe->FU0, FV3, FPipe->getCurrentTime(), dt);
	FF1 = FF;
	// Second EulerStep call, to get a different prediction.
	EulerStep(FU12, FV3, FPipe->getCurrentTime(), dt);
	FPipe->FU1 = (FU12 + FPipe->FU0) / 2.;
	FF2 = FF;
	// The inter-cell fluxes are the average of FF1 and FF2.
	FF = (FF1 + FF2) / 2.;
}

void TMUSCL::Solve() {
	IntegrateWithoutUpdating();
	FPipe->UpdateStateVector();
}

void use_muscl(const Pipe_ptr & pipe, const RiemannSolver_pt & rs,
	const Limiter_pt & limiter) {
	unique_ptr<TMUSCL> method(new TMUSCL());
	method->Connect(pipe);
	method->setRiemannSolver(rs);
	method->setLimiter(limiter);
	pipe->setMethod(move(method));
}
