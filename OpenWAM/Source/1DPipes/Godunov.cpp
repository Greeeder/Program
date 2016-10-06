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
 * \file Godunov.cpp
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
 * This file defines a Godunov finite-volume integrator for
 * one-dimensional pipes.
 */

#include "Godunov.hpp"
#include "BoundaryCondition.hpp"


TGodunov::TGodunov() {
	FMaxCourant = 0.9;
	FName = "Godunov";
}

TGodunov::TGodunov(const Pipe_ptr & pipe) {
	FMaxCourant = 0.9;
	FName = "Godunov";
	Connect(pipe);
	setPTU(1E5, 300, 0.);
}

TGodunov::~TGodunov() {}

void TGodunov::ComputeFlux(const RowArray & U, RowArray & W, const RowVector & Gamma, const RowVector & Gamma1) {
	/*
	 * Computes the direct flow vector.
	 * They are needed to compute the real fluxes between cells when solving
	 * the Riemann probleam at each cell interface.
	 */
	int m = U.rows();
	int n = U.cols();
	RowVector U1U0 = U.row(1) / (U.row(0) + 1E-10);

	W.row(0) = U.row(1);
	W.row(1) = U.row(2) * Gamma1 - (Gamma - 3.) * U.row(1) * U1U0 / 2.;
	W.row(2) = Gamma * U.row(2) * U1U0 - Gamma1 * U.row(1) * U1U0.square() / 2.;
	W.bottomRows(m - 3) = U.bottomRows(m - 3).rowwise() * U1U0;

	for(int i = 0; i < n; i++) {
		if(DoubEqZero(U(1, i))) {
			W(0, i) = 0.;
			W(2, i) = 0.;
			W.block(3, i, m - 3, 1).setZero();
		}
	}
}

void TGodunov::ComputeMaxTimeStep() {
	// The maximum time-step is limited due to the CFL condition.
	double c = (FPipe->getSpeedNB().abs() + FPipe->getSpeedOfSoundNB()).maxCoeff();
	FMaxTimeStep = FMaxCourant * FPipe->Fdx(0) / c;
}

void TGodunov::ComputeSource1(const RowArray & U, RowArray & V, const RowVector & A, const RowVector & Gamma1) {
	// Computes the source term due to area change in the pipe.
	auto n = A.cols() - 1;
	V.setZero();
	V.row(1) = (U.row(2) - U.row(1) * U.row(1) / U.row(0) / 2.)
		* (FPipe->FGamma - 1.) * (A.tail(n) - A.head(n));
}

void TGodunov::ComputeSource2(const RowArray & U, RowArray & V, const RowVector & D, const RowVector & f, const RowVector & q, const double dx) {
	// Fills the source terms due to friction and heat transfer.
	V.setZero();
	V.row(1) = - f * U.row(1) * U.row(1).abs() * __cons::Pi * D * dx / U.row(0);
	V.row(2) = q;
}

void TGodunov::Connect(const Pipe_ptr & pipe) {
	TPipeMethod::Connect(pipe);
	int m = pipe->FU0.rows();
	int n = pipe->FMx.cols();
	pipe->FU0.setZero(m, n);
	pipe->FU1.setZero(m, n);
	FW.setZero(m, n);
	FF.setZero(m, n + 1);
	auto s = pipe->WorkingFluid[0]->GetComposition().size();
	FFY.setZero(s, n + 1);
	FV1.setZero(m, n);
	FV2.setZero(m, n);
	pipe->Fx_sv = pipe->FMx;
	n -= 1;
	Fa_l.setZero(n);
	Fa_r.setZero(n);
	Fe_l.setZero(n);
	Fe_r.setZero(n);
	Fp_l.setZero(n);
	Fp_r.setZero(n);
	Frho_l.setZero(n);
	Frho_r.setZero(n);
	FU_l.setZero(m, n);
	FU_r.setZero(m, n);
	Fu_l.setZero(n);
	Fu_r.setZero(n);
	FW_l.setZero(m, n);
	FW_r.setZero(m, n);
}

RowArray TGodunov::getFluxes() const {
	return FF;
}

void TGodunov::setPTU(double p, double T, double u) {
	auto n_nodes = FPipe->Fx.size() - 1;
	FPipe->FU0.setZero(3, n_nodes);
	FPipe->FU1.setZero(3, n_nodes);
	double R = FPipe->WorkingFluid[0]->getR();
	double Cv = FPipe->WorkingFluid[0]->FunCv(T);
	double Cp = Cv + R;
	double g = Cp / Cv;
	double mu = FPipe->WorkingFluid[0]->FunVisc(T);
	FPipe->FR.setConstant(1, n_nodes, R);
	FPipe->FGamma.setConstant(1, n_nodes, g);
	FPipe->Fcv.setConstant(1, n_nodes, Cv);
	FPipe->Fcp.setConstant(1, n_nodes, Cp);
	RowVector rho = p / (FPipe->FR * T);
	FPipe->FViscosity.setConstant(1, n_nodes, mu);
	FPipe->FU0.row(0) = rho;
	FPipe->FU0.row(1) = rho * u;
	FPipe->FU0.row(2) = rho * FPipe->Fcv * T + rho * u * u / 2.;
	FPipe->Ftemperature.setConstant(1, n_nodes, T);
	FPipe->Fpressure.setConstant(1, n_nodes, p);
	FPipe->Fspeed.setConstant(1, n_nodes, u);
	FPipe->Fhi.setZero(1, n_nodes);
	FPipe->FRe.setZero(1, n_nodes);
	FPipe->FTWPipe.setZero(1, n_nodes);
	FPipe->FNin = n_nodes;
	UpdateFlowVariables();
}

void TGodunov::setPTU(const RowVector& p, const RowVector& T, const RowVector& u) {
	auto n_nodes = FPipe->Fx.size() - 1;
	if((p.size() == n_nodes) & (T.size() == n_nodes) & (u.size() == n_nodes)) {
		FPipe->FU0.setZero(3, n_nodes);
		FPipe->FU1.setZero(3, n_nodes);
		double R = FPipe->WorkingFluid[0]->getR();
		double Cv = FPipe->WorkingFluid[0]->FunCv(T(0));
		double Cp = Cv + R;
		double g = Cp / Cv;
		double mu = FPipe->WorkingFluid[0]->FunVisc(T(0));
		FPipe->FR.setConstant(1, n_nodes, R);
		FPipe->FGamma.setConstant(1, n_nodes, g);
		FPipe->Fcv.setConstant(1, n_nodes, Cv);
		FPipe->Fcp.setConstant(1, n_nodes, Cp);
		FPipe->FViscosity.setConstant(1, n_nodes, mu);
		FPipe->FU0.row(0) = p / FPipe->FR / T;
		FPipe->FU0.row(1) = FPipe->FU0.row(0) * u;
		FPipe->FU0.row(2) = FPipe->FU0.row(0) * (FPipe->Fcv * T + u.square() / 2.);
		FPipe->Fhi.setZero(1, n_nodes);
		FPipe->FRe.setZero(1, n_nodes);
		FPipe->FTWPipe.setZero(1, n_nodes);
		FPipe->FNin = n_nodes;
		UpdateFlowVariables(); //TODO
	}
}

void TGodunov::setRiemannSolver(RiemannSolver_pt rs) {
	Frs = rs;
}


void TGodunov::IntegrateWithoutUpdating() {
	auto n = FPipe->FU0.cols();
	double dt = FPipe->getTimeStep();
	FPipe->FIsIntegrated = false;
	/*
	 * The system is solved using Godunov's reconstruction
	 * and a single explicit Euler step.
	 * First, the direct fluxes have to be computed.
	 */
	ComputeFlux(FPipe->FU0, FW, FPipe->FGamma, FPipe->FGamma1);
	// Vectors with the state at both sides of each interface are updated.
	UpdateLeftRightVars();
	/*
	 * Solves the Riemann problem at each cell interface, generating FF.
	 * FF contains the fluxes at the cell interfaces.
	 */
	SolveCentralCells();
	// Computes some extra source terms.
	ComputeSource1(FPipe->FU0, FV1, FPipe->FArea, FPipe->FGamma1);
	ComputeSource2(FPipe->FU0, FV2, FPipe->FDcell, FPipe->FFric, FPipe->FQint, FPipe->FXref);
	FV1 += FPipe->FSource->ComputeSource(FPipe->FCurrentTime, dt);
	// Fill the value of the fluxes at the left and right ends of the duct.
	FF.col(0) = FPipe->FLeftBC->Flux(FPipe->FCurrentTime, dt);
	FF.col(n) = FPipe->FRightBC->Flux(FPipe->FCurrentTime, dt);
	/*
	 * With all the flux and source terms, the updated state vector can be computed.
	 * It is computed using an explicit Euler step.
	 */
	for(auto i = 0; i < 3; i++) {
		FPipe->FU1.row(i) = (FF.row(i).head(n) * FPipe->FArea.head(n)
							 - FF.row(i).tail(n) * FPipe->FArea.tail(n)
							 + FV1.row(i) + FV2.row(i)) / FPipe->FVolume * dt
							+ FPipe->FU0.row(i);
	}
}

void TGodunov::UpdateComposition(){
	// Updates the composition of each cell.
	auto n = FPipe->WorkingFluid[0]->GetComposition().size();

	if (n > 1){
		auto m = n - 1;

		/*
		 * The composition is convected, so the mass flow have to be checked
		 * at each cell interface.
		 */
		for (auto i = 1; i < FPipe->FNin; i++){
			if (FF(0, i) > 0){
				FFY.col(i).head(m) = FF(0, i) * FPipe->FArea(i) * FPipe->WorkingFluid[i - 1]->GetComposition().head(m);
			}
			else{
				FFY.col(i).head(m) = FF(0, i) * FPipe->FArea(i) * FPipe->WorkingFluid[i]->GetComposition().head(m);
			}
		}
		// The convection is also checked at the boundary conditions.
		if (FF.row(0).head(1)(0) < 0)
			FFY.col(0).head(m) = FF(0, 0) * FPipe->FArea(0) * FPipe->WorkingFluid[0]->GetComposition().head(m);
		else
			FFY.col(0).head(m) = FF(0, 0) * FPipe->FArea(0) * FPipe->FLeftBC->FluidBC()->GetComposition().head(m);

		if (FF.row(0).tail(1)(0) > 0.)
			FFY.col(FPipe->FNin).head(m) = FF.row(0).tail(1)(0) * FPipe->FArea.tail(1)(0) * FPipe->WorkingFluid.back()->GetComposition().head(m);
		else
			FFY.col(FPipe->FNin).head(m) = FF.row(0).tail(1)(0) * FPipe->FArea.tail(1)(0) * FPipe->FRightBC->FluidBC()->GetComposition().head(m);

		// With the composition fluxes, the new composition can be computed.
		RowVector Ynew = RowVector::Zero(n);
		RowVector Yflux = RowVector::Zero(n);
		for (auto i = 0; i < FPipe->FNin; i++){
			Ynew.head(m) = FPipe->FU0(0, i) * FPipe->WorkingFluid[i]->GetComposition().head(m);
			Yflux.head(m) = -(FFY.col(i + 1).head(m) - FFY.col(i).head(m)) / FPipe->FVolume(i) * FPipe->getTimeStep();
			Ynew.head(m) = (Ynew.head(m) + Yflux.head(m)) / FPipe->FU1(0, i);
			Ynew.tail(1) = 1 - Ynew.head(m).sum();
			FPipe->WorkingFluid[i]->SetComposition(Ynew);
		}
	}

}

void TGodunov::UpdateFlowVariables() {
	/*
	 * Updates variables such as the density or the temperature using the
	 * state vector.
	 */
	FPipe->Frho = FPipe->FU0.row(0);
	FPipe->Fspeed = FPipe->FU0.row(1) / FPipe->FU0.row(0);
// 	for (auto i = 0; i < FPipe->Frho.size(); i++){
// 		FPipe->Ftemperature(i) = FPipe->WorkingFluid[i]->getTemperature(FPipe->FU0(2, i) / FPipe->Frho(i) - pow2(FPipe->Fspeed(i)) / 2.);
// 		FPipe->Fpressure(i) = FPipe->Frho(i) * FPipe->WorkingFluid[i]->FunR() * FPipe->Ftemperature(i);
// 	}
	FPipe->Ftemperature =
		(FPipe->FU0.row(2) / FPipe->Frho - pow2(FPipe->Fspeed) / 2.) / FPipe->Fcv;
	FPipe->Fpressure = FPipe->Frho * FPipe->FR * FPipe->Ftemperature;
	//double T = FPipe->WorkingFluid[0]->getTemperature(FPipe->FU0(2, 0) / FPipe->Frho(0));
	//FPipe->Fpressure = (FPipe->FU0.row(2) - FPipe->FU0.row(1) * FPipe->Fspeed / 2.) * (FPipe->FGamma - 1.);
	//FPipe->Ftemperature = FPipe->Fpressure / FPipe->Frho / FPipe->FR;
	UpdateGasProperties();

	FPipe->Fa = (FPipe->FGamma * FPipe->FR * FPipe->Ftemperature).sqrt();
	//for(auto i = 0; i < FPipe->FA_A.size(); i++) {
	//	FPipe->FA_A(i) = (FPipe->Fa(i) / __cons::ARef) / pow(FPipe->Fpressure(i) / __cons::PRef / 1E5,
	//					 FPipe->FGamma1(i) / 2. / FPipe->FGamma(i));
	//}
	FPipe->Fbeta = (FPipe->Fa - (FPipe->FGamma - 1.) / 2. * FPipe->Fspeed) / __cons::ARef;
	FPipe->Flambda = (FPipe->Fa + (FPipe->FGamma - 1.) / 2. * FPipe->Fspeed) / __cons::ARef;
	
	FPipe->FTotalTemperature = FPipe->Ftemperature + 0.5 * pow2(FPipe->Fspeed) / FPipe->Fcp;
	FPipe->FTotalPressure = FPipe->Fpressure
		* pow(FPipe->FTotalTemperature / FPipe->Ftemperature, FPipe->FGamma / FPipe->FGamma1);
	FPipe->FRe = FPipe->Frho * FPipe->Fspeed.abs() * FPipe->FDcell / FPipe->FViscosity;
	ComputeMaxTimeStep();
}

void TGodunov::UpdateLeftRightVars() {
	// Update some vectors with the state at both sides of each interface.
	auto n = FPipe->FU0.cols() - 1;
	Fa_l = FPipe->Fa.tail(n);
	Fa_r = FPipe->Fa.head(n);
	Fe_l = FPipe->FU0.row(2).tail(n) / FPipe->FU0.row(0).tail(n);
	Fe_r = FPipe->FU0.row(2).head(n) / FPipe->FU0.row(0).head(n);
	Fp_l = FPipe->Fpressure.tail(n);
	Fp_r = FPipe->Fpressure.head(n);
	Frho_l = FPipe->Frho.tail(n);
	Frho_r = FPipe->Frho.head(n);
	Fu_l = FPipe->Fspeed.tail(n);
	Fu_r = FPipe->Fspeed.head(n);
	FU_l = FPipe->FU0.rightCols(n);
	FU_r = FPipe->FU0.leftCols(n);
	FW_l = FW.rightCols(n);
	FW_r = FW.leftCols(n);
}

void TGodunov::UpdateGasProperties() {
	/* TODO */
	for (auto i = 0; i < FPipe->WorkingFluid.size(); i++){
		FPipe->FR(i) = FPipe->WorkingFluid[i]->getR();
		FPipe->Fcv(i) = FPipe->WorkingFluid[i]->FunCv(FPipe->Ftemperature(i));
		FPipe->FViscosity(i) = FPipe->WorkingFluid[i]->FunVisc(FPipe->Ftemperature(i));
	}

	FPipe->Fcp = FPipe->Fcv + FPipe->FR;
	FPipe->FGamma = FPipe->Fcp / FPipe->Fcv;
	FPipe->FGamma1 = FPipe->FGamma - 1.;
}

void TGodunov::Solve() {
	IntegrateWithoutUpdating();
	FPipe->UpdateStateVector();
}

void TGodunov::SolveCentralCells() {
	/*
	 * Solves the Riemann problem at each cell interface, generating FF.
	 * FF contains the fluxes at the cell interfaces.
	 */
	auto n = FPipe->FU0.cols() - 1;
	FF.middleCols(1, n) = Frs(FU_r, FU_l, FW_r, FW_l, Frho_r, Frho_l, Fp_r, Fp_l, Fe_r, Fe_l, Fu_r, Fu_l, Fa_r, Fa_l);
}

void use_godunov(const Pipe_ptr & pipe, const RiemannSolver_pt & rs) {
	unique_ptr<TGodunov> method(new TGodunov());
	method->Connect(pipe);
	method->setRiemannSolver(rs);
	pipe->setMethod(move(method));
}

