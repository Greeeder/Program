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
 * @file LaxWendroff.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
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
 * This file defines a Lax Wendroff finite-differences integrator for
 * one-dimensional pipes.
 */

#include "LaxWendroff.hpp"
#include "BoundaryCondition.hpp"

TLaxWendroff::TLaxWendroff() {
	FMaxCourant = 1.;
	FName = "Lax-Wendroff";
}

TLaxWendroff::TLaxWendroff(const Pipe_ptr & pipe) {
	FMaxCourant = 1.;
	FName = "Lax-Wendroff";
	Connect(pipe);
	setPTU(1E5, 300, 0.);
}

TLaxWendroff::~TLaxWendroff() {}

void TLaxWendroff::ComputeFlux(const RowArray & U, RowArray & W, const RowVector & Gamma, const RowVector & Gamma1) {
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

void TLaxWendroff::ComputeMaxTimeStep() {
	double c = (FPipe->getSpeedNB().abs() + FPipe->getSpeedOfSoundNB()).maxCoeff();
	FMaxTimeStep = FMaxCourant * FPipe->Fdx(0) / c;
}

void TLaxWendroff::ComputeSource2(const RowArray & U, RowArray & V, const RowVector & A, const RowVector &q,
	const RowVector &f) {

	if(DoubEqZero(FPipe->FCoefAjusFric)){
		V.row(1).setZero();
	}
	else{
		V.setZero();
		V.row(1) = U.row(0) * f* (U.row(1) / U.row(0)).pow(2) * U.row(1).cwiseSign() / (A * __cons::Pi).sqrt()
			* FPipe->FCoefAjusFric;
	}
	if (DoubEqZero(FPipe->FCoefAjusTC)) {  /* q=0 */
		V.row(2).setZero();
	}
	else {
		V.row(2) = q / FPipe->FXref;
	}


}

void TLaxWendroff::ComputeSource1(const RowArray & U, RowArray & V, const RowVector & A, const RowVector & Gamma1) {
	V.setZero();
	V.row(1) = -(U.row(2) - U.row(1).square() / (U.row(0) + 1E-10) / 2.) * Gamma1 / A;
	for(int i = 0; i < U.cols(); i++) {
		if(DoubEqZero(U(1, i))) {
			V(1, i) = -U(2, i) * Gamma1(i) / A(i);
		}
	}
}

void TLaxWendroff::Connect(const Pipe_ptr & pipe) {
	TPipeMethod::Connect(pipe);
	int m = pipe->FU0.rows();
	int n = pipe->FU0.cols();
	FW.setZero(m, n);
	FV1.setZero(m, n);
	FV2.setZero(m, n);
	n -= 1;
	pipe->Fx_sv = pipe->Fx;
	FArea_12 = (pipe->FArea.rightCols(n) + pipe->FArea.leftCols(n)) / 2.;
	FDerLinArea_12 = (FArea_12.tail(n - 1) - FArea_12.head(n - 1)) / pipe->FXref;
	Fhi12.setZero(n);
	Frho12.setZero(n);
	FRe12.setZero(n);
	FTWPipe12.setZero(n);
	FGamma12.setZero(n);
	FR12.setZero(n);
	FGamma1_12.setZero(n);
	Fx1.setZero(m, n);
	Fx2.setZero(m, n);
	Fx3.setZero(m, n);
	Fx4.setZero(m, n);
	FU_12.setZero(m, n);
	FW_12.setZero(m, n);
	FV1_12.setZero(m, n);
	FV2_12.setZero(m, n);
	n -= 1;
	Fx1_12.setZero(m, n);
	Fx2_12.setZero(m, n);
	Fx3_12.setZero(m, n);
	Fx4_12.setZero(m, n);
}

RowArray TLaxWendroff::getFluxes() const {
	return FW;
}

void TLaxWendroff::IntegrateWithoutUpdating() {
	double dt = FPipe->getTimeStep();
	FPipe->FIsIntegrated = false;
	auto n = FPipe->getNin() - 1;
	SolveCentralNodes();
	ComputePipeEndCharacteristics(dt);
	FPipe->FLeftBC->computePTUFromCharacteristics(FPipe->FCurrentTime, dt);
	FPipe->FRightBC->computePTUFromCharacteristics(FPipe->FCurrentTime, dt);
	setPTU1(FPipe->FLeftBC->getPressure(), FPipe->FLeftBC->getTemperature(), FPipe->FLeftBC->getSpeed(), 0);
	setPTU1(FPipe->FRightBC->getPressure(), FPipe->FRightBC->getTemperature(), FPipe->FRightBC->getSpeed(), n);
}

void TLaxWendroff::Solve() {
	IntegrateWithoutUpdating();
	FPipe->UpdateStateVector();
}

void TLaxWendroff::SolveCentralNodes() {
	int m = FPipe->FU0.rows();
	int n = FPipe->FU0.cols() - 1;
	int p = FUY.rows();
	double dt = FPipe->getTimeStep();
	double dtdx = dt / FPipe->FXref;
	double dt2 = dt / 2.;

	ComputeFlux(FPipe->FU0, FW, FPipe->FGamma, FPipe->FGamma1);
	ComputeSource1(FPipe->FU0, FV1, FPipe->FArea, FPipe->FGamma1);
	ComputeSource2(FPipe->FU0, FV2, FPipe->FArea, FPipe->FQint, FPipe->FFric);

	double x1, x2, x3, x4;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			x1 = FPipe->FU0(i, j) + FPipe->FU0(i, j + 1);
			x2 = -dtdx * (FW(i, j + 1) - FW(i, j));
			x3 = (FV1(i, j + 1) + FV1(i, j)) * FPipe->FDerLinArea(j) * (-dt2);
			x4 = -dt2 * (FV2(i, j + 1) + FV2(i, j));
			FU_12(i, j) = (x1 + x2 + x3 + x4) / 2.;
		}
	}
	//for(int i = 0; i < n; i++) {
	//	FTWPipe12(i) = (FPipe->FTWPipe(0, i) + FPipe->FTWPipe(0, i + 1)) / 2.;
	//}
	FQint_12 = (FPipe->FQint.head(n) + FPipe->FQint.tail(n)) / 2.;
	FFric_12 = (FPipe->FFric.head(n) + FPipe->FFric.tail(n)) / 2.;
	FGamma12 = (FPipe->FGamma.head(n) + FPipe->FGamma.tail(n)) / 2.;
	FGamma1_12 = (FPipe->FGamma1.head(n) + FPipe->FGamma1.tail(n)) / 2.;

	ComputeFlux(FU_12, FW_12, FGamma12, FGamma1_12);
	ComputeSource1(FU_12, FV1_12, FArea_12, FGamma1_12);
	ComputeSource2(FU_12, FV2_12, FArea_12, FQint_12, FFric_12);

	for(int i = 0; i < m; i++) {
		for(int j = 1; j < n; j++) {
			x1 = FPipe->FU0(i, j);
			x2 = -dtdx * (FW_12(i, j) - FW_12(i, j - 1));
			x3 = (FV1_12(i, j) + FV1_12(i, j - 1)) * FDerLinArea_12(j - 1) * (-dt2);
			x4 = -dt2 * (FV2_12(i, j) + FV2_12(i, j - 1));
			FPipe->FU1(i, j) = (x1 + x2 + x3 + x4);
		}
	}

	/*! Solve species tracking */
	if (p > 1){

		/*! Fist step */
		FWY = FPipe->FU0.row(1) / FPipe->FU0.row(0) * FUY;

		for (int i = 0; i < p - 1; i++){
			FUY_12.row(i) = (FUY.row(i).head(n) + FUY.row(i).tail(n) - dtdx * (FWY.row(i).tail(n)
				 - FWY.row(i).head(n))) / 2.;
		}

		FWY = FU_12.row(1) / FU_12.row(0) * FUY_12;

		/*! Second step */
		for (int i = 0; i < p - 1; i++){
			FUY.row(i).segment(1, n - 1) = FUY.row(i).segment(1, n - 1) - dtdx * (FWY.row(i).tail(n - 1)
				 - FWY.row(i).head(n - 1));
			FY.row(i).segment(1, n - 1) = FUY.row(i).segment(1, n - 1) / FPipe->FU1.row(0).segment(1, n - 1);
		}

		FY.row(p - 1) = 1 - FUY.topRows(p - 1).colwise().sum();
		FUY.row(p - 1) = FY.row(p - 1) * FPipe->FU1.row(0);
	}
}

void TLaxWendroff::SolveCentralNodesSpecies(){

//	for (int i = 1; i < FPipe->WorkingFluid.size()-1; i++){
//		for (int j = 0; j < FPipe->WorkingFluid[i]->GetComposition().size(); j++){
//
//		}
//	}
}

void TLaxWendroff::setPTU(double p, double T, double u) {
	auto n_nodes = FPipe->Fx.size();
	FPipe->FU0.setZero(3, n_nodes);
	FPipe->FU1.setZero(3, n_nodes);
	FPipe->FR.setConstant(1, n_nodes, 287.);
	FPipe->FGamma.setConstant(1, n_nodes, 1.4);
	FPipe->Fcv.setConstant(1, n_nodes, 287. / (1.4 - 1.));
	FPipe->Fcp.setConstant(1, n_nodes, 287. * 1.4 / (1.4 - 1.));
	RowVector rhoA = p / (FPipe->FR * T) * FPipe->FArea;
	FPipe->FGamma.setConstant(1.4);
	FPipe->Fcp = FPipe->FR * FPipe->FGamma / (FPipe->FGamma - 1.);
	FPipe->Fcv = FPipe->Fcp / FPipe->FGamma;
	FPipe->FU0.row(0) = rhoA;
	FPipe->FU0.row(1) = rhoA * u;
	FPipe->FU0.row(2) = rhoA * FPipe->Fcv * T + rhoA * u * u / 2.;
	FPipe->Fhi.setZero(1, n_nodes);
	FPipe->FRe.setZero(1, n_nodes);
	FPipe->FTWPipe.setZero(1, n_nodes);
	FPipe->FNin = n_nodes;
	UpdateFlowVariables();
}

void TLaxWendroff::setPTU(double p, double T, double u, unsigned int i) {
	auto n_nodes = FPipe->Fx.size();
	double rhoA = p / (FPipe->FR(i) * T) * FPipe->FArea(i);
	/* TODO FGamma(T)! */
	FPipe->FU0(0, i) = rhoA;
	FPipe->FU0(1, i) = rhoA * u;
	FPipe->FU0(2, i) = rhoA * FPipe->Fcv(i) * T + rhoA * u * u / 2.;
// 	UpdateFlowVariables();
}

void TLaxWendroff::setPTU1(double p, double T, double u, unsigned int i) {
	auto n_nodes = FPipe->Fx.size();
	double rhoA = p / (FPipe->FR(i) * T) * FPipe->FArea(i);
	/* TODO FGamma(T)! */
	FPipe->FU1(0, i) = rhoA;
	FPipe->FU1(1, i) = rhoA * u;
	FPipe->FU1(2, i) = rhoA * FPipe->Fcv(i) * T + rhoA * u * u / 2.;
// 	UpdateFlowVariables();
}


void TLaxWendroff::setPTU(const RowVector& p, const RowVector& T, const RowVector& u) {
	auto n_nodes = FPipe->Fx.size();
	if((p.size() == n_nodes) & (T.size() == n_nodes) & (u.size() == n_nodes)) {
		FPipe->FU0.setZero(3, n_nodes);
		FPipe->FU1.setZero(3, n_nodes);
		FPipe->FR.setConstant(1, n_nodes, 287.);
		FPipe->FGamma.setConstant(1, n_nodes, 1.4);
		FPipe->Fcv.setConstant(1, n_nodes, 287. / (1.4 - 1.));
		FPipe->Fcp.setConstant(1, n_nodes, 287. * 1.4 / (1.4 - 1.));
		FPipe->FU0.row(0) = p / FPipe->FR / T * FPipe->FArea;
		FPipe->FU0.row(1) = FPipe->FU0.row(0) * u;
		FPipe->FU0.row(2) = FPipe->FU0.row(0) * (FPipe->Fcv * T + u.square() / 2.);
		FPipe->Fhi.setZero(1, n_nodes);
		FPipe->FRe.setZero(1, n_nodes);
		FPipe->FTWPipe.setZero(1, n_nodes);
		FPipe->FNin = n_nodes;
		UpdateFlowVariables(); //TODO
	}
}

void TLaxWendroff::UpdateFlowVariables() {
	/* TODO: BCs */
	FPipe->Frho = FPipe->FU0.row(0) / FPipe->FArea;
	FPipe->Fspeed = FPipe->FU0.row(1) / FPipe->FU0.row(0);
	FPipe->Fpressure = (FPipe->FU0.row(2) - FPipe->FU0.row(1) * FPipe->Fspeed / 2.) * (FPipe->FGamma - 1.)
					   / FPipe->FArea;
	FPipe->Ftemperature = FPipe->Fpressure / FPipe->Frho / FPipe->FR;
	FPipe->Fa = (FPipe->FGamma * FPipe->FR * FPipe->Ftemperature).sqrt();
	for(auto i = 0; i < FPipe->FA_A.size(); i++) {
		FPipe->FA_A(i) = (FPipe->Fa(i) / __cons::ARef) / pow(FPipe->Fpressure(i) / __cons::PRef / 1E5,
						 FPipe->FGamma1(i) / 2. / FPipe->FGamma(i));
	}
	FPipe->Fbeta = (FPipe->Fa - (FPipe->FGamma - 1.) / 2. * FPipe->Fspeed) / __cons::ARef;
	FPipe->Flambda = (FPipe->Fa + (FPipe->FGamma - 1.) / 2. * FPipe->Fspeed) / __cons::ARef;
	UpdateGasProperties();
	FPipe->FTotalTemperature = FPipe->Ftemperature + 0.5 * pow2(FPipe->Fspeed) / FPipe->Fcp;
	FPipe->FTotalPressure = FPipe->Fpressure
		* pow(FPipe->FTotalTemperature / FPipe->Ftemperature, FPipe->FGamma / FPipe->FGamma1);
	ComputeMaxTimeStep();
}

void TLaxWendroff::UpdateGasProperties() {
	/* TODO */
	FPipe->FGamma1 = FPipe->FGamma - 1.;
}

void use_laxwendroff(const Pipe_ptr & pipe) {
	unique_ptr<TLaxWendroff> method(new TLaxWendroff());
	method->Connect(pipe);
	pipe->setMethod(move(method));
}

