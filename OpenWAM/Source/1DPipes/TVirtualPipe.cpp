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
 * @file TVirtualPipe.cpp
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
 * This file defines a pipe connection boundary condition and helper functions.
 */

#include "TVirtualPipe.hpp"
#include "Godunov.hpp"

TVirtualPipe::TVirtualPipe() {
	FCells.resize(2);
	FFluxSign.resize(2);
	FFluxSign[0].setOnes(3);
	FFluxSign[1].setOnes(3);
	//FPipes.resize(2);
	FFlowObj.resize(2);
	FCurrentTime = -1.;
	FDischCoef = 1.;

	FValve = make_unique<TCDFijo>();
	
}

TVirtualPipe::TVirtualPipe(const TFlowObject_ptr & fobj_0, nmPipeEnd pipe_end_0, const TFlowObject_ptr & fobj_1, 
	nmPipeEnd pipe_end_1) : TVirtualPipe(){
	double dx = 0.;
	RowVector x(3);
	RowVector D(3);
	FCurrentTime = -1.;
	FFlowObj[0] = fobj_0.get();
	FFlowObj[1] = fobj_1.get();

	if (pipe_end_0 == nmLeft) {
		FCells[0] = 0;
		D.setConstant(__geom::Circle_diameter(FFlowObj[0]->getArea()(0)));
	}
	else {
		FCells[0] = FFlowObj[0]->getNin() - 1;
		D.setConstant(__geom::Circle_diameter(FFlowObj[0]->getArea().tail(1)(0)));
	}
	if (pipe_end_1 == nmLeft) {
		FCells[1] = 0;
		D.setConstant(min(D(0), __geom::Circle_diameter(FFlowObj[1]->getArea()(0))));
	}
	else {
		FCells[1] = FFlowObj[1]->getNin() - 1;
		D.setConstant(min(D(0), __geom::Circle_diameter(FFlowObj[1]->getArea().tail(1)(0))));
	}
	if (pipe_end_0 == nmLeft) {
		FFluxSign[0](0) = -1;
		FFluxSign[0](2) = -1;
	}
	if (pipe_end_1 == nmRight) {
		FFluxSign[1](0) = -1;
		FFluxSign[1](2) = -1;
	}
	dx = FFlowObj[0]->getDeltaX();
	dx = min(dx, FFlowObj[1]->getDeltaX());
	x.setZero();
	x(1) = dx;
	x(2) = 2 * dx;
	setGeometry(x, dx, D);
	FGeometricalArea = FArea(1);
	Fpressure.setOnes(2);
	Fspeed.setZero(2);
	Ftemperature.setOnes(2);
	WorkingFluid.resize(2);
	WorkingFluid[0] = make_shared<TFluid>(FFlowObj[0]->getFluidComponents());
	WorkingFluid[0]->SetComposition(FFlowObj[0]->getFluid(FCells[0])->GetComposition());
	WorkingFluid[1] = make_shared<TFluid>(FFlowObj[1]->getFluidComponents());
	WorkingFluid[1]->SetComposition(FFlowObj[1]->getFluid(FCells[1])->GetComposition());
}

TVirtualPipe::~TVirtualPipe() {}

ColVector TVirtualPipe::Flow(double t, unsigned int pipe_number) {
	std::lock_guard<std::mutex> lock(vp_mtx);
	if(needsUpdate(t)) {
		UpdateStateVector(t);
		FMethod->ComputeFlux(FU0, FMethod->FW, FGamma, FGamma1);
		FMethod->UpdateLeftRightVars();
		FMethod->SolveCentralCells();
		setArea();
	}
	return FMethod->getFluxes().col(1) * FFluxSign[pipe_number] * FArea(1);
}

ColVector TVirtualPipe::Flow()
{
	std::lock_guard<std::mutex> lock(vp_mtx);
	UpdateStateVector();
	FMethod->ComputeFlux(FU0, FMethod->FW, FGamma, FGamma1);
	FMethod->UpdateLeftRightVars();
	FMethod->SolveCentralCells();
	setArea();
	return FMethod->getFluxes().col(1);
}

double TVirtualPipe::CalculateDCin(double t)
{
	FValve->GetCDin(t);
	FDischCoef = FValve->getCDTubVol();
	return FDischCoef;
}

double TVirtualPipe::CalculateDCout(double t)
{
	FValve->GetCDout(t);
	FDischCoef = FValve->getCDVolTub();
	return FDischCoef;
}

double TVirtualPipe::getVirtualArea() const {
	return FArea(1);
}

bool TVirtualPipe::needsUpdate(double t) const {
	auto t_pipe_min = min(FFlowObj[0]->getCurrentTime(), FFlowObj[1]->getCurrentTime());
	if(FCurrentTime < t_pipe_min) {
		return true;
	} else {
		return false;
	}
}


void TVirtualPipe::setArea(double A) {
	FArea(1) = A;
}

void TVirtualPipe::setArea() {
	FArea(1) = FGeometricalArea;
}

void TVirtualPipe::UpdateStateVector(double time) {
	FCurrentTime = time;
	UpdateStateVector();
}

void TVirtualPipe::UpdateStateVector()
{
	Fpressure(0) = FFlowObj[0]->getPressure(FCells[0]);
	Fpressure(1) = FFlowObj[1]->getPressure(FCells[1]);
	Fspeed(0) = FFlowObj[0]->getSpeed(FCells[0]) * FFluxSign[0](0);
	Fspeed(1) = FFlowObj[1]->getSpeed(FCells[1]) * FFluxSign[1](0);
	Ftemperature(0) = FFlowObj[0]->getTemperature(FCells[0]);
	Ftemperature(1) = FFlowObj[1]->getTemperature(FCells[1]);
	WorkingFluid[0]->SetComposition(FFlowObj[0]->getFluid(FCells[0])->GetComposition());
	WorkingFluid[1]->SetComposition(FFlowObj[1]->getFluid(FCells[1])->GetComposition());
	FMethod->setPTU(Fpressure, Ftemperature, Fspeed);
}
