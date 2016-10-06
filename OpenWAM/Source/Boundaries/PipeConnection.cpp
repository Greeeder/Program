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
 * \file PipeConnection.cpp
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
 * This file defines a pipe connection boundary condition and helper functions.
 */

#include "PipeConnection.hpp"
#include "Godunov.hpp"

TPipeConnection::TPipeConnection() {}

TPipeConnection::TPipeConnection(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end_0,
	const VirtualPipe_ptr & virtual_pipe,
	unsigned int pipe_number) : TBoundaryCondition(pipe_0, pipe_end_0) {
	FVirtualPipe = virtual_pipe;
	FPipeNumber = pipe_number;
}

TPipeConnection::~TPipeConnection() {}

ColVector TPipeConnection::Flux(double t, double dt) {
	FFlux_0 = FVirtualPipe->Flow(t, FPipeNumber) / getArea();
	if (FVirtualPipe->getVirtualArea() < getArea()) {
		FFlux_0.row(1) = FFlux_0.row(1)
			+ FVirtualPipe->getPressure()(FPipeNumber)
			* (1. - FVirtualPipe->getVirtualArea() / getArea());
	}
	return FFlux_0;
}

TFluid_ptr TPipeConnection::FluidBC(){
	if (FPipeNumber == 0)
		return FVirtualPipe->getWorkingFluid().back();
	else
		return FVirtualPipe->getWorkingFluid().front();
}

double TPipeConnection::getAeff() const {
	return FVirtualPipe->getVirtualArea();
}

void TPipeConnection::setAeff(double Aeff) {
	if (Aeff < getArea()) {
		FVirtualPipe->setArea(Aeff);
	}
}

void TPipeConnection::updateBCCharacteristics(double t, double dt) {
	/* TODO */
}

void attach_pipes(const Pipe_ptr pipe_0, nmPipeEnd pipe_end_0,
	const Pipe_ptr pipe_1, nmPipeEnd pipe_end_1) {
	VirtualPipe_ptr virtual_pipe = make_shared<TVirtualPipe>(pipe_0,
		pipe_end_0, pipe_1, pipe_end_1);
	unique_ptr<TGodunov> method(new TGodunov());
	method->Connect(virtual_pipe);
	method->setRiemannSolver(&KT);
	virtual_pipe->FMethod = std::move(method);
	virtual_pipe->FMethod->setPTU(1E5, 300, 0);
	BoundaryCondition_ptr first(new TPipeConnection(pipe_0,
		pipe_end_0, virtual_pipe, 0));
	BoundaryCondition_ptr second(new TPipeConnection(pipe_1,
		pipe_end_1, virtual_pipe, 1));
	if(pipe_end_0 == nmLeft) {
		dynamic_cast<TPipeConnection*>(first.get())->FArea =
			pipe_0->getArea()(0);
		pipe_0->setLeftBC(move(first));
	} else {
		dynamic_cast<TPipeConnection*>(first.get())->FArea =
			pipe_0->getArea().tail(1)(0);
		pipe_0->setRightBC(move(first));
	}
	if(pipe_end_1 == nmLeft) {
		dynamic_cast<TPipeConnection*>(second.get())->FArea =
			pipe_1->getArea()(0);
		pipe_1->setLeftBC(move(second));
	} else {
		dynamic_cast<TPipeConnection*>(second.get())->FArea =
			pipe_1->getArea().tail(1)(0);
		pipe_1->setRightBC(move(second));
	}
}
