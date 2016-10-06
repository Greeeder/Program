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
 * @file ClosedEnd.cpp
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
 * This file defines a closed end boundary condition and helper functions.
 */

#include "ClosedEnd.hpp"

TClosedEnd::TClosedEnd() :
	TBoundaryCondition() {
}

TClosedEnd::TClosedEnd(const Pipe_ptr & pipe_0, nmPipeEnd pipe_end) :
	TBoundaryCondition(pipe_0, pipe_end) {
}

ColVector TClosedEnd::Flux(double t, double dt) {
	FCurrentTime = FPipe_0->getCurrentTime();
	FFlux_0(1) = FPipe_0->getPressure(FCell_0);
	return FFlux_0;
}

TFluid_ptr TClosedEnd::FluidBC(){
	return FPipe_0->getWorkingFluid().front();
}

void TClosedEnd::updateBCCharacteristics(double t, double dt) {
	if(FPipeEnd_0 == nmRight) {
		setBeta(getLambda());
	} else {
		setLambda(getBeta());
	}
}

void close_pipe_end(const Pipe_ptr & pipe, nmPipeEnd pipe_end) {
	BoundaryCondition_ptr wall(new TClosedEnd(pipe, pipe_end));
	if(pipe_end == nmLeft) {
		pipe->setLeftBC(move(wall));
	} else {
		pipe->setRightBC(move(wall));
	}
}
