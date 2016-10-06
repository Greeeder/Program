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
 * \file Connect0Dto1D.cpp
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
 * This file defines a boundary condition with constant pressure
 * and temperature.
 */

#include "Connect0Dto1D.hpp"
#include "TBasicPlenum.h"

TConnect0Dto1D::TConnect0Dto1D() {}

TConnect0Dto1D::TConnect0Dto1D(const TFlowObject_ptr & pipe, nmPipeEnd pipe_end,
	const VirtualPipe_ptr & virtual_pipe) : TBoundaryCondition(pipe, pipe_end){
	FVirtualPipe = virtual_pipe;
	OneD = dynamic_cast<TPipe*>(FVirtualPipe->getFlowObject(0));
	ZeroD = dynamic_cast<TBasicPlenum*>(FVirtualPipe->getFlowObject(1));
	FPipeNumber = 0.;
	T_buffer = ZeroD->getTemperature();
	WorkFluid = ZeroD->getFluid();
	InletFlow = false;
}

TConnect0Dto1D::TConnect0Dto1D(const TFlowObject_ptr & zerod, int con_index,
	const VirtualPipe_ptr & virtual_pipe) : TBoundaryCondition(zerod, con_index){
	FVirtualPipe = virtual_pipe;
	OneD = dynamic_cast<TPipe*>(FVirtualPipe->getFlowObject(0));
	ZeroD = dynamic_cast<TBasicPlenum*>(FVirtualPipe->getFlowObject(1));
	T_buffer = ZeroD->getTemperature();
	WorkFluid = ZeroD->getFluid();
	FPipeNumber = 1.;
}

ColVector TConnect0Dto1D::Flux(double t, double dt) {

	bool update = FVirtualPipe->needsUpdate(t);
	auto entry = FVirtualPipe->getCell(1);
	double rho;
	double T_BC = ZeroD->getTemperature();
	double p_BC = ZeroD->getPressure();
	auto signo = 1;
	auto damping = 0.0005;
	int condition = 0;
	auto flow = 0.;

	if (FPipeEnd_0 == nmLeft && FPipeNumber == 0) signo = -1;

	auto v_pipe = signo * OneD->getSpeed(FCell_0);

	if (update) {
		
		auto v_BC = ZeroD->getSpeed(entry);
		auto T_pipe = OneD->getTemperature(FCell_0);
		auto T_t_pipe = OneD->getTotalTemperature(FCell_0);
		auto p_pipe = OneD->getPressure(FCell_0);
		auto p_t_pipe = OneD->getTotalPressure(FCell_0);
		auto cp = OneD->getcp(FCell_0);
		auto g = OneD->getGamma(FCell_0);
		auto R = OneD->getR(FCell_0);
		auto p_obj = p_BC;

		flow = FVirtualPipe->Flow()(0);
		if (flow > 1.e-5) {
			condition = 1;
		}
		else if(flow < -1.e-5) {
			condition = -1;
		}

		if (condition == -1) {
			if (InletFlow) {
				//if (v_BC > 0) v_BC = 0;
				InletFlow = false;
			}
			auto cd = FVirtualPipe->CalculateDCout(t);
			T_buffer = T_buffer + damping * (T_BC - T_buffer);
			if (cd > 1.e-5) {
				if (p_t_pipe >= p_BC) {

					auto p_star = p_t_pipe / pow(2 / (g + 1.), g / (1. - g));
					auto v = [&](double p) {
						if (p < p_star) {
							p = p_star;
						}
						return sqrt(2. * cp * T_t_pipe
							* (1. - pow(p / p_t_pipe, (g - 1) / g)));
					};
					auto rho = [&](double p) {
						if (p < p_star) {
							p = p_star;
						}
						return p_t_pipe / R / T_t_pipe * pow(p / p_t_pipe, 1. / g);
					};
					auto loop = [&](double p) {
						return rho(p_BC) * v(p_BC) * cd - rho(p) * v(p);
					};
					p_obj = zbrent(loop, p_BC, p_t_pipe, 1E-6);
				}
				else if (p_BC > p_pipe) {

					auto v = [&](double p1, double p2) {
						return sqrt(2. * cp * T_BC
							* (1. - pow(p2 / p1, (g - 1) / g)));
					};
					auto rho = [&](double p1, double p2) {
						return p1 / R / T_BC * pow(p2 / p1, 1. / g);
					};
					auto p_star = [&](double p) {
						return p / pow(2 / (g + 1.), g / (1. - g));
					};
					auto loop = [&](double p) {
						auto p2_1 = p_star(p_BC);
						auto p2_2 = p_star(p);
						if (p2_1 < p_pipe) {
							p2_1 = p_pipe;
						}
						if (p2_2 < p_pipe) {
							p2_2 = p_pipe;
						}
						return rho(p_BC, p2_1) * v(p_BC, p2_1) * cd
							- rho(p, p2_2) * v(p, p2_2);
					};
					p_obj = zbrent(loop, p_pipe, p_BC, 1E-6);
				}
				else {
					p_obj = p_BC;
				}
				ZeroD->setPtTtU(entry, p_obj, T_buffer, v_BC);
			}
			else {
				p_obj = p_pipe;
				v_BC = -signo * OneD->getSpeed(FCell_0);
				ZeroD->setPTU(entry, p_obj, T_pipe, v_BC);
			}
		} else if (condition == 1){
			if (!InletFlow) {
				//if (v_BC < 0) v_BC = 0;
				InletFlow = true;
			}
			auto cd = FVirtualPipe->CalculateDCin(t);
			if (cd > 1.e-5) {
				if (p_t_pipe >= p_BC) {

					T_buffer = T_buffer + damping * (T_t_pipe - pow2(v_BC) / 2. / cp - T_buffer);
					auto p_star = p_t_pipe / pow(2 / (g + 1.), g / (1. - g));
					auto v = [&](double p) {
						if (p < p_star) {
							p = p_star;
						}
						return sqrt(2. * cp * T_t_pipe
							* (1. - pow(p / p_t_pipe, (g - 1) / g)));
					};
					auto rho = [&](double p) {
						if (p < p_star) {
							p = p_star;
						}
						return p_t_pipe / R / T_t_pipe * pow(p / p_t_pipe, 1. / g);
					};
					auto loop = [&](double p) {
						return rho(p_BC) * v(p_BC) * cd - rho(p) * v(p);
					};
					p_obj = zbrent(loop, p_BC, p_t_pipe, 1E-6);
				}
				else if (p_BC > p_pipe) {

					T_buffer = T_buffer + damping * (T_pipe - T_buffer);
					auto v = [&](double p) {
						return sqrt(2. * cp * T_buffer
							* (1. - pow(p_pipe / p, (g - 1) / g)));
					};
					auto rho = [&](double p) {
						return p / R / T_buffer * pow(p_pipe / p, 1. / g);
					};
					auto loop = [&](double p) {
						return rho(p_BC) * v(p_BC) * cd - rho(p) * v(p);
					};
					p_obj = zbrent(loop, p_pipe, p_BC, 1E-6);
				}
				else {
					p_obj = p_BC;
				}
			}
			else {
				p_obj = p_pipe;
				v_BC = -signo * OneD->getSpeed(FCell_0);
			}
			ZeroD->setPTU(entry, p_obj, T_buffer, v_BC);
		}
		else {
			//T_buffer = T_pipe;
			p_obj = p_pipe;
			v_BC = - signo * OneD->getSpeed(FCell_0);
			ZeroD->setPTU(entry, p_obj, T_pipe, v_BC);
		}
	}
	FFlux_0 = FVirtualPipe->Flow(t, FPipeNumber) / getArea();
	if (FVirtualPipe->getVirtualArea() < getArea()) {

		FFlux_0.row(1) = FFlux_0.row(1)
			+ FVirtualPipe->getPressure()(FPipeNumber)
			* (1. - FVirtualPipe->getVirtualArea() / getArea());

	}
	if (update) {
		rho = ZeroD->getPressure(entry) / WorkFluid->FunR() / ZeroD->getTemperature(entry);
		auto v_BC = FFlux_0(0) / rho;
		if (FPipeEnd_0 == nmLeft) {
			v_BC = -v_BC;
		}
		auto p = FVirtualPipe->getPressure(Uint(1));
		auto T = FVirtualPipe->getTemperature(Uint(1));
		ZeroD->setPTU(FVirtualPipe->getCell(1), p, T, v_BC);

		//if (v_BC != 0) {
		//	if (flow / v_BC < 0) {
		//		auto p_pipe = OneD->getPressure(FCell_0);
		//		auto T_pipe = OneD->getTemperature(FCell_0);
		//		auto p_plm = FVirtualPipe->getPressure(Uint(1));
		//		auto T_plm = FVirtualPipe->getTemperature(Uint(1));
		//		auto v_pipe = -signo * OneD->getSpeed(FCell_0);
		//		ZeroD->setPTU(entry, p_pipe, T_pipe, v_pipe);
		//		FVirtualPipe->Flow();
		//		FFlux_0.setZero();
		//		FFlux_0(1) = p_pipe;

		//		ZeroD->setPTU(entry, p_plm, T_plm, v_BC);
		//	}
		//}

		if (ZeroD->getSpeed(entry) > 0.) {
			WorkFluid = FVirtualPipe->getFluid(0);
		}
		else {
			WorkFluid = FVirtualPipe->getFluid(1);
		}
	}


	return FFlux_0;
}

TFluid_ptr TConnect0Dto1D::FluidBC(){

	if (FFlux_0(0, 0) > 0){
		return FVirtualPipe->getFluid(0);
	}
	else{
		return FVirtualPipe->getFluid(1);
	}
}

void TConnect0Dto1D::setPT(double p, double T) {
	//FCCPipe->setPTU(p, T, 0.);
}

void attach_to_0Dto1Dconnection(const int id, const TFlowObject_ptr& pipe, nmPipeEnd pipe_end,
	const TFlowObject_ptr& zeroD, TTipoValvula_ptr valve) {

	TBasicPlenum* plenum = dynamic_cast<TBasicPlenum*>(zeroD.get());
	TPipe* duct = dynamic_cast<TPipe*>(pipe.get());
	plenum->AppendNewEntry();

	nmPipeEnd zd_end = nmLeft;
	VirtualPipe_ptr virtual_pipe = make_shared<TVirtualPipe>(pipe,
		pipe_end, zeroD, zd_end);

	virtual_pipe->setName("Boundary_" + to_string(id));

	unique_ptr<TGodunov> method(new TGodunov());
	virtual_pipe->WorkingFluid.resize(virtual_pipe->FNin - 1);
	if (pipe_end == nmLeft){
		virtual_pipe->WorkingFluid.front() = pipe->getFluid(0);
	}
	else{
		virtual_pipe->WorkingFluid.front() = pipe->getFluid(pipe->getNin() - 1);
	}
	virtual_pipe->WorkingFluid.back() = zeroD->getFluid(0);

	virtual_pipe->setValve(move(valve));

	//set_numerical_method(virtual_pipe, "Godunov");
	method->Connect(virtual_pipe);
	method->setRiemannSolver(&HLLC);
	virtual_pipe->FMethod = std::move(method);

	virtual_pipe->FMethod->setPTU(1E5, 300, 0);

	BoundaryCondition_ptr first(new TConnect0Dto1D(pipe, pipe_end, virtual_pipe));

	BoundaryCondition_ptr second(new TConnect0Dto1D(zeroD, virtual_pipe->getCell(1) - 1, virtual_pipe));

	if (pipe_end == nmLeft) {
		dynamic_cast<TConnect0Dto1D*>(first.get())->FArea =
			pipe->getArea()(0);
		dynamic_cast<TConnect0Dto1D*>(second.get())->FArea =
			pipe->getArea()(0);
		duct->setLeftBC(move(first));
	}
	else {
		dynamic_cast<TConnect0Dto1D*>(first.get())->FArea =
			pipe->getArea().tail(1)(0);
		dynamic_cast<TConnect0Dto1D*>(second.get())->FArea =
			pipe->getArea().tail(1)(0);
		duct->setRightBC(move(first));
	}
	plenum->setBoundary(move(second));

}
