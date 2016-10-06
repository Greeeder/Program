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
 * \file PressureBC.cpp
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

#include "PressureBC.hpp"

TPressureBC::TPressureBC() {
	T_buffer = __cons::TRef;
	FCD = 1.;
}

TPressureBC::TPressureBC(const Pipe_ptr & pipe,
	nmPipeEnd pipe_end, const ConstantConditionsPipe_ptr & pipe_cc,
	const VirtualPipe_ptr & virtual_pipe):
	TConstantConditionsBC(pipe, pipe_end, pipe_cc, virtual_pipe) {
		T_buffer = FCCPipe->getTemperature()(0);
		FCD = 1.;
	}

ColVector TPressureBC::Flux(double t, double dt) {

	int condition = 0;
	auto flow = 0.;
	auto signo = 1;
	auto damping = 0.005;
	if (FPipeEnd_0 == nmLeft && FPipeNumber == 0) 
		signo = -1;

	FCCPipe->setTimeStep(dt);
	FCCPipe->Solve(t);
	if (FVirtualPipe->needsUpdate(t)) {
		auto p_BC = FCCPipe->getPressure()(0);
		double T_BC = FCCPipe->getTemperature()(0);
		auto v_BC = FCCPipe->getSpeed()(0);
		auto T_pipe = FPipe_0->getTemperature(FCell_0);
		auto T_t_pipe = FPipe_0->getTotalTemperature(FCell_0);
		auto p_pipe = FPipe_0->getPressure(FCell_0);
		auto v_pipe = FPipe_0->getSpeed(FCell_0);
//		auto rho_pipe = FPipe_0->getDensity(FCell_0);
		auto p_t_pipe = FPipe_0->getTotalPressure(FCell_0);
		auto cp = FPipe_0->getcp(FCell_0);
		auto g = FPipe_0->getGamma(FCell_0);
		auto R = FPipe_0->getR(FCell_0);
		auto p_obj = p_BC;

		flow = FVirtualPipe->Flow()(0);
		if (flow > 1.e-5) {
			condition = 1;
		}
		else if (flow < -1.e-5) {
			condition = -1;
		}

		if (condition == -1) {
			auto cd = getCD();
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
				FCCPipe->setPtTtU(p_obj, T_BC, v_BC);
			}
			else {
				p_obj = p_pipe;
				v_BC = -signo * v_pipe;
				FCCPipe->setPTU(p_obj, T_pipe, v_BC);
			}
			
		} else if (condition == 1){
			auto cd = getCD();
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
				v_BC = -signo * v_pipe;
			}
			FCCPipe->setPTU(p_obj, T_buffer, v_BC);
		}
		else {
			p_obj = p_pipe;
			v_BC = - signo * v_pipe;
		}
		TPipeConnection::Flux(t, dt);
		auto v_0 = v_BC;
		v_BC = FFlux_0(0) / FCCPipe->getDensity(Uint(0));
		if (FPipeEnd_0 == nmLeft) {
			v_BC = -v_BC;
		}

		FCCPipe->setPTU(p_BC, T_BC, v_BC);

		if (v_BC > 0.) {
			WorkFluid = FVirtualPipe->getFluid(0);
		}
		else {
			WorkFluid = FVirtualPipe->getFluid(1);
		}
	}
	return FFlux_0;
}

double TPressureBC::getCD() const {
	return FCD;
}

void TPressureBC::setCD(double CD) {
	FCD = CD;
}


void TPressureBC::setPT(double p, double T) {
	FCCPipe->setPTU(p, T, FCCPipe->getSpeed()(0));
}

void attach_to_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
	double p, double T, TFluid_ptr fluid) {
	double D = 0.;
	if (pipe_end == nmLeft) {
		D = pipe->getD()(0);
	} else {
		D = pipe->getD().tail(1)(0);
	}
	ConstantConditionsPipe_ptr constant_pipe =
		create_constant_conditions_pipe(p, T, 0, D, 1E10, 1, fluid);

	nmPipeEnd vp_end = nmLeft;
	VirtualPipe_ptr virtual_pipe =
		make_shared<TVirtualPipe>(pipe, pipe_end, constant_pipe, vp_end);
	unique_ptr<TGodunov> method(new TGodunov());
	virtual_pipe->WorkingFluid.resize(virtual_pipe->FNin - 1);
	for (auto i = 0; i < virtual_pipe->FNin - 1; i++){
		virtual_pipe->WorkingFluid[i] = make_shared<TFluid>(fluid.get());
	}
	method->Connect(virtual_pipe);
	method->setRiemannSolver(&HLLC);
	virtual_pipe->FMethod = std::move(method);

	virtual_pipe->FMethod->setPTU(1E5, 300, 0);
	PressureBC_ptr first(new TPressureBC(pipe, pipe_end,
		constant_pipe, virtual_pipe));
	if(pipe_end == nmLeft) {
		dynamic_cast<TPressureBC*>(first.get())->FArea =
			pipe->getArea()(0);
		pipe->setLeftBC(move(first));
	} else {
		dynamic_cast<TPressureBC*>(first.get())->FArea =
			pipe->getArea().tail(1)(0);
		pipe->setRightBC(move(first));
	}
}
