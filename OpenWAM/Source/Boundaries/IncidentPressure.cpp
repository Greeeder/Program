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
 * \file IncidentPressure.cpp
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
 * This file defines a boundary condition that imposes an incident pressure
 * wave.
 */

#include "IncidentPressure.hpp"

TIncidentPressureBC::TIncidentPressureBC() {
	FTimeVector.setZero(1);
	FIncidentPressure.setOnes(1) * 101325.;
	FA_A.setOnes(1);
}

TIncidentPressureBC::TIncidentPressureBC(const Pipe_ptr& pipe,
	nmPipeEnd pipe_end, const ConstantConditionsPipe_ptr& pipe_cc,
	const VirtualPipe_ptr& virtual_pipe):
	TConstantConditionsBC(pipe, pipe_end, pipe_cc, virtual_pipe) {
	FTimeVector.setZero(1);
	FIncidentPressure.setOnes(1) * 101325.;
	FA_A.setOnes(1);
}

void TIncidentPressureBC::setPA_A(const RowVector& t, const RowVector& p,
	const RowVector& A_A) {
	FTimeVector = t;
	FIncidentPressure = p;
	FA_A = A_A;
}

ColVector TIncidentPressureBC::Flux(double t, double dt) {
	double g = FPipe_0->getGamma(FCell_0);
	double p = periodic_linear_interp(FTimeVector, FIncidentPressure, t);
	double A = periodic_linear_interp(FTimeVector, FA_A, t);
	double lambda = 1.;
	double beta = 1.;
	beta = (2. * pow(p / __units::BarToPa(__cons::PRef),
		(g - 1.) / (2. * g)) - 1.) * A;
	if (FPipeEnd_0 == nmLeft) {
		lambda = FPipe_0->getBeta(Uint(0));
	} else {
		lambda = FPipe_0->getLambda(FCell_0);
	}
	double a = __cons::ARef * (lambda + beta) / 2.;
	double u = __cons::ARef * (lambda - beta) / (g - 1.);
	double T = pow2(a) / g / FPipe_0->getR(FCell_0);
	double pres = __units::BarToPa(__cons::PRef)
		* pow(a / __cons::ARef / A, 2. * g / (g - 1.));
	FCCPipe->setPTU(pres, T, u);
	return TConstantConditionsBC::Flux(t, dt);
}

void attach_to_incident_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
	const RowVector& t, const RowVector& p, const RowVector& A_A,
	TFluid_ptr fluid) {
	double D = 0.;
	if (pipe_end == nmLeft) {
		D = pipe->getD()(0);
	} else {
		D = pipe->getD().tail(1)(0);
	}
	ConstantConditionsPipe_ptr constant_pipe =
		create_constant_conditions_pipe(__units::BarToPa(__cons::PRef),
			__units::degCToK(__cons::TRef), 0, D, 1E10, 1, fluid);

	nmPipeEnd vp_end = nmLeft;
	VirtualPipe_ptr virtual_pipe =
		make_shared<TVirtualPipe>(pipe, pipe_end, constant_pipe, vp_end);
	unique_ptr<TGodunov> method(new TGodunov());
	virtual_pipe->WorkingFluid.resize(virtual_pipe->FNin - 1);
	for (auto i = 0; i < virtual_pipe->FNin - 1; i++){
		virtual_pipe->WorkingFluid[i] = make_shared<TFluid>(fluid.get());
	}
	method->Connect(virtual_pipe);
	method->setRiemannSolver(&KT);
	virtual_pipe->FMethod = std::move(method);

	virtual_pipe->FMethod->setPTU(1E5, 300, 0);
	IncidentPressureBC_ptr first(new TIncidentPressureBC(pipe, pipe_end,
		constant_pipe, virtual_pipe));
	first->setPA_A(t, p, A_A);
	if(pipe_end == nmLeft) {
		dynamic_cast<TIncidentPressureBC*>(first.get())->FArea =
			pipe->getArea()(0);
		pipe->setLeftBC(move(first));
	} else {
		dynamic_cast<TIncidentPressureBC*>(first.get())->FArea =
			pipe->getArea().tail(1)(0);
		pipe->setRightBC(move(first));
	}
}

void load_incident_pressure_data(const string& file_name,
	char separator, const string& t_label, const string& p_label,
	const string& A_A_label, const string& t_unit, const string& p_unit,
	RowVector* t, RowVector* p, RowVector* A_A) {
	auto data = load_table(file_name, separator);
	*t = data[t_label];
	*p = data[p_label];
	*A_A = data[A_A_label];
	if (t->size() < 2) {
		std::string error_msg = "Not enough time points in file "
			+ file_name + ", variable key " + t_label;
		throw Exception(error_msg);
	}
	if (p->size() < 2) {
		std::string error_msg = "Not enough pressure points in file "
			+ file_name + ", variable key " + p_label;
		throw Exception(error_msg);
	}
	if (A_A->size() < 2) {
		std::string error_msg = "Not enough entropy level points in file "
			+ file_name + ", variable key " + A_A_label;
		throw Exception(error_msg);
	}
	auto to_s = [&] (double x) {return to_seconds(x, t_unit);};
	auto to_Pa = [&] (double x) {return __units::BarToPa(to_bar(x, p_unit));};
	t->unaryExpr(to_s);
	p->unaryExpr(to_Pa);
}