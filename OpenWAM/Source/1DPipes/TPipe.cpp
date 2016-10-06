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
 * @file TPipe.cpp
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
 * This file defines a basic one-dimensional pipe.
 */

#include "TPipe.hpp"
#include "AllPipeMethods.hpp"
#include "AllExtraSourceTerms.hpp"
#include "BoundaryCondition.hpp"

PipeInsResults::PipeInsResults()
{
	Distance = 0.;
	MassFlow = false;
	Pressure = false;
	Temperature = false;
	TotalPressure = false;
	TotalTemperature = false;
	Velocity = false;
	LeftPressure = false;
	RightPressure = false;
	R = false;
	Cp = false;
	Cv = false;
	Gamma = false;
}

PipeOutput::PipeOutput() {
	Location = 0.;
	LeftPressure.avg = false;
	LeftPressure.ins = false;
	LeftPressure.ind = 0;
	MassFlow.avg = false;
	MassFlow.ins = false;
	MassFlow.ind = 0;
	Pressure.avg = false;
	Pressure.ins = false;
	Pressure.ind = 0;
	RightPressure.avg = false;
	RightPressure.ins = false;
	RightPressure.ind = 0;
	Temperature.avg = false;
	Temperature.ins = false;
	Temperature.ind = 0;
	TotalPressure.avg = false;
	TotalPressure.ins = false;
	TotalPressure.ind = 0;
	TotalTemperature.avg = false;
	TotalTemperature.ins = false;
	TotalTemperature.ind = 0;
	Velocity.avg = false;
	Velocity.ins = false;
	Velocity.ind = 0;
}

TPipe::TPipe() {
	FCurrentTime = 0.;
	FXref = 1;
	FTimeStep = 1;
	FIsIntegrated = true;
	FPipeFriction = TPipeFriction_ptr(new TPipeFriction());

	Is_plot = false;
	Is_cycleout = false;
}

TPipe::~TPipe() {}


double TPipe::getMaxTimeStep() {
	return FMethod->getMaxTimeStep();
}

RowVector TPipe::getArea() const {
	return FArea;
}

double TPipe::getArea(Uint i) const {
	if (i < FArea.size() - 1) {
		return FArea(i);
	} else {
		return FArea.tail(1)(0);
	}
}

RowVector TPipe::getA_A() const {
	std::lock_guard<std::mutex> lock(mtx);
	return (Fa / __cons::ARef) / pow(Fpressure / __units::BarToPa(__cons::PRef),
		FGamma1 / 2. / FGamma);
}

double TPipe::getA_A(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i > FNin) {
		i = FNin;
	}
	return (Fa(i) / __cons::ARef) / pow(Fpressure(i) / __units::BarToPa(__cons::PRef),
		FGamma1(i) / 2. / FGamma(i));
}

double TPipe::getA_A(double x) const {
// 	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	auto g = getGamma(x);
	return (getSpeedOfSound(x) / __cons::ARef) / pow(getPressure(x) / __units::BarToPa(__cons::PRef),
		(g - 1.) / 2. / g);
}

RowVector TPipe::getBeta() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Fbeta;
}

double TPipe::getBeta(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Fbeta(i);
	} else {
		return Fbeta(FNin - 1);
	}
}

double TPipe::getBeta(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Fbeta, x);
}

RowVector TPipe::getcp() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Fcp;
}

double TPipe::getcp(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Fcp(i);
	} else {
		return Fcp(FNin - 1);
	}
}

double TPipe::getcp(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Fcp, x);
}

RowVector TPipe::getD() const {
	return FD;
}

double TPipe::getD(Uint i) const {
	return FD(i);
}

double TPipe::getDeltaX() const {
	return FXref;
}

RowVector TPipe::getDensity() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Frho;
}

double TPipe::getDensity(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Frho(i);
	} else {
		return Frho(FNin - 1);
	}
}

double TPipe::getDensity(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Frho, x);
}

TGenericSource* TPipe::getExtraSources() const {
	return FSource.get();
}

TFluid_ptr TPipe::getFluid(Uint i) const {
	return WorkingFluid[i];
}

TComponentArray_ptr TPipe::getFluidComponents() const {
	return WorkingFluid[0]->GetComponents();
}


RowVector TPipe::getGamma() const {
	std::lock_guard<std::mutex> lock(mtx);
	return FGamma;
}

double TPipe::getGamma(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (i < FNin)
	{
		return FGamma(i);
	}
	else
	{
		return FGamma(FNin - 1);
	}
}

double TPipe::getGamma(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, FGamma, x);
}

RowVector TPipe::getLambda() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Flambda;
}

double TPipe::getLambda(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Flambda(i);
	} else {
		return Flambda(FNin - 1);
	}
}

double TPipe::getLambda(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Flambda, x);
}

TBoundaryCondition* TPipe::getLeftBC() const {
	return FLeftBC.get();
}

int TPipe::getLeftNode() const {
	return FNodeLeft;
}

double TPipe::getLeftPressure(double x) const {
	auto p = getPressure(x);
	auto g = getGamma(x);
	auto a = pow(p / __units::BarToPa(__cons::PRef), (g - 1.) / (2. * g));
	auto b = 1. - (g - 1.) / 2. * getSpeed(x) / getSpeedOfSound(x);
	return __units::BarToPa(__cons::PRef) * pow(0.5 * (1. + a * b),
		2. * g / (g - 1.));
}

RowVector TPipe::getMArea() const {
	return FMArea;
}

double TPipe::getMassFlow(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx(0)) {
		return FLeftBC->getMassFlow();
	} else if (x > Fx.tail(1)(0)) {
		return FRightBC->getMassFlow();
	} else {
		return linear_interp(Fx, FMethod->getFluxes().row(0) * FArea, x);
	}
}

Uint TPipe::getNin() const {
	return FNin;
}

Uint TPipe::getPipeNumber() const {
	return FPipeNumber;
}

RowVector TPipe::getPressure() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Fpressure;
}

double TPipe::getPressure(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Fpressure(i);
	} else {
		return Fpressure(FNin - 1);
	}
}

double TPipe::getPressure(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Fpressure, x);
}

RowVector TPipe::getR() const {
	std::lock_guard<std::mutex> lock(mtx);
	return FR;
}

double TPipe::getR(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return FR(i);
	} else {
		return FR(FNin - 1);
	}
}

double TPipe::getR(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, FR, x);
}

TBoundaryCondition* TPipe::getRightBC() const {
	return FRightBC.get();
}

int TPipe::getRightNode() const {
	return FNodeRight;
}

double TPipe::getRightPressure(double x) const {
	auto p = getPressure(x);
	auto g = getGamma(x);
	auto a = pow(p / __units::BarToPa(__cons::PRef), (g - 1.) / (2. * g));
	auto b = 1. + (g - 1.) / 2. * getSpeed(x) / getSpeedOfSound(x);
	return __units::BarToPa(__cons::PRef) * pow(0.5 * (1. + a * b),
		2. * g / (g - 1.));
}

RowVector TPipe::getSpeed() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Fspeed;
}

double TPipe::getSpeed(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Fspeed(i);
	} else {
		return Fspeed(FNin - 1);
	}
}

double TPipe::getSpeed(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Fspeed, x);
}

RowVector TPipe::getSpeedNB() const {
	return Fspeed;
}


RowVector TPipe::getSpeedOfSound() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Fa;
}

double TPipe::getSpeedOfSound(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Fa(i);
	} else {
		return Fa(FNin - 1);
	}
}

double TPipe::getSpeedOfSound(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Fa, x);
}

RowVector TPipe::getSpeedOfSoundNB() const {
	return Fa;
}

RowArray TPipe::getStateVector() const
{
	std::lock_guard<std::mutex> lock(mtx);
	return FU0;
}

ColArray TPipe::getStateVector(Uint i) const
{
	std::lock_guard<std::mutex> lock(mtx);
	return FU0.col(i);
}

RowVector TPipe::getTemperature() const {
	std::lock_guard<std::mutex> lock(mtx);
	return Ftemperature;
}

double TPipe::getTemperature(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return Ftemperature(i);
	} else {
		return Ftemperature(FNin - 1);
	}
}

double TPipe::getTemperature(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, Ftemperature, x);
}

RowVector TPipe::getTotalPressure() const {
	std::lock_guard<std::mutex> lock(mtx);
	return FTotalPressure;
}

double TPipe::getTotalPressure(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return FTotalPressure(i);
	} else {
		return FTotalPressure(FNin - 1);
	}
}

double TPipe::getTotalPressure(double x) const {
	std::lock_guard<std::mutex> lock(mtx);
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, FTotalPressure, x);
}

RowVector TPipe::getTotalTemperature() const {
	std::lock_guard<std::mutex> lock(mtx);
	return FTotalTemperature;
}

double TPipe::getTotalTemperature(Uint i) const {
	std::lock_guard<std::mutex> lock(mtx);
	if(i < FNin) {
		return FTotalTemperature(i);
	} else {
		return FTotalTemperature(FNin - 1);
	}
}

double TPipe::getTotalTemperature(double x) const {
	if (x < Fx_sv(0)) x = Fx_sv(0);
	return linear_interp(Fx_sv, FTotalTemperature, x);
}

RowVector TPipe::getVolume() const {
	return FVolume;
}

double TPipe::getVolume(Uint i) const {
	if(i < FNin) {
		return FVolume(i);
	} else {
		return FVolume(FNin - 1);
	}
}

RowVector TPipe::getX() const {
	return Fx;
}

void TPipe::HeaderForCycleOutput(std::vector<std::string> & output) const{

	double pos = 0.;
	string base = "Pipe/" + FName + "/";
	string base2 = base;

	for (const auto& out : FPipeOutput) {
		pos = out.Location;
		base2 = base + std::to_string(pos) + "/";
		if (out.LeftPressure.avg) {
			output.push_back(base2 + "Pressure[kg/s]/Left");
		}
		if (out.MassFlow.avg) {
			output.push_back(base2 + "MassFlow[kg/s]");
		}
		if (out.Pressure.avg) {
			output.push_back(base2 + "Pressure[bar]/Static");
		}
		if (out.RightPressure.avg) {
			output.push_back(base2 + "Pressure[kg/s]/Right");
		}
		if (out.Temperature.avg) {
			output.push_back(base2 + "Temperature[degC]/Static");
		}
		if (out.TotalPressure.avg) {
			output.push_back(base2 + "Pressure[bar]/Total");
		}
		if (out.TotalTemperature.avg) {
			output.push_back(base2 + "Temperature[degC]/Total");
		}
		if (out.Velocity.avg) {
			output.push_back(base2 + "Velocity[m/s]");
		}
	}
}

void TPipe::HeaderForPlots(std::vector<std::string> & output) const{
	double pos = 0.;
	string base = "Pipe/" + FName + "/";
	string base2 = base;

	for (auto& out : FPipeOutput) {
		pos = out.Location;
		base2 = base + std::to_string(pos) + "/";
		if (out.LeftPressure.ins) {
			out.LeftPressure.ind = output.size();
			output.push_back(base2 + "Pressure[bar]/Left");
		}
		if (out.MassFlow.ins) {
			out.MassFlow.ind = output.size();
			output.push_back(base2 + "MassFlow[kg/s]");
		}
		if (out.Pressure.ins) {
			out.Pressure.ind = output.size();
			output.push_back(base2 + "Pressure[bar]/Static");
		}
		if (out.RightPressure.ins) {
			out.RightPressure.ind = output.size();
			output.push_back(base2 + "Pressure[bar]/Right");
		}
		if (out.Temperature.ins) {
			out.Temperature.ind = output.size();
			output.push_back(base2 + "Temperature[degC]/Static");
		}
		if (out.TotalPressure.ins) {
			out.TotalPressure.ind = output.size();
			output.push_back(base2 + "Pressure[bar]/Total");
		}
		if (out.TotalTemperature.ins) {
			out.TotalTemperature.ind = output.size();
			output.push_back(base2 + "Temperature[degC]/Total");
		}
		if (out.Velocity.ins) {
			out.Velocity.ind = output.size();
			output.push_back(base2 + "Velocity[m/s]");
		}
	}
}

void TPipe::IntegrateWithoutUpdating() {
	FFric = FPipeFriction->FrictionCoefficient(FRe);
	FQint = FHeatTransfer->Heat(Ftemperature, FWallHT->getTwallint().array(), FRe, WorkingFluid);
	FMethod->IntegrateWithoutUpdating();
}

void TPipe::IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const
{
	double pos = 0.;
	Map<ArrayXf> dt(Output.at(0).data(), Output[0].size());
	double dt_sum = dt.sum();
	for (const auto& out : FPipeOutput) {
		if (out.LeftPressure.avg) {
			Map<ArrayXf> p(Output.at(out.LeftPressure.ind).data(), Output.at(out.LeftPressure.ind).size());
			Cicle.push_back((p * dt).sum() / dt_sum);
		}
		if (out.MassFlow.avg) {
			Map<ArrayXf> m(Output.at(out.MassFlow.ind).data(), Output.at(out.MassFlow.ind).size());
			Cicle.push_back((m * dt).sum() / dt_sum);
		}
		if (out.Pressure.avg) {
			Map<ArrayXf> p(Output.at(out.Pressure.ind).data(), Output.at(out.Pressure.ind).size());
			Cicle.push_back((p * dt).sum() / dt_sum);
		}
		if (out.RightPressure.avg) {
			Map<ArrayXf> p(Output.at(out.RightPressure.ind).data(), Output.at(out.RightPressure.ind).size());
			Cicle.push_back((p * dt).sum() / dt_sum);
		}
		if (out.Temperature.avg) {
			Map<ArrayXf> T(Output.at(out.Temperature.ind).data(), Output.at(out.Temperature.ind).size());
			Map<ArrayXf> m(Output.at(out.MassFlow.ind).data(), Output.at(out.MassFlow.ind).size());
			Cicle.push_back((T * m * dt).sum() / (m * dt).sum());
		}
		if (out.TotalPressure.avg) {
			Map<ArrayXf> p(Output.at(out.TotalPressure.ind).data(), Output.at(out.TotalPressure.ind).size());
			Cicle.push_back((p * dt).sum() / dt_sum);
		}
		if (out.TotalTemperature.avg) {
			Map<ArrayXf> T(Output.at(out.TotalTemperature.ind).data(), Output.at(out.TotalTemperature.ind).size());
			Map<ArrayXf> m(Output.at(out.MassFlow.ind).data(), Output.at(out.MassFlow.ind).size());
			Cicle.push_back((T * m * dt).sum() / (m * dt).sum());
		}
		if (out.Velocity.avg) {
			Map<ArrayXf> v(Output.at(out.Velocity.ind).data(), Output.at(out.Velocity.ind).size());
			Cicle.push_back((v * dt).sum() / dt_sum);
		}
	}
}

void TPipe::ReadInsResults(const xml_node & node) {
    auto node_res = GetNodeChild(node, "Pip_InsOutput");
	for(auto node_point = GetNodeChild(node_res, "Opi_Point"); node_point;
		node_point = node_point.next_sibling("Opi_Point")) {
		FInsOutput = true;
		InsResults.push_back(PipeInsResults());
		InsResults.back().Distance = GetXMLLength(node_point, "Distance");
		for(auto parameter = node_point.attribute("parameter");
			parameter; parameter = parameter.next_attribute()) {
			std::string par = parameter.value();
			if (par == "Cp") {
				InsResults.back().Cp = true;
			} else if(par == "Cv") {
				InsResults.back().Cv = true;
			} else if(par == "GasConstant") {
				InsResults.back().R = true;
			} else if(par == "Gamma") {
				InsResults.back().Gamma = true;
			} else if(par == "GasTemperature") {
				InsResults.back().Temperature = true;
			} else if(par == "LeftPressure") {
				InsResults.back().LeftPressure = true;
			} else if(par == "MassFlow") {
				InsResults.back().MassFlow =true;
			} else if(par == "Pressure") {
				InsResults.back().Pressure = true;
			} else if(par == "RightPressure") {
				InsResults.back().RightPressure = true;
			} else if(par == "TotalPressure") {
				InsResults.back().TotalPressure = true;
			} else if(par == "TotalTemperature") {
				InsResults.back().TotalTemperature = true;
			} else if(par == "Velocity") {
				InsResults.back().Velocity = true;
			} else {
				std::cout << "WARNING: Unknown variable: " << par <<
					std::endl;
			}
		}
	}
}

void TPipe::ReadOutput(const xml_node& node) {

	string param;
	bool Avg;
	auto node_res = node.child("Pip_OutputResults");
	for (const auto& node_point: node_res.children("Output")) {
		FPipeOutput.push_back(PipeOutput());
		if (node_point.attribute("Location"))
			FPipeOutput.back().Location = GetAttributeAsDouble(node_point, "Location");
		else
			FPipeOutput.back().Location = 0.;

		for (const auto& node_par: node_point.children("Parameter")) {
			Is_plot = true;
			param = node_par.attribute("Name").as_string();
			Avg = node_par.attribute("Average").as_bool();
			if (Avg)
				Is_cycleout = true;
			if (param == "LeftPressure") {
				FPipeOutput.back().LeftPressure.ins = true;
				FPipeOutput.back().LeftPressure.avg = Avg;
			}
			else if (param == "MassFlow") {
				FPipeOutput.back().MassFlow.ins = true;
				FPipeOutput.back().MassFlow.avg = Avg;
			}
			else if (param == "Pressure") {
				FPipeOutput.back().Pressure.ins = true;
				FPipeOutput.back().Pressure.avg = Avg;
			}
			else if (param == "RightPressure") {
				FPipeOutput.back().RightPressure.ins = true;
				FPipeOutput.back().RightPressure.avg = Avg;
			}
			else if (param == "Temperature") {
				FPipeOutput.back().Temperature.ins = true;
				FPipeOutput.back().Temperature.avg = Avg;
				if (Avg)
					FPipeOutput.back().MassFlow.ins = true;
			}
			else if (param == "TotalPressure") {
				FPipeOutput.back().TotalPressure.ins = true;
				FPipeOutput.back().TotalPressure.avg = Avg;
			}
			else if (param == "TotalTemperature") {
				FPipeOutput.back().TotalTemperature.ins = true;
				FPipeOutput.back().TotalTemperature.avg = Avg;
				if (Avg)
					FPipeOutput.back().MassFlow.ins = true;
			}
			else if (param == "Velocity") {
				FPipeOutput.back().Velocity.ins = true;
				FPipeOutput.back().Velocity.avg = Avg;
			}

		}

	}
}

void TPipe::saveState(xml_node* base) {
	auto pipe_node = base->append_child("Bop_NewPipe");
	pipe_node.append_attribute("Pipe_ID") = getPipeNumber();
	pipe_node.append_attribute("NodeL_ID") = getLeftNode(); 
	pipe_node.append_attribute("NodeR_ID") = getRightNode();
	pipe_node.append_attribute("ParallelPipes") = 1;
	pipe_node.append_attribute("Name") = getName().c_str();
	auto gas_props = pipe_node.append_child("GasProperties");
	auto units = gas_props.append_child("Units");
	units.append_attribute("Pressure") = "Pa";
	units.append_attribute("Temperature") = "K";
	units.append_attribute("Speed") = "m/s";
	units.append_attribute("Length") = "m";
	for (Uint i = 0; i < getNin(); i++) {
		auto node = gas_props.append_child("GasProps_Point");
		node.append_attribute("Position") = Fx_sv(i);
		node.append_attribute("Pressure") = getPressure(i);
		node.append_attribute("Temperature") = getTemperature(i);
		node.append_attribute("Velocity") = getSpeed(i);
	}
	auto composition = gas_props.append_child("Composition");
	auto method = pipe_node.append_child("Pip_NumericalMethod");
	method.append_attribute("MeshSize") = getDeltaX();
	method.append_attribute("Courant") = FMethod->getCourant();
	method.append_attribute("Scheme") = FMethod->getName().c_str();
	auto geo = pipe_node.append_child("Pipe_Geometry");
	geo.append_attribute("Diameter") = getD(Uint(0)); 
	units = geo.append_child("Units");
	units.append_attribute("Length") = "m";
	for (auto i = 1; i <= FNumStretch; i++) {
		auto node = geo.append_child("Geo_Stretch");
		node.append_attribute("Stretch_ID") = i;
		node.append_attribute("Length") = FLStretch(i);
		node.append_attribute("Diameter") = FDStretch(i);
	}
}

void TPipe::setBCs(BoundaryCondition_ptr leftBC, BoundaryCondition_ptr rightBC) {
	setLeftBC(move(leftBC));
	setRightBC(move(rightBC));
}

void TPipe::setCourant(double C) {
	FMethod->setCourant(C);
}

void TPipe::setExtraSource(Source_ptr source) {
	FSource = move(source);
}

void TPipe::setGeometry(const RowVector& x, double dx, const RowVector& D) {
	double n_nodes = round((x.tail(1)(0) - x(0)) / dx) + 1;
// 	if (n_nodes < 3)
// 	{
// 		n_nodes = 3;
// 	}
	dx = (x.tail(1)(0) - x(0)) / (n_nodes - 1);
	Fdx.setConstant(n_nodes - 1, dx);
	FXref = Fdx(0);
	Fx.setLinSpaced(n_nodes, x(0), x.tail(1)(0));
	FMx = (Fx.head(n_nodes - 1) + Fx.tail(n_nodes - 1)) / 2.;
	FD = linear_interp(x, D, Fx);
	FDcell = (FD.head(n_nodes - 1) + FD.tail(n_nodes - 1)) / 2.;
	FArea = __cons::Pi * FD * FD / 4.;
	FVolume = dx * __cons::Pi / 12 * (FD.head(n_nodes - 1).pow(2) + FD.tail(n_nodes - 1).pow(2)
		+ FD.head(n_nodes - 1) * FD.tail(n_nodes - 1));
	FMArea = FVolume / FXref;
	FDerLinArea = (FArea.tail(n_nodes - 1) - FArea.head(n_nodes - 1)) / Fdx;

	FU0.setZero(3, n_nodes);
	FNin = n_nodes;
}

void TPipe::setLeftBC(BoundaryCondition_ptr leftBC) {
	FLeftBC = move(leftBC);
}

void TPipe::setMethod(PipeMethod_ptr method) {
	FMethod = move(method);
}

void TPipe::setRightBC(BoundaryCondition_ptr rightBC) {
	FRightBC = move(rightBC);
}

void TPipe::setPtTtU(double p_t, double T_t, double u) {
	FMethod->setPtTtU(p_t, T_t, u);
}

void TPipe::setPtTtU(const RowVector & p_t, const RowVector & T_t,
	const RowVector & u) {
	FMethod->setPtTtU(p_t, T_t, u);
}


void TPipe::setPTU(double p, double T, double u) {
	FMethod->setPTU(p, T, u);
}

void TPipe::setPTU(const RowVector& p, const RowVector& T, const RowVector& u) {
	FMethod->setPTU(p, T, u);
}

void TPipe::setPTU(const RowVector& p, const RowVector& T, const RowVector& u,
	const RowVector& x) {
	auto p_ = linear_interp(x, p, Fx_sv);
	auto T_ = linear_interp(x, T, Fx_sv);
	auto u_ = linear_interp(x, u, Fx_sv);
	FMethod->setPTU(p_, T_, u_);
}

void TPipe::Solve() {
	IntegrateWithoutUpdating();
	UpdateStateVector();
}

void TPipe::Solve(double t) {
	while (FCurrentTime <= t) {
		double dt = getMaxTimeStep();
		FFric = FPipeFriction->FrictionCoefficient(FRe);
		FQint = FHeatTransfer->Heat(Ftemperature.array(), FWallHT->getTwallint().array(), FRe.array(), WorkingFluid);
		FWallHT->AddInternalHeat(FQint * FTimeStep);
		setTimeStep(dt);
		Solve();
	}
}

void TPipe::StoreOutput(std::vector<std::vector<float>>& Output) const{
	double pos = 0.;
	for (const auto& out : FPipeOutput) {
		pos = out.Location;
		if (out.LeftPressure.ins) {
			Output.at(out.LeftPressure.ind).push_back(__units::PaToBar(getLeftPressure(pos)));
		}
		if (out.MassFlow.ins) {
			Output.at(out.MassFlow.ind).push_back(getMassFlow(pos));
		}
		if (out.Pressure.ins) {
			Output.at(out.Pressure.ind).push_back(__units::PaToBar(getPressure(pos)));
		}
		if (out.RightPressure.ins) {
			Output.at(out.RightPressure.ind).push_back(__units::PaToBar(getRightPressure(pos)));
		}
		if (out.Temperature.ins) {
			Output.at(out.Temperature.ind).push_back(__units::KTodegC(getTemperature(pos)));
		}
		if (out.TotalPressure.ins) {
			Output.at(out.TotalPressure.ind).push_back(__units::PaToBar(getTotalPressure(pos)));
		}
		if (out.TotalTemperature.ins) {
			Output.at(out.TotalTemperature.ind).push_back(__units::KTodegC(getTotalTemperature(pos)));
		}
		if (out.Velocity.ins) {
			Output.at(out.Velocity.ind).push_back(getSpeed(pos));
		}
	}
}

void TPipe::UpdateStateVector() {
	std::lock_guard<std::mutex> lock(mtx);
	FMethod->UpdateComposition();
	FU0 = FU1;
	FCurrentTime += FTimeStep;
	FMethod->UpdateFlowVariables();
	FIsIntegrated = true;
}

void TPipe::WriteInsHeader(std::stringstream& output) const {
	std::string base_label;
	std::string label;
	std::ostringstream TextDist;
	TextDist.precision(8);

	for (auto result: InsResults) {
		TextDist.str("");
		TextDist.clear();
		TextDist << result.Distance;
		base_label = "\t" + PutLabel(5004) + "/"
			+ std::to_string(getPipeNumber()) + "/" + TextDist.str() + "m/";
		if (result.Cp == true) {
			label = base_label + PutLabel(4031) + PutLabel(927);
			output << label.c_str();
		}
		if (result.Cv == true) {
			label = base_label + PutLabel(4032) + PutLabel(927);
			output << label.c_str();
		}
		if (result.Gamma == true) {
			label = base_label + PutLabel(4017) + PutLabel(901);
			output << label.c_str();
		}
		if (result.LeftPressure == true) {
			label = base_label + PutLabel(4006) + PutLabel(908)
				+ "/" + PutLabel(3030);
			output << label.c_str();
		}
		if (result.MassFlow == true) {
			label = base_label + PutLabel(4008) + PutLabel(904);
			output << label.c_str();
		}
		if (result.Pressure == true) {
			label = base_label + PutLabel(4006) + PutLabel(908)
				+ "/" + PutLabel(3042);
			output << label.c_str();
		}
		if (result.R == true) {
			label = base_label + PutLabel(4033) + PutLabel(927);
			output << label.c_str();
		}
		if (result.RightPressure == true) {
			label = base_label + PutLabel(4006) + PutLabel(908)
				+ "/" + PutLabel(3031);
			output << label.c_str();
		}
		if (result.Temperature == true) {
			label = base_label + PutLabel(4005) + PutLabel(910)
				+ "/" + PutLabel(3042);
			output << label.c_str();
		}
		if (result.TotalPressure == true) {
			label = base_label + PutLabel(4006) + PutLabel(908)
				+ "/" + PutLabel(3026);
			output << label.c_str();
		}
		if (result.TotalTemperature == true) {
			label = base_label + PutLabel(4005) + PutLabel(910)
				+ "/" + PutLabel(3026);
			output << label.c_str();
		}
		if (result.Velocity == true) {
			label = base_label + PutLabel(4007) + PutLabel(909);
			output << label.c_str();
		}
	}
}

void TPipe::WriteInsHeader(std::vector<std::string>& output) const {
	std::string base_label;
	std::string label;
	std::ostringstream TextDist;
	TextDist.precision(8);

	for (auto result : InsResults) {
		TextDist.str("");
		TextDist.clear();
		TextDist << result.Distance;
		base_label = "\t" + PutLabel(5004) + "/"
			+ std::to_string(getPipeNumber()) + "/" + TextDist.str() + "m/";
		if (result.MassFlow == true) {
			output.push_back(base_label + PutLabel(4008) + PutLabel(904));
		}
		if (result.Pressure == true) {
			output.push_back(base_label + PutLabel(4006) + PutLabel(908));
		}
		if (result.Temperature == true) {
			output.push_back(base_label + PutLabel(4005) + PutLabel(910));
		}
		if (result.TotalPressure == true) {
			output.push_back(base_label + PutLabel(3026) + PutLabel(4006) + PutLabel(908));
		}
		if (result.TotalTemperature == true) {
			output.push_back(base_label + PutLabel(3026) + PutLabel(4005) + PutLabel(910));
		}
		if (result.Velocity == true) {
			output.push_back(base_label + PutLabel(4007) + PutLabel(909));
		}
	}
}

void TPipe::WriteInsResults(std::stringstream& output) const {
	double pos = 0.;
	for (auto result: InsResults) {
		pos = result.Distance;
		if (result.Cp == true) {
			output << "\t" << getcp(pos);
		}
		if (result.Cv == true) {
			output << "\t" << getcp(pos) / getGamma(pos);
		}
		if (result.Gamma == true) {
			output << "\t" << getGamma(pos);
		}
		if (result.LeftPressure == true) {
			output << "\t" << __units::PaToBar(getLeftPressure(pos));
		}
		if (result.MassFlow == true) {
			output << "\t" << getMassFlow(pos);
		}
		if (result.Pressure == true) {
			output << "\t" << __units::PaToBar(getPressure(pos));
		}
		if (result.R == true) {
			output << "\t" << getR(pos);
		}
		if (result.RightPressure == true) {
			output << "\t" << __units::PaToBar(getRightPressure(pos));
		}
		if (result.Temperature == true) {
			output << "\t" << __units::KTodegC(getTemperature(pos));
		}
		if (result.TotalPressure == true) {
			output << "\t" << __units::PaToBar(getTotalPressure(pos));
		}
		if (result.TotalTemperature == true) {
			output << "\t" << __units::KTodegC(getTotalTemperature(pos));
		}
		if (result.Velocity == true) {
			output << "\t" << getSpeed(pos);
		}
	}
}

Pipe_ptr create_pipe(const xml_node& node,
	const TComponentArray_ptr& fluid,
	const std::map<string, TSolid*> &MDB)
{
	auto pipe = make_shared<TPipe>();
	pipe->setName(node.attribute("Name").as_string());
	if (node.attribute("Name")) {
		pipe->setName(node.attribute("Name").as_string());
	}
	else {
		pipe->setName(node.attribute("Pipe_ID").as_string());
	}

	pipe->FPipeNumber = IDtoInt(node.attribute("Pipe_ID").as_string()) - 1;

	pipe->FNodeLeft = IDtoInt(node.attribute("NodeL_ID").as_string());
	pipe->FNodeRight = IDtoInt(node.attribute("NodeR_ID").as_string());

	pipe->FNIdenticalPipes = GetAttributeAsInt(node, "ParallelPipes");

	auto node_geometry = GetNodeChild(node, "Pip_Geometry");
	auto node_num_method = GetNodeChild(node, "Pip_NumericalMethod");
	std::string numerical_method = node_num_method.attribute("Scheme").as_string();
	
	pipe->FMeshSize = GetAttributeAsDouble(node_num_method, "MeshSize");

	pipe->FNumStretch = CountNodes(node_geometry, "Geo_Stretch");

	pipe->FDStretch.resize(pipe->FNumStretch + 1);
	pipe->FLStretch.resize(pipe->FNumStretch + 1);

	pipe->FDStretch[0] = GetAttributeAsDouble(node_geometry, "Diameter");
	pipe->FLStretch[0] = 0.;

	auto id = 0;
	for(auto node_stretch = GetNodeChild(node_geometry, "Geo_Stretch");
		node_stretch; node_stretch = node_stretch.next_sibling("Geo_Stretch")) {
		id = GetAttributeAsInt(node_stretch, "Stretch_ID");
		if(id > pipe->FNumStretch) {
			std::cout << "ERROR: The stretches of the pipe are not correctly ordered in pipe: "
				<< pipe->FPipeNumber << std::endl;
			throw Exception("");
		}
		pipe->FLStretch[id] = GetAttributeAsDouble(node_stretch, "Length");
		pipe->FDStretch[id] = GetAttributeAsDouble(node_stretch, "Diameter");
	}

	ArrayXd LStr(pipe->FLStretch.size());
	LStr(0) = pipe->FLStretch(0);
	for(auto i = 1; i < pipe->FLStretch.size(); ++i) {
		LStr(i) = pipe->FLStretch(i) + LStr(i - 1);
	}

	pipe->setGeometry(LStr, pipe->FMeshSize, pipe->FDStretch);

	auto node_prop = GetNodeChild(node, "GasProperties");

	RowVector Y;
	RowVector Pini(1);
	RowVector Tini(1);
	RowVector Vmean(1);
	RowVector xini(1);
	Pini.setConstant(__units::BarToPa(__cons::PRef));
	Tini.setConstant(__cons::TRef);
	Vmean.setZero();
	xini.setZero();
	Y.setZero(fluid->size());
	if (node_prop.child("GasProps_Point")) {
		ReadGasProperties(node_prop, &Pini, &Tini, &Vmean, &Y, fluid, &xini);
	} else {
		ReadGasProperties(node_prop, Pini[0], Tini[0], Vmean[0], Y, fluid);
	}

	pipe->WorkingFluid.resize(pipe->FD.size() - 1);
	for (auto i = 0; i < pipe->WorkingFluid.size(); i++) {
		pipe->WorkingFluid[i] = make_shared<TFluid>(fluid);
		pipe->WorkingFluid[i]->SetComposition(Y);
	}

	set_numerical_method(pipe, numerical_method);
	set_null_source(pipe);

	pipe->setPTU(Pini, Tini, Vmean, xini);
	pipe->setCourant(GetAttributeAsDouble(node_num_method, "Courant"));
	pipe->ReadInsResults(node); 
	pipe->ReadOutput(node);


	// HEAT TRANSFER -------
	xml_node node_ht = GetNodeChild(node, "Pip_HeatTransfer");
	double IntHeatMult = GetAttributeAsDouble(node_ht, "IntMultiplier");
	string Type = node_ht.attribute("HT_Type").as_string();

	if (Type == "IntakePipe"){
		pipe->FHeatTransfer = make_unique<TPipeIntHeatIntake>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "ExhaustPipe"){
		pipe->FHeatTransfer = make_unique<TPipeIntHeatExhaust>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "IntakePort"){
		pipe->FHeatTransfer = make_unique<TPortIntHeatIntake>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "ExhaustPort"){
		pipe->FHeatTransfer = make_unique<TPortIntHeatExhaust>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else{
		cout << "Internal heat transfer in pipe not correctly defined" << endl;
	}

	string WallCalc = node_ht.attribute("WallCalculation").as_string();
	pipe->FWallHT = make_shared<TPipeHeatT>();
	pipe->FWallHT->ReadHeatTransferData(node, pipe->FNin, pipe->FXref, MDB);
	

	if (WallCalc != "Constant"){

		pipe->FWallHT->BuildMatrix(pipe->FDcell);
		xml_node node_eheat = GetNodeChild(node_ht, "Pht_External");
		double velocity = GetAttributeAsDouble(node_eheat, "Velocity");
		double extMult = GetAttributeAsDouble(node_eheat, "ExtMultiplier");
		double emis = GetAttributeAsDouble(node_eheat, "Emissivity");
		string Type2 = node_eheat.attribute("Type").as_string();
		TPipeExtHeat *extheatobj;
		if (Type2 == "AirCooled"){
			extheatobj = new TPipeExtHeatA(pipe->FNin, velocity, extMult, pipe->FWallHT->getExtDiameter().array(),
				pipe->FXref, emis);
		}
		else if (Type2 == "WaterCooled"){
			double Tw = GetAttributeAsDouble(node_eheat, "WaterTemp");
			extheatobj = new TPipeExtHeatW(pipe->FNin, velocity, extMult, pipe->FWallHT->getExtDiameter().array(),
				pipe->FXref, emis, Tw);
		}
		else if (Type2 == "Port"){

		}
		pipe->FWallHT->AddExtHeatObject(extheatobj);
	}

	xml_node node_friction = GetNodeChild(node, "Pip_Friction");
	string FricType = node_friction.attribute("Type").as_string();
	if (FricType == "Colebrook"){
		pipe->FPipeFriction = make_unique<TPipeFrictionColebrook>();
	}
	else{
		pipe->FPipeFriction = make_unique<TPipeFriction>();
	}
	pipe->FPipeFriction->ReadFrictionData(node);
	pipe->FPipeFriction->setRelativeRugosity(pipe->FDcell);
	pipe->FFric.setZero(pipe->FNin);

	return pipe;
}
