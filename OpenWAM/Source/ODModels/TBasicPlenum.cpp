/**
 * @file TBasicPlenum.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 30 de mar. de 2016
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
 * This file include the different methods to solve a basic 0D element.
 */
#include "TBasicPlenum.h"

TBasicPlenum::TBasicPlenum() {
	// TODO Auto-generated constructor stub
	//FMass_in->resize(0);
	//FEnth_in->resize(0);
	FCurrentTime = 0.;
}

TBasicPlenum::~TBasicPlenum() {
	// TODO Auto-generated destructor stub
}

void TBasicPlenum::AppendNewEntry(){
	stPlenumEntry entry;
	entry.p = FZeroD->getPressure();
	entry.T = FZeroD->getTemperature();
	entry.u = 0.;
	PlenumEntry.push_back(entry);
	FMass_in.push_back(0.);
	FEnth_in.push_back(0.);
	TFluid_ptr fluid;
	FFluid_in.push_back(fluid);

}

TComponentArray_ptr TBasicPlenum::getFluidComponents() const {
	return FWorkingFluid->GetComponents();
}

double TBasicPlenum::getMaxTimeStep()
{
	auto pmin = 1.e7;
	auto pmax = 0.;
	auto p = 1.;
	auto enthalpy = 0.;

	auto dt = 1.;

	//for (int i = 0; i < Boundaries.size(); i++) {
	//	Boundaries[i]->Flux(FCurrentTime, FTimeStep);
	//	enthalpy += Boundaries[i]->getEnthalpyFlow();
	//	p = Boundaries[i]->getVirtualPipe()->getPressure(Uint(0));
	//	if (p < pmin)
	//		pmin = p;
	//	p = Boundaries[i]->getVirtualPipe()->getTotalPressure(Uint(0));
	//	if (p > pmax)
	//		pmax = p;
	//}
	//if (enthalpy == 0.)
	//	return 1.;

	//auto rvcv = FVolume * FZeroD->getCv() / FZeroD->getR();
	//auto pdep = FZeroD->getPressure();
	//if(enthalpy < 0)
	//	dt = rvcv * (pmin - pdep) / enthalpy;
	//else
	//	dt = rvcv * (pmax - pdep) / enthalpy;

	//if (dt < 0) dt = 1e-6;
	return dt;
}

void TBasicPlenum::IntegrateWithoutUpdating(){

	for (int i = 0; i < Boundaries.size(); i++){
		Boundaries[i]->Flux(FCurrentTime, FTimeStep);

		FMass_in.at(i) = Boundaries[i]->getMassFlow() * FTimeStep;
		FEnth_in.at(i) = Boundaries[i]->getEnthalpyFlow() * FTimeStep / FMass_in.at(i);
		FFluid_in.at(i) = Boundaries[i]->getFluid();
	}

	//! Mass balance
	for (int i = 0; i < Boundaries.size(); i++){
		FWorkingFluid->AppendFluid(FFluid_in.at(i), FMass_in.at(i), FZeroD->getMass());
		FZeroD->AppendMass(FMass_in.at(i));
	}

	//! Energy balance
	FZeroD->SolveNewPressure(&FMass_in, &FEnth_in, make_shared<TFluidArray>(FFluid_in), FWorkingFluid);

}

void TBasicPlenum::set0DModel(T0DModel_ptr zerod){
	FZeroD = move(zerod);
}

void TBasicPlenum::setInitialConditions(double p, double T){

	double mcil = p * FVolume / __R::Air / T;
	FZeroD->Initialize(p, T, mcil, FVolume);

}

void TBasicPlenum::setFluid(const TComponentArray_ptr& com, RowVector Y){

	FWorkingFluid = make_shared<TFluid>(com);
	FWorkingFluid->SetComposition(Y);
}

//void TBasicPlenum::setTimeStep(double dt){
//	FTimeStep = dt;
//	FCurrentTime += dt;
//}


void TBasicPlenum::Solve(){

	IntegrateWithoutUpdating();
	UpdateStateVector();

}


void TBasicPlenum::Solve(double t){

	FCurrentTime = t;
	FTimeStep = FCurrentTime - FPreviousTime;

	IntegrateWithoutUpdating();
	UpdateStateVector();

}

void TBasicPlenum::setPtTtU(int i, double p, double T, double u){

	PlenumEntry[i].T = T - pow2(u) / 2 / FWorkingFluid->FunCp(T);
	PlenumEntry[i].p = p * pow(PlenumEntry[i].T / T, __Gamma::G9(FWorkingFluid->FunGamma(T)));
	PlenumEntry[i].u = u;
}

void TBasicPlenum::setPTU(int i, double p, double T, double u){
	PlenumEntry[i].p = p;
	PlenumEntry[i].T = T;
	PlenumEntry[i].u = u;
}


void TBasicPlenum::UpdateStateVector(){

	FZeroD->UpdateStateVector();
	FCurrentTime += FTimeStep;
}

void TBasicPlenum::HeaderForCycleOutput(std::vector<std::string>& output) const
{
	string base = "Pipe/" + FName + "/";
	if (FPlenumOutput.Mass.avg) {
		output.push_back(base + "Mass[kg]");
	}
	if (FPlenumOutput.Pressure.avg) {
		output.push_back(base + "Pressure[bar]");
	}
	if (FPlenumOutput.Temperature.avg) {
		output.push_back(base + "Temperature[degC]");
	}
	if (FPlenumOutput.Volume.avg) {
		output.push_back(base + "Volume[m**3]");
	}
}

void TBasicPlenum::HeaderForPlots(std::vector<std::string>& output) const
{
	string base = "Plenum/" + FName + "/";
	if (FPlenumOutput.Mass.ins) {
		FPlenumOutput.Mass.ind = output.size();
		output.push_back(base + "Mass[kg]");
	}
	if (FPlenumOutput.Pressure.ins) {
		FPlenumOutput.Pressure.ind = output.size();
		output.push_back(base + "Pressure[bar]");
	}
	if (FPlenumOutput.Temperature.ins) {
		FPlenumOutput.Temperature.ind = output.size();
		output.push_back(base + "Temperature[degC]");
	}
	if (FPlenumOutput.Volume.ins) {
		FPlenumOutput.Volume.ind = output.size();
		output.push_back(base + "Volume[m**3]");
	}
}

void TBasicPlenum::IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const
{
	double pos = 0.;
	Map<ArrayXf> dt(Output.at(0).data(), Output[0].size());
	double dt_sum = dt.sum();

	if (FPlenumOutput.Mass.avg) {
		Map<ArrayXf> m(Output.at(FPlenumOutput.Mass.ind).data(), Output.at(FPlenumOutput.Mass.ind).size());
		Cicle.push_back((m * dt).sum() / dt_sum);
	}
	if (FPlenumOutput.Pressure.avg) {
		Map<ArrayXf> m(Output.at(FPlenumOutput.Pressure.ind).data(), Output.at(FPlenumOutput.Pressure.ind).size());
		Cicle.push_back((m * dt).sum() / dt_sum);
	}
	if (FPlenumOutput.Temperature.avg) {
		Map<ArrayXf> m(Output.at(FPlenumOutput.Temperature.ind).data(), Output.at(FPlenumOutput.Temperature.ind).size());
		Cicle.push_back((m * dt).sum() / dt_sum);
	}
	if (FPlenumOutput.Volume.avg) {
		Map<ArrayXf> m(Output.at(FPlenumOutput.Volume.ind).data(), Output.at(FPlenumOutput.Volume.ind).size());
		Cicle.push_back((m * dt).sum() / dt_sum);
	}
}

void TBasicPlenum::StoreOutput(std::vector<std::vector<float>>& output) const
{
	if (FPlenumOutput.Mass.ins) {
		output.at(FPlenumOutput.Mass.ind).push_back(FZeroD->getMass());
	}
	if (FPlenumOutput.Pressure.ins) {
		output.at(FPlenumOutput.Pressure.ind).push_back(__units::PaToBar(FZeroD->getPressure()));
	}
	if (FPlenumOutput.Temperature.ins) {
		output.at(FPlenumOutput.Temperature.ind).push_back(__units::KTodegC(FZeroD->getTemperature()));
	}
	if (FPlenumOutput.Volume.ins) {
		output.at(FPlenumOutput.Volume.ind).push_back(FVolume);
	}
}

void TBasicPlenum::ReadOutput(const xml_node & node)
{
	string param;
	bool Avg;
	auto node_out = GetNodeChild(node, "Output");
	FPlenumOutput = PlenumOutput();
	for (auto node_par = GetNodeChild(node_out, "Parameter"); node_par;
		node_par = node_par.next_sibling("Parameter")) {
		Is_plot = true;
		param = node_par.attribute("Name").as_string();
		Avg = node_par.attribute("Average").as_bool();
		if (Avg)
			Is_cycleout = true;
		if (param == "Mass") {
			FPlenumOutput.Mass.ins = true;
			FPlenumOutput.Mass.avg = Avg;
		}
		else if (param == "Pressure") {
			FPlenumOutput.Pressure.ins = true;
			FPlenumOutput.Pressure.avg = Avg;
		}
		else if (param == "Temperature") {
			FPlenumOutput.Temperature.ins = true;
			FPlenumOutput.Temperature.avg = Avg;
		}
		else if(param == "Volume"){
			FPlenumOutput.Volume.ins = true;
			FPlenumOutput.Volume.avg = Avg;
		}
	}

}

void TBasicPlenum::ReadAverageOutputXML(xml_node node_cyl) {



}

void TBasicPlenum::HeaderAverageOutput(stringstream& medoutput) {



}

void TBasicPlenum::PrintAverageOutput(stringstream& medoutput) {


}

void TBasicPlenum::IntegrateAverageOutput(double TActual) {



}

void TBasicPlenum::SetAverageOutput() {



}

void TBasicPlenum::ReadInstantaneousOutputXML(xml_node node_cyl) {



}


void TBasicPlenum::HeaderInstantaneousOutput(stringstream& insoutput) {



}


void TBasicPlenum::PrintInstantaneuosOutput(stringstream& insoutput) {



}


void TBasicPlenum::SetInstantaneousOutput() {



}

BasicPlenum_ptr create_basicplenum(const xml_node& node,
	const TComponentArray_ptr& fluid){

	auto plenum = make_shared<TBasicPlenum>();

	plenum->FBasicPlenum_ID = IDtoInt(node.attribute("ID").as_string()) - 1;

	if (node.attribute("Name")){
		plenum->FName = node.attribute("Name").as_string();
	}
	else{
		plenum->FName = node.attribute("ID").as_string();
	}

	plenum->FVolume = GetAttributeAsDouble(node, "Volume");

	auto node_prop = GetNodeChild(node, "GasProperties");

	double Pini, Tini, Vmean;
	RowVector Y;
	Y.setZero(fluid->size());
	ReadGasProperties(node_prop, Pini, Tini, Vmean, Y, fluid);

	plenum->FWorkingFluid = make_shared<TFluid>(fluid);
	plenum->FWorkingFluid->SetComposition(Y);

	unique_ptr<T0DModel> ZD(new T0DModel());
	double Mini = Pini * plenum->FVolume / plenum->FWorkingFluid->FunR() / Tini;
	ZD->Initialize(Pini, Tini, Mini, plenum->FVolume);
	plenum->FZeroD = std::move(ZD);

	return plenum;

}

PlenumOutput::PlenumOutput()
{
	Mass.avg = false;				//!< true if store average mass
	Mass.ins = false;				//!< true if store instantaneous mass
	Mass.ind = 0;				//!< Column index where instantaneous mass is stored
	Pressure.avg = false;			//!< true if store average pressure
	Pressure.ins = false;			//!< true if store instantanous pressure
	Pressure.ind = 0;			//!< Column index where instantaneous pressure is stored
	Temperature.avg = false;		//!< true if store averaged temperature
	Temperature.ins = false;		//!< true if store instantaneous temperature
	Temperature.ind = 0;		//!< Column index where instantaneous temperature is stored
	Volume.avg = false;			//!< true if store averaged volume
	Volume.ins = false;			//!< true if store instantaneous volume
	Volume.ind = 0;			//!< Column index where instantaneous volume is stored
}
