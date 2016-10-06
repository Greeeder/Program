/**
 * @file TCylinder.cpp
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
 * This file include the different methods to solve the thermodynamics in the cylinder.
 */
#include "TCylinder.h"
#include "HeatTransferCombMECSwirl.h"

TCylinder::TCylinder() {

	FCycle = 0;

	FInj = false;
	FComb = false;
	FFirstCycle = true;

}

TCylinder::~TCylinder() {
	// TODO Auto-generated destructor stub
}

void TCylinder::IntegrateWithoutUpdating() {

	FCrankMech->setTime(FTimeStep);
	//FCrankMech->setCurrentAngle(FCrankMech->getPreviousAngle() + FAngleStep);

	for (int i = 0; i < Boundaries.size(); i++) {
		Boundaries[i]->Flux(FCurrentTime, FTimeStep);

		FMass_in.at(i) = Boundaries[i]->getMassFlow() * FTimeStep;
		FEnth_in.at(i) = Boundaries[i]->getEnthalpyFlow() * FTimeStep / FMass_in.at(i);
		FFluid_in.at(i) = Boundaries[i]->getFluid();
	}

	//! Mass balance
	for (int i = 0; i < Boundaries.size(); i++) {
		FWorkingFluid->AppendFluid(FFluid_in.at(i), FMass_in.at(i), FZeroD->getMass());
		FZeroD->AppendMass(FMass_in.at(i));
	}

	double dfql = 0.;

	if (!FFirstCycle) {
		//! Fuel injected
		double MFuelIns = FInjector->getFuelRate(__units::RadToDeg(FCombAngle), FCurrentTime) * FTimeStep;
		if (MFuelIns > 0) {
			FInj = true;
			FComb = true;
			FWorkingFluid->AppendFluid(FInjector->getFluid(), MFuelIns, FZeroD->getMass());
			FZeroD->AppendMass(MFuelIns);
		}
		else {
			FInj = false;
		}

		if (FComb) {
			//! Heat released
			double MFuel = FInjector->getFuelInjected();
			FHRL1 = FCombustion->getHRL(__units::RadToDeg(FCombAngle));
			dfql = (FHRL1 - FHRL0) * MFuel * FFuel->FunHV() * FTimeStep;
		}
	}

	//! Energy balance

	HeatTransfer::stInicializarPorAngulo init;
	init.ang = FCrankMech->getCurrentAngle();
	init.areacil = FCrankMech->getAreaCyl();
	init.ciclo_cerrado = true;
	init.deltat = FTimeStep;
	init.VCIL = FCrankMech->getVolume();

	double vol = FCrankMech->getVolume(FCrankMech->getCurrentAngle() + FAngleStep, FZeroD->getPressure());

	dynamic_cast<HeatTransfer *>(FZeroD->getHeatTransfer())->InicializarPorAngulo(init);

	FZeroD->SolveNewPressure(vol, dfql, &FMass_in, &FEnth_in, make_shared<TFluidArray>(FFluid_in), FWorkingFluid);

}

void TCylinder::setAngle(double angle) {

	int cycle = floor((angle + FPhaseDiference) / FAngleCycle);
	double ang0 = angle + FPhaseDiference - cycle * FAngleCycle;
	if (cycle > FCycle) {
		FCrankMech->setPreviousAngle(FCrankMech->getCurrentAngle() - FAngleCycle);
	}
	if (ang0 >= FIVCAngle && FCrankMech->getPreviousAngle() < FIVCAngle) {
		FClosedLoop = true;

		HeatTransfer::stInicializarAlRCA init;
		init.CW1 = 0.;
		init.CW2 = 0.;
		init.PCA = FZeroD->getPressure();
		init.TCA = FZeroD->getTemperature();
		init.VCA = FCrankMech->getVolume();
		init.TCIL = 500;
		init.TCUL = 500;
		init.TPIS = 500;
		dynamic_cast<HeatTransfer *>(FZeroD->getHeatTransfer())->InicializarAlRCA(init);
	}
	if (ang0 >= FEVOAngle && FCrankMech->getPreviousAngle() < FEVOAngle) {
		FClosedLoop = false;
	}

	FCrankMech->setCurrentAngle(ang0);

}

void TCylinder::setCombustion(TCombustion_ptr Comb) {
	//if (Comb->getMode() == "Wiebe"){
	//	FCombustion = make_unique<TCombustion>(new TWiebe(dynamic_cast<TWiebe*>(Comb)));
	//}
	//else if (Comb->getMode() == "MultiWiebe"){
	//	FCombustion = make_unique<TCombustion>(new TMultiWiebe(dynamic_cast<TMultiWiebe*>(Comb)));
	//}
	//else if (Comb->getMode() == "HRL"){
	//	FCombustion = make_unique<TCombustion>(new THRL(dynamic_cast<THRL*>(Comb)));
	//}
	//else if (Comb->getMode() == "CombustionModel"){
	//	FCombustion = make_unique<TCombustion>(new TCombustionModel());
	//}
	FCombustion = move(Comb);
}

void TCylinder::setCrankMechanism(TCrankMechanism_ptr CrMch) {
	FCrankMech = move(CrMch);
}

void TCylinder::setDeltaAngle(double da)
{
	FAngleStep = da;
}

void TCylinder::setInitialConditions(double p, double T, string EngineType) {
	double vcil = FCrankMech->getVolumeStatic(FCrankMech->getCurrentAngle());
	double v = FCrankMech->getVolumeStatic(__cons::Pi);
	double pcil = 1;
	double mcil = 0;
	double Tcil = 300;
	if (EngineType == "4Strokes") {
		if (FCrankMech->getCurrentAngle() < __cons::Pi || FCrankMech->getCurrentAngle() > __cons::Pi * 3) {
			mcil = p * v / __R::Air / T;
			pcil = p * pow(v / vcil, __Gamma::G);
			Tcil = pcil * vcil / __R::Air / mcil;
		}
		else {
			pcil = p;
			Tcil = T;
			mcil = pcil * vcil / __R::Air / Tcil;
		}
		FZeroD->Initialize(pcil, Tcil, mcil, vcil);
	}
	else if (EngineType == "2Strokes") {

	}
	FCrankMech->setVolume(vcil);
	FCrankMech->setVolume0(vcil);
}

void TCylinder::setEngineSpeed(double n) {
	FCrankMech->setRotationalSpeed(n);
}

void TCylinder::setFluid(const TComponentArray_ptr& com, RowVector Y) {

	FWorkingFluid = make_shared<TFluid>(com);
	FWorkingFluid->SetComposition(Y);
}

void TCylinder::setTimeStep(double dt) {
	FTimeStep = dt;
}


void TCylinder::Solve() {

	IntegrateWithoutUpdating();
	UpdateStateVector();
}


void TCylinder::Solve(double t) {

	FCurrentTime = t;
	FTimeStep = FCurrentTime - FPreviousTime;

}


void TCylinder::UpdateStateVector() {

	FNetTorque = (FZeroD->getPressureStep() - __Ambient::p_Pa) * FCrankMech->getProjectedArea(FCrankMech->getCurrentAngle());
	FWork = FZeroD->getdWi();

	FZeroD->UpdateStateVector();
	FHRL0 = FHRL1;
	FCombAngle += FAngleStep;
	if (FCombAngle > FAngleCycle / 2.) {
		FCombAngle -= FAngleCycle;
		FInj = false;
		FComb = false;
		FFirstCycle = false;
	}

	FCrankMech->UpdateStateData();

	FCurrentTime += FTimeStep;
}

void TCylinder::HeaderForCycleOutput(std::vector<std::string>& output) const
{
	string base = "Cylinder/" + FName + "/";

	if (FCylinderOutput.Mass.avg) {
		output.push_back(base + "TrapedMass[mg]");
	}
	if (FCylinderOutput.IntakeMass.avg) {
		output.push_back(base + "IntakeMass[mg]");
	}
	if (FCylinderOutput.ExhaustMass.avg) {
		output.push_back(base + "ExhaustMass[mg]");
	}
	if (FCylinderOutput.FuelMass.avg) {
		output.push_back(base + "FuelMass[mg]");
	}
	if (FCylinderOutput.BlowBy.avg) {
		output.push_back(base + "Blow-By[mg]");
	}
	if (FCylinderOutput.HeatPower.avg) {
		output.push_back(base + "HeatExhanged[J]");
	}
	if (FCylinderOutput.HeatReleased.avg) {
		output.push_back(base + "HeatReleased[J]");
	}
}

void TCylinder::HeaderForPlots(std::vector<std::string>& output) const
{
	string base = "Cylinder/" + FName + "/";
	if (FCylinderOutput.Angle.ins) {
		FCylinderOutput.Angle.ind = output.size();
		output.push_back(base + "Angle[deg]");
	}
	if (FCylinderOutput.Pressure.ins) {
		FCylinderOutput.Pressure.ind = output.size();
		output.push_back(base + "Pressure[bar]");
	}
	if (FCylinderOutput.Temperature.ins) {
		FCylinderOutput.Temperature.ind = output.size();
		output.push_back(base + "Temperature[degC]");
	}
	if (FCylinderOutput.Mass.ins) {
		FCylinderOutput.Mass.ind = output.size();
		output.push_back(base + "Mass[mg]");
	}
	if (FCylinderOutput.Volume.ins) {
		FCylinderOutput.Volume.ind = output.size();
		output.push_back(base + "Volume[l]");
	}
	if (FCylinderOutput.IntakeMass.ins) {
		FCylinderOutput.IntakeMass.ind = output.size();
		output.push_back(base + "IntakeMass[g/s]");
	}
	if (FCylinderOutput.ExhaustMass.ins) {
		FCylinderOutput.ExhaustMass.ind = output.size();
		output.push_back(base + "ExhaustMass[g/s]");
	}
	if (FCylinderOutput.FuelMass.ins) {
		FCylinderOutput.FuelMass.ind = output.size();
		output.push_back(base + "FuelMassRate[g/s]");
	}
	if (FCylinderOutput.BlowBy.ins) {
		FCylinderOutput.BlowBy.ind = output.size();
		output.push_back(base + "Blow-By[g/s]");
	}
	if (FCylinderOutput.HeatPower.ins) {
		FCylinderOutput.HeatPower.ind = output.size();
		output.push_back(base + "HeatPowerExhanged[W]");
	}
	if (FCylinderOutput.HeatReleased.ins) {
		FCylinderOutput.Pressure.ind = output.size();
		output.push_back(base + "HeatPowerReleased[W]");
	}

}

void TCylinder::IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const
{
	Map<ArrayXf> dt(Output.at(0).data(), Output[0].size());
	double dt_sum = dt.sum();

	if (FCylinderOutput.Mass.avg) {
		//!< \todo Get the trapped mass
		Map<ArrayXf> m(Output.at(FCylinderOutput.Mass.ind).data(), Output.at(FCylinderOutput.Mass.ind).size());
		Cicle.push_back((m * dt).sum() / dt_sum);
	}
	if (FCylinderOutput.IntakeMass.avg) {
		//!< \todo Units?
		Map<ArrayXf> mi(Output.at(FCylinderOutput.IntakeMass.ind).data(), Output.at(FCylinderOutput.IntakeMass.ind).size());
		Cicle.push_back((mi * dt).sum());
	}
	if (FCylinderOutput.ExhaustMass.avg) {
		Map<ArrayXf> me(Output.at(FCylinderOutput.ExhaustMass.ind).data(), Output.at(FCylinderOutput.ExhaustMass.ind).size());
		Cicle.push_back((me * dt).sum());
	}
	if (FCylinderOutput.FuelMass.avg) {
		Map<ArrayXf> fm(Output.at(FCylinderOutput.FuelMass.ind).data(), Output.at(FCylinderOutput.FuelMass.ind).size());
		Cicle.push_back((fm * dt).sum());
	}
	if (FCylinderOutput.BlowBy.avg) {
		Map<ArrayXf> bb(Output.at(FCylinderOutput.BlowBy.ind).data(), Output.at(FCylinderOutput.BlowBy.ind).size());
		Cicle.push_back((bb * dt).sum());
	}
	if (FCylinderOutput.HeatPower.avg) {
		Map<ArrayXf> hp(Output.at(FCylinderOutput.HeatPower.ind).data(), Output.at(FCylinderOutput.HeatPower.ind).size());
		Cicle.push_back((hp * dt).sum());
	}
	if (FCylinderOutput.HeatReleased.avg) {
		Map<ArrayXf> hr(Output.at(FCylinderOutput.HeatReleased.ind).data(), Output.at(FCylinderOutput.HeatReleased.ind).size());
		Cicle.push_back((hr * dt).sum());
	}

}

void TCylinder::StoreOutput(std::vector<std::vector<float>>& output) const
{
	if (FCylinderOutput.Angle.ins) {
		output.at(FCylinderOutput.Mass.ind).push_back(__units::DegToRad(FCrankMech->getCurrentAngle()));
	}
	if (FCylinderOutput.Pressure.ins) {
		output.at(FCylinderOutput.Pressure.ind).push_back(__units::PaToBar(FZeroD->getPressure()));
	}
	if (FCylinderOutput.Temperature.ins) {
		output.at(FCylinderOutput.Temperature.ind).push_back(__units::KTodegC(FZeroD->getTemperature()));
	}
	if (FCylinderOutput.Mass.ins) {
		output.at(FCylinderOutput.Mass.ind).push_back(__units::From_kilo(FZeroD->getMass()));
	}
	if (FCylinderOutput.Volume.ins) {
		output.at(FCylinderOutput.Volume.ind).push_back(__units::From_kilo(FCrankMech->getVolume()));
	}
	if (FCylinderOutput.IntakeMass.ins) {
		double m = IntakeValve.at(0)->getMassFlow();
		for (int i = 1; i < IntakeValve.size(); i++) {
			m += IntakeValve.at(0)->getMassFlow();
		}
		output.at(FCylinderOutput.IntakeMass.ind).push_back(m);
	}
	if (FCylinderOutput.ExhaustMass.ins) {
		double m = ExhaustValve.at(0)->getMassFlow();
		for (int i = 1; i < ExhaustValve.size(); i++) {
			m += ExhaustValve.at(0)->getMassFlow();
		}
		output.at(FCylinderOutput.ExhaustMass.ind).push_back(0);
	}
	if (FCylinderOutput.FuelMass.ins) {
		output.at(FCylinderOutput.FuelMass.ind).push_back(FInjector->getFuelRate());
	}
	if (FCylinderOutput.BlowBy.ins) {
		output.at(FCylinderOutput.ExhaustMass.ind).push_back(0);
	}
	if (FCylinderOutput.HeatPower.ins) {
		output.at(FCylinderOutput.HeatPower.ind).push_back(0);
	}
	if (FCylinderOutput.HeatReleased.ins) {
		output.at(FCylinderOutput.HeatReleased.ind).push_back(0);
	}

}

void TCylinder::ReadOutput(const xml_node & node)
{

	string param;
	bool Avg;
	FCylinderOutput = CylinderOutput();
	auto node_out = GetNodeChild(node, "Output");
	if (node_out)
		FCylinderOutput.Angle.ins = true;
	for (auto node_par = GetNodeChild(node_out, "Parameter"); node_par;
		node_par = node_par.next_sibling("Parameter")) {
		Is_plot = true;
		param = node_par.attribute("Name").as_string();
		Avg = node_par.attribute("Average").as_bool();
		if (Avg)
			Is_cycleout = true;
		if (param == "Pressure") {
			FCylinderOutput.Pressure.ins = true;
		}
		else if (param == "Temperature") {
			FCylinderOutput.Temperature.ins = true;
		}
		else if (param == "Mass") {
			FCylinderOutput.Mass.ins = true;
			FCylinderOutput.Mass.avg = Avg;
		}
		else if (param == "Volume") {
			FCylinderOutput.Volume.ins = true;
		}
		else if (param == "IntakeMass") {
			FCylinderOutput.IntakeMass.ins = true;
			FCylinderOutput.IntakeMass.avg = Avg;
		}
		else if (param == "ExhaustMass") {
			FCylinderOutput.ExhaustMass.ins = true;
			FCylinderOutput.ExhaustMass.avg = Avg;
		}
		else if (param == "FuelMass") {
			FCylinderOutput.FuelMass.ins = true;
			FCylinderOutput.FuelMass.avg = Avg;
		}
		else if (param == "BlowBy") {
			FCylinderOutput.BlowBy.ins = true;
			FCylinderOutput.BlowBy.avg = Avg;
		}
		else if (param == "HeatTransfer") {
			FCylinderOutput.HeatPower.ins = true;
			FCylinderOutput.HeatPower.avg = Avg;
		}
		else if (param == "HeatReleased") {
			FCylinderOutput.HeatReleased.ins = true;
			FCylinderOutput.HeatReleased.avg = Avg;
		}
	}

}

void TCylinder::ReadAverageOutputXML(xml_node node_cyl) {


	FAVGOutput.TrabajoNeto = false;
	FAVGOutput.TrabajoNetoMED = 0.;
	FAVGOutput.TrabajoNetoSUM = 0.;
	FAVGOutput.PresionMediaNeta = false;
	FAVGOutput.PresionMediaNetaMED = 0.;
	FAVGOutput.TrabajoBombeo = false;
	FAVGOutput.TrabajoBombeoMED = 0.;
	FAVGOutput.TrabajoBombeoSUM = 0.;
	FAVGOutput.PresionMediaBombeo = false;
	FAVGOutput.PresionMediaBombeoMED = 0.;
	FAVGOutput.CalorCombustion = false;
	FAVGOutput.CalorCombustionMED = 0.;
	FAVGOutput.CalorCombustionSUM = 0.;
	FAVGOutput.CalorCilindro = false;
	FAVGOutput.CalorCilindroMED = 0.;
	FAVGOutput.CalorCilindroSUM = 0.;
	FAVGOutput.CalorCulata = false;
	FAVGOutput.CalorCulataMED = 0.;
	FAVGOutput.CalorCulataSUM = 0.;
	FAVGOutput.CalorPiston = false;
	FAVGOutput.CalorPistonMED = 0.;
	FAVGOutput.CalorPistonSUM = 0.;
	FAVGOutput.PresionMediaIndicada = false;
	FAVGOutput.PresionMediaIndicadaMED = 0.;
	FAVGOutput.MasaAtrapada = false;
	FAVGOutput.MasaAtrapadaMED = 0.;
	FAVGOutput.TemperaturaCilindroInterna = false;
	FAVGOutput.TemperaturaCilindroInternaMED = 0.;
	FAVGOutput.TemperaturaCilindroInternaSUM = 0.;
	FAVGOutput.TemperaturaCilindroMedia = false;
	FAVGOutput.TemperaturaCilindroMediaMED = 0.;
	FAVGOutput.TemperaturaCilindroMediaSUM = 0.;
	FAVGOutput.TemperaturaCilindroExterna = false;
	FAVGOutput.TemperaturaCilindroExternaMED = 0.;
	FAVGOutput.TemperaturaCilindroExternaSUM = 0.;
	FAVGOutput.TemperaturaPistonInterna = false;
	FAVGOutput.TemperaturaPistonInternaMED = 0.;
	FAVGOutput.TemperaturaPistonInternaSUM = 0.;
	FAVGOutput.TemperaturaPistonMedia = false;
	FAVGOutput.TemperaturaPistonMediaMED = 0.;
	FAVGOutput.TemperaturaPistonMediaSUM = 0.;
	FAVGOutput.TemperaturaPistonExterna = false;
	FAVGOutput.TemperaturaPistonExternaMED = 0.;
	FAVGOutput.TemperaturaPistonExternaSUM = 0.;
	FAVGOutput.TemperaturaCulataInterna = false;
	FAVGOutput.TemperaturaCulataInternaMED = 0.;
	FAVGOutput.TemperaturaCulataInternaSUM = 0.;
	FAVGOutput.TemperaturaCulataMedia = false;
	FAVGOutput.TemperaturaCulataMediaMED = 0.;
	FAVGOutput.TemperaturaCulataMediaSUM = 0.;
	FAVGOutput.TemperaturaCulataExterna = false;
	FAVGOutput.TemperaturaCulataExternaMED = 0.;
	FAVGOutput.TemperaturaCulataExternaSUM = 0.;
	FAVGOutput.NITMedio = false;
	FAVGOutput.NITMedioMED = 0.;
	FAVGOutput.NITMedioSUM = 0.;
	FAVGOutput.AFRMedio = false;
	FAVGOutput.AFRMedioMED = 0.;
	FAVGOutput.MasaBlowBy = false;
	FAVGOutput.MasaBlowByMED = 0.;
	FAVGOutput.MasaBlowBySUM = 0.;
	FAVGOutput.MasaAdmision = false;
	FAVGOutput.MasaAdmisionMED = 0.;
	FAVGOutput.MasaEscape = false;
	FAVGOutput.MasaEscapeMED = 0.;
	FAVGOutput.TemperaturaMedia = false;
	FAVGOutput.TemperaturaMediaMED = 0.;
	FAVGOutput.TemperaturaMediaSUM = 0.;
	FAVGOutput.Swirl = false;
	FAVGOutput.SwirlMED = 0.;
	FAVGOutput.RendVolumetrico = false;
	FAVGOutput.RendVolumetricoMED = 0.;
	FAVGOutput.DensidadReferenciaSUM = 0.;
	FAVGOutput.MasaCortocircuito = false;
	FAVGOutput.MasaCortocircuitoMED = 0.;
	FAVGOutput.MasaCortocircuitoSUM = 0.;

	FAVGOutput.Tiempo0 = 0.;
	FAVGOutput.TiempoSUM = 0.;

	FPrintAVG = false;

	xml_node node_avgout = GetNodeChild(node_cyl, "Cyl_AvgOutput");

	if (node_avgout) {
		FPrintAVG = true;
		for (xml_attribute parameter = node_avgout.attribute("Parameter"); parameter; parameter = parameter.next_attribute()) {

			if (parameter.value() == "NetWork") {
				FAVGOutput.TrabajoNeto = true;
			}
			else if (parameter.value() == "NMEP") {
				FAVGOutput.PresionMediaNeta = true;
			}
			else if (parameter.value() == "PumpingWork") {
				FAVGOutput.TrabajoBombeo = true;
			}
			else if (parameter.value() == "PMEP") {
				FAVGOutput.PresionMediaBombeo = true;
			}
			else if (parameter.value() == "CombustionHeat") {
				FAVGOutput.CalorCombustion = true;
			}
			else if (parameter.value() == "LinerHeat") {
				FAVGOutput.CalorCilindro = true;
			}
			else if (parameter.value() == "CylHeadHeat") {
				FAVGOutput.CalorCulata = true;
			}
			else if (parameter.value() == "PistonHeat") {
				FAVGOutput.CalorPiston = true;
			}
			else if (parameter.value() == "IMEP") {
				FAVGOutput.PresionMediaIndicada = true;
			}
			else if (parameter.value() == "TrappedMass") {
				FAVGOutput.MasaAtrapada = true;
			}
			else if (parameter.value() == "IntCylTemperature") {
				FAVGOutput.TemperaturaCilindroInterna = true;
			}
			else if (parameter.value() == "MedCylTemperature") {
				FAVGOutput.TemperaturaCilindroMedia = true;
			}
			else if (parameter.value() == "ExtCylTemperature") {
				FAVGOutput.TemperaturaCilindroExterna = true;
			}
			else if (parameter.value() == "IntPistonTemperature") {
				FAVGOutput.TemperaturaPistonInterna = true;
			}
			else if (parameter.value() == "MedPistonTemperature") {
				FAVGOutput.TemperaturaPistonMedia = true;
			}
			else if (parameter.value() == "ExtPistonTemperature") {
				FAVGOutput.TemperaturaPistonExterna = true;
			}
			else if (parameter.value() == "IntCylHeadTemperature") {
				FAVGOutput.TemperaturaCulataInterna = true;
			}
			else if (parameter.value() == "MedCylHeadTemperature") {
				FAVGOutput.TemperaturaCulataMedia = true;
			}
			else if (parameter.value() == "ExtCylHeadTemperature") {
				FAVGOutput.TemperaturaCulataExterna = true;
			}
			else if (parameter.value() == "NIT") {
				FAVGOutput.NITMedio = true;
			}
			else if (parameter.value() == "AFR") {
				FAVGOutput.AFRMedio = true;
			}
			else if (parameter.value() == "BlowByMass") {
				FAVGOutput.MasaBlowBy = true;
			}
			else if (parameter.value() == "IntakeMass") {
				FAVGOutput.MasaAdmision = true;
			}
			else if (parameter.value() == "ExhaustMass") {
				FAVGOutput.MasaEscape = true;
			}
			else if (parameter.value() == "ShortCircuitMass") {
				FAVGOutput.MasaCortocircuito = true;
			}
			else if (parameter.value() == "Temperature") {
				FAVGOutput.TemperaturaMedia = true;
			}
			else if (parameter.value() == "Swirl") {
				FAVGOutput.Swirl = true;
			}
			else if (parameter.value() == "VolumetricEfficiency") {
				FAVGOutput.RendVolumetrico = true;
			}
			else {
				std::cout << "Resultados medios en cilindro " << FCylinder_ID << " no implementado " << std::endl;
			}
		}
	}

}

void TCylinder::HeaderAverageOutput(stringstream& medoutput) {

	std::string Label;

	if (FPrintAVG) {
		if (FAVGOutput.TrabajoNeto) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4021) + PutLabel(
				907) + "/" + PutLabel(3001);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.PresionMediaNeta) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4006) + PutLabel(
				908) + "/" + PutLabel(3002);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TrabajoBombeo) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4021) + PutLabel(
				907) + "/" + PutLabel(3003);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.PresionMediaBombeo) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4006) + PutLabel(
				908) + "/" + PutLabel(3004);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.CalorCombustion) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				907) + "/" + PutLabel(3005);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.CalorCilindro) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				907) + "/" + PutLabel(3006);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.CalorCulata) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				907) + "/" + PutLabel(3007);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.CalorPiston) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				907) + "/" + PutLabel(3008);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.PresionMediaIndicada) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4006) + PutLabel(
				908) + "/" + PutLabel(3009);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.MasaAtrapada) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3010);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCilindroInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3011);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCilindroMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3012);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCilindroExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3013);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaPistonInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3011);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaPistonMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3012);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaPistonExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3013);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCulataInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3011);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCulataMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3012);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaCulataExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3013);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.NITMedio) {
			//! \todo NIT in cylinders?

			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3014) + PutLabel(
			//		903) + "/" + PutLabel(3018) + "/" + std::to_string(i + 1);
			//	medoutput << Label.c_str();
			//}
			//Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + PutLabel(3014) + PutLabel(903) + "/" + PutLabel(
			//	3018) + "/" + PutLabel(3020);
			//medoutput << Label.c_str();
		}
		if (FAVGOutput.AFRMedio) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4015) + PutLabel(901);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.MasaBlowBy) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3016);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.MasaAdmision) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3017);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.MasaEscape) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3018);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.MasaCortocircuito) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3019);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.TemperaturaMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3020);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.Swirl) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3021) + PutLabel(901);
			medoutput << Label.c_str();
		}
		if (FAVGOutput.RendVolumetrico) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4011) + PutLabel(
				910) + "/" + PutLabel(3022);
			medoutput << Label.c_str();
		}
	}

}

void TCylinder::PrintAverageOutput(stringstream& medoutput) {

	if (FPrintAVG) {
		if (FAVGOutput.TrabajoNeto)
			medoutput << "\t" << FAVGOutput.TrabajoNetoMED;
		if (FAVGOutput.PresionMediaNeta)
			medoutput << "\t" << FAVGOutput.PresionMediaNetaMED;
		if (FAVGOutput.TrabajoBombeo)
			medoutput << "\t" << FAVGOutput.TrabajoBombeoMED;
		if (FAVGOutput.PresionMediaBombeo)
			medoutput << "\t" << FAVGOutput.PresionMediaBombeoMED;
		if (FAVGOutput.CalorCombustion)
			medoutput << "\t" << FAVGOutput.CalorCombustionMED;
		if (FAVGOutput.CalorCilindro)
			medoutput << "\t" << FAVGOutput.CalorCilindroMED;
		if (FAVGOutput.CalorCulata)
			medoutput << "\t" << FAVGOutput.CalorCulataMED;
		if (FAVGOutput.CalorPiston)
			medoutput << "\t" << FAVGOutput.CalorPistonMED;
		if (FAVGOutput.PresionMediaIndicada)
			medoutput << "\t" << FAVGOutput.PresionMediaIndicadaMED;
		if (FAVGOutput.MasaAtrapada)
			medoutput << "\t" << FAVGOutput.MasaAtrapadaMED;
		if (FAVGOutput.TemperaturaCilindroInterna)
			medoutput << "\t" << FAVGOutput.TemperaturaCilindroInternaMED;
		if (FAVGOutput.TemperaturaCilindroMedia)
			medoutput << "\t" << FAVGOutput.TemperaturaCilindroMediaMED;
		if (FAVGOutput.TemperaturaCilindroExterna)
			medoutput << "\t" << FAVGOutput.TemperaturaCilindroExternaMED;
		if (FAVGOutput.TemperaturaPistonInterna)
			medoutput << "\t" << FAVGOutput.TemperaturaPistonInternaMED;
		if (FAVGOutput.TemperaturaPistonMedia)
			medoutput << "\t" << FAVGOutput.TemperaturaPistonMediaMED;
		if (FAVGOutput.TemperaturaPistonExterna)
			medoutput << "\t" << FAVGOutput.TemperaturaPistonExternaMED;
		if (FAVGOutput.TemperaturaCulataInterna)
			medoutput << "\t" << FAVGOutput.TemperaturaCulataInternaMED;
		if (FAVGOutput.TemperaturaCulataMedia)
			medoutput << "\t" << FAVGOutput.TemperaturaCulataMediaMED;
		if (FAVGOutput.TemperaturaCulataExterna)
			medoutput << "\t" << FAVGOutput.TemperaturaCulataExternaMED;
		if (FAVGOutput.NITMedio) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	medoutput << "\t" << FAVGOutput.NITMED[i];
			//}
			//medoutput << "\t" << FAVGOutput.NITMedioMED;
		}
		if (FAVGOutput.AFRMedio)
			medoutput << "\t" << FAVGOutput.AFRMedioMED;
		if (FAVGOutput.MasaBlowBy)
			medoutput << "\t" << FAVGOutput.MasaBlowByMED;
		if (FAVGOutput.MasaAdmision)
			medoutput << "\t" << FAVGOutput.MasaAdmisionMED;
		if (FAVGOutput.MasaEscape)
			medoutput << "\t" << FAVGOutput.MasaEscapeMED;
		if (FAVGOutput.MasaCortocircuito)
			medoutput << "\t" << FAVGOutput.MasaCortocircuitoMED;
		if (FAVGOutput.TemperaturaMedia)
			medoutput << "\t" << FAVGOutput.TemperaturaMediaMED;
		if (FAVGOutput.Swirl)
			medoutput << "\t" << FAVGOutput.SwirlMED;
		if (FAVGOutput.RendVolumetrico)
			medoutput << "\t" << FAVGOutput.RendVolumetricoMED;
	}

}

void TCylinder::IntegrateAverageOutput(double TActual) {


	if (FPrintAVG) {
		double DeltaT = TActual - FAVGOutput.Tiempo0;
		// double DeltaAngulo=360.*FMotor->PutRegimen(360.*FMotor->getRegimen()/60.*DeltaT);

		FAVGOutput.TrabajoNetoSUM += FZeroD->getdWi();
		if (FCrankMech->getCurrentAngle() > 180. && FCrankMech->getCurrentAngle() < 540.) {
			FAVGOutput.TrabajoBombeoSUM += FZeroD->getdWi();
		}

		if (FAVGOutput.CalorCombustion && DeltaT > 0.)
			//FAVGOutput.CalorCombustionSUM += FCalor.Liberado;
			if (FAVGOutput.CalorCilindro)
				//FAVGOutput.CalorCilindroSUM += FCalor.TransCilindro * DeltaT;
				if (FAVGOutput.CalorCulata)
					//FAVGOutput.CalorCulataSUM += FCalor.TransCulata * DeltaT;
					if (FAVGOutput.CalorPiston)
						//FAVGOutput.CalorPistonSUM += FCalor.TransPiston * DeltaT;

						if (FAVGOutput.TemperaturaCilindroInterna)
							//FAVGOutput.TemperaturaCilindroInternaSUM += FTempPared[0].Cylinder * DeltaT;
							if (FAVGOutput.TemperaturaCilindroMedia)
								//FAVGOutput.TemperaturaCilindroMediaSUM += FTempPared[1].Cylinder * DeltaT;
								if (FAVGOutput.TemperaturaCilindroExterna)
									//FAVGOutput.TemperaturaCilindroExternaSUM += FTempPared[2].Cylinder * DeltaT;

									if (FAVGOutput.TemperaturaCulataInterna)
										//FAVGOutput.TemperaturaCulataInternaSUM += FTempPared[0].Culata * DeltaT;
										if (FAVGOutput.TemperaturaCulataMedia)
											//FAVGOutput.TemperaturaCulataMediaSUM += FTempPared[1].Culata * DeltaT;
											if (FAVGOutput.TemperaturaCulataExterna)
												//FAVGOutput.TemperaturaCulataExternaSUM += FTempPared[2].Culata * DeltaT;

												if (FAVGOutput.TemperaturaPistonInterna)
													//FAVGOutput.TemperaturaPistonInternaSUM += FTempPared[0].Piston * DeltaT;
													if (FAVGOutput.TemperaturaPistonMedia)
														//FAVGOutput.TemperaturaPistonMediaSUM += FTempPared[1].Piston * DeltaT;
														if (FAVGOutput.TemperaturaPistonExterna)
															//FAVGOutput.TemperaturaPistonExternaSUM += FTempPared[2].Piston * DeltaT;

															if (FAVGOutput.NITMedio)
																//FAVGOutput.NITMedioSUM += FNIT * DeltaT;
																if (FAVGOutput.MasaBlowBy)
																	//FAVGOutput.MasaBlowBySUM += FMasaBlowBy;
																	if (FAVGOutput.MasaCortocircuito)
																		//FAVGOutput.MasaCortocircuitoSUM += FMasaCortocircuito;
																		if (FAVGOutput.TemperaturaMedia)
																			FAVGOutput.TemperaturaMediaSUM += FZeroD->getTemperature();
		if (FAVGOutput.RendVolumetrico)
			//FAVGOutput.DensidadReferenciaSUM += FDensidadReferencia * FDeltaT;

			FAVGOutput.TiempoSUM += DeltaT;
		FAVGOutput.Tiempo0 = TActual;
	}
}

void TCylinder::SetAverageOutput() {

	if (FPrintAVG) {
		if (FAVGOutput.TrabajoNeto || FAVGOutput.PresionMediaNeta || FAVGOutput.PresionMediaIndicada) {
			FAVGOutput.TrabajoNetoMED = FAVGOutput.TrabajoNetoSUM;
			FAVGOutput.TrabajoNetoSUM = 0.;
		}
		if (FAVGOutput.PresionMediaNeta || FAVGOutput.PresionMediaIndicada) {
			FAVGOutput.PresionMediaNetaMED = __units::PaToBar(FAVGOutput.TrabajoNetoMED /
				FCrankMech->getVolumeDisplaced());
		}
		if (FAVGOutput.TrabajoBombeo || FAVGOutput.PresionMediaBombeo
			|| FAVGOutput.PresionMediaIndicada) {
			FAVGOutput.TrabajoBombeoMED = FAVGOutput.TrabajoBombeoSUM;
			FAVGOutput.TrabajoBombeoSUM = 0.;
		}
		if (FAVGOutput.PresionMediaBombeo || FAVGOutput.PresionMediaIndicada) {
			FAVGOutput.PresionMediaBombeoMED = __units::PaToBar(-FAVGOutput.TrabajoBombeoMED /
				FCrankMech->getVolumeDisplaced());
		}
		if (FAVGOutput.CalorCombustion) {
			FAVGOutput.CalorCombustionMED = FAVGOutput.CalorCombustionSUM;
			FAVGOutput.CalorCombustionSUM = 0.;
		}
		if (FAVGOutput.CalorCilindro) {
			FAVGOutput.CalorCilindroMED = FAVGOutput.CalorCilindroSUM;
			FAVGOutput.CalorCilindroSUM = 0.;
		}
		if (FAVGOutput.CalorCulata) {
			FAVGOutput.CalorCulataMED = FAVGOutput.CalorCulataSUM;
			FAVGOutput.CalorCulataSUM = 0.;
		}
		if (FAVGOutput.CalorPiston) {
			FAVGOutput.CalorPistonMED = FAVGOutput.CalorPistonSUM;
			FAVGOutput.CalorPistonSUM = 0.;
		}
		if (FAVGOutput.PresionMediaIndicada) {
			FAVGOutput.PresionMediaIndicadaMED = FAVGOutput.PresionMediaNetaMED +
				FAVGOutput.PresionMediaBombeoMED;
		}
		if (FAVGOutput.MasaAtrapada) {
			//FAVGOutput.MasaAtrapadaMED = FMasaAtrapada;
		}
		if (FAVGOutput.TemperaturaCilindroInterna) {
			FAVGOutput.TemperaturaCilindroInternaMED = FAVGOutput.TemperaturaCilindroInternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCilindroInternaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaCilindroMedia) {
			FAVGOutput.TemperaturaCilindroMediaMED = FAVGOutput.TemperaturaCilindroMediaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCilindroMediaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaCilindroExterna) {
			FAVGOutput.TemperaturaCilindroExternaMED = FAVGOutput.TemperaturaCilindroExternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCilindroExternaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaCulataInterna) {
			FAVGOutput.TemperaturaCulataInternaMED = FAVGOutput.TemperaturaCulataInternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCulataInternaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaCulataMedia) {
			FAVGOutput.TemperaturaCulataMediaMED = FAVGOutput.TemperaturaCulataMediaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCulataMediaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaCulataExterna) {
			FAVGOutput.TemperaturaCulataExternaMED = FAVGOutput.TemperaturaCulataExternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaCulataExternaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaPistonInterna) {
			FAVGOutput.TemperaturaPistonInternaMED = FAVGOutput.TemperaturaPistonInternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaPistonInternaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaPistonMedia) {
			FAVGOutput.TemperaturaPistonMediaMED = FAVGOutput.TemperaturaPistonMediaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaPistonMediaSUM = 0.;
		}
		if (FAVGOutput.TemperaturaPistonExterna) {
			FAVGOutput.TemperaturaPistonExternaMED = FAVGOutput.TemperaturaPistonExternaSUM /
				FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaPistonExternaSUM = 0.;
		}
		if (FAVGOutput.AFRMedio) {
			//FAVGOutput.AFRMedioMED = FAFR;
		}
		if (FAVGOutput.NITMedio) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FAVGOutput.NITMED[i] = FValvEsc[i].NITSUM / FAVGOutput.TiempoSUM;
			//	FValvEsc[i].NITSUM = 0.;
			//}
			//FAVGOutput.NITMedioMED = FAVGOutput.NITMedioSUM / FAVGOutput.TiempoSUM;
			//FAVGOutput.NITMedioSUM = 0.;
		}

		if (FAVGOutput.MasaBlowBy) {
			FAVGOutput.MasaBlowByMED = FAVGOutput.MasaBlowBySUM;
			FAVGOutput.MasaBlowBySUM = 0.;
		}
		if (FAVGOutput.MasaAdmision) {
			//FAVGOutput.MasaAdmisionMED = FMasaPorAdmision;
		}
		if (FAVGOutput.MasaEscape) {
			//FAVGOutput.MasaEscapeMED = FMasaPorEscape;
		}
		if (FAVGOutput.MasaCortocircuito) {
			FAVGOutput.MasaCortocircuitoMED = FAVGOutput.MasaCortocircuitoSUM;
			FAVGOutput.MasaCortocircuitoSUM = 0.;
		}
		if (FAVGOutput.TemperaturaMedia) {
			FAVGOutput.TemperaturaMediaMED = FAVGOutput.TemperaturaMediaSUM / FAVGOutput.TiempoSUM;
			FAVGOutput.TemperaturaMediaSUM = 0.;
		}
		if (FAVGOutput.Swirl) {
			//FAVGOutput.SwirlMED = FSwirlSUM / FAVGOutput.TiempoSUM;
			//FSwirlSUM = 0.;
		}
		if (FAVGOutput.RendVolumetrico) {
			//double DensidadReferencia = FAVGOutput.DensidadReferenciaSUM / FAVGOutput.TiempoSUM;
			//FAVGOutput.RendVolumetricoMED = FMasaAtrapada / DensidadReferencia / FCrankMech->getVolumeDisplaced();
			//FAVGOutput.DensidadReferenciaSUM = 0.;
		}

		FAVGOutput.TiempoSUM = 0;
	}

}

void TCylinder::ReadInstantaneousOutputXML(xml_node node_cyl) {

	int nvars, var;

	FINSOutput.Pressure = false;
	FINSOutput.PresionINS = 0.;
	FINSOutput.Temperature = false;
	FINSOutput.TemperaturaINS = 0.;
	FINSOutput.MomentoAngularEsc = false;
	FINSOutput.MomentoAngularTotalEscINS = 0.;
	FINSOutput.MomentoAngularAdm = false;
	FINSOutput.MomentoAngularTotalAdmINS = 0.;
	FINSOutput.GastoEsc = false;
	FINSOutput.GastoTotalEscINS = 0.;
	FINSOutput.GastoAdm = false;
	FINSOutput.GastoTotalAdmINS = 0.;
	FINSOutput.MachEsc = false;
	FINSOutput.MachAdm = false;
	FINSOutput.SeccionEfectivaAdm = false;
	FINSOutput.SeccionEfectivaTotalAdmINS = 0.;
	FINSOutput.SeccionEfectivaEsc = false;
	FINSOutput.SeccionEfectivaTotalEscINS = 0.;
	FINSOutput.Masa = false;
	FINSOutput.MasaINS = 0.;
	FINSOutput.Volumen = false;
	FINSOutput.VolumenINS = 0.;
	FINSOutput.CoeficienteWoschni = false;
	FINSOutput.CoeficienteWoschniINS = 0.;
	FINSOutput.MasaCombustible = false;
	FINSOutput.MasaCombustibleINS = 0.;
	FINSOutput.FQL = false;
	FINSOutput.FQLINS = 0.;
	FINSOutput.TemperaturaCilindroInterna = false;
	FINSOutput.TemperaturaCilindroInternaINS = 0.;
	FINSOutput.TemperaturaCilindroMedia = false;
	FINSOutput.TemperaturaCilindroMediaINS = 0.;
	FINSOutput.TemperaturaCilindroExterna = false;
	FINSOutput.TemperaturaCilindroExternaINS = 0.;
	FINSOutput.TemperaturaPistonInterna = false;
	FINSOutput.TemperaturaPistonInternaINS = 0.;
	FINSOutput.TemperaturaPistonMedia = false;
	FINSOutput.TemperaturaPistonMediaINS = 0.;
	FINSOutput.TemperaturaPistonExterna = false;
	FINSOutput.TemperaturaPistonExternaINS = 0.;
	FINSOutput.TemperaturaCulataInterna = false;
	FINSOutput.TemperaturaCulataInternaINS = 0.;
	FINSOutput.TemperaturaCulataMedia = false;
	FINSOutput.TemperaturaCulataMediaINS = 0.;
	FINSOutput.TemperaturaCulataExterna = false;
	FINSOutput.TemperaturaCulataExternaINS = 0.;
	FINSOutput.NIT = false;
	FINSOutput.GastoCortocircuito = false;
	FINSOutput.GastoCortocircuitoINS = 0.;
	FINSOutput.ParInstantaneo = false;
	FINSOutput.ParInstantaneoINS = 0.;
	FINSOutput.GastoBlowBy = false;
	FINSOutput.GastoBlowByINS = 0.;
	FINSOutput.FraccionMasica = false;
	FINSOutput.FraccionINS = new double[FWorkingFluid->size()];
	for (int i = 0; i < FWorkingFluid->size(); i++) {
		FINSOutput.FraccionINS[i] = 0.;
	}
	FINSOutput.Gamma = false;
	FINSOutput.GammaINS = 0.;

	FINSOutput.HeatHead = false;
	FINSOutput.HeatHeadINS = 0;
	FINSOutput.HeatCyl = false;
	FINSOutput.HeatCylINS = 0;
	FINSOutput.HeatPis = false;
	FINSOutput.HeatPisINS = 0;

	FPrintINS = false;

	xml_node node_insout = GetNodeChild(node_cyl, "Cyl_InsOutput");
	if (node_insout) {
		FPrintINS = true;
		for (xml_attribute parameter = node_insout.attribute("Parameter"); parameter; parameter = parameter.next_attribute()) {

			if (parameter.value() == "Pressure") {
				FINSOutput.Pressure = true;
			}
			else if (parameter.value() == "Temperature") {
				FINSOutput.Temperature = true;
			}
			else if (parameter.value() == "ExhaustAngularMomentum") {
				FINSOutput.MomentoAngularEsc = true;
			}
			else if (parameter.value() == "IntakeAngulaMomentum") {
				FINSOutput.MomentoAngularAdm = true;
			}
			else if (parameter.value() == "ExhaustMassFlow") {
				FINSOutput.GastoEsc = true;
			}
			else if (parameter.value() == "IntakeMassFlow") {
				FINSOutput.GastoAdm = true;
			}
			else if (parameter.value() == "ExhaustMach") {
				FINSOutput.MachEsc = true;
			}
			else if (parameter.value() == "IntakeMach") {
				FINSOutput.MachAdm = true;
			}
			else if (parameter.value() == "IntakeEffectiveSection") {
				FINSOutput.SeccionEfectivaAdm = true;
			}
			else if (parameter.value() == "ExhaustEffectiveSection") {
				FINSOutput.SeccionEfectivaEsc = true;
			}
			else if (parameter.value() == "Mass") {
				FINSOutput.Masa = true;
			}
			else if (parameter.value() == "Volume") {
				FINSOutput.Volumen = true;
			}
			else if (parameter.value() == "FuelMass") {
				FINSOutput.MasaCombustible = true;
			}
			else if (parameter.value() == "HRL") {
				FINSOutput.FQL = true;
			}
			else if (parameter.value() == "WoschniCoef") {
				FINSOutput.CoeficienteWoschni = true;
			}
			else if (parameter.value() == "IntCylinderTemperature") {
				FINSOutput.TemperaturaCilindroInterna = true;
			}
			else if (parameter.value() == "MedCylinderTeperature") {
				FINSOutput.TemperaturaCilindroMedia = true;
			}
			else if (parameter.value() == "ExtCylinderTemperature") {
				FINSOutput.TemperaturaCilindroExterna = true;
			}
			else if (parameter.value() == "IntPistonTemperature") {
				FINSOutput.TemperaturaPistonInterna = true;
			}
			else if (parameter.value() == "MedPistonTemperature") {
				FINSOutput.TemperaturaPistonMedia = true;
			}
			else if (parameter.value() == "ExtPistonTemperature") {
				FINSOutput.TemperaturaPistonExterna = true;
			}
			else if (parameter.value() == "IntCylHeadTemperature") {
				FINSOutput.TemperaturaCulataInterna = true;
			}
			else if (parameter.value() == "MedCylHeadTemperature") {
				FINSOutput.TemperaturaCulataMedia = true;
			}
			else if (parameter.value() == "ExtCylHeadTemperature") {
				FINSOutput.TemperaturaCulataExterna = true;
			}
			else if (parameter.value() == "NIT") {
				FINSOutput.NIT = true;
			}
			else if (parameter.value() == "Torque") {
				FINSOutput.ParInstantaneo = true;
			}
			else if (parameter.value() == "ShortCircuitMassFlow") {
				FINSOutput.GastoCortocircuito = true;
			}
			else if (parameter.value() == "BlowByMassFlow") {
				FINSOutput.GastoBlowBy = true;
			}
			else if (parameter.value() == "MassFraction") {
				FINSOutput.FraccionMasica = true;
			}
			else if (parameter.value() == "SpecificHeatRatio") {
				FINSOutput.Gamma = true;
			}
			else if (parameter.value() == "CylHeadHeat") {
				FINSOutput.HeatHead = true;
			}
			else if (parameter.value() == "CylinderHeat") {
				FINSOutput.HeatCyl = true;
			}
			else if (parameter.value() == "PistonHeat") {
				FINSOutput.HeatPis = true;
			}
			else {
				cout << "Instantaneous parameter " << parameter << " is not correct for cylinder " << FCylinder_ID << endl;
			}
		}
	}

}

void TCylinder::setAngle()
{
	FCrankMech->setPreviousAngle(FCrankMech->getCurrentAngle());

	double angle = FPhaseDiference;
	if (FPhaseDiference == FAngleCycle)
		angle = 0.;

	FCrankMech->setCurrentAngle(angle);

	FCombAngle = FCrankMech->getCurrentAngle();
	if (FCombAngle > FAngleCycle / 2.) {
		FCombAngle -= FAngleCycle;
		FInj = false;
		FComb = false;
	}

}


void TCylinder::HeaderInstantaneousOutput(stringstream& insoutput) {

	if (FPrintINS) {

		char *label1, *label2, *label3;
		std::string Label;
		string kkkk;
		// cadena Label7;

		if (FINSOutput.Pressure) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4006) + PutLabel(908);
			insoutput << Label.c_str();
		}
		if (FINSOutput.Temperature) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3041);
			insoutput << Label.c_str();
		}
		if (FINSOutput.MomentoAngularEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4013) + PutLabel(
			//		901) + "/" + PutLabel(3018) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4013) + PutLabel(
			//		901) + "/" + PutLabel(3018) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.MomentoAngularAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4013) + PutLabel(
			//		901) + "/" + PutLabel(3017) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4013) + PutLabel(
			//		901) + "/" + PutLabel(3017) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.GastoEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
			//		904) + "/" + PutLabel(3018) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
			//		904) + "/" + PutLabel(3018) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.GastoAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
			//		904) + "/" + PutLabel(3017) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
			//		904) + "/" + PutLabel(3017) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.MachEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4014) + PutLabel(
			//		901) + "/" + PutLabel(3018) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.MachAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4014) + PutLabel(
			//		901) + "/" + PutLabel(3017) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.SeccionEfectivaEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3029) + PutLabel(
			//		4012) + PutLabel(915) + "/" + PutLabel(3018) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3029) + PutLabel(
			//		4012) + PutLabel(915) + "/" + PutLabel(3018) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.SeccionEfectivaAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3029) + PutLabel(
			//		4012) + PutLabel(915) + "/" + PutLabel(3017) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3029) + PutLabel(
			//		4012) + PutLabel(915) + "/" + PutLabel(3017) + "/" + PutLabel(3026);
			//	insoutput << Label.c_str();
			//}
		}
		if (FINSOutput.Masa) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3041);
			insoutput << Label.c_str();
		}
		if (FINSOutput.Volumen) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4003) + PutLabel(912);
			insoutput << Label.c_str();
		}
		if (FINSOutput.CoeficienteWoschni) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4015) + PutLabel(911);
			insoutput << Label.c_str();
		}
		if (FINSOutput.MasaCombustible) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4004) + PutLabel(
				913) + "/" + PutLabel(3027);
			insoutput << Label.c_str();
		}
		if (FINSOutput.FQL) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4018) + PutLabel(901);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCilindroInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3011);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCilindroMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3012);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCilindroExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3006) + "/" + PutLabel(3013);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaPistonInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3011);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaPistonMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3012);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaPistonExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3008) + "/" + PutLabel(3013);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCulataInterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3011);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCulataMedia) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3012);
			insoutput << Label.c_str();
		}
		if (FINSOutput.TemperaturaCulataExterna) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4005) + PutLabel(
				910) + "/" + PutLabel(3007) + "/" + PutLabel(3013);
			insoutput << Label.c_str();
		}
		if (FINSOutput.NIT) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3014) + PutLabel(
			//		903) + "/" + PutLabel(3018) + "/" + std::to_string(i);
			//	insoutput << Label.c_str();
			//}
			//Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(3014) + PutLabel(
			//	903) + "/" + PutLabel(3018) + "/" + PutLabel(3026);
			//insoutput << Label.c_str();
		}
		if (FINSOutput.ParInstantaneo) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4016) + PutLabel(917);
			insoutput << Label.c_str();
		}
		if (FINSOutput.GastoCortocircuito) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
				904) + "/" + PutLabel(3019);
			insoutput << Label.c_str();
		}
		if (FINSOutput.GastoBlowBy) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4008) + PutLabel(
				904) + "/" + PutLabel(3019);
			insoutput << Label.c_str();
		}
		if (FINSOutput.FraccionMasica) {
			for (int i = 0; i < FWorkingFluid->size(); i++) {
				Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4023) + PutLabel(
					901) + "/" + FWorkingFluid->GetComponent(i)->getName();
				insoutput << Label.c_str();
			}
		}
		if (FINSOutput.Gamma) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4017) + PutLabel(901);
			insoutput << Label.c_str();
		}
		if (FINSOutput.HeatHead) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				903) + "/" + PutLabel(3007);
			insoutput << Label.c_str();
		}
		if (FINSOutput.HeatCyl) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				903) + "/" + PutLabel(3006);
			insoutput << Label.c_str();
		}
		if (FINSOutput.HeatPis) {
			Label = "\t" + PutLabel(5001) + "/" + std::to_string(FCylinder_ID) + "/" + PutLabel(4010) + PutLabel(
				903) + "/" + PutLabel(3008);
			insoutput << Label.c_str();
		}
	}

}


void TCylinder::PrintInstantaneuosOutput(stringstream& insoutput) {

	if (FPrintINS) {

		if (FINSOutput.Pressure)
			insoutput << "\t" << FINSOutput.PresionINS;
		if (FINSOutput.Temperature)
			insoutput << "\t" << FINSOutput.TemperaturaINS;
		if (FINSOutput.MomentoAngularEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	insoutput << "\t" << FINSOutput.MomentoAngularEscINS[i];
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	insoutput << "\t" << FINSOutput.MomentoAngularTotalEscINS;
			//}
		}
		if (FINSOutput.MomentoAngularAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	insoutput << "\t" << FINSOutput.MomentoAngularAdmINS[i];
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	insoutput << "\t" << FINSOutput.MomentoAngularTotalAdmINS;
			//}
		}
		if (FINSOutput.GastoEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	insoutput << "\t" << FINSOutput.GastoEscINS[i];
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	insoutput << "\t" << FINSOutput.GastoTotalEscINS;
			//}
		}
		if (FINSOutput.GastoAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	insoutput << "\t" << FINSOutput.GastoAdmINS[i];
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	insoutput << "\t" << FINSOutput.GastoTotalAdmINS;
			//}
		}
		if (FINSOutput.MachEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	insoutput << "\t" << FINSOutput.MachEscINS[i];
			//}
		}
		if (FINSOutput.MachAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	insoutput << "\t" << FINSOutput.MachAdmINS[i];
			//}
		}
		if (FINSOutput.SeccionEfectivaEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	insoutput << "\t" << FINSOutput.SeccionEfectivaEscINS[i];
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	insoutput << "\t" << FINSOutput.SeccionEfectivaTotalEscINS;
			//}
		}
		if (FINSOutput.SeccionEfectivaAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	insoutput << "\t" << FINSOutput.SeccionEfectivaAdmINS[i];
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	insoutput << "\t" << FINSOutput.SeccionEfectivaTotalAdmINS;
			//}
		}
		if (FINSOutput.Masa)
			insoutput << "\t" << FINSOutput.MasaINS;
		if (FINSOutput.Volumen)
			insoutput << "\t" << FINSOutput.VolumenINS;
		if (FINSOutput.CoeficienteWoschni)
			insoutput << "\t" << FINSOutput.CoeficienteWoschniINS;
		if (FINSOutput.MasaCombustible)
			insoutput << "\t" << FINSOutput.MasaCombustibleINS;
		if (FINSOutput.FQL)
			insoutput << "\t" << FINSOutput.FQLINS;
		if (FINSOutput.TemperaturaCilindroInterna)
			insoutput << "\t" << FINSOutput.TemperaturaCilindroInternaINS;
		if (FINSOutput.TemperaturaCilindroMedia)
			insoutput << "\t" << FINSOutput.TemperaturaCilindroMediaINS;
		if (FINSOutput.TemperaturaCilindroExterna)
			insoutput << "\t" << FINSOutput.TemperaturaCilindroExternaINS;
		if (FINSOutput.TemperaturaPistonInterna)
			insoutput << "\t" << FINSOutput.TemperaturaPistonInternaINS;
		if (FINSOutput.TemperaturaPistonMedia)
			insoutput << "\t" << FINSOutput.TemperaturaPistonMediaINS;
		if (FINSOutput.TemperaturaPistonExterna)
			insoutput << "\t" << FINSOutput.TemperaturaPistonExternaINS;
		if (FINSOutput.TemperaturaCulataInterna)
			insoutput << "\t" << FINSOutput.TemperaturaCulataInternaINS;
		if (FINSOutput.TemperaturaCulataMedia)
			insoutput << "\t" << FINSOutput.TemperaturaCulataMediaINS;
		if (FINSOutput.TemperaturaCulataExterna)
			insoutput << "\t" << FINSOutput.TemperaturaCulataExternaINS;
		if (FINSOutput.NIT) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	insoutput << "\t" << FINSOutput.NITINS[i];
			//}
			//insoutput << "\t" << FINSOutput.NITTotalINS;
		}
		if (FINSOutput.ParInstantaneo)
			insoutput << "\t" << FINSOutput.ParInstantaneoINS;
		if (FINSOutput.GastoCortocircuito)
			insoutput << "\t" << FINSOutput.GastoCortocircuitoINS;
		if (FINSOutput.GastoBlowBy)
			insoutput << "\t" << FINSOutput.GastoBlowByINS;
		if (FINSOutput.FraccionMasica) {
			for (int i = 0; i < FWorkingFluid->size(); i++) {
				insoutput << "\t" << FINSOutput.FraccionINS[i];
			}
		}
		if (FINSOutput.Gamma)
			insoutput << "\t" << FINSOutput.GammaINS;

		if (FINSOutput.HeatHead)
			insoutput << "\t" << FINSOutput.HeatHeadINS;
		if (FINSOutput.HeatCyl)
			insoutput << "\t" << FINSOutput.HeatCylINS;
		if (FINSOutput.HeatPis)
			insoutput << "\t" << FINSOutput.HeatPisINS;

	}

}


void TCylinder::SetInstantaneousOutput() {

	if (FPrintINS) {

		double gastoesc = 0., gastoadm = 0.;
		double secefectotaladm = 0., secefectotalesc = 0.;

		if (FINSOutput.Pressure)
			FINSOutput.PresionINS = FZeroD->getPressure();
		if (FINSOutput.Temperature)
			FINSOutput.TemperaturaINS = FZeroD->getTemperature();
		if (FINSOutput.MomentoAngularEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FINSOutput.MomentoAngularEscINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getMomento();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	FINSOutput.MomentoAngularTotalEscINS = FMomentoAngularEsc;
			//}
		}
		if (FINSOutput.MomentoAngularAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	FINSOutput.MomentoAngularAdmINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getMomento();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	FINSOutput.MomentoAngularTotalAdmINS = FMomentoAngularAdm;
			//}
		}
		if (FINSOutput.GastoEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FINSOutput.GastoEscINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getMassflow();
			//	gastoesc += dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getMassflow();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	FINSOutput.GastoTotalEscINS = gastoesc;
			//}
		}
		if (FINSOutput.GastoAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	FINSOutput.GastoAdmINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getMassflow();
			//	gastoadm += dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getMassflow();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	FINSOutput.GastoTotalAdmINS = gastoadm;
			//}
		}
		if (FINSOutput.MachEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FINSOutput.MachEscINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getMach();
			//}
		}
		if (FINSOutput.MachAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	FINSOutput.MachAdmINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getMach();
			//}
		}
		if (FINSOutput.SeccionEfectivaEsc) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FINSOutput.SeccionEfectivaEscINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getSeccionEficaz();
			//	secefectotalesc += dynamic_cast<TCCCilindro*>(FCCValvulaEsc[i])->getSeccionEficaz();
			//}
			//if (FNumeroUnionesEsc > 1) {
			//	FINSOutput.SeccionEfectivaTotalEscINS = secefectotalesc;
			//}
		}
		if (FINSOutput.SeccionEfectivaAdm) {
			//for (int i = 0; i < FNumeroUnionesAdm; i++) {
			//	FINSOutput.SeccionEfectivaAdmINS[i] = dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getSeccionEficaz();
			//	secefectotaladm += dynamic_cast<TCCCilindro*>(FCCValvulaAdm[i])->getSeccionEficaz();
			//}
			//if (FNumeroUnionesAdm > 1) {
			//	FINSOutput.SeccionEfectivaTotalAdmINS = secefectotaladm;
			//}
		}
		if (FINSOutput.Masa)
			FINSOutput.MasaINS = FZeroD->getMass();
		if (FINSOutput.Volumen)
			FINSOutput.VolumenINS = FCrankMech->getVolume();
		if (FINSOutput.CoeficienteWoschni)
			FINSOutput.CoeficienteWoschniINS = 0;
		if (FINSOutput.MasaCombustible)
			FINSOutput.MasaCombustibleINS = FInjector->getFuelInjected();
		if (FINSOutput.FQL)
			FINSOutput.FQLINS = FCombustion->getHRL();
		if (FINSOutput.TemperaturaCilindroInterna)
			FINSOutput.TemperaturaCilindroInternaINS = 0;
		if (FINSOutput.TemperaturaCilindroMedia)
			FINSOutput.TemperaturaCilindroMediaINS = 0;
		if (FINSOutput.TemperaturaCilindroExterna)
			FINSOutput.TemperaturaCilindroExternaINS = 0;
		if (FINSOutput.TemperaturaPistonInterna)
			FINSOutput.TemperaturaPistonInternaINS = 0;
		if (FINSOutput.TemperaturaPistonMedia)
			FINSOutput.TemperaturaPistonMediaINS = 0;
		if (FINSOutput.TemperaturaPistonExterna)
			FINSOutput.TemperaturaPistonExternaINS = 0;
		if (FINSOutput.TemperaturaCulataInterna)
			FINSOutput.TemperaturaCulataInternaINS = 0;
		if (FINSOutput.TemperaturaCulataMedia)
			FINSOutput.TemperaturaCulataMediaINS = 0;
		if (FINSOutput.TemperaturaCulataExterna)
			FINSOutput.TemperaturaCulataExternaINS = 0;
		if (FINSOutput.NIT) {
			//for (int i = 0; i < FNumeroUnionesEsc; i++) {
			//	FINSOutput.NITINS[i] = FValvEsc[i].NIT;
			//}
			//FINSOutput.NITTotalINS = FNIT;
		}
		if (FINSOutput.ParInstantaneo)
			FINSOutput.ParInstantaneoINS = 0;
		if (FINSOutput.GastoCortocircuito)
			FINSOutput.GastoCortocircuitoINS = 0;
		if (FINSOutput.GastoBlowBy)
			FINSOutput.GastoBlowByINS = 0;
		if (FINSOutput.FraccionMasica) {
			for (int i = 0; i < FWorkingFluid->size(); i++) {
				FINSOutput.FraccionINS[i] = FWorkingFluid->GetY(i);
			}
		}
		if (FINSOutput.Gamma)
			FINSOutput.GammaINS = FWorkingFluid->FunGamma(FZeroD->getTemperature());
		if (FINSOutput.HeatHead)
			FINSOutput.HeatHeadINS = 0;
		if (FINSOutput.HeatCyl)
			FINSOutput.HeatCylINS = 0;
		if (FINSOutput.HeatPis)
			FINSOutput.HeatPisINS = 0;

	}

}

TCylinder_ptr create_cylinder(xml_node node_engine, xml_node node_cylinder, TComponentArray_ptr WorkFluid) {

	auto Cylinder = make_shared<TCylinder>();

	int ID = IDtoInt(node_cylinder.attribute("Cyl_ID").as_string()) - 1;

	string EngineType = node_engine.attribute("Type").as_string();
	Cylinder->FAngleCycle = 4 * __cons::Pi;
	if (EngineType == "2Strokes") {
		Cylinder->FAngleCycle = __cons::Pi_x_2;
	}

	//!< GEOMETRY

	double speed = GetAttributeAsDouble(node_engine, "Speed");

	xml_node node_geo = GetNodeChild(node_engine, "Eng_Geometry");
	TCrankMechanism_ptr CrankMech = make_unique<TCrankMechanism>(speed);
	CrankMech->ReadInputDataXML(node_geo, EngineType);

	Cylinder->setCrankMechanism(move(CrankMech));

	//! INITAL ANGLE
	int ncyl = CountNodes(node_engine, "Cylinder");

	Cylinder->FPhaseDiference = Cylinder->FAngleCycle / (double)ncyl;
	Cylinder->FOrder = GetAttributeAsInt(node_cylinder, "Order") - 1;
	Cylinder->FPhaseDiference *= Cylinder->FOrder;
	Cylinder->FPhaseDiference = Cylinder->FAngleCycle - Cylinder->FPhaseDiference;
	Cylinder->setAngle();

	//!< INITAL DATA
	xml_node node_gas = GetNodeChild(node_engine, "GasProperties");
	double PressureIC;
	double TemperatureIC;
	double V;
	RowVector Y;
	Y.setZero(WorkFluid->size());
	ReadGasProperties(node_gas, PressureIC, TemperatureIC, V, Y, WorkFluid);
	Cylinder->setFluid(WorkFluid, Y);

	T0DModel_ptr zerod = make_unique<T0DModel>();
	zerod->setFluid(Cylinder->FWorkingFluid);

	//!< HEAT TRANSFER
	BasicHeatTransfer *ht;
	xml_node node_ht = GetNodeChild(node_engine, "HeatTransfer");
	string type = node_ht.attribute("Type").as_string();
	if (type == "CI-Swirl") {
		ht = new HeatTransferCombMECSwirl();
	}
	ht->ReadData(node_engine);

	zerod->setHeatTransfer(ht);

	Cylinder->set0DModel(move(zerod));

	Cylinder->setInitialConditions(PressureIC, TemperatureIC, EngineType);

	//!< COMBUSTION

	TCombustion_ptr Combustion;

	xml_node node_comb = GetNodeChild(node_engine, "Combustion");
	xml_node node_cbtype;
	if (node_comb.child("CombustionModel")) {

	}
	else if (node_comb.child("Wiebe")) {
		Combustion = make_unique<TWiebe>();
		node_cbtype = GetNodeChild(node_comb, "Wiebe");
	}
	else if (node_comb.child("MultiWiebe")) {
		Combustion = make_unique<TMultiWiebe>();
		node_cbtype = GetNodeChild(node_comb, "MultiWiebe");
	}
	else if (node_comb.child("HRL")) {
		Combustion = make_unique<THRL>();
		node_cbtype = GetNodeChild(node_comb, "HRL");
	}
	Combustion->ReadCombustionData(node_cbtype);

	Cylinder->setCombustion(move(Combustion));

	//!< INJECTOR

	xml_node node_injsys = GetNodeChild(node_engine, "InjectionSystem");
	xml_node node_injector = GetNodeChild(node_injsys, "Injector");

	Cylinder->FInjector = createInjector(node_injector, WorkFluid);

	//!< OUTPUT

	Cylinder->ReadAverageOutputXML(node_cylinder);
	Cylinder->ReadInstantaneousOutputXML(node_cylinder);

	Cylinder->ReadOutput(node_cylinder);

	return Cylinder;

}
