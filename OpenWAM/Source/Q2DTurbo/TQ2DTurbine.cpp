/*--------------------------------------------------------------------------------*\
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
 * @file TQ2DTurbine.cpp
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
 * This file defines a quasi-two-dimensional turbine class.
 */

#include "TQ2DTurbine.h"

#include "TCondicionContorno.h"
#include "TCCDeposito.h"
#include "TRotorTurbina.h"
#include "TEstatorTurbina.h"
#include "TTubo.h"

TQ2DTurbine::TQ2DTurbine() {
	FPresRef = 1.;
	FNumeroTurbina = 1;
}

double TQ2DTurbine::GetEfficiency() {
	return FTurbineEfficiency.mean();
}

double TQ2DTurbine::GetEfficiency() const {
	return FTurbineEfficiency.mean();
}

vector<Pipe_ptr> TQ2DTurbine::getInletPipes() const {
	return std::vector<Pipe_ptr>(1,
		dynamic_cast<TQ2DAcousticTurbine*>(FAcTurb)->FInlet);
}

Pipe_ptr TQ2DTurbine::getOutletPipe() const {
	return dynamic_cast<TQ2DAcousticTurbine*>(FAcTurb)->FOutlet;
}

vector<Volute_ptr> TQ2DTurbine::getVolutes() const {
	return std::vector<Volute_ptr>(1,
		dynamic_cast<TQ2DAcousticTurbine*>(FAcTurb)->FVolute);
}

void TQ2DTurbine::ReadXML(const xml_node & node,
	const TComponentArray_ptr& fluid, const std::map<string, TSolid*> &MDB) {
	try {
		std::string turbine_type = node.attribute("Type").as_string();
		auto p_0 = 0.;
		auto T_0 = 0.;
		FName = node.attribute("Name").as_string();
		auto node_pipes = GetNodeChild(node, "BlockOfPipes");
		auto node_pipe = node_pipes.child("Bop_NewPipe");
		auto inlet = create_pipe(node_pipe, fluid, MDB);
		node_pipe = node_pipe.next_sibling();
		auto volute = create_volute(node_pipe, fluid, MDB);
		node_pipe = node_pipe.next_sibling();
		auto outlet = create_pipe(node_pipe, fluid, MDB);
		FPressure = (inlet->getPressure((unsigned int)0)
			+ outlet->getPressure((unsigned int)0)) / 2.;
		FTemperature = (inlet->getTemperature((unsigned int)0)
			+ outlet->getTemperature((unsigned int)0)) / 2.;
		attach_pipes(inlet, nmRight, volute, nmLeft);

		auto FluidBC = make_shared<TFluid>(fluid);
		FluidBC->SetComposition(outlet->getFluid(0)->GetComposition());

		auto node_plenum = GetNodeChild(node, "Plm_Turbine");
		std::string tipoturb = node_plenum.attribute("Turbine_type").as_string();
		if(tipoturb == "FixedDC") {
			FTipoTurbina = nmFixedTurbine;
		} else if(tipoturb == "VariableDC") {
			FTipoTurbina = nmVariableGeometry;
		} else if(tipoturb == "Map") {
			FTipoTurbina = nmTurbineMap;
		}
		FDiametroRodete = GetXMLLength(node_plenum, "WheelDiameter");
		if(FTipoTurbina == nmTurbineMap) {
			FDiametroRodeteOut = GetXMLLength(node_plenum, "WheelDiameterOut");
			FDiametroTuerca = GetXMLLength(node_plenum, "NutDiameter");
			FDiametroTurbinaIn = GetXMLLength(node_plenum, "InletDiameter");
			FCriticalAngle = GetXMLAngle(node_plenum, "CriticalAngle");
			FBladeHeight = GetXMLLength(node_plenum, "BladeHeight");
			FBeta = GetXMLAngle(node_plenum, "Beta");
			FBeta = __units::DegToRad(FBeta);

			FMapa = new TTurbineMap();
			FMapa->LoadTurbineMapXML(node_plenum);

			if(FMapa->IsFixed()) {
				FVaneless = GetAttributeAsBool(node_plenum, "Vaneless");
				if(!FVaneless) {
					FZ0 = GetAttributeAsInt(node_plenum, "Z0");
					VGTtoAlpha2 = GetAttributeAsDouble(node_plenum, "VGTtoAlpha2");
				}
			} else {
				VGTtoAlpha1 = GetAttributeAsDouble(node_plenum, "VGTtoAlpha1");
				VGTtoAlpha2 = GetAttributeAsDouble(node_plenum, "VGTtoAlpha2");
				FLTE = GetAttributeAsDouble(node_plenum, "LTE");
				FR2geom = GetAttributeAsDouble(node_plenum, "R2geom");
				FZ0 = GetAttributeAsInt(node_plenum, "Z0");
			}

			for(auto node_plenum_ctrl = GetNodeChild(node_plenum, "Actuator"); node_plenum_ctrl;
				node_plenum_ctrl = node_plenum_ctrl.next_sibling("Actuator")) {

				std::string ctrl = node_plenum_ctrl.attribute("Parameter").value();
				if(ctrl == "Rack") {
					FRackIsControlled = true;
					FNumControlObject = GetAttributeAsInt(node_plenum_ctrl, "CtrlID");
				}
			}

			if(!FRackIsControlled)
				FRack = GetAttributeAsDouble(node_plenum, "Rack");
			FCalRendTurbina = nmRendMapa;

			attach_to_pressure_BC(outlet, nmLeft, FPressure, FTemperature, FluidBC);
			if (turbine_type == "Q2DTurbine") {
				FQ2DAcTurb = make_shared<TQ2DAcousticTurbine>(inlet, volute, outlet);
				set_lateral_nozzles(volute);
			} else if (turbine_type == "1DTurbine") {
				attach_to_pressure_BC(volute, nmRight, FPressure, FTemperature, FluidBC);
				FQ2DAcTurb = make_shared<T1DAcousticTurbine>(inlet, volute, outlet);
			}

			FQ2DAcTurb->setVoluteOutletConditions(FPressure, FTemperature);
			FQ2DAcTurb->setRotorInletConditions(FPressure, FTemperature);
			
			FVolumen = FBladeHeight * (
				__geom::Circle_area(volute->getX().tail(1)(0) / __cons::Pi)
				- __geom::Circle_area(FDiametroRodete)) * 0.4;
			FDensidad = FPressure / FTemperature / __R::Air;
			FMasa = FVolumen * FDensidad;
			FPressure = __units::PaToBar(FPressure);
			FTemperature = __units::KTodegC(FTemperature);
			FAcTurb = FQ2DAcTurb.get();
			FMembers.push_back(FQ2DAcTurb);
		} else {
			std::string rdturb = node_plenum.attribute("Eff_type").value();
			if(rdturb == "Watson") {
				FCalRendTurbina = nmWatson;
			} else if(rdturb == "Polinomy") {
				FCalRendTurbina = nmPolinomio;
				FRcoptima = GetAttributeAsDouble(node_plenum, "BSRopt");
				FRcmaxima = GetAttributeAsDouble(node_plenum, "BSRmax");
				FRendmaximo = GetAttributeAsDouble(node_plenum, "Effmax");
			} else if(rdturb == "ExtCalculated") {
				FCalRendTurbina = nmCalcExtRD;
			} else {
				std::cout << "ERROR: Unknown method to calculate turbine efficiency " << std::endl;
			}
		}
		FAjustRendTurb = GetAttributeAsDouble(node_plenum, "FitEfficiency");

		if(node.child("Trb_InsOutput"))
			ReadInsResults(node);
		if(node.child("Trb_AvgOutput"))
			ReadAverageResultsTurbXML(node);

	FBSR.setZero(FQ2DAcTurb->VoluteOutletp().size());
	FTurbineEfficiency.setZero(FQ2DAcTurb->VoluteOutletp().size());
	FRMezcla = __R::Air;
	FCvMezcla = __Gamma::Cv;
	FCpMezcla = __Gamma::Cp;
	FGamma = __Gamma::G;
	FFraccionMasicaEspecie = new double[2];
	FFraccionMasicaEspecie[0] = 1.;
	FFraccionMasicaEspecie[1] = 0.;
	
	} catch(exception & N) {
		std::cout << "ERROR: TQ2DTurbine::ReadXML in turbine " << FNumeroTurbina << std::endl;
		std::cout << "Error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TQ2DTurbine::setInletBC(BoundaryCondition_ptr bc) {
	FQ2DAcTurb->FInlet->setLeftBC(move(bc));
}

void TQ2DTurbine::setOutletBC(BoundaryCondition_ptr bc) {
	FQ2DAcTurb->FOutlet->setRightBC(move(bc));
}

void TQ2DTurbine::Solve() {
	UpdateWorkingPoint(FTimeStep + FCurrentTime);
	FQ2DAcTurb->setVoluteOutletConditions(__units::BarToPa(FPressure),
		__units::degCToK(FTemperature));
	FQ2DAcTurb->setRotorInletConditions(__units::BarToPa(FPressure),
		__units::degCToK(FTemperature));
	TIntegrableGroup::Solve();
}

void TQ2DTurbine::Solve(double t) {
	double new_t = 0.;
	while ((t - FCurrentTime) > 0.) {
		FTimeStep = getMaxTimeStep();
		Solve();
	}
}

void TQ2DTurbine::UpdateWorkingPoint(double t) {

	double delta_t = 0.;
	double incrRelCin = 0.;
	double dd = 0.;
	double b = 0.;
	double c = 0.;
	double d = 0.;
	
	try {

		delta_t = t - FTimeTurbina;
		FTimeTurbina = t;

		if(delta_t > 0) {
			FTrabajoIsenInstTotal = 0.;

			FInletMassFlow = FQ2DAcTurb->VoluteOutletMassFlow();
			FInletTemperature = FQ2DAcTurb->VoluteOutletT();
			FInletPressure = FQ2DAcTurb->VoluteOutletp();
			FOutletPressure = FQ2DAcTurb->OutletNozzlep();
			FOutletMassFlow = FQ2DAcTurb->RotorInletMassFlow();
// TODO: c_p, R...
			FRMezcla = 287;
			Fcp_med.setConstant(FInletMassFlow.size(), __Gamma::Cp);
			FInletTotalTemperature = FQ2DAcTurb->VoluteOutletT0();
			FInletTotalPressure = FQ2DAcTurb->VoluteOutletp0();
			FInletTotalPressure.setConstant(__units::BarToPa(FQ2DAcTurb->P30()));
			FInletTotalTemperature.setConstant(FQ2DAcTurb->T30());
			FOutletPressure = __units::BarToPa(FQ2DAcTurb->P4());
// TODO: HTM
			FReducedMassFlow = FInletMassFlow * FInletTotalTemperature.sqrt()
				/ FInletTotalPressure;
			FReducedSpeed = FRegimen / FInletTotalTemperature.sqrt();
			FExpansionRatio = FInletTotalPressure / FOutletPressure;
			for (auto i = 0; i < FExpansionRatio.size(); i++) {
				if (FExpansionRatio(i) < 1.) {
					FExpansionRatio(i) = 1.;
				}
			}
			FOutletIsentropicTemperature = FInletTotalTemperature * pow(1. / FExpansionRatio,
				FRMezcla / Fcp_med);
			if(FTipoTurbina == nmTurbineMap) {
				if(FRackIsControlled) {
					FRack = FRackController->Output(FTime);
				}
				auto MassAcum = FInletMassFlow.sum();
				double RotorEF = 0;
				if(FThereIsHTM)
					FMaxEfficiency = FMapa->MaxEfficiency(FRack);
// 				for(auto i = 0; i < FReducedSpeed.size(); i++) {
// 					if(FMapa->IsExtrapolated()) {
// 						FMapa->InterpExtendedMap(FReducedSpeed(i) / 60., FExpansionRatio(i),
// 							FBSR(i), FRack);
// 					} else {
// 						FMapa->CurrentEffectiveSection(FReducedSpeed(i) / 60.,  FExpansionRatio(i), FRack,
// 							__units::degCToK(FTemperature) / FInletTotalTemperature(i));
// 					}
				if(FMapa->IsExtrapolated()) {
					FMapa->InterpExtendedMap(FReducedSpeed(0) / 60., FExpansionRatio(0),
						FBSR(0), FRack);
				} else {
					FMapa->CurrentEffectiveSection(FReducedSpeed(0) / 60.,  FExpansionRatio(0), FRack,
						__units::degCToK(FTemperature) / FInletTotalTemperature(0));
				}
				for(auto i = 0; i < FReducedSpeed.size(); i++) {

					FQ2DAcTurb->setVoluteNozzleEffectiveArea(FMapa->StatorEF() / FInletMassFlow.size(), i);

					if(FInletMassFlow(i) > 0) {
						RotorEF += FMapa->RotorEF() * FInletMassFlow(i);
					}
					FTurbineEfficiency(i) = FMapa->EffTurb();
				}
				if(MassAcum > 0)
					RotorEF /= MassAcum;
				FQ2DAcTurb->setRotorOutletEffectiveArea(RotorEF);
			}

			FBSR = __units::RPMToRPS(FRegimen) * __cons::Pi * FDiametroRodete / sqrt(2 *
				Fcp_med * (FInletTotalTemperature - FOutletIsentropicTemperature + 1E-6));
			for (auto i = 0; i < FBSR.size(); i++)
			{
				if (FBSR(i) < 0.) {
					FBSR(i) = 0.;
				}
			}
			if (FCalRendTurbina == nmWatson) {
				FTurbineEfficiency = 0.004022 + 1.55766 * FBSR - 0.511626 * pow2(FBSR) - 0.121795
					* pow3(FBSR) - 0.445804 * pow4(FBSR);
				for (auto i = 0; i < FBSR.size(); i++) {
					if(FBSR(i) <= 0 || FBSR(i) >= 1.19) {
						FTurbineEfficiency(i) = 0;
					}
				}
			} else if (FCalRendTurbina == nmPolinomio) {
				dd = 2. * pow3(FRcoptima) * pow2(FRcmaxima) - pow2(FRcoptima) * pow3(FRcmaxima) - pow4(FRcoptima) * FRcmaxima;
				b = FRendmaximo * (3. * pow2(FRcmaxima) * pow2(FRcoptima) - 2 * pow3(FRcmaxima) * FRcoptima) / dd;
				c = FRendmaximo * (pow3(FRcmaxima) - 3. * FRcmaxima * pow2(FRcoptima)) / dd;
				d = FRendmaximo * (2. * FRcmaxima * FRcoptima - pow2(FRcmaxima)) / dd;
				FTurbineEfficiency = b * FBSR + c * pow2(FBSR) + d * pow3(FBSR);
				for (auto i = 0; i < FBSR.size(); i++) {
					if (FBSR(i) >= FRcmaxima || FBSR(i) <= 0) {
						FTurbineEfficiency(i) = 0.;
					}
				}
			} else if (FCalRendTurbina == nmCalcExtRD) {
				for (auto i = 0; i < FBSR.size(); i++) {
					if (FDatosTGV[FNumeroTurbinaTGV].Rendimiento[i] < 0)
					{
						FTurbineEfficiency(i) = 0.;
					} else {
						FTurbineEfficiency(i) = FDatosTGV[FNumeroTurbinaTGV].Rendimiento[i];
					}
				}
			} else if(FCalRendTurbina == nmRendMapa) {
				// Already computed.
			} else {
				std::cout << "ERROR: Calculo del rendimiento de la turbina desconocido" << std::endl;
					throw Exception("");
			}
			FTurbineEfficiency *= FAjustRendTurb;

			FTrabajoFluido = 0;
			incrRelCin = 0;

			FIsentropicWork = (FInletMassFlow * Fcp_med
				* (FInletTotalTemperature - FOutletIsentropicTemperature)) * delta_t;
			FTrabajoIsenInstTotal = FIsentropicWork.sum();
// TODO: Work, power...
// 				// Para el calculo de la potencia del paso.
// 				// FRendInstantaneo+=FRendTurbina[i]*FGastoEntrada[i];
// 				// FGastoEntradaTotal+=FGastoEntrada[i];
// 				incrRelCin = FRendTurbina[i] * FRelacionCinematica[i] * FTrabajoIsen;
// 
// 				FTrabajoTotal += FTrabajoIsen;
// 				FRelacionCinAcum[i] += incrRelCin;
// 				FRelacionCinGlobalAcum += incrRelCin;
// 				FPonderacionRelacionCinematica[i] += FRendTurbina[i] * FTrabajoIsen;
			FTrabajoReal = (FIsentropicWork * FTurbineEfficiency).sum();
			FTrabajoRealPaso += FTrabajoReal;
			FTrabajoFluido = FTrabajoReal
				- (FInletMassFlow * Fcp_med * (FInletTotalTemperature - FInletTemperature)).sum() * delta_t;
			if (delta_t != 0.) {
				FPotencia = FTrabajoFluido / delta_t;
			}
			else {
				FPotencia = 0.;
			}
			FDeltaPaso += delta_t;
		} else {
			FTrabajoIsenInstTotal = 0.;
			FTrabajoFluido = 0.;
		}
	UpdatePlenum(t);
	} catch(exception &N) {
		std::cout << "ERROR: TQ2DTurbine::CalculaCondicionTurbina en la turbina: " << FNumeroTurbina << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TQ2DTurbine::UpdatePlenum(double t) {
	double H = 0.; // Entalpia de entrada
	double Energia = 0.;
	double MasaEntrante, FraccionMasicaAcum = 0.;
	auto DeltaT = t - FTime;
	double g = 0., v = 0., a = 0., m = 0.;
	int SignoFlujo = 1;
	double MassIn = 0.;
	double MassOut = 0.;
	double EnthalpyIn = 0.;
	double EnthalpyOut = 0.;
	if(FThereIsHTM)
		FHeatPower = FHTM->TurbineHeatFlow();
	double Heat = FHeatPower * DeltaT;

	try {
		FMasa0 = FMasa;
		MasaEntrante = 0.;
		H = 0.;
		MassIn = FInletMassFlow.sum() * DeltaT;
		MassOut = FOutletMassFlow * DeltaT;
		FMasa += MassIn - MassOut;
		EnthalpyIn = (FInletMassFlow * FCpMezcla
			* FQ2DAcTurb->VoluteOutletT0()).sum() * DeltaT;
		EnthalpyOut = MassOut * FCpMezcla * __units::degCToK(FTemperature);
		FTemperature = (__units::degCToK(FTemperature) * FMasa0 * FCvMezcla
			+ (EnthalpyIn - EnthalpyOut - FTrabajoReal - Heat))
			/ (FMasa * FCvMezcla);
		FPressure = __units::PaToBar(FMasa / FVolumen * FRMezcla * FTemperature);
		FPresionIsen = pow(FPressure / FPresRef, __Gamma::G5(FGamma));
		FTemperature = __units::KTodegC(FTemperature);

		if(!(DoubEqZero(DeltaT))) {
			if(FCalculoEspecies == nmCalculoCompleto) {

				FRMezcla = CalculoCompletoRMezcla(FFraccionMasicaEspecie[0], FFraccionMasicaEspecie[1], FFraccionMasicaEspecie[2], 0,
												  FCalculoGamma, nmMEP);
				FCpMezcla = CalculoCompletoCpMezcla(FFraccionMasicaEspecie[0], FFraccionMasicaEspecie[1], FFraccionMasicaEspecie[2], 0,
													__units::degCToK(FTemperature), FCalculoGamma, nmMEP);
				FGamma = CalculoCompletoGamma(FRMezcla, FCpMezcla, FCalculoGamma);

			} else if(FCalculoEspecies == nmCalculoSimple) {

				FRMezcla = CalculoSimpleRMezcla(FFraccionMasicaEspecie[0], FFraccionMasicaEspecie[1], FCalculoGamma, nmMEP);
				FCvMezcla = CalculoSimpleCvMezcla(__units::degCToK(FTemperature), FFraccionMasicaEspecie[0], FFraccionMasicaEspecie[1],
												  FCalculoGamma, nmMEP);
				FGamma = CalculoSimpleGamma(FRMezcla, FCvMezcla, FCalculoGamma);

			}

		}
		FTime = t;
	} catch(exception & N) {
		std::cout << "ERROR: TQ2DTurbine::UpdatePlenum en la turbina " << FNumeroTurbina << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TQ2DTurbine::UpdateStateVector() {
	TIntegrableGroup::UpdateStateVector();
}

void TQ2DTurbine::AcumulaMedias(double Tiempo) {
// TODO
}

void TQ2DTurbine::AsignaEntradaSalidaCC() {
// TODO
}

void TQ2DTurbine::CabeceraResultadosInstantTurb(stringstream& insoutput) {
	WriteInsHeader(insoutput);
}

void TQ2DTurbine::CabeceraResultadosMedTurb(stringstream& medoutput) {
// TODO
}

void TQ2DTurbine::CalculaCondicionTurbina(double TimeCalculo) {
// TODO
}

void TQ2DTurbine::CalculaResultadosMediosTurb() {
// TODO
}

void TQ2DTurbine::ImprimeResultadosInstantTurb(stringstream& insoutput) {
	WriteInsResults(insoutput);
}

void TQ2DTurbine::ImprimeResultadosMediosPantalla() {
// TODO
}

void TQ2DTurbine::ImprimeResultadosMedTurb(stringstream& medoutput) {
// TODO
}

void TQ2DTurbine::IniciaMedias() {
// TODO
}

void TQ2DTurbine::LeeResultadosInstantTurb(const char* FileWAM, fpos_t& filepos) {
// TODO
}

void TQ2DTurbine::LeeResultadosInstantTurbXML(xml_node node) {
	ReadInsResults(node);
}

void TQ2DTurbine::ReadAverageResultsTurb(const char* FileWAM, fpos_t& filepos) {
// TODO
}

void TQ2DTurbine::ReadAverageResultsTurbXML(xml_node node) {
// TODO
}

void TQ2DTurbine::ReadInsResults(const xml_node& node)
{
	try {
		auto node_ins = GetNodeChild(node, "Trb_InsOutput");
		for(auto parameter = node_ins.attribute("Parameter");
			parameter; parameter = parameter.next_attribute()) {
			std::string par = parameter.value();
			if(par == "Power") {
				FResInstantTurbina.Potencia = true;
				FInsOutput = true;
			} else if(par == "Efficiency") {
				FResInstantTurbina.Rendimiento = true;
				FInsOutput = true;
			} else if(par == "BSR") {
				FResInstantTurbina.RelaCinematica = true;
				FInsOutput = true;
			} else if(par == "CorrectedMassFlow") {
				FResInstantTurbina.GastoCorregido = true;
				FInsOutput = true;
			} else if(par == "CorrectedSpeed") {
				FResInstantTurbina.RegimenCorregido = true;
				FInsOutput = true;
			} else if(par == "ExpansionRatio") {
				FResInstantTurbina.RelacionExpansion = true;
				FInsOutput = true;
			} else if(par == "VoluteOutletPressure") {
				FResInstantTurbina.VoluteOutletPressure = true;
				FInsOutput = true;
			} else if(par == "PlenumPressure") {
				FResInstantTurbina.PlenumPressure = true;
				FInsOutput = true;
			} else {
				std::cout << "Instantaneous results in turbine " << FNumeroTurbina
					<< " are not implemented: " << par << std::endl;
			}
		}
	} catch(exception & N) {
		std::cout << "ERROR: TQ2DTurbine::ReadInsResults in turbine " << FNumeroTurbina << std::endl;
		std::cout << "Error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}


void TQ2DTurbine::ResultadosInstantTurb() {
// TODO
}

void TQ2DTurbine::WriteInsHeader(stringstream& insoutput) const {
	try {
		std::string Label;

		if(FResInstantTurbina.Potencia) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4009) + PutLabel(903);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.Rendimiento) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4011) + PutLabel(901);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.RelaCinematica) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4027) + PutLabel(
				901) + "/" + std::to_string(0);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.GastoCorregido) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4024) + PutLabel(
				905) + "/" + std::to_string(0);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.RegimenCorregido) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4025) + PutLabel(
				906) + "/" + std::to_string(0);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.RelacionExpansion) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/" + PutLabel(4026) + PutLabel(
				901) + "/" + std::to_string(0);
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.VoluteOutletPressure) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/Pressure(bar)/Volute";
			insoutput << Label.c_str();
		}
		if(FResInstantTurbina.PlenumPressure) {
			Label = "\t" + PutLabel(5009) + "/" + std::to_string(FNumeroTurbina) + "/Pressure(bar)/Plenum";
			insoutput << Label.c_str();
		}
		// fclose(fich);
	} catch(exception &N) {
		std::cout << "ERROR: TQ2DTurbine::WriteInsHeader in turbine " << FNumeroTurbina << std::endl;
		std::cout << "Error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TQ2DTurbine::WriteInsResults(std::stringstream & output) const {
		try {
		if(FResInstantTurbina.Potencia)
			output << "\t" << FPotencia;
		if(FResInstantTurbina.Rendimiento)
			output << "\t"
				<< (FTurbineEfficiency * FInletMassFlow).sum()
				/ (FInletMassFlow.sum() + 1E-10);
			// Adiabatic [-]
		if(FResInstantTurbina.RelaCinematica) {
			output << "\t"
				<< (FBSR * FInletMassFlow).sum() / (FInletMassFlow.sum() + 1E-10);
			// [-]
		}
		if(FResInstantTurbina.GastoCorregido) {
			output << "\t" << FReducedMassFlow.sum() * 1E6;
			// [kg * sqrt(K) / (s * MPa)]
		}
		if(FResInstantTurbina.RegimenCorregido) {
			output << "\t"
				<< (FReducedSpeed * FInletMassFlow).sum()
				/ (FInletMassFlow.sum() + 1E-10);
			// [rpm / sqrt(K)]
		}
		if(FResInstantTurbina.RelacionExpansion) {
			output << "\t" << FExpansionRatio.mean();
			// Total-to-static [-]
		}
		if(FResInstantTurbina.VoluteOutletPressure) {
			output << "\t" << __units::PaToBar(FQ2DAcTurb->VoluteOutletp().mean());
			// [bar]
		}
		if(FResInstantTurbina.PlenumPressure) {
			output << "\t" << FPressure;
			// [bar]
		}
		// fclose(fich);
	} catch(exception &N) {
		std::cout << "ERROR: TQ2DTurbine::WriteInsResults in turbine "
			<< FNumeroTurbina << std::endl;
		std::cout << "Error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}
