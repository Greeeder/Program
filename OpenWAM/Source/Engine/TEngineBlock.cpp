/**
 * @file TEngineBlock.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 de abr. de 2016
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
 * Includes all properties and methos to solve the physics of the engine block and its components.
 */
#include "TEngineBlock.h"
#include "TCrankMechanism.h"
#include "AllCombustionType.h"

TEngineBlock::TEngineBlock() {
	// TODO Auto-generated constructor stub

}

TEngineBlock::~TEngineBlock() {
	// TODO Auto-generated destructor stub
}

void TEngineBlock::ReadEngineBlockXML(xml_node node_openwam, TComponentArray_ptr WorkFluid){

	xml_node node_engine = GetNodeChild(node_openwam, "NewEngineBlock");

	FEngineType = node_engine.attribute("Type").as_string();
	if (FEngineType == "4Strokes"){
		FAngleCycle = 4 * __cons::Pi;
	}
	else if (FEngineType == "2Strokes"){
		FAngleCycle = __cons::Pi_x_2;
	}
	else{
		cout << "ERROR: Engine type " << FEngineType << " in not define, please select 4Strokes or 2Strokes" << endl;
	}

	FSpeed = GetAttributeAsDouble(node_engine, "Speed");
	FAngle = 0.;

	for (xml_node node_cyl = node_engine.child("Cylinder"); node_cyl; node_cyl = node_cyl.next_sibling("Cylinder")){
		Cylinder.push_back(create_cylinder(node_engine, node_cyl, WorkFluid));
		FMembers.push_back(Cylinder.back());
	}

	FInjectionSystem = make_unique<TInjectionSystem>();
	for (auto cyl : Cylinder){
		FInjectionSystem->addInjector(cyl->getInjector());
	}

	FTotalVolume = CountNodes(node_engine, "Cylinder") * Cylinder[0]->getVolumeDisplaced();

	for (xml_node node_ctrl = node_engine.child("Actuator"); node_ctrl; node_ctrl = node_ctrl.next_sibling("Actuator")) {

		std::string CtrlParam = node_ctrl.attribute("Parameter").value();
		if (CtrlParam == "RPM") {
			FRPMControllerID = IDtoInt(node_ctrl.attribute("CtrlID").as_string());
			FRPMControlled = true;
		}
	}
}

void TEngineBlock::ReadAVGOutputXML(xml_node node_engine) {


	FAverageOutput.ParNeto = false;
	FAverageOutput.ParNetoSUM = 0.;
	FAverageOutput.ParNetoMED = 0.;
	FAverageOutput.PMN = false;
	FAverageOutput.PMNMED = 0.;
	FAverageOutput.ParEfectivo = false;
	FAverageOutput.ParEfectivoSUM = 0.;
	FAverageOutput.ParEfectivoMED = 0.;
	FAverageOutput.ParEfectivoCiclo = false;
	FAverageOutput.ParEfectivoCicloMED = 0.;
	FAverageOutput.PME = false;
	FAverageOutput.PMEMED = 0.;
	FAverageOutput.Potencia = false;
	FAverageOutput.PotenciaMED = 0.;
	FAverageOutput.PotenciaCiclo = false;
	FAverageOutput.PotenciaCicloMED = 0.;
	FAverageOutput.MasaAdmision = false;
	FAverageOutput.MasaAdmisionMED = 0.;
	FAverageOutput.MasaAdmisionSUM = 0.;
	FAverageOutput.MasaFuel = false;
	FAverageOutput.MasaFuelMED = 0.;
	FAverageOutput.MasaFuelSUM = 0.;
	FAverageOutput.RegimenGiro = false;
	FAverageOutput.RegimenGiroSUM = 0.;
	FAverageOutput.RegimenGiroMED = 0.;
	FAverageOutput.RendimientoVolumetrico = false;
	FAverageOutput.RendimientoVolumetricoMED = 0.;
	FAverageOutput.RendimientoVolumetricoAtm = false;
	FAverageOutput.RendimientoVolumetricoAtmMED = 0.;
	FAverageOutput.ParPerdidasMecanicas = false;
	FAverageOutput.ParPerdidasMecanicasSUM = 0.;
	FAverageOutput.ParPerdidasMecanicasMED = 0.;
	FAverageOutput.ParResistente = false;
	FAverageOutput.ParResistenteSUM = 0.;
	FAverageOutput.ParResistenteMED = 0.;
	FAverageOutput.VelocidadVehiculo = false;
	FAverageOutput.VelocidadVehiculoSUM = 0.;
	FAverageOutput.VelocidadVehiculoMED = 0.;
	FAverageOutput.DensidadReferenciaSUM = 0.;
	FAverageOutput.DensidadReferenciaMED = 0.;
	FAverageOutput.MasaTuboReferenciaSUM = 0.;
	FAverageOutput.MasaTuboReferenciaMED = 0.;
	FAverageOutput.GastoTuboReferenciaSUM = 0.;
	FAverageOutput.GastoTuboReferenciaMED = 0.;
	FAverageOutput.MasaAtrapada = false;
	FAverageOutput.MasaAtrapadaMED = 0.;
	FAverageOutput.TrabajoNeto = false;
	FAverageOutput.TrabajoNetoSUM = 0.;
	FAverageOutput.TrabajoNetoMED = 0.;
	FAverageOutput.TrabajoBombeo = false;
	FAverageOutput.TrabajoBombeoSUM = 0.;
	FAverageOutput.TrabajoBombeoMED = 0.;
	FAverageOutput.PMNCiclo = false;
	FAverageOutput.PMNCicloMED = 0.;
	FAverageOutput.PME = false;
	FAverageOutput.PMECicloMED = 0.;
	FAverageOutput.PMBCiclo = false;
	FAverageOutput.PMBCicloMED = 0.;
	FAverageOutput.PMICiclo = false;
	FAverageOutput.PMICicloMED = 0.;
	FAverageOutput.RendEfectivo = false;
	FAverageOutput.RendEfectivoMED = 0.;
	FAverageOutput.RendIndicado = false;
	FAverageOutput.RendIndicadoMED = 0.;
	FAverageOutput.ConsumoEspecifico = false;
	FAverageOutput.ConsumoEspecificoMED = 0.;
	FAverageOutput.Dosado = false;
	FAverageOutput.DosadoMED = 0.;
	FAverageOutput.AFR = false;
	FAverageOutput.AFRMED = 0.;
	FAverageOutput.Swirl = false;
	FAverageOutput.SwirlMED = 0.;

	FAverageOutput.TiempoSUM = 0.;
	FAverageOutput.Tiempo0 = 0.;

	xml_node node_avg = GetNodeChild(node_engine, "Eng_AvgOutput");
	for (xml_attribute parameter = node_avg.attribute("Parameter"); parameter; parameter.next_attribute()) {
		if (parameter.value() == "NetTorque") {
			FAverageOutput.ParNeto = true;
		}
		else if (parameter.value() == "EffectiveTorque") {
			FAverageOutput.ParEfectivo = true;
		}
		else if (parameter.value() == "EffectiveTorque_cycle") {
			FAverageOutput.ParEfectivoCiclo = true;
		}
		else if (parameter.value() == "MechLossesTorque") {
			FAverageOutput.ParPerdidasMecanicas = true;
		}
		else if (parameter.value() == "NetWork") {
			FAverageOutput.TrabajoNeto = true;
		}
		else if (parameter.value() == "PumpingWork") {
			FAverageOutput.TrabajoBombeo = true;
		}
		else if (parameter.value() == "NMEP") {
			FAverageOutput.PMN = true;
		}
		else if (parameter.value() == "BMEP") {
			FAverageOutput.PME = true;
		}
		else if (parameter.value() == "NMEP_cycle") {
			FAverageOutput.PMNCiclo = true;
		}
		else if (parameter.value() == "BMEP_cycle") {
			FAverageOutput.PMECiclo = true;
		}
		else if (parameter.value() == "IMEP_cycle") {
			FAverageOutput.PMICiclo = true;
		}
		else if (parameter.value() == "PMEP_cycle") {
			FAverageOutput.PMBCiclo = true;
		}
		else if (parameter.value() == "Power") {
			FAverageOutput.Potencia = true;
		}
		else if (parameter.value() == "Power_cycle") {
			FAverageOutput.PotenciaCiclo = true;
		}
		else if (parameter.value() == "IntakeMass") {
			FAverageOutput.MasaAdmision = true;
		}
		else if (parameter.value() == "FuelMass") {
			FAverageOutput.MasaFuel = true;
		}
		else if (parameter.value() == "TrappedMass") {
			FAverageOutput.MasaAtrapada = true;
		}
		else if (parameter.value() == "Speed") {
			FAverageOutput.RegimenGiro = true;
		}
		else if (parameter.value() == "VolumetricEfficiency") {
			FAverageOutput.RendimientoVolumetrico = true;
		}
		else if (parameter.value() == "VolumetricEfficiency_amb") {
			FAverageOutput.RendimientoVolumetricoAtm = true;
		}
		else if (parameter.value() == "EffectiveEfficiency") {
			FAverageOutput.RendEfectivo = true;
		}
		else if (parameter.value() == "IndicatedEfficiency") {
			FAverageOutput.RendIndicado = true;
		}
		else if (parameter.value() == "BSFC") {
			FAverageOutput.ConsumoEspecifico = true;
		}
		else if (parameter.value() == "ResistantTorque") {
			FAverageOutput.ParResistente = true;
		}
		else if (parameter.value() == "VehicleSpeed") {
			FAverageOutput.VelocidadVehiculo = true;
		}
		else if (parameter.value() == "FueltoAirRatio") {
			FAverageOutput.Dosado = true;
		}
		else if (parameter.value() == "AFR") {
			FAverageOutput.AFR = true;
		}
		else if (parameter.value() == "Swirl") {
			FAverageOutput.Swirl = true;
		}
		else {
			std::cout << "Average output in the engine (" << parameter.value() << ") not defined" << std::endl;
		}
	}

}

void TEngineBlock::HeaderAVGOutput(stringstream& medoutput) {

	std::string Label;

	if (FAverageOutput.ParNeto) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4016) + PutLabel(917) + "/" + PutLabel(3001);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.ParEfectivo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4016) + PutLabel(917) + "/" + PutLabel(3029) + "/" + PutLabel(3032);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.ParEfectivoCiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4016) + PutLabel(917) + "/" + PutLabel(3029) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.ParPerdidasMecanicas) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4016) + PutLabel(917) + "/" + PutLabel(3033);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.TrabajoNeto) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4021) + PutLabel(907) + "/" + PutLabel(3001);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.TrabajoBombeo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4021) + PutLabel(907) + "/" + PutLabel(3003);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PMN) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3002) + "/" + PutLabel(3032);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PME) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3034) + "/" + PutLabel(3032);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PMNCiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3002) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PMECiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3034) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PMICiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3009) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.PMBCiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4006) + PutLabel(908) + "/" + PutLabel(3004) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.Potencia) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4009) + PutLabel(903) + "/" + PutLabel(3032);
	}
	if (FAverageOutput.PotenciaCiclo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4009) + PutLabel(903) + "/" + PutLabel(4020);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.MasaAdmision) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4004) + PutLabel(913) + "/" + PutLabel(3017);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.MasaFuel) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4004) + PutLabel(913) + "/" + PutLabel(3027);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.MasaAtrapada) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4004) + PutLabel(913) + "/" + PutLabel(3010);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.RegimenGiro) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4022) + PutLabel(918);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.RendimientoVolumetrico) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4011) + PutLabel(901) + "/" + PutLabel(3022) + "/" + PutLabel(3017);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.RendimientoVolumetricoAtm) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4011) + PutLabel(901) + "/" + PutLabel(3022) + "/" + PutLabel(3035);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.RendEfectivo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4011) + PutLabel(901) + "/" + PutLabel(3029);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.RendIndicado) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(4011) + PutLabel(901) + "/" + PutLabel(3036);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.ConsumoEspecifico) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3037) + PutLabel(924);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.ParResistente) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3038) + PutLabel(917);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.VelocidadVehiculo) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3039) + PutLabel(925);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.Dosado) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3040) + PutLabel(901);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.AFR) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3015) + PutLabel(901);
		medoutput << Label.c_str();
	}
	if (FAverageOutput.Swirl) {
		Label = "\t" + PutLabel(5002) + "/" + PutLabel(3021) + PutLabel(901);
		medoutput << Label.c_str();
	}

}

void TEngineBlock::IntegrateAVGOutput(double TActual, int CilindroActual) {

	double partotalins = 0.;

	double DeltaT = TActual - FAverageOutput.Tiempo0;
	// double DeltaAngulo=360.*FRegimen/60.*DeltaT;

	if (FAverageOutput.ParNeto || FAverageOutput.PMN || FAverageOutput.ParEfectivo || FAverageOutput.PME
		|| FAverageOutput.Potencia) {
		for (int i = 0; i < Cylinder.size(); i++) {
			partotalins += Cylinder[i]->getNetTorque() * DeltaT;
		}
		FAverageOutput.ParNetoSUM += partotalins;
	}

	//if (FAverageOutput.ParEfectivo || FAverageOutput.PME || FAverageOutput.Potencia) {
	//	FAverageOutput.ParEfectivoSUM += partotalins - FParPerdidasMecanicas * DeltaT;
	//}

	//if (FCilindro[0]->getAnguloActual() > FCilindro[0]->getDistribucion().CA
	//	&& FCilindro[0]->getAnguloActual() <= FCilindro[0]->getDistribucion().CA + (FCilindro[0]->getAnguloActual() -
	//	FCilindro[0]->getAnguloAnterior())) {
	//	if (FPrimeravezAcumulaFuel) {
	//		for (int i = 0; i < FGeom.NCilin; i++) {
	//			FAverageOutput.MasaFuelSUM += FCilindro[i]->getMasaFuel();
	//		}
	//		FPrimeravezAcumulaFuel = false;
	//	}
	//}

	//if (FAverageOutput.RendimientoVolumetrico || FAverageOutput.Dosado) {
	//	FAverageOutput.DensidadReferenciaSUM += FTuboRendVol->GetDensidad(FNodoMedio) * DeltaT;
	//	FAverageOutput.GastoTuboReferenciaSUM += (FTuboRendVol->GetVelocidad(FNodoMedio) * __cons::ARef *
	//		FTuboRendVol->GetArea(FNodoMedio) * FTuboRendVol->GetDensidad(FNodoMedio)) * DeltaT;

	//}

	if (FAverageOutput.RendimientoVolumetrico || FAverageOutput.Potencia || FAverageOutput.RendimientoVolumetricoAtm) {
		FAverageOutput.RegimenGiroSUM += FSpeed * DeltaT;
	}

	//if (FAverageOutput.ParPerdidasMecanicas) {
	//	FAverageOutput.ParPerdidasMecanicasSUM += FParPerdidasMecanicas * DeltaT;
	//}
	//if (FAverageOutput.ParResistente) {
	//	FAverageOutput.ParResistenteSUM += FParResistente * DeltaT;
	//}
	//if (FAverageOutput.VelocidadVehiculo) {
	//	FAverageOutput.VelocidadVehiculoSUM += FVelocidadVehiculo * DeltaT;
	//}
	//if (FAverageOutput.Dosado) {
	//	FAverageOutput.MasaTuboReferenciaSUM += (FTuboRendVol->GetVelocidad(FNodoMedio) * __cons::ARef * FTuboRendVol->GetArea(
	//		FNodoMedio) * FTuboRendVol->GetDensidad(FNodoMedio)) / FGeom.NCilin
	//		/ (FRegimen / 120) * DeltaT;
	//}

	/* for(int i=0;i<FGeom.NCilin;i++){
	FAverageOutput.TrabajoNetoSUM+=FCilindro[i]->getPreMed()*1e5*(FCilindro[i]->getVolumen()-FCilindro[i]->getVolumen0());
	if(FCilindro[i]->getAnguloActual()>180. && FCilindro[i]->getAnguloActual()<540.){
	FAverageOutput.TrabajoBombeoSUM+=FCilindro[i]->getPreMed()*1e5*(FCilindro[i]->getVolumen()-FCilindro[i]->getVolumen0());
	}
	} */

	FAverageOutput.TrabajoNetoSUM += Cylinder[CilindroActual - 1]->getWork();

	//if (FCilindro[CilindroActual - 1]->getAnguloActual() > 180. && FCilindro[CilindroActual - 1]->getAnguloActual() < 540.) {
	//	FAverageOutput.TrabajoBombeoSUM += __units::BarToPa(FCilindro[CilindroActual - 1]->getPreMed())
	//		* (FCilindro[CilindroActual - 1]->getVolumen() - FCilindro[CilindroActual - 1]->getVolumen0());
	//}

	FAverageOutput.TiempoSUM += DeltaT;
	FAverageOutput.Tiempo0 = TActual;

}


void TEngineBlock::ComputeAVGOutput() {

	double DensidadAtm = 0.;
	double MasaAtrapadaSUM = 0.;
	double FraccionAireFrescoSUM = 0., AFRSUM = 0.;
	double swirltotal = 0.;

	if (FAverageOutput.RegimenGiro || FAverageOutput.Potencia || FAverageOutput.RendimientoVolumetricoAtm
		|| FAverageOutput.RendimientoVolumetrico) {
		FAverageOutput.RegimenGiroMED = FAverageOutput.RegimenGiroSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.RegimenGiroSUM = 0.;
	}
	if (FAverageOutput.ParNeto || FAverageOutput.PMN) {
		FAverageOutput.ParNetoMED = FAverageOutput.ParNetoSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.ParNetoSUM = 0.;
	}
	if (FAverageOutput.PMN) {
		FAverageOutput.PMNMED = __units::PaToBar(FAverageOutput.ParNetoMED * 4. * __cons::Pi) / FTotalVolume;
	}
	if (FAverageOutput.ParEfectivo || FAverageOutput.PME || FAverageOutput.Potencia) {
		FAverageOutput.ParEfectivoMED = FAverageOutput.ParEfectivoSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.ParEfectivoSUM = 0.;
	}
	if (FAverageOutput.PME) {
		FAverageOutput.PMEMED = __units::PaToBar(FAverageOutput.ParEfectivoMED * 4. * __cons::Pi) / FTotalVolume;
	}
	if (FAverageOutput.Potencia) {
		FAverageOutput.PotenciaMED = __units::To_kilo(FAverageOutput.ParEfectivoMED * __cons::Pi_x_2 * __units::RPMToRPS(
			FAverageOutput.RegimenGiroMED));
	}
	if (FAverageOutput.MasaAdmision) {
		//for (int i = 0; i < FGeom.NCilin; i++) {
		//	FAverageOutput.MasaAdmisionSUM += FCilindro[i]->getMasaPorAdmision();
		//}
		//FAverageOutput.MasaAdmisionMED = FAverageOutput.MasaAdmisionSUM / FGeom.NCilin;
		//FAverageOutput.MasaAdmisionSUM = 0.;
	}

	//FAverageOutput.MasaFuelMED = FAverageOutput.MasaFuelSUM / FGeom.NCilin;
	//FAverageOutput.MasaFuelSUM = 0.;
	//FPrimeravezAcumulaFuel = true;

	//if (FAverageOutput.MasaAtrapada || FAverageOutput.RendimientoVolumetrico || FAverageOutput.RendimientoVolumetricoAtm
	//	|| FAverageOutput.AFR) {
	//	for (int i = 0; i < FGeom.NCilin; i++) {
	//		MasaAtrapadaSUM += FCilindro[i]->getMasaAtrapada();
	//		FraccionAireFrescoSUM += FCilindro[i]->getFraccionAireFresco();
	//	}
	//	FAverageOutput.MasaAtrapadaMED = MasaAtrapadaSUM / FGeom.NCilin;
	//	FAverageOutput.FraccionAireFrescoMED = FraccionAireFrescoSUM / FGeom.NCilin;
	//}
	if (FAverageOutput.RendimientoVolumetrico || FAverageOutput.RendimientoVolumetricoAtm) {
		FAverageOutput.GastoTuboReferenciaMED = FAverageOutput.GastoTuboReferenciaSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.GastoTuboReferenciaSUM = 0.;
	}
	if (FAverageOutput.RendimientoVolumetrico) {
		FAverageOutput.DensidadReferenciaMED = FAverageOutput.DensidadReferenciaSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.RendimientoVolumetricoMED = fabs(FAverageOutput.GastoTuboReferenciaMED) /
			FAverageOutput.DensidadReferenciaMED / FTotalVolume
			/ (FAverageOutput.RegimenGiroMED / 60. * (360. / FAngleCycle));
		FAverageOutput.DensidadReferenciaSUM = 0;
	}
	if (FAverageOutput.ParPerdidasMecanicas) {
		FAverageOutput.ParPerdidasMecanicasMED = FAverageOutput.ParPerdidasMecanicasSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.ParPerdidasMecanicasSUM = 0.;
	}
	if (FAverageOutput.ParResistente) {
		FAverageOutput.ParResistenteMED = FAverageOutput.ParResistenteSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.ParResistenteSUM = 0.;
	}
	if (FAverageOutput.VelocidadVehiculo) {
		FAverageOutput.VelocidadVehiculoMED = FAverageOutput.VelocidadVehiculoSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.VelocidadVehiculoSUM = 0.;
	}
	if (FAverageOutput.RendimientoVolumetricoAtm) {
		DensidadAtm = __Ambient::p_bar / (__R::Air * __Ambient::T_K);
		FAverageOutput.RendimientoVolumetricoAtmMED = fabs(FAverageOutput.GastoTuboReferenciaMED) / DensidadAtm /
			FTotalVolume / (FAverageOutput.RegimenGiroMED / 120.);
	}
	if (FAverageOutput.Dosado) {
		FAverageOutput.MasaTuboReferenciaMED = FAverageOutput.MasaTuboReferenciaSUM / FAverageOutput.TiempoSUM;
		FAverageOutput.DosadoMED = FAverageOutput.MasaFuelMED / fabs(FAverageOutput.MasaTuboReferenciaMED);
		FAverageOutput.MasaTuboReferenciaSUM = 0.;
	}
	if (FAverageOutput.AFR) {
		//for (int i = 0; i < FGeom.NCilin; i++) {
		//	AFRSUM += FCilindro[i]->getAFR();
		//}
		//FAverageOutput.AFRMED = AFRSUM / FGeom.NCilin;
		//AFRSUM = 0.;
		//// FAverageOutput.AFRMED=FAverageOutput.MasaAtrapadaMED*FAverageOutput.FraccionAireFrescoMED/FAverageOutput.MasaFuelMED;
	}
	if (FAverageOutput.Swirl) {
		//for (int i = 0; i < FGeom.NCilin; i++) {
		//	swirltotal += FCilindro[i]->getSwirlSUM() / FAverageOutput.TiempoSUM;
		//	/* Valor medio de Swirl para el cilindro i */
		//}
		//FAverageOutput.SwirlMED = swirltotal / FGeom.NCilin;
	}

	FAverageOutput.TrabajoNetoMED = FAverageOutput.TrabajoNetoSUM;
	FAverageOutput.TrabajoNetoSUM = 0.;

	/* for(int i=0;i<FGeom.NCilin;i++){
	FAverageOutput.TrabajoBombeoSUM+=FCilindro[i]->getTrabajoBombeo();
	} */

	FAverageOutput.TrabajoBombeoMED = FAverageOutput.TrabajoBombeoSUM;
	FAverageOutput.TrabajoBombeoSUM = 0.;

	FAverageOutput.PMNCicloMED = __units::PaToBar(FAverageOutput.TrabajoNetoMED / FTotalVolume);

	FAverageOutput.PMBCicloMED = __units::PaToBar(FAverageOutput.TrabajoBombeoMED / FTotalVolume);

	FAverageOutput.PMICicloMED = FAverageOutput.PMNCicloMED - FAverageOutput.PMBCicloMED;

	//FPMPMMotor = FPerMec.Coef0 + FPerMec.Coef1 * FAverageOutput.RegimenGiroMED / 60. - FPerMec.Coef2 *
	//	FAverageOutput.RegimenGiroMED * FAverageOutput.RegimenGiroMED / 3600
	//	+ FPerMec.Coef3 * FAverageOutput.PMICicloMED;
	//FAverageOutput.PMECicloMED = FAverageOutput.PMNCicloMED - FPMPMMotor;

	FAverageOutput.ParEfectivoCicloMED = __units::BarToPa(FAverageOutput.PMECicloMED) * FTotalVolume /
		(__units::DegToRad(FAngleCycle));

	FAverageOutput.PotenciaCicloMED = __units::To_kilo(__cons::Pi_x_2 * FAverageOutput.ParEfectivoCicloMED *
		__units::RPMToRPS(FAverageOutput.RegimenGiroMED));

	if (FAverageOutput.MasaFuelMED == 0) {
		FAverageOutput.RendIndicadoMED = 0.;
		FAverageOutput.RendEfectivoMED = 0.;
		FAverageOutput.ConsumoEspecificoMED = 0.;
	}
	else {
		//FAverageOutput.RendIndicadoMED = (FAverageOutput.TrabajoNetoMED - FAverageOutput.TrabajoBombeoMED) /
		//	(FGeom.NCilin * FAverageOutput.MasaFuelMED * FPoderCalorifico);
		//FAverageOutput.RendEfectivoMED = __units::BarToPa(FAverageOutput.PMECicloMED) * FGeom.CilindradaTotal /
		//	(FGeom.NCilin * FAverageOutput.MasaFuelMED * FPoderCalorifico);
		//FAverageOutput.ConsumoEspecificoMED = 3.6e9 / (FAverageOutput.RendEfectivoMED * FPoderCalorifico);
	}

	FAverageOutput.TiempoSUM = 0;

}

void TEngineBlock::PrintAVGOutput(stringstream& medoutput) {

	// FILE *fich=fopen(FileSALIDA,"a");

	if (FAverageOutput.ParNeto)
		medoutput << "\t" << FAverageOutput.ParNetoMED;
	if (FAverageOutput.ParEfectivo)
		medoutput << "\t" << FAverageOutput.ParEfectivoMED;
	if (FAverageOutput.ParEfectivoCiclo)
		medoutput << "\t" << FAverageOutput.ParEfectivoCicloMED;
	if (FAverageOutput.ParPerdidasMecanicas)
		medoutput << "\t" << FAverageOutput.ParPerdidasMecanicasMED;
	if (FAverageOutput.TrabajoNeto)
		medoutput << "\t" << FAverageOutput.TrabajoNetoMED;
	if (FAverageOutput.TrabajoBombeo)
		medoutput << "\t" << FAverageOutput.TrabajoBombeoMED;
	if (FAverageOutput.PMN)
		medoutput << "\t" << FAverageOutput.PMNMED;
	if (FAverageOutput.PME)
		medoutput << "\t" << FAverageOutput.PMEMED;
	if (FAverageOutput.PMNCiclo)
		medoutput << "\t" << FAverageOutput.PMNCicloMED;
	if (FAverageOutput.PMECiclo)
		medoutput << "\t" << FAverageOutput.PMECicloMED;
	if (FAverageOutput.PMICiclo)
		medoutput << "\t" << FAverageOutput.PMICicloMED;
	if (FAverageOutput.PMBCiclo)
		medoutput << "\t" << FAverageOutput.PMBCicloMED;
	if (FAverageOutput.Potencia)
		medoutput << "\t" << FAverageOutput.PotenciaMED;
	if (FAverageOutput.PotenciaCiclo)
		medoutput << "\t" << FAverageOutput.PotenciaCicloMED;
	if (FAverageOutput.MasaAdmision)
		medoutput << "\t" << FAverageOutput.MasaAdmisionMED;
	if (FAverageOutput.MasaFuel)
		medoutput << "\t" << FAverageOutput.MasaFuelMED;
	if (FAverageOutput.MasaAtrapada)
		medoutput << "\t" << FAverageOutput.MasaAtrapadaMED;
	if (FAverageOutput.RegimenGiro)
		medoutput << "\t" << FAverageOutput.RegimenGiroMED;
	if (FAverageOutput.RendimientoVolumetrico)
		medoutput << "\t" << FAverageOutput.RendimientoVolumetricoMED;
	if (FAverageOutput.RendimientoVolumetricoAtm)
		medoutput << "\t" << FAverageOutput.RendimientoVolumetricoAtmMED;
	if (FAverageOutput.RendEfectivo)
		medoutput << "\t" << FAverageOutput.RendEfectivoMED;
	if (FAverageOutput.RendIndicado)
		medoutput << "\t" << FAverageOutput.RendIndicadoMED;
	if (FAverageOutput.ConsumoEspecifico)
		medoutput << "\t" << FAverageOutput.ConsumoEspecificoMED;
	if (FAverageOutput.ParResistente)
		medoutput << "\t" << FAverageOutput.ParResistenteMED;
	if (FAverageOutput.VelocidadVehiculo)
		medoutput << "\t" << FAverageOutput.VelocidadVehiculoMED;
	if (FAverageOutput.Dosado)
		medoutput << "\t" << FAverageOutput.DosadoMED;
	if (FAverageOutput.AFR)
		medoutput << "\t" << FAverageOutput.AFRMED;
	if (FAverageOutput.Swirl)
		medoutput << "\t" << FAverageOutput.SwirlMED;

	// fclose(fich);
}

void TEngineBlock::SetRPMController(vector<shared_ptr<TController>> Controller) {
	if (FRPMControlled)
		FRPMController = Controller[FRPMControllerID - 1];
}

void TEngineBlock::Solve(){

	FAngleStep = __units::RPMToRad_s(FSpeed) * FTimeStep;

	UpdateWorkingPoint();

	TIntegrableGroup::Solve();

	FAngle += FAngleStep;

}

void TEngineBlock::UpdateWorkingPoint(){
	
	for (auto cil : Cylinder){
		cil->setEngineSpeed(FSpeed);
		cil->setAngle(FAngle);
		cil->setDeltaAngle(FAngleStep);
		cil->setTimeStep(FTimeStep);
	}
}

void TEngineBlock::ReadOutput(const xml_node & node)
{
	string param;
	bool Avg;
	auto node_out = GetNodeChild(node, "Output");
	//FPlenumOutput = PlenumOutput();
	for (auto node_par = GetNodeChild(node_out, "Parameter"); node_par;
		node_par = node_par.next_sibling("Parameter")) {
		Is_plot = true;
		param = node_par.attribute("Name").as_string();
		Avg = node_par.attribute("Average").as_bool();
		if (Avg)
			Is_cycleout = true;
		if (param == "Mass") {}
	}
}

void TEngineBlock::HeaderForCycleOutput(std::vector<std::string>& output) const
{
}

void TEngineBlock::HeaderForPlots(std::vector<std::string>& output) const
{
}

void TEngineBlock::IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const
{
}

void TEngineBlock::StoreOutput(std::vector<std::vector<float>>& output) const
{
}

EngineBlockOuput::EngineBlockOuput()
{
}
