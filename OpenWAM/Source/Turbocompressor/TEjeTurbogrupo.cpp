/* --------------------------------------------------------------------------------*\
==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
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


 \*-------------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------
#pragma hdrstop

#include "TEjeTurbogrupo.h"

#ifdef __BORLANDC__
#include <vcl.h>
#endif
#include "TCompresor.h"
#include "TCompTubDep.h"
#include "TTurbina.h"
#include "fluids.h"
#include "TFluidOil.h"

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TEjeTurbogrupo::TEjeTurbogrupo(int i, int ncilin) {

	FNumeroEje = i + 1;
	FNumeroCompresor = NULL;
	FCompresor = NULL;
	FNumeroTurbina = NULL;
	FTurbina = NULL;
	FNumCilindros = ncilin;

	FResMediosEje.Regimen = false;
	FNumCiclo = 0;

	FRPMControlled = false;
	FTime = 0;

	FThereIsHTM = false;

//	FHTM = NULL;
	FHTM = NULL;

//	FMechLosses = NULL;
	FMechLosses = NULL;
	FMechPower = 0;

	FConvTemp.Time0 = 0;

	FConvTemp.T1u = 0;
	FConvTemp.T1d = 0;

	FConvTemp.T2u = 0;
	FConvTemp.T2d = 0;

	FConvTemp.T3u = 0;
	FConvTemp.T3d = 0;

	FConvTemp.MassC = 0;
	FConvTemp.MassT = 0;
	FConvTemp.HeatC1 = 0;
	FConvTemp.HeatC2 = 0;
	FConvTemp.TimeSum = 0;
	FConvTemp.MechP = 0.;
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

TEjeTurbogrupo::~TEjeTurbogrupo() {

	if(FNumeroCompresor != NULL)
		delete[] FNumeroCompresor;
	if(FCompresor != NULL)
		delete[] FCompresor;

	if(FNumeroTurbina != NULL)
		delete[] FNumeroTurbina;
	if(FTurbina != NULL)
		delete[] FTurbina;

	//if(FHTM != NULL)
	//	delete FHTM;
	if(FHTM != NULL)
		delete FHTM;

	if(FMechLosses != NULL)
		delete FMechLosses;

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ReadTurbochargerAxis(const char *FileWAM, fpos_t &filepos, TCompresor **Compressor,
		TTurbina **Turbine) {
	try {
		int variacion = 0, htm = 0;

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		// Lectura del regimen inicial y si se mantiene constante durante el calculo.
		fscanf(fich, "%lf %d", &FRegimenEje, &variacion);
		switch(variacion) {
		case 0:
			FVariacionRegimen = nmVariable;
			break;
		case 1:
			FVariacionRegimen = nmFixed;
			break;
		default:
			std::cout << "ERROR: Reading turbocharger speed variation in axis: " << FNumeroEje << std::endl;
			throw Exception("");
		}
		if(FVariacionRegimen == nmVariable) {
			fscanf(fich, "%lf ", &FMomentoInercia);
		}

		fscanf(fich, "%d ", &FNumCompresoresAcoplados);
		FCompresor = new TCompresor*[FNumCompresoresAcoplados];
		FNumeroCompresor = new int[FNumCompresoresAcoplados];
		for(int i = 0; i < FNumCompresoresAcoplados; ++i) {
			// Identifica que compresores estan acoplados a este eje.
			fscanf(fich, "%d ", &FNumeroCompresor[i]);
			FCompresor[i] = Compressor[FNumeroCompresor[i] - 1];
			FCompresor[i]->IniciaMedias();
		}

		fscanf(fich, "%d ", &FNumTurbinasAcopladas);
		FTurbina = new TTurbina*[FNumTurbinasAcopladas];
		FNumeroTurbina = new int[FNumTurbinasAcopladas];
		for(int i = 0; i < FNumTurbinasAcopladas; ++i) {
			// Identifica que turbinas estan acopladas a este eje.
			fscanf(fich, "%d ", &FNumeroTurbina[i]);
			// numidturbina es un numero que necesita WAMer y que a nosotros no nos interesa.
			FTurbina[i] = Turbine[FNumeroTurbina[i] - 1];
			FTurbina[i]->PutRegimen(FRegimenEje);
		}

		fscanf(fich, "%d ", &FControllerID);
		if(FControllerID > 0)
			FRPMControlled = true;

		// basura para wamer
		int numide = 0;
		fscanf(fich, "%d ", &numide);

		fscanf(fich, "%d ", &htm);
		if(htm == 1) {

			if(FNumTurbinasAcopladas != 1 || FNumCompresoresAcoplados != 1) {
				std::cout << "ERROR: Turbocharger heat trasfer model is not adapted for more than one turbine" << std::endl;
				std::cout << "       or more than one compressor in the same turbocharger" << std::endl;
			}
			FThereIsHTM = true;
			// Geometrical parameters journal bearing
			fscanf(fich, "%lf %lf %lf ", &FDShaft, &FJournalBLengh, &FHD);
			// Geometrical parameters thrust bearing
			fscanf(fich, "%lf %lf", &FTthrustBRmin, &FTthrustBRmax);
			// Geometrical parameters inlet ports
			fscanf(fich, "%lf %lf", &FDoil, &FDwater);
			// Fitting coefficients
			fscanf(fich, "%lf %lf %lf %lf %lf", &FJournalB_K, &Fk_m, &Fk_tb, &FCAC, &FCAT);
			// Wheel areas
			fscanf(fich, "%lf %lf", &FCWArea, &FTWArea);
			double DT = 0., LT = 0.;
			fscanf(fich, "%lf %lf", &DT, &LT);
			double DC = 0., LC = 0.;
			fscanf(fich, "%lf %lf", &DC, &LC);
			double DH = 0., LH = 0.;
			fscanf(fich, "%lf %lf", &DH, &LH);

			// Oil properties.
			fscanf(fich, "%lf %lf %lf ", &FMoil, &FToil, &FPoil);
			double K1 = 0., K2 = 0., K3 = 0.;
			// Oil Voeguel parameters.
			fscanf(fich, "%lf %lf %lf ", &K1, &K2, &K3);
			Oil = new TFluidOil();
			dVector K;
			K.push_back(K1);
			K.push_back(K2);
			K.push_back(K3);
			Oil->GetComponent(0)->setCoeficients("Viscosity", K);

			// Water properties in turbocharger.
			fscanf(fich, "%lf %lf ", &FTwater, &FMwater);
//			FWater = new stHTMwater();

			FHTM = new TTC_HTM2(Oil);

			int Conv = 0;
			fscanf(fich, "%d ", &Conv);
			FConvTemp.FastConvergence = false;
			if(Conv == 1) {
				FConvTemp.FastConvergence = true;
				fscanf(fich, "%lf %lf ", &FConvTemp.TimeStep, &FConvTemp.MaxTime);
				FConvTemp.NewTime = FConvTemp.TimeStep;
				if(FConvTemp.MaxTime == 0.)
					FConvTemp.FastConvergence = false;
			}

			FMechLosses = new TurboBearings(Oil, FJournalBLengh, FDShaft / 2, FHD, FJournalB_K, FCAC, FCAT, FCWArea, FTWArea, Fk_m,
											FTthrustBRmin, FTthrustBRmax, Fk_tb);

			//TODO or not
			//FHTM->Read_HTM(fich);

			//FHTM->TurbochargerData(FDShaft, FHD, FDoil, FDwater, DT, LT, DC, LC, DH, LH);

		}

		fgetpos(fich, &filepos);
		fclose(fich);

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadTurbochargerAxis in the boundary condition: " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TEjeTurbogrupo::ReadTurbochargerAxisXML(xml_node node_tch, TCompresor **Compressor, TTurbina **Turbine) {
	try {
		int htm;
		int i;

		FRegimenEje = GetXMLRotationalSpeed(node_tch, "Speed");
		if(GetAttributeAsBool(node_tch, "FixedSpeed")) {
			FVariacionRegimen = nmFixed;
		} else {
			FVariacionRegimen = nmVariable;
			FMomentoInercia = GetXMLInertia(node_tch, "Inertia");
		}
		FNumCompresoresAcoplados = CountNodes(node_tch, "Tch_Compressor");
		FCompresor = new TCompresor*[FNumCompresoresAcoplados];
		FNumeroCompresor = new int[FNumCompresoresAcoplados];
		i = 0;
		for(xml_node node_cmp = GetNodeChild(node_tch, "Tch_Compressor"); node_cmp;
			node_cmp = node_cmp.next_sibling("Tch_Compressor")) {
			FNumeroCompresor[i] = GetAttributeAsInt(node_cmp, "Compressor_ID");
			FCompresor[i] = Compressor[FNumeroCompresor[i] - 1];
			FCompresor[i]->IniciaMedias();
		}
		FNumTurbinasAcopladas = CountNodes(node_tch, "Tch_Turbine");
		FTurbina = new TTurbina*[FNumTurbinasAcopladas];
		FNumeroTurbina = new int[FNumTurbinasAcopladas];
		i = 0;
		for(xml_node node_cmp = GetNodeChild(node_tch, "Tch_Turbine"); node_cmp;
			node_cmp = node_cmp.next_sibling("Tch_Turbine")) {
			FNumeroTurbina[i] = GetAttributeAsInt(node_cmp, "Turbine_ID");
			FTurbina[i] = Turbine[FNumeroTurbina[i] - 1];
			FTurbina[i]->PutRegimen(FRegimenEje);
		}

		std::string Parameter;

		for(xml_node node_act = node_tch.child("Actuator"); node_act; node_act = node_act.next_sibling("Actuator")) {
			Parameter = node_act.attribute("Parameter").value();
			if(Parameter == "Speed") {
				FControllerID = GetAttributeAsInt(node_act, "CtrlID");
				FRPMControlled = true;
			}
		}

		if(node_tch.child("Tch_HTML")) {

			xml_node node_htm = GetNodeChild(node_tch, "Tch_HTML");
			if(FNumTurbinasAcopladas != 1 || FNumCompresoresAcoplados != 1) {
				std::cout << "ERROR: Turbocharger heat transfer model is not adapted for more than one turbine" << std::endl;
				std::cout << "       or more than one compressor in the same shaft" << std::endl;
			}
			FThereIsHTM = true;

			xml_node node_bear = GetNodeChild(node_htm, "Htm_Bearings");
			FDShaft = GetXMLLength(node_bear, "ShaftDiameter");
			FJournalBLengh = GetXMLLength(node_bear, "JournalBearingLength");
			FHD = GetXMLLength(node_bear, "Clearence");
			FTthrustBRmin = GetXMLLength(node_bear, "ThrustBearingMinR");
			FTthrustBRmax = GetXMLLength(node_bear, "ThrustBearingMaxR");

			xml_node node_flu = GetNodeChild(node_htm, "Htm_Fluids");
			xml_node node_oil = GetNodeChild(node_flu, "Flu_Oil");

			FDoil = GetXMLLength(node_oil, "PortDiameter");
			FMoil = GetXMLMassFlow(node_oil, "MassFlow");
			FToil = GetXMLTemperature(node_oil, "Temperature");
			FPoil = GetXMLPressure(node_oil, "Pressure");

			double K1 = 0., K2 = 0., K3 = 0.;
			K1 = GetAttributeAsDouble(node_oil, "VoeguelK1");
			K2 = GetAttributeAsDouble(node_oil, "VoeguelK2");
			K3 = GetAttributeAsDouble(node_oil, "VoeguelK3");



//			FOil = new stHTMoil();
			Oil = new TFluidOil();

			dVector K;
			K.push_back(GetAttributeAsDouble(node_oil, "VoeguelK1"));
			K.push_back(GetAttributeAsDouble(node_oil, "VoeguelK2"));
			K.push_back(GetAttributeAsDouble(node_oil, "VoeguelK3"));
			Oil->GetComponent(0)->setCoeficients("Viscosity", K);

			//FOil->mu_c1 = K1;
			//FOil->mu_c2 = K2;
			//FOil->mu_c3 = K3;

			if(node_flu.child("Flu_Water")) {
				xml_node node_water = GetNodeChild(node_flu, "Flu_Water");

				FDwater = GetXMLLength(node_water, "PortDiameter");
				FMwater = GetXMLMassFlow(node_water, "MassFlow");
				FTwater = GetXMLTemperature(node_water, "Temperature");

//				FWater = new stHTMwater();
			}

			xml_node node_ml = GetNodeChild(node_htm, "Htm_MechLosses");

			FJournalB_K = GetAttributeAsDouble(node_ml, "JournalB_K");
			Fk_m = GetAttributeAsDouble(node_ml, "K_m");
			Fk_tb = GetAttributeAsDouble(node_ml, "K_tb");
			FCAC = GetAttributeAsDouble(node_ml, "CAC");
			FCAT = GetAttributeAsDouble(node_ml, "CAT");

			xml_node node_geom = GetNodeChild(node_htm, "Htm_Geometry");

			FCWArea = GetXMLArea(node_geom, "CompWheelArea");
			FTWArea = GetXMLArea(node_geom, "TurbWheelArea");
			double DT = 0., LT = 0.;
			DT = GetXMLLength(node_geom, "TurbExtDiameter");
			LT = GetXMLLength(node_geom, "TurbExtLength");
			double DC = 0., LC = 0.;
			DC = GetXMLLength(node_geom, "CompExtDiameter");
			LC = GetXMLLength(node_geom, "CompExtLength");
			double DH = 0., LH = 0.;
			DH = GetXMLLength(node_geom, "HousExtDiameter");
			LH = GetXMLLength(node_geom, "HousExtLength");

			//FMechLosses = new TurboBearings(FOil, FJournalBLengh, FDShaft / 2, FHD, FJournalB_K, FCAC, FCAT, FCWArea, FTWArea, Fk_m,
			//								FTthrustBRmin, FTthrustBRmax, Fk_tb);

			FMechLosses = new TurboBearings(Oil, FJournalBLengh, FDShaft / 2, FHD, FJournalB_K, FCAC, FCAT, FCWArea, FTWArea, Fk_m,
											FTthrustBRmin, FTthrustBRmax, Fk_tb);

			//FHTM = new TTC_HTM(FOil);

			FHTM = new TTC_HTM2(Oil);

			xml_node node_ht = GetNodeChild(node_htm, "Htm_HeatTransfer");

			xml_node node_set = GetNodeChild(node_ht, "Htr_Settings");
			FConvTemp.FastConvergence = GetAttributeAsBool(node_set, "FastConvergence");
			if(FConvTemp.FastConvergence) {
				FConvTemp.TimeStep = GetXMLTime(node_set, "TimeStep");
				FConvTemp.MaxTime = GetXMLTime(node_set, "MaxTime");
				FConvTemp.NewTime = FConvTemp.TimeStep;
			}

			//FHTM->Read_HTMXML(node_ht);

			FHTM->Read_HTMXML(node_ht);

			//FHTM->TurbochargerData(FDShaft, FHD, FDoil, FDwater, DT, LT, DC, LC, DH, LH);

			FHTM->TurbochargerData(FDoil, FDwater, DT, LT, DC, LC, DH, LH);

		}

		if(node_tch.child("Tch_AvgOutput"))
			ReadAverageResultsEjeXML(node_tch);
		if(node_tch.child("Tch_InsOutput"))
			ReadInstantaneousResultsEjeXML(node_tch);

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadTurbochargerAxis in the boundary condition: " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::CalculaEjesTurbogrupo(double Theta, nmTipoModelado SimulationType, double Time,
		double CrankAngle) {
	try {

		// Calculo del nuevo regimen del turbogrupo.

		double MechWork = 0;
		//  HAY QUE PASAR EL VALOR DE TAMB
		double DeltaTime = Time - FTime;
		FTime = Time;

		if(FThereIsHTM) {
			double p1 = __units::BarToPa(FCompresor[0]->AcousticC()->P1());
			double p2 = __units::BarToPa(FCompresor[0]->AcousticC()->P2());
			double p3 = __units::BarToPa(FTurbina[0]->AcousticT()->P3());
			/* The mechanical losses model needs the stator outlet pressure */
			double p_so = __units::BarToPa(FTurbina[0]->getPressure());
			double p4 = __units::BarToPa(FTurbina[0]->AcousticT()->P4());
			FMechPower = FMechLosses->P_oil(__units::degCToK(FToil), __units::RPMToRad_s(FRegimenEje), p1, p2, p_so, p4, FMoil);
			MechWork = FMechPower * DeltaTime;
			if(FSumTrabajoTurbinas > MechWork) {
				FMechEff = 1 - MechWork / FSumTrabajoTurbinas;
			} else {
				FMechEff = 0;
			}
			double T1 = FCompresor[0]->AcousticC()->T1();
			double T2 = FCompresor[0]->AcousticC()->T2d();
			double HC = FHTM->CompressorHeatFlow();

			double cr = FCompresor[0]->AcousticC()->CompRatio();
			double efc = FCompresor[0]->getEfficiency();

			double T3 = FTurbina[0]->AcousticT()->T3();
			double er = FTurbina[0]->AcousticT()->ExpRatio();
			double eft = FTurbina[0]->GetEfficiency();

			FHTM->PutEffTurbMax(FTurbina[0]->GetMaxEff());

			if(FConvTemp.FastConvergence) {

				FConvTemp.T1u += T1 * fabs(FCompresor[0]->AcousticC()->MassFlowIn()) * DeltaTime;
				FConvTemp.T1d += fabs(FCompresor[0]->AcousticC()->MassFlowIn()) * DeltaTime;
				FConvTemp.T2u += T2 * fabs(FCompresor[0]->AcousticC()->MassFlowOut()) * DeltaTime;
				FConvTemp.T2d += fabs(FCompresor[0]->AcousticC()->MassFlowOut()) * DeltaTime;
				FConvTemp.T3u += T3 * fabs(FTurbina[0]->AcousticT()->MassIn()) * DeltaTime;
				FConvTemp.T3d += fabs(FTurbina[0]->AcousticT()->MassIn()) * DeltaTime;
				FConvTemp.MassC += FCompresor[0]->AcousticC()->MassFlowIn() * DeltaTime;
				FConvTemp.MassT += FTurbina[0]->AcousticT()->MassIn() * DeltaTime;
				FConvTemp.MechP += FMechPower * DeltaTime;
				FConvTemp.HeatC1 += HC * DeltaTime;

				FConvTemp.TimeSum += DeltaTime;

				if(Time > FConvTemp.NewTime) {
					double T1c = FConvTemp.T1u / FConvTemp.T1d;
					double T2c = FConvTemp.T2u / FConvTemp.T2d;
					double T3c = FConvTemp.T3u / FConvTemp.T3d;
					double MassC = FConvTemp.MassC / FConvTemp.TimeSum;
					double MassT = FConvTemp.MassT / FConvTemp.TimeSum;
					double HeatC1 = FConvTemp.HeatC1 / FConvTemp.TimeSum;
					double MPc = FConvTemp.MechP / FConvTemp.TimeSum;

					FHTM->setKFromGas(T3c, MassT);

					double T2cavg = T2c;
					if (MassC > 1e-7){
						T2cavg = FHTM->getT2avg(T2c, MassC, HeatC1);
						if (T2cavg < T1c) T2cavg = T1c;
						else if (T2cavg > T3c) T2cavg = T3c;
					}
					FHTM->setKFromAir(T2cavg, MassC);
					FHTM->setKFromOil(__units::degCToK(FToil), FMoil, MPc);
					FHTM->setKFromCoolant(__units::degCToK(FTwater), FMwater);
					FHTM->setKFromAmbient(__Ambient::T_K, __Ambient::p_Pa);
					FHTM->setKRadMatrix();
					FHTM->setKRadiatioFromAmbient(__Ambient::T_K);

					FHTM->SolveImplicit(FConvTemp.TimeSum, 0);

					FConvTemp.T1u = 0;
					FConvTemp.T1d = 0;
					FConvTemp.T2u = 0;
					FConvTemp.T2d = 0;
					FConvTemp.T3u = 0;
					FConvTemp.T3d = 0;
					FConvTemp.MassC = 0;
					FConvTemp.MassT = 0;
					FConvTemp.HeatC1 = 0;
					FConvTemp.MechP = 0;

					FConvTemp.NewTime += FConvTemp.TimeStep;

					FConvTemp.TimeSum = 0;

					if(Time > FConvTemp.MaxTime)
						FConvTemp.FastConvergence = false;
				}
			}

			FHTM->setHeatFromGas(DeltaTime, T3, FTurbina[0]->AcousticT()->MassIn());
			double T2avg = T2;
			if (FCompresor[0]->getMassflow() > 1e-10){
				T2avg = FHTM->getT2avg(T2, FCompresor[0]->getMassflow(), HC);

				if (T2avg < T1) T2avg = T1;
				else if (T2avg > T3) T2avg = T3;
			}

			

			FHTM->setHeatFromAir(DeltaTime, T2avg, FCompresor[0]->getMassflow());
			FHTM->setHeatFromOil(DeltaTime, __units::degCToK(FToil), FMoil, FMechPower);
			FHTM->setHeatFromCoolant(DeltaTime, __units::degCToK(FTwater), FMwater);
			FHTM->setHeatFromAmbient(DeltaTime, __Ambient::T_K, __Ambient::p_Pa);
			FHTM->setKRadMatrix();
			FHTM->setRadiatioFromAmbient(DeltaTime, __Ambient::T_K);

			FHTM->SolveExplicit(DeltaTime, 0);

			//FCompresor[0]->AcousticC()->PutHeatPower(FHTM->Comp_Heat_Flow());
			//FCompresor[0]->AcousticC()->PutHeatPowerIn(FHTM->Comp_Heat_Flow_In());

			FCompresor[0]->AcousticC()->PutHeatPower(FHTM->CompressorHeatFlow());
			FCompresor[0]->AcousticC()->PutHeatPowerIn(0.);
		}

		FSumTrabajoCompresores = 0.;
		FSumTrabajoTurbinas = 0.;

		for(int i = 0; i < FNumCompresoresAcoplados; i++) {
			FCompresor[i]->CalculoPotenciaPaso();
			FSumTrabajoCompresores += (FCompresor[i]->getPotenciaPaso() * DeltaTime);
		}

		for(int i = 0; i < FNumTurbinasAcopladas; i++) {
			FTurbina[i]->CalculoPotenciaPaso();
			FSumTrabajoTurbinas += (FTurbina[i]->getPotenciaPaso() * DeltaTime);
		}

		if(FRPMControlled) {
			FRegimenEje = FController->Output(FTime);
		} else {
			if(FVariacionRegimen == nmVariable) {
				if((Theta > 2880. && SimulationType != nmEstacionario) || (SimulationType == nmEstacionario && Theta > 50.
						&& FNumCilindros == 0)
				   || (SimulationType == nmEstacionario && Theta > 2880 && FNumCilindros != 0)) {

					FDeltaReg = __units::Rad_sToRPM(__units::Rad_sToRPM((FSumTrabajoTurbinas - FSumTrabajoCompresores - MechWork) /
													(FMomentoInercia * FRegimenEje)));
					FRegimenEje += FDeltaReg;

					if(FRegimenEje < 10000)
						FRegimenEje = 10000;
				}
			}
		}

		// Actualizacion del regimen de los compresores acoplados al eje y su interpolacion.
		for(int i = 0; i < FNumCompresoresAcoplados; ++i) {
			if(FCompresor[i]->getModeloCompresor() == nmCompPlenums || FCompresor[i]->getModeloCompresor() == nmCompPipes) {
				FCompresor[i]->InterpolaValoresMapa(FRegimenEje);
			} else if(FCompresor[i]->getModeloCompresor() == nmCompOriginal) {
				if(FVariacionRegimen == nmVariable
				   || (dynamic_cast<TCompTubDep*>(FCompresor[i]))->getEntradaCompresor() != nmAtmosphere) {
					FCompresor[i]->InterpolaValoresMapa(FRegimenEje);
				}
			}
		}
		// Actualizacion del regimen de las turbinas acopladas al eje.
		for(int j = 0; j < FNumTurbinasAcopladas; ++j) {
			FTurbina[j]->PutRegimen(FRegimenEje);
		}

		// Acumulacion de valores medios.

		if(FResMediosEje.Regimen) {
			FResMediosEje.RegimenSUM += FRegimenEje * DeltaTime;
		}
		FResMediosEje.TiempoSUM += DeltaTime;

		// printf("%lf \n",FRegimenEje);

		// Salida de resultados por pantalla.
		// printf("%lf %lf\n",CrankAngle,Theta);

		if(CrankAngle - FAngle0 <= 0. && Theta >= 750.) {
			FNumCiclo++;
			printf("\n");
			printf("*****************************************************\n");
			printf("***TURBOCHARGER AVERAGE VALUES %3d **CYCLE N. %3d ***\n", FNumeroEje, FNumCiclo);
			printf("*****************************************************\n");
			printf("\n");
			for(int j = 0; j < FNumTurbinasAcopladas; ++j) {
				FTurbina[j]->ImprimeResultadosMediosPantalla();
			}
			for(int j = 0; j < FNumCompresoresAcoplados; ++j) {
				FCompresor[j]->CalculaMedias();
				printf("COMPRESSOR WORK %d       = %6.3lf Julios \n", FNumeroCompresor[j], FCompresor[j]->getTrabCiclo());
				printf("COMPRESSOR EFFICIENCY %d = %6.3lf \n", FNumeroCompresor[j], FCompresor[j]->getRendMed());
				printf("COMPRESSOR MASS FLOW %d  = %6.3lf g/s\n", FNumeroCompresor[j], FCompresor[j]->getGastoMed() * 1000);
				printf("COMPRESSOR RATIO %d      = %6.3lf \n", FNumeroCompresor[j], FCompresor[j]->getRCMed());
			}
			printf("TURBOCHARGER SPEED      = %6.3lf r.p.m.\n", FResMediosEje.RegimenSUM / FResMediosEje.TiempoSUM);
			printf("*****************************************************\n\n");
		}
		FAngle0 = CrankAngle;

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::CalculaEjesTurbogrupo in the boundary condition: " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ReadAverageResultsEje(const char* FileWAM, fpos_t & filepos) {
	try {
		int nvars = 0, var = 0;

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		fscanf(fich, "%d ", &nvars);
		for(int i = 0; i < nvars; i++) {
			fscanf(fich, "%d ", &var);
			switch(var) {
			case 0:
				FResMediosEje.Regimen = true;
				break;
			default:
				std::cout << "Resultados medios en Axis " << FNumeroEje << " no implementados " << std::endl;
			}
		}

		fgetpos(fich, &filepos);
		fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadAverageResultsEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TEjeTurbogrupo::ReadAverageResultsEjeXML(xml_node node_shaft) {
	try {
		int nvars, var;

		FResMediosEje.Regimen = false;

		xml_node node_avg = node_shaft.child("Tch_AvgOutput");

		for(xml_attribute parameter = node_avg.attribute("Parameter"); parameter; parameter.next_attribute()) {
			if(parameter.value() == "Speed") {
				FResMediosEje.Regimen = true;
			} else {
				std::cout << "Resultados medios en Axis " << FNumeroEje << " no implementados " << std::endl;
			}
		}

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadAverageResultsEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::CabeceraResultadosMedEje(stringstream & medoutput) {
	try {
		// FILE *fich=fopen(FileSALIDA,"a");
		std::string Label;

		if(FResMediosEje.Regimen) {
			Label = "\t" + PutLabel(5007) + "/" + std::to_string(FNumeroEje) + "/" + PutLabel(4022);
			medoutput << Label.c_str();
		}

		// fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::CabeceraResultadosMedEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ImprimeResultadosMedEje(stringstream & medoutput) {
	try {
		// FILE *fich=fopen(FileSALIDA,"a");

		if(FResMediosEje.Regimen)
			medoutput << "\t" << FResMediosEje.RegimenMED;

		// fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ImprimeResultadosMedEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::IniciaMedias() {
	try {

		FResMediosEje.RegimenSUM = 0.;
		FResMediosEje.TiempoSUM = 0.;
		FResMediosEje.Tiempo0 = 0.;

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::IniciaMedias en el eje: " << FNumeroEje << std::endl;
		// std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ResultadosMediosEje() {
	try {

		if(FResMediosEje.Regimen) {
			FResMediosEje.RegimenMED = FResMediosEje.RegimenSUM / FResMediosEje.TiempoSUM;
			FResMediosEje.RegimenSUM = 0.;
		}
		FResMediosEje.TiempoSUM = 0;

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ResultadosMediosEje en el eje: " << FNumeroEje << std::endl;
		// std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::AcumulaResultadosMediosEje(double Actual) {
	try {
		/* Lo que se hace en esta funcion se realiza dentro del calculo del eje, para asi poder
		 llevar a cabo la salida de resultados medios por pantalla. */
		double Delta = Actual - FResMediosEje.Tiempo0;

		if(FResMediosEje.Regimen) {
			FResMediosEje.RegimenSUM += FRegimenEje * Delta;
		}
		FResMediosEje.TiempoSUM += Delta;
		FResMediosEje.Tiempo0 = Delta;

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::AcumulaResultadosMediosEje en el eje: " << FNumeroEje << std::endl;
		// std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ReadInstantaneousResultsEje(const char* FileWAM, fpos_t & filepos) {
	try {
		int nvars = 0, var = 0;

		FILE *fich = fopen(FileWAM, "r");
		fsetpos(fich, &filepos);

		FResInstantEje.Regimen = false;
		FResInstantEje.MechPower = false;
		FResInstantEje.MechEff = false;
		FResInstantEje.NodeTemp = false;
		FResInstantEje.HeatFlow = false;

		fscanf(fich, "%d ", &nvars);
		for(int i = 0; i < nvars; i++) {
			fscanf(fich, "%d ", &var);
			switch(var) {
			case 0:
				FResInstantEje.Regimen = true;
				break;
			case 1:
				FResInstantEje.MechPower = true;
				break;
			case 2:
				FResInstantEje.MechEff = true;
				break;
			case 3:
				FResInstantEje.NodeTemp = true;
				break;
			case 4:
				FResInstantEje.HeatFlow = true;
				break;
			default:
				std::cout << "Instantaneous results in axis " << FNumeroEje << " are not implemented " << std::endl;
			}
		}

		fgetpos(fich, &filepos);
		fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadInstantaneousResultsEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

void TEjeTurbogrupo::ReadInstantaneousResultsEjeXML(xml_node node_shaft) {
	try {

		FResInstantEje.Regimen = false;
		FResInstantEje.MechPower = false;
		FResInstantEje.MechEff = false;
		FResInstantEje.NodeTemp = false;
		FResInstantEje.HeatFlow = false;

		xml_node node_ins = node_shaft.child("Tch_InsOutput");

		for(xml_attribute parameter = node_ins.attribute("Parameter"); parameter; parameter.next_attribute()) {
			if(parameter.value() == "Speed") {
				FResInstantEje.Regimen = true;
			} else if(parameter.value() == "MechanicalLosses") {
				FResInstantEje.MechPower = true;
			} else if(parameter.value() == "MechanicalEfficiency") {
				FResInstantEje.MechEff = true;
			} else if(parameter.value() == "NodeTemperature") {
				FResInstantEje.NodeTemp = true;
			} else if(parameter.value() == "HeatFlow") {
				FResInstantEje.HeatFlow = true;
			} else {
				std::cout << "Instantaneous results in axis " << FNumeroEje << " are not implemented " << std::endl;
			}
		}

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ReadInstantaneousResultsEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::HeaderInstantaneousResultsEje(stringstream & insoutput) {
	try {
		// FILE *fich=fopen(FileSALIDA,"a");

		std::string Label;

		if(FResInstantEje.Regimen) {
			Label = "\t" + PutLabel(5007) + "/" + std::to_string(FNumeroEje) + "/" + PutLabel(4022) + PutLabel(918);
			insoutput << Label.c_str();
		}
		if(FResInstantEje.MechPower) {
			Label = "\t" + PutLabel(5007) + "/" + std::to_string(FNumeroEje) + "/" + PutLabel(4029) + PutLabel(4009) + PutLabel(
						903);
			insoutput << Label.c_str();
		}
		if(FResInstantEje.MechEff) {
			Label = "\t" + PutLabel(5007) + "/" + std::to_string(FNumeroEje) + "/" + PutLabel(4029) + PutLabel(4011) + PutLabel(
						901);
			insoutput << Label.c_str();
		}
		if(FThereIsHTM) {
			if(FResInstantEje.NodeTemp) {
				FHTM->HeaderInsTemperatures(insoutput, FNumeroEje);
			}
			if(FResInstantEje.HeatFlow) {
				FHTM->HeaderInsHeatFlow(insoutput, FNumeroEje);
			}
		}

		// fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::HeaderInstantaneousResultsEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ImprimeResultadosInstantaneosEje(stringstream & insoutput) {
	try {
		// FILE *fich=fopen(FileSALIDA,"a");

		if(FResInstantEje.Regimen)
			insoutput << "\t" << FResInstantEje.RegimenINS;
		if(FResInstantEje.MechPower)
			insoutput << "\t" << FResInstantEje.MechPowerINS;
		if(FResInstantEje.MechEff)
			insoutput << "\t" << FResInstantEje.MechEffINS;
		if(FThereIsHTM) {
			if(FResInstantEje.NodeTemp) {
				FHTM->PrintInsTemperatures(insoutput);
			}
			if(FResInstantEje.HeatFlow) {
				FHTM->PrintInsHeatFlow(insoutput);
			}
		}

		// fclose(fich);
	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ImprimeResultadosInstantaneosEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::ResultadosInstantEje() {
	try {
		if(FResInstantEje.Regimen)
			FResInstantEje.RegimenINS = FRegimenEje;
		if(FResInstantEje.MechPower)
			FResInstantEje.MechPowerINS = FMechPower;
		if(FResInstantEje.MechEff)
			FResInstantEje.MechEffINS = FMechEff;

	} catch(exception & N) {
		std::cout << "ERROR: TEjeTurbogrupo::ResultadosInstantEje en el eje " << FNumeroEje << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int TEjeTurbogrupo::GetNumeroCompresor(int i) {

	return FNumeroCompresor[i];

}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

void TEjeTurbogrupo::AsignaRPMController(TController * *Controller) {
	if(FRPMControlled)
		FController = Controller[FControllerID - 1];
}

void TEjeTurbogrupo::InitializeTurbocharger(double Tamb) {

	if(FThereIsHTM) {

		//FHTM->AsignTCMechLosses(FMechLosses);

		FHTM->AsignTCMechLosses(FMechLosses);

		FHTM->setOilConditions(FMoil, __units::degCToK(FToil));
		FHTM->setCoolantConditions(FMwater, __units::degCToK(FTwater));

		FHTM->setDiameters("AIR", FCompresor[0]->AcousticC()->Din(), FCompresor[0]->AcousticC()->Dout());
		FHTM->setDiameters("GAS", FTurbina[0]->AcousticT()->DInTot(), FTurbina[0]->AcousticT()->DOut());


		//FHTM->TurbochargerWorkingPoint(FRegimenEje, 0.9, FMoil, __units::degCToK(FToil), FPoil, __units::degCToK(FTwater),
		//							   FMwater);
		//FHTM->CompressorData(FCompresor[0]->GetMap()->getPresionRef(), FCompresor[0]->GetMap()->getTempRef(),
		//					 FCompresor[0]->GetMap()->getTempMeasure(), FCompresor[0]->AcousticC()->Din(),
		//					 FCompresor[0]->AcousticC()->Din());
		//FHTM->TurbineData(1, 300, FTurbina[0]->getMap()->getTempMeasure(), FTurbina[0]->AcousticT()->DInTot(),
		//				  FTurbina[0]->AcousticT()->DOut());

		//FTurbina[0]->AsignTCHTM(FHTM);

		FTurbina[0]->AsignTCHTM(FHTM);

		//FCompresor[0]->AsignTCHTM(FHTM);

		FCompresor[0]->AsignTCHTM(FHTM);

		//  TEMPERATURA AMBIENTE

		double CMT = FCompresor[0]->GetMap()->getTempMeasure();
		double TMT = FTurbina[0]->getMap()->getTempMeasure();

		FTurbina[0]->PreprocessMap(CMT);

		FCompresor[0]->PreprocessMap(TMT, sqrt(FCWArea * __cons::_4_Pi));

		FCompresor[0]->InterpolaValoresMapa(FRegimenEje);

		FHTM->InitializeTemp(FTurbina[0]->AcousticT()->T3(), FCompresor[0]->AcousticC()->T2(),
							 __units::degCToK(FToil), __units::degCToK(FTwater));

	}

}

void TEjeTurbogrupo::Put_Conditions(int ID, double val) {
	switch(ID) {
	case 2:
		FToil = val;
		break;
	case 3:
		if(val == 0.) {
			cout << "ERROR: Mass flow is equal to 0. Have you added the corresponding input to WiringCMT?" << endl;
			FMoil = 1e-12;
		} else {
			FMoil = val;
		}

		break;
	case 4:
		FTwater = val;
		break;
	case 5:
		FMwater = val;
		break;
	case 6:
		FTamb = val;
		break;
	case 7:
		FHTM->Put_ExternalVelocity(val);
		break;
	case 8:
		FHTM->Put_ExternalHeat(val);
		break;
	case 9:
		FRegimenEje = val;
		break;
	case 10:
		if(val < 0) {
			FVariacionRegimen = nmFixed;
		} else {
			FVariacionRegimen = nmVariable;
		}
		break;
	}
}

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

#pragma package(smart_init)
