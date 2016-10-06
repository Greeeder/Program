/* --------------------------------------------------------------------------------*\
 |==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
 |----------------------------------------------------------------------------------
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

#include "TExternalController.h"

TExternalController::TExternalController(int i) :
	TController(nmExtCtrl, i) {

	fID = i + 1;
}

TExternalController::~TExternalController() {
}

void TExternalController::AssignExternalValue(double Value) {

	fOutput = Value;

}

double TExternalController::Output(double Time) {

	return fOutput;
}

void TExternalController::LeeController(const char *FileWAM, fpos_t &filepos) {

	FILE *fich = fopen(FileWAM, "r");
	fsetpos(fich, &filepos);

	fscanf(fich, "%lf ", &fOutput);

	fgetpos(fich, &filepos);
	fclose(fich);

}

void TExternalController::LeeControllerXML(xml_node node_ctrl) {

	xml_node node_ec = GetNodeChild(node_ctrl, "Ctr_ExtCtr");

	fOutput = GetAttributeAsDouble(node_ec, "Output");
}

void TExternalController::AsignaObjetos(TSensor **Sensor, TController **Controller) {

}

void TExternalController::LeeResultadosMedControlador(const char *FileWAM, fpos_t &filepos) {

	int nvars = 0, var = 0;

	FILE *fich = fopen(FileWAM, "r");
	fsetpos(fich, &filepos);

	fscanf(fich, "%d ", &nvars);
	for(int i = 0; i < nvars; i++) {
		fscanf(fich, "%d ", &var);
		switch(var) {
		case 0:
			FResMediosCtrl.Output = true;
			break;
		default:
			std::cout << "Resultados medios en Controlador " << fID << " no implementados " << std::endl;
		}
	}

	fgetpos(fich, &filepos);
	fclose(fich);

}

void TExternalController::LeeResultadosMedControladorXML(xml_node node_ctrl) {

	xml_node node_avg = GetNodeChild(node_ctrl, "Ctrl:AvgOutput");
	for(xml_attribute parameter = node_avg.attribute("Parameter"); parameter; parameter.next_attribute()) {
		if(parameter.value() == "Output") {
			FResMediosCtrl.Output = true;
		} else {
			std::cout << "Resultados medios en Controlador " << fID << " no implementados " << std::endl;
		}
	}
}

void TExternalController::LeeResultadosInsControlador(const char *FileWAM, fpos_t &filepos) {

	int nvars = 0, var = 0;

	FILE *fich = fopen(FileWAM, "r");
	fsetpos(fich, &filepos);

	fscanf(fich, "%d ", &nvars);
	for(int i = 0; i < nvars; i++) {
		fscanf(fich, "%d ", &var);
		switch(var) {
		case 0:
			FResInstantCtrl.Output = true;
			break;
		default:
			std::cout << "Resultados instantaneos en Controlador " << fID << " no implementados " << std::endl;
		}
	}

	fgetpos(fich, &filepos);
	fclose(fich);

}

void TExternalController::LeeResultadosInsControladorXML(xml_node node_ctrl) {

	xml_node node_ins = GetNodeChild(node_ctrl, "Ctrl:InsOutput");
	for(xml_attribute parameter = node_ins.attribute("Parameter"); parameter; parameter.next_attribute()) {
		if(parameter.value() == "Output") {
			FResInstantCtrl.Output = true;
		} else {
			std::cout << "Resultados instantaneos en Controlador " << fID << " no implementados " << std::endl;
		}
	}
}

void TExternalController::CabeceraResultadosMedControlador(stringstream& medoutput) {

	std::string Label;

	if(FResMediosCtrl.Output) {
		Label = "\t" + PutLabel(705) + std::to_string(fID) + PutLabel(901);
		medoutput << Label.c_str();
	}

}

void TExternalController::CabeceraResultadosInsControlador(stringstream& insoutput) {

	std::string Label;

	if(FResInstantCtrl.Output) {
		Label = "\t" + PutLabel(705) + std::to_string(fID) + PutLabel(901);
		insoutput << Label.c_str();
	}

}

void TExternalController::ImprimeResultadosMedControlador(stringstream& medoutput) {

	std::string Label;

	if(FResMediosCtrl.Output) {
		medoutput << "\t" << FResMediosCtrl.OutputMED;
	}

}

void TExternalController::ImprimeResultadosInsControlador(stringstream& insoutput) {

	std::string Label;

	if(FResInstantCtrl.Output) {
		insoutput << "\t" << FResInstantCtrl.OutputINS;
	}

}

void TExternalController::IniciaMedias() {

	FResMediosCtrl.OutputSUM = 0.;
	FResMediosCtrl.TiempoSUM = 0.;
	FResMediosCtrl.Tiempo0 = 0.;

}

void TExternalController::ResultadosMediosController() {

	if(FResMediosCtrl.Output) {
		FResMediosCtrl.OutputMED = FResMediosCtrl.OutputSUM / FResMediosCtrl.TiempoSUM;
		FResMediosCtrl.OutputSUM = 0.;
	}

	FResMediosCtrl.TiempoSUM = 0;

}

void TExternalController::AcumulaResultadosMediosController(double Actual) {

	/* Lo que se hace en esta funci�n se realiza dentro del calculo del eje, para as� poder
	 llevar a cabo la salida de resultados medios por pantalla. */
	double Delta = Actual - FResMediosCtrl.Tiempo0;

	if(FResMediosCtrl.Output) {
		FResMediosCtrl.OutputSUM += fOutput * Delta;
	}

	FResMediosCtrl.TiempoSUM += Delta;
	FResMediosCtrl.Tiempo0 = Actual;

}

void TExternalController::ResultadosInstantController() {

	if(FResInstantCtrl.Output)
		FResInstantCtrl.OutputINS = fOutput;

}

// ---------------------------------------------------------------------------

#pragma package(smart_init)
