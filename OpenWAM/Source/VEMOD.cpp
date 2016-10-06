// ---------------------------------------------------------------------------
//#define VERBOSE

#include <iostream>
#include "VEMOD.h"

#include "TOpenWAM.h"

static std::vector<TOpenWAM*> Library;

namespace __Ambient{
	double p_Pa = 100000;					//!< The ambient pressure [Pa]
	double T_K = 300;						//!< The ambient temperature [K]
	double p_bar = 1;						//!< The ambient pressure [Pa]
	double T_degC = 27;						//!< The ambient temperature [K]
	double HR = 50;							//!< The humidity [%]
	RowVector Y_amb = RowVector::Ones(1);	//!< The ambient mass fraction [-]
}

Uint verbosity = 0;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

EXTERNC void LoadModel(const int icase, char* FileName){

	if (Library[icase - 1] == NULL) {
		Library[icase - 1] = new TOpenWAM();
	} else {
		Library[icase - 1] = NULL;
		Library[icase - 1] = new TOpenWAM();
	}

	Library[icase - 1]->ReadInputDataXML(FileName);

	Library[icase - 1]->InitialParameterTCM();
	
	Library[icase - 1]->ConnectFlowElements();
	
	Library[icase - 1]->ConnectControlElements();
	
	Library[icase - 1]->InitializeParameters();
	
	Library[icase - 1]->InitializeOutput();

}
// (JS) prueba commit 22-07
EXTERNC void RunModel(const int icase, const double t){
	do {
	
		Library[icase - 1]->DetermineTimeStep(t);
	
		Library[icase - 1]->NewEngineCycle();
	
		Library[icase - 1]->CalculateFlowCommon();
	
	} while(!Library[icase - 1]->GetIs_EndStep());
}

EXTERNC double getSensor(const int icase, const int sensor_id){
	return 0.;
}

EXTERNC void setActuator(const int icase, const int actuator_id, const double value){

}





//void __stdcall INITCMODELF(int icase, int itch) {
//	if (Library[icase - 1][itch - 1] == NULL) {
//		Library[icase - 1][itch - 1] = new TOpenWAM();
//	} else {
//		Library[icase - 1][itch - 1] = NULL;
//		Library[icase - 1][itch - 1] = new TOpenWAM();
//	}
//}
//
//void __stdcall LOADCMODELF(int icase, int itch, char *FileName) {
//
//	char file[256];
//	int i = 1;
//
//	Library[icase - 1][itch - 1]->ReadInputDataXML(FileName);
//}
//
//void __stdcall ALLOCATEMODELF(int icase, int itch) {
//
//	Library[icase - 1][itch - 1]->InitialParameterTCM();
//
//	Library[icase - 1][itch - 1]->ConnectFlowElements();
//
//	Library[icase - 1][itch - 1]->ConnectControlElements();
//
//	Library[icase - 1][itch - 1]->InitializeParameters();
//
//	Library[icase - 1][itch - 1]->InitializeOutput();
//
//}
//
//// void __stdcall UPDATEBOUNDARY(int i, double U0, double U1, double T0, double T1, double P0,
//// double P1, double t) {
////
//// Library->UpdateExternalBoundary(i, U0, U1, T0, T1, P0, P1, t);
//// }
////
//// void __stdcall INITIATEBOUNDARY(int i, double D0, double D1, double dX) {
////
//// Library->InitiateExternalBoundary(i, D0, D1, dX);
//// }
//
//void __stdcall UPDATEBOUNDARY(int icase, int itch, int i, double U0, double T0, double P0, double t) {
//
//	Library[icase - 1][itch - 1]->UpdateExternalBoundary(i, U0, T0, P0, t);
//}
//
//void __stdcall UPDATECONTROLLER(int icase, int itch, int i, double Val, double t) {
//
//	Library[icase - 1][itch - 1]->UpdateControllers(i, Val, t);
//
//}
//
//void __stdcall INITIATEBOUNDARY(int icase, int itch, int i, double A) {
//
//	Library[icase - 1][itch - 1]->InitiateExternalBoundary(i, A);
//}
//
//void __stdcall RUNSTEP(int icase, int itch, double t) {
//
//	do {
//
//		Library[icase - 1][itch - 1]->DetermineTimeStep(t);
//
//		Library[icase - 1][itch - 1]->NewEngineCycle();
//
//		Library[icase - 1][itch - 1]->CalculateFlowCommon();
//
//	} while(!Library[icase - 1][itch - 1]->GetIs_EndStep());
//}
//
//void __stdcall LOADNEWDATA(int icase, int itch, int i, double* p, double* T, double* u) {
//
//	Library[icase - 1][itch - 1]->LoadNewData(i, p, T, u);
//
//#ifdef VERBOSE
//	cout << i << "\t" << *p << "\t" << *T << "\t" << *u << endl;
//	if(isnan(*p) || isnan(*T) || isnan(*u)) {
//		cerr << "ERROR: model error" << endl;
//	}
//#endif
//
//}
//
//void __stdcall CLOSEMODEL(int icase, int itch) {
//
//	Library[icase - 1][itch - 1]->GeneralOutput();
//
//	delete Library[icase - 1][itch - 1];
//}
//
//void __stdcall GETPLOT(int icase, int itch, int ID, double* val) {
//
//	*val = Library[icase - 1][itch - 1]->Get_Output(ID);
//
//#ifdef VERBOSE
//	cout << ID << "\t" << *val << endl;
//	if(isnan(*val)) {
//		cerr << "ERROR: model error" << endl;
//	}
//#endif
//
//}
//
//void __stdcall UPDATEHTMDATA(int icase, int itch, int TC, int ID, double val) {
//
//	Library[icase - 1][itch - 1]->UpdateHTMData(TC, ID, val);
//
//}

