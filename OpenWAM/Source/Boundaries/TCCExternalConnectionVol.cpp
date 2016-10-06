// ---------------------------------------------------------------------------

#pragma hdrstop

#include "TCCExternalConnectionVol.h"
#include "TTubo.h"

// ---------------------------------------------------------------------------


TCCExternalConnectionVol::TCCExternalConnectionVol() {
	FTuboExtremo = NULL;
	FTime0 = 0;

	FTimeSum = 0;

	FT_BoundarySum = 0;
	FU_BoundarySum = 0;
	FP_BoundarySum = 0;
}

TCCExternalConnectionVol::TCCExternalConnectionVol(nmTypeBC TipoCC, int numCC, nmTipoCalculoEspecies SpeciesModel,
		int numeroespecies, nmCalculoGamma GammaCalculation, bool ThereIsEGR) :
	TCondicionContorno(TipoCC, numCC, SpeciesModel, numeroespecies, GammaCalculation, ThereIsEGR) {

	FTuboExtremo = NULL;
	FTime0 = 0;

	FTimeSum = 0;

	FT_BoundarySum = 0;
	FU_BoundarySum = 0;
	FP_BoundarySum = 0;

	// FComposicion = NULL;

}

TCCExternalConnectionVol::~TCCExternalConnectionVol() {
}

void TCCExternalConnectionVol::UpdateCurrentExternalProperties(double U0, double T0, double P0, double t) {

	FUExt = U0;
	FTExt = T0;
	FPExt = P0;

	FCurrentTime = t;
}

void TCCExternalConnectionVol::AsignGeometricalData(double A) {

	FAExt = A;
}

void TCCExternalConnectionVol::ExternalCharacteristics(double Time) {

	// Calculo Entropia
	double x = 0., Px = 0., Tx = 0., Ux = 0., Ax = 0.;
	// double Deltat = Time - FTime0;

	double g = __Gamma::G;
	double gg = __Gamma::G_5;

	double A0 = Sqrt(g * __R::Air * FTExt);

	FA_AExt = A0 / pow(__units::PaToBar(FPExt), gg) / __cons::ARef;

	FK_CExt = (A0 + (g - 1) / 2 * FUExt) / __cons::ARef;

	FA_Dep = sqrt(pow2(A0) + (g - 1) / 2 * pow2(FUExt)) / __cons::ARef;

}

void TCCExternalConnectionVol::CalculaCondicionContorno(double Time) {

	double flujo = 0.;

	ExternalCharacteristics(Time);

	flujo = (FK_CExt / FA_AExt) / (*FCC / FTuboExtremo[0].Entropia);

	if(flujo < 0.999995) {

		double A = FTuboExtremo[0].Entropia / (FTuboExtremo[0].Entropia + FA_AExt) * (*FCC + FK_CExt);
		FT_Boundary = pow2(A * __cons::ARef) / __R::Air / __Gamma::G;
		FU_Boundary = (A - *FCC) / __Gamma::G_3;
		*FCD = A + __Gamma::G_3 * FU_Boundary;
		FU_Boundary *= __cons::ARef;
		FP_Boundary = pow(A / FTuboExtremo[0].Entropia, __Gamma::G_4);
	} else if(flujo > 1.000005) {
		//double A = FA_AExt / (FTuboExtremo[0].Entropia + FA_AExt) * (*FCC + FK_CExt);
		double xyx = FA_AExt / FTuboExtremo[0].Entropia;
		double a = pow2(__Gamma::G_3 * xyx) + __Gamma::G_3;
		double b = __Gamma::G_1 * *FCC * pow2(xyx);
		double c = pow2(xyx * *FCC) - pow2(FA_Dep);
		FU_Boundary = (-b + sqrt(pow(b, 2.) - a * 4. * c)) / (2. * a);
		FTuboExtremo[0].Entropia = FA_AExt;
		double A = sqrt(pow(FA_Dep, 2.) - __Gamma::G_3 * pow(FU_Boundary, 2.));
		*FCD = A + __Gamma::G_3 * FU_Boundary;
		*FCC = A - __Gamma::G_3 * FU_Boundary;
		FT_Boundary = pow2(A * __cons::ARef) / __R::Air / __Gamma::G;
		FU_Boundary *= __cons::ARef;
		FP_Boundary = pow(A / FTuboExtremo[0].Entropia, __Gamma::G_4);

	} else {
		double A = *FCC;
		FT_Boundary = pow2(A * __cons::ARef) / __R::Air / __Gamma::G;
		FU_Boundary = 0;
		*FCD = *FCC;
		FP_Boundary = pow(A / FTuboExtremo[0].Entropia, __Gamma::G_4);
	}
	double deltat = Time - FTime0;

	FT_BoundarySum = FT_Boundary * deltat;
	FU_BoundarySum = FU_Boundary * deltat;
	FP_BoundarySum = FP_Boundary * deltat;

	FTimeSum = deltat;

	FTime0 = Time;

}

void TCCExternalConnectionVol::ReadBoundaryData(const char *FileWAM, fpos_t &filepos, int NumberOfPipes, TTubo **Pipe,
		int nDPF, TDPF **DPF) {

	int i = 0;

	FTuboExtremo = new stTuboExtremo[1];
	FTuboExtremo[0].Pipe = NULL;

	while(FNumeroTubosCC < 1 && i < NumberOfPipes) {
		if(Pipe[i]->getNodoIzq() == FNumeroCC) {
			FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
			FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmLeft;
			FCC = &(FTuboExtremo[FNumeroTubosCC].Beta);
			FCD = &(FTuboExtremo[FNumeroTubosCC].Landa);
			FNodoFin = 0;
			FIndiceCC = 0;
			FNumeroTubosCC++;
		}
		if(Pipe[i]->getNodoDer() == FNumeroCC) {
			FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
			FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmRight;
			FCC = &(FTuboExtremo[FNumeroTubosCC].Landa);
			FCD = &(FTuboExtremo[FNumeroTubosCC].Beta);
			FNodoFin = FTuboExtremo[FNumeroTubosCC].Pipe->getNin() - 1;
			FIndiceCC = 1;
			FNumeroTubosCC++;
		}
		i++;
	}

	// Inicializacion del transporte de especies quimicas.
	FFraccionMasicaEspecie = new double[FNumeroEspecies - FIntEGR];
	for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
		FFraccionMasicaEspecie[i] = FTuboExtremo[0].Pipe->GetFraccionMasicaInicial(i);
	}

	FILE *fich = fopen(FileWAM, "r");
	fsetpos(fich, &filepos);

	fscanf(fich, "%d ", &FID);

	fgetpos(fich, &filepos);
	fclose(fich);

}

void TCCExternalConnectionVol::ReadBoundaryDataXML(xml_node node_connect, int NumberOfPipes, TTubo **Pipe, int nDPF,
		TDPF **DPF) {

	int i = 0;

	FTuboExtremo = new stTuboExtremo[1];
	FTuboExtremo[0].Pipe = NULL;

	while(FNumeroTubosCC < 1 && i < NumberOfPipes) {
		if(Pipe[i]->getNodoIzq() == FNumeroCC) {
			FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
			FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmLeft;
			FCC = &(FTuboExtremo[FNumeroTubosCC].Beta);
			FCD = &(FTuboExtremo[FNumeroTubosCC].Landa);
			FNodoFin = 0;
			FIndiceCC = 0;
			FNumeroTubosCC++;
		}
		if(Pipe[i]->getNodoDer() == FNumeroCC) {
			FTuboExtremo[FNumeroTubosCC].Pipe = Pipe[i];
			FTuboExtremo[FNumeroTubosCC].TipoExtremo = nmRight;
			FCC = &(FTuboExtremo[FNumeroTubosCC].Landa);
			FCD = &(FTuboExtremo[FNumeroTubosCC].Beta);
			FNodoFin = FTuboExtremo[FNumeroTubosCC].Pipe->getNin() - 1;
			FIndiceCC = 1;
			FNumeroTubosCC++;
		}
		i++;
	}

	// Inicializacion del transporte de especies quimicas.
	FFraccionMasicaEspecie = new double[FNumeroEspecies - FIntEGR];
	for(int i = 0; i < FNumeroEspecies - FIntEGR; i++) {
		FFraccionMasicaEspecie[i] = FTuboExtremo[0].Pipe->GetFraccionMasicaInicial(i);
	}

	xml_node node_ext = GetNodeChild(node_connect, "Con_ExternalConnection");
	FID = GetAttributeAsDouble(node_ext, "Connection_ID");

}

void TCCExternalConnectionVol::LoadNewData(double* p, double* T, double* u) {

	if(FTimeSum > 0) {
		*p = FP_BoundarySum / FTimeSum;
		*T = FT_BoundarySum / FTimeSum;
		*u = FU_BoundarySum / FTimeSum;
	} else {
		*p = FP_Boundary;
		*T = FT_Boundary;
		*u = FU_Boundary;
	}

	FT_BoundarySum = 0;
	FU_BoundarySum = 0;
	FP_BoundarySum = 0;

	FTimeSum = 0;
}

#pragma package(smart_init)
