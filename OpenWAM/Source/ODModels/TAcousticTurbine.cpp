// ---------------------------------------------------------------------------

#pragma hdrstop

#include "TAcousticTurbine.h"

// ---------------------------------------------------------------------------

TAcousticTurbine::TAcousticTurbine(iVector InletPipeID, iVector VoluteID, int OutletPipeID) {

	FInletPipeID.resize(InletPipeID.size());
	FVoluteID.resize(InletPipeID.size());
	FInletPipe.resize(InletPipeID.size());
	FVolute.resize(InletPipeID.size());

	for(int i = 0; i < InletPipeID.size(); i++) {
		FInletPipeID[i] = InletPipeID[i];
		FVoluteID[i] = VoluteID[i];
	}
	FOutletPipeID = OutletPipeID;
}

TAcousticTurbine::TAcousticTurbine() {
}

TAcousticTurbine::~TAcousticTurbine() {
}

double TAcousticTurbine::T3(int i) const {

	double a = FInletPipe[i]->GetAsonido(0) * __cons::ARef;
	double g = FInletPipe[i]->GetGamma(0);
	double R = FInletPipe[i]->GetRMezcla(0);

	return a * a / g / R;
}

double TAcousticTurbine::T3() const {

	double T3Sum = 0;
	double MasSum = 0;
	double ret_val = 0;
	for(int i = 0; i < FInletPipe.size(); i++) {
		MasSum += MassIn(i);
	}
	if(MasSum == 0) {
		for(int i = 0; i < FInletPipe.size(); i++) {
			T3Sum += T3(i);
		}
		T3Sum /= FInletPipe.size();
	} else {
		for(int i = 0; i < FInletPipe.size(); i++) {
			T3Sum += T3(i) * MassIn(i);
		}
		T3Sum /= MasSum;
	}
	return T3Sum;
}

double TAcousticTurbine::P3(int i) const {

	return FInletPipe[i]->GetPresion(0);
}

double TAcousticTurbine::P3() const {

	double P3Sum = P3(0);
	for(int i = 1; i < FInletPipe.size(); i++) {
		P3Sum += P3(i);
	}
	P3Sum /= FInletPipe.size();

	return P3Sum;
}

double TAcousticTurbine::P30(int i) const {

	double p = FInletPipe[i]->GetPresion(0);
	double a = FInletPipe[i]->GetAsonido(0) * __cons::ARef;
	double v = FInletPipe[i]->GetVelocidad(0) * __cons::ARef;
	double g = FInletPipe[i]->GetGamma(0);

	double p3 = p * pow(1 + (g - 1) / 2 * pow2(v / a), g / (g - 1));

	return p3;
}

double TAcousticTurbine::T30(int i) const {
	double p = FInletPipe[i]->GetPresion(0);
	double a = FInletPipe[i]->GetAsonido(0) * __cons::ARef;
	double v = FInletPipe[i]->GetVelocidad(0) * __cons::ARef;
	double g = FInletPipe[i]->GetGamma(0);
	double R = FInletPipe[i]->GetRMezcla(0);

	double T0 = (pow2(a) + (g - 1) / 2 * pow2(v)) / g / R;

	return T0;
}

double TAcousticTurbine::DiabEfficiency(int i) const {
	double g = FInletPipe[i]->GetGamma(0);
	double eff = (T30(i) - T4()) / (T30(i) * (1 - pow(1 / ExpRatio(i), (g - 1) / g)));

	return eff;
}

double TAcousticTurbine::ExpRatio(int i) const {

	return P30(i) / P4();

}

double TAcousticTurbine::ExpRatio() const {

	double ERSum = ExpRatio(0);
	for(int i = 1; i < FInletPipe.size(); i++) {
		ERSum += ExpRatio(i);
	}
	ERSum /= FInletPipe.size();

	return ERSum;

}

double TAcousticTurbine::P4() const {

	int n = FOutletPipe->getNin() - 1;

	return FOutletPipe->GetPresion(n);

}

double TAcousticTurbine::T4() const {

	int n = FOutletPipe->getNin() - 1;

	return pow2(FOutletPipe->GetAsonido(n) * __cons::ARef) / FOutletPipe->GetGamma(n) / FOutletPipe->GetRMezcla(0);

}

double TAcousticTurbine::MassIn(int i) const {

	return FInletPipe[i]->GetDensidad(0) * FInletPipe[i]->GetVelocidad(0) * __cons::ARef * SIn(i);
}

double TAcousticTurbine::MassIn() const {

	double MasSum = MassIn(0);
	for(int i = 1; i < FInletPipe.size(); i++) {
		MasSum += MassIn(i);
	}

	return MasSum;
}

double TAcousticTurbine::MassOut() const {

	int n = FOutletPipe->getNin() - 1;

	return FOutletPipe->GetDensidad(n) * FOutletPipe->GetVelocidad(n) * __cons::ARef * SOut();
}

double TAcousticTurbine::DInTot() const {

	double DInSum = 0.0;
	for(int i = 0; i < FInletPipe.size(); i++) {
		DInSum += pow2(DIn(i));
	}

	return Sqrt(DInSum);
}

RowVector TAcousticTurbine::VoluteOutletp() const {

	RowVector p(FInletPipe.size());
	for(int i = 0; i < FInletPipe.size(); i++) {
		p(i) = __units::BarToPa(FInletPipe[i]->GetPresion(FInletPipe[i]->getNin() - 1));
	}

	return p;
}


#pragma package(smart_init)
