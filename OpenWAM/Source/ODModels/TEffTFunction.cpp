// ---------------------------------------------------------------------------

#pragma hdrstop

#include "TEffTFunction.h"
// ---------------------------------------------------------------------------

TEffTFunction::TEffTFunction() {

	z.resize(3, 1);

}

TEffTFunction::TEffTFunction(double _K1, double _A0, double _beta, double _K3, double _g, double _Dwheel) :
	K1(_K1), A0(_A0), beta(_beta), K3(_K3), g(_g), Dwheel(_Dwheel) {

	z.resize(3, 1);

}

TEffTFunction::~TEffTFunction() {
}

double TEffTFunction::val(double blade, double Aeff, double ER) {

	if(blade <= 0.)
		return 1e-6;
	else {
		K3 = pow2(blade) * (1 - pow(1 / ER, (g - 1) / g));
		return -K1 * pow2(blade) + ((2 * (Aeff / A0) * (z[2] * sin(-z[0] * blade + z[1]) + (sqrt(K1 / 2) * tan(beta)))) * (pow(
										1 - (K3 / pow2(blade)), 1 / (g - 1))) * blade);
	}
}

void TEffTFunction::fefft(dVector & x, dVector & a, double&y, dVector & dyda) {

	double C1 = a[1] - x[0] * a[0];
	double C2 = 1 / (g - 1);
	double C3 = pow(1 - K3 / pow2(x[0]), C2);
	double C4 = 2 * x[1] / A0;

	y = -K1 * pow2(x[0]) + ((C4 * (a[2] * (sin(C1) + sqrt(K1 / 2) * tan(beta)))) * C3 * x[0]);

	dyda[0] = -C4 * pow2(x[0]) * a[2] * cos(C1) * C3;

	dyda[1] = C4 * x[0] * a[2] * cos(C1) * C3;

	dyda[2] = C4 * x[0] * (sin(C1) + tan(beta) * sqrt(K1 / 2));

}

void TEffTFunction::fit(dVector Speed, dVector ER, dVector EFF, dVector MR, TAeffFunction * AeffFun) {

	x.resize(Speed.size());
	y.resize(Speed.size());
	sig.resize(Speed.size());

	double aux1 = sqrt(g * __R::Air);
	double aux2 = 1 / g;
	double aux3 = g - 1;
	double aux4 = 2 / aux3;
	double aux5 = aux3 / g;
	double aux6 = sqrt(aux4 * g * __R::Air);

	for(unsigned int i = 0; i < Speed.size(); i++) {
		x[i].resize(2);

		x[i][0] = Speed[i] * __cons::Pi * Dwheel / (aux6 * sqrt(1 - pow(1 / ER[i], aux5)));
		x[i][1] = 1e-6 * MR[i] * aux1 / (g * pow(1 / ER[i], aux2) * sqrt(aux4 * (1 - pow(1 / ER[i], aux5))));
		//x[i][1] = AeffFun->val(x[i][0], 2 / (ER[i] + 1), EFF[i]);
		y[i] = EFF[i];
		sig[i] = 1;

	}

	CorrEffT = new stCorrEffT(K1, A0, beta, K3, g);

	z[0] = 1;
	z[1] = 1.;
	z[2] = 10;

	LMM = new FitmrqMult(x, y, sig, z, CorrEffT, 1e-6);

	for(int i = 0; i < 2; i++) {
		LMM->max_limit(i, 1.5);
		LMM->min_limit(i, 0.5);
	}
	//LMM->hold(0, 1.);
	//LMM->hold(1, 1.);

	LMM->max_limit(2, 30.);
	LMM->min_limit(2, 0.);

	LMM->fit();

	z[0] = LMM->a[0];
	z[1] = LMM->a[1];
	z[2] = LMM->a[2];

	FILE *fich = fopen("TurbineExtrp.txt", "a");
	fprintf(fich, "Chi square: %g\n", LMM->chisq);
	fprintf(fich, "BSR\EffectiveSection\tRealEfficiency\tFittedEfficiency\n");
	for (unsigned int i = 0; i < Speed.size(); i++) {
		fprintf(fich, "%g\t", x[i][0]);
		fprintf(fich, "%g\t", x[i][1]);
		fprintf(fich, "%g\t", EFF[i]);
		fprintf(fich, "%g\n", val(x[i][0], x[i][1], ER[i]));
	}
	fclose(fich);

}

void TEffTFunction::FixParameter(int i, double val) {

	//LMM->free(0);
	//LMM->free(1);

	for (int i = 0; i < 2; i++) {
		LMM->max_limit(i, 1.5);
		LMM->min_limit(i, 0.5);
	}

	LMM->hold(i, val);
}

void TEffTFunction::fit2() {

	LMM->fit();

	z[0] = LMM->a[0];
	z[1] = LMM->a[1];
	z[2] = LMM->a[2];

}

TEffTFunction_v2::TEffTFunction_v2() {

}

TEffTFunction_v2::TEffTFunction_v2(double _K1, double _A0, double _A3, double _beta, stBladeMechanism *BM, double C_1,
								   double C_2, double _g, double _Dwheel) :
	K1(_K1), A0(_A0), beta(_beta), g(_g), Dwheel(_Dwheel), C1(C_1), C2(C_2), A3(_A3) {

	BladeMech = BM;
	Fun_Eff = new stFun_Eff(K1, A0, beta, A3, g, BladeMech, C1, C2, Dwheel);
}

TEffTFunction_v2::~TEffTFunction_v2() {

}

double TEffTFunction_v2::val(double bsr, double Aeff, double Nred, double VGT) {

	return Fun_Eff->val(bsr, Aeff, Nred, VGT);
}

void TEffTFunction_v2::fefft(dVector& x, dVector &a, double &y, dVector &dyda) {

	/*K3 = pow2(x[2] * __cons::Pi * Dwheel) / 2 / Cp;
	 double C = C1 + C2 * x[3];
	 double z = Fun_z(x[0], x[2], x[3], a);
	 double K2 = x[1] * (K4 * C * z + K5);
	 y = -K1 * pow2(x[0]) + K2 * (pow(1 - (K3 / pow2(x[0])), 1 / (g - 1))) * x[0];

	 double Cte = x[1] * K4 * C;

	 dyda[0] = -Cte * x[2] * x[0];

	 dyda[1] = -Cte * x[0];

	 dyda[2] = Cte * x[2];

	 dyda[3] = Cte * pow2(x[3]);

	 dyda[4] = Cte * x[3];

	 dyda[5] = Cte;*/

}

void TEffTFunction_v2::fit(dVector BSR, dVector N, dVector EFF, dVector VGT, dVector Aeff) {

	dMatrix x;
	x.resize(BSR.size());
	dVector y;
	y.resize(BSR.size());
	dVector sig;
	sig.resize(BSR.size());
	dVector dy, aa;
	aa.resize(6, 1.0);
	dy.resize(6, 1.0);

	for(unsigned int i = 0; i < BSR.size(); i++) {
		x[i].resize(4);

		x[i][0] = BSR[i];
		x[i][1] = Aeff[i];
		x[i][2] = N[i];
		x[i][3] = VGT[i];
		y[i] = EFF[i];
		sig[i] = 1;
	}

	CorrEffT = new stCorrEffT2(Fun_Eff);

	aa[0] = 2.80E-05;
	aa[1] = 2.83E+00;
	aa[2] = 0.000211318;
	aa[3] = -0.000220089;
	aa[4] = 0.039615966;
	aa[5] = 1.735936821;

	Fun_Eff->Asign_a(aa);

	LMM = new FitmrqMult(x, y, sig, aa, CorrEffT, 1e-6);

	LMM->max_limit(0, 9);
	LMM->min_limit(0, 0);

	LMM->max_limit(1, 9);
	LMM->min_limit(1, 0);

	LMM->max_limit(2, 9);
	LMM->min_limit(2, 0);

	if(fabs(BladeMech->sVGTtoAlpha1) < 1e-10) {
		LMM->hold(3, 0.);
	} else {
		LMM->max_limit(3, 9);
		LMM->min_limit(3, -9);
	}

	if(fabs(BladeMech->sVGTtoAlpha1) < 1e-10) {
		LMM->hold(4, 0.);
	} else {
		LMM->max_limit(4, 9);
		LMM->min_limit(4, -9);
	}

	LMM->max_limit(5, 9);
	LMM->min_limit(5, -9);

	double temp = Fun_Eff->val(BSR[0], Aeff[0], N[0], VGT[0]);
	double error = temp - EFF[0];

	LMM->fit();

	z = LMM->a;

	Fun_Eff->Asign_a(LMM->a);

	FILE *fich = fopen("TurbineExtrp.txt", "a");
	fprintf(fich, "Chi square: %g\n", LMM->chisq);
	fprintf(fich, "BSR\tEffectiveSection\tSpeed\tPosition\tRealEfficiency\tFittedEfficiency\n");
	for (unsigned int i = 0; i < BSR.size(); i++) {
		fprintf(fich, "%g\t", BSR[i]);
		fprintf(fich, "%g\t", Aeff[i]);
		fprintf(fich, "%g\t", N[i]);
		fprintf(fich, "%g\t", VGT[i]);
		fprintf(fich, "%g\t", EFF[i]);
		fprintf(fich, "%g\n", Fun_Eff->val(BSR[i], Aeff[i], N[i], VGT[i]));
	}
	fclose(fich);

}

void TEffTFunction_v2::fit2() {

}

inline double TEffTFunction_v2::Fun_z(double bsr, double n, double vgt, dVector& z_) {
	return -(z_[0] * n + z_[1]) * bsr + (z_[2] * n + z_[3] * pow2(vgt) + z_[4] * vgt + z_[5]);
}

#pragma package(smart_init)
