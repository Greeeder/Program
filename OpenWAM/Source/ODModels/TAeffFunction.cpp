// ---------------------------------------------------------------------------

#pragma hdrstop

#include "TAeffFunction.h"
// ---------------------------------------------------------------------------

TAeffFunction::TAeffFunction() {
}

TAeffFunction::TAeffFunction(double _AR, double _AN, double _diam, double _g, double _R, double _rendMEAN,
							 double _Dwheel) :
	AR(_AR), AN(_AN), diam(_diam), g(_g), R(_R), rendMEAN(_rendMEAN), Dwheel(_Dwheel) {
}

TAeffFunction::~TAeffFunction() {

}

double TAeffFunction::val(double blade, double REmod, double RendTS) {

	double C1 = (1 + (b + pow2(blade) * diam) / rendMEAN);
	if(C1 < 0)
		C1 = 1e-10;
	else
		C1 = sqrt(C1);

	// double C1 = sqrt(1 + (b + pow2(blade) * diam) / rendMEAN);
	double C4 = AR * a;
	double C6 = sqrt(pow2((AR * REmod * c * d) / (AN * (RendTS * (pow(REmod * d, (g - 1) / g) - 1) + 1))) + 1);

	return C4 * C1 / C6;
}

void TAeffFunction::faeff(dVector& x, dVector &a, double &y, dVector &dyda) {

	double C1 = sqrt(1 + (a[1] + pow2(x[0]) * diam) / rendMEAN);
	double C2 = pow2(AR * x[1] * a[2] * a[3]);
	double C31 = x[1] * a[3];
	double C32 = (g - 1) / g;
	double C33 = x[2] * (pow(C31, C32) - 1) + 1;
	double C3 = pow2(AN * C33);
	double C4 = AR * a[0];
	double C5 = C2 / C3 + 1;
	double C6 = sqrt(C5);
	double C7 = pow3(C6);

	y = C4 * C1 / C6;

	dyda[0] = (AR * C1) / C6;

	dyda[1] = C4 / (2 * rendMEAN * C6 * C1);

	dyda[2] = -(C4 * C2 / a[2] * C1) / (C7 * C3);

	dyda[3] = -(C4 * (2 * C2 / a[3] / C3 - (2 * x[1] * C2 * x[2] * pow(C31,
											C32 - 1) * (g - 1)) / (g * C3 * C33)) * C1) / (2 * C7);
}

void TAeffFunction::fit(dVector Speed, dVector ER, dVector EFF, dVector MR) {

	dMatrix x;
	x.resize(Speed.size());
	dVector y;
	y.resize(Speed.size());
	dVector sig;
	sig.resize(Speed.size());
	dVector aa, dy;
	aa.resize(4, 1);
	dy.resize(4, 1);

	double aux1 = sqrt(g * R);
	double aux2 = 1 / g;
	double aux3 = g - 1;
	double aux4 = 2 / aux3;
	double aux5 = aux3 / g;
	double aux6 = sqrt(aux4 * g * R);

	for(unsigned int i = 0; i < Speed.size(); i++) {
		x[i].resize(3);
		x[i][0] = Speed[i] * __cons::Pi * Dwheel / (aux6 * sqrt(1 - pow(1 / ER[i], aux5)));
		x[i][1] = 2 / (ER[i] + 1);
		x[i][2] = EFF[i];
		y[i] = 1e-6 * MR[i] * aux1 / (g * pow(1 / ER[i], aux2) * sqrt(aux4 * (1 - pow(1 / ER[i], aux5))));
		sig[i] = 0.0002;

	}

	CorrAeff = new stCorrAeff(AR, AN, diam, g, rendMEAN);

	// CorrAeff->Fun(x[0], aa, y[0], dy);

	LMM = new FitmrqMult(x, y, sig, aa, CorrAeff, 1e-10);

	for(int i = 0; i < 4; i++) {
		LMM->max_limit(i, 3.);
		LMM->min_limit(i, 0.);
	}

	LMM->fit();

	a = LMM->a[0];
	b = LMM->a[1];
	c = LMM->a[2];
	d = LMM->a[3];

	FILE *fich = fopen("TurbineExtrp.txt", "a");
	fprintf(fich, "Chi square: %g\n", LMM->chisq);
	fprintf(fich, "BSR\tRE_mod\tEfficiency\tRealEffectiveSection\tFittedEffectiveSection\n");
	for (unsigned int i = 0; i < Speed.size(); i++) {
		fprintf(fich, "%g\t", x[i][0]);
		fprintf(fich, "%g\t", x[i][1]);
		fprintf(fich, "%g\t", x[i][2]);
		fprintf(fich, "%g\t", y[i]);
		fprintf(fich, "%g\n", val(x[i][0], x[i][1], x[i][2]));
	}
	fclose(fich);

}

TAeffFunction_v2::TAeffFunction_v2() {
}

TAeffFunction_v2::TAeffFunction_v2(double _AR, double _diam, double _g, double _rendMEAN, stBladeMechanism *BM) :
	AR(_AR), diam(_diam), g(_g), rendMEAN(_rendMEAN) {
	BladeMech = BM;
	g1 = (g - 1) / g;
	Fun_Aeff = new stFun_Aeff(AR, diam, g, rendMEAN, BladeMech);
}

TAeffFunction_v2::~TAeffFunction_v2() {

}

double TAeffFunction_v2::val(double blade, double RE, double RendTS, double Pos) {

	return Fun_Aeff->val(blade, RE, RendTS, Pos);
}

void TAeffFunction_v2::fit(dVector BSR, dVector Pos, dVector EFF, dVector ER, dVector AEFF) {

	dMatrix x;
	x.resize(BSR.size());
	dVector y;
	y.resize(BSR.size());
	dVector sig;
	sig.resize(BSR.size());
	dVector aa, dy;
	aa.resize(7, 1);
	dy.resize(7, 1);

	for(unsigned int i = 0; i < BSR.size(); i++) {
		x[i].resize(4);
		x[i][0] = BSR[i];
		x[i][1] = Pos[i];
		x[i][2] = EFF[i];
		x[i][3] = ER[i];
		y[i] = AEFF[i];
		sig[i] = 0.0002;

	}

	CorrAeff = new stCorrAeff2(Fun_Aeff);

	//Initial values
	aa[0] = 0.47;
	aa[1] = 0;
	aa[2] = 1.037364431;
	aa[3] = -0.00817246;
	aa[4] = 1.310154002;
	aa[5] = 0;
	aa[6] = 0.530250916;

	aa[0] = 0.4;
	aa[1] = 0;
	aa[2] = 1;
	aa[3] = 0;
	aa[4] = 1;
	aa[5] = 0;
	aa[6] = 0.5;

	Fun_Aeff->Asign_a(aa);

	LMM = new FitmrqMult(x, y, sig, aa, CorrAeff, 1e-10);

	LMM->max_limit(0, 1.0);
	LMM->min_limit(0, 0.1);

	if(fabs(BladeMech->sVGTtoAlpha1) < 1e-10) {
		LMM->hold(1, 0.0);
	} else {
		LMM->max_limit(1, 0.02763);
		LMM->min_limit(1, 0.);
	}

	LMM->max_limit(2, 2.75);
	LMM->min_limit(2, 0);

	if(fabs(BladeMech->sVGTtoAlpha1) < 1e-10) {
		LMM->hold(3, 0.);
	} else {
		LMM->max_limit(3, 0);
		LMM->min_limit(3, -0.0308);
	}

	LMM->max_limit(4, 5);
	LMM->min_limit(4, 0.28);

	if(fabs(BladeMech->sVGTtoAlpha1) < 1e-10) {
		LMM->hold(5, 0.);
	} else {
		LMM->max_limit(5, 0);
		LMM->min_limit(5, -0.0075);
	}

	LMM->max_limit(6, 2);
	LMM->min_limit(6, 0.6963);

	LMM->fit();

	a = LMM->a[0];
	b1 = LMM->a[1];
	b2 = LMM->a[2];
	c1 = LMM->a[3];
	c2 = LMM->a[4];
	d1 = LMM->a[5];
	d2 = LMM->a[6];

	Fun_Aeff->Asign_a(LMM->a);

	FILE *fich = fopen("TurbineExtrp.txt", "a");
	fprintf(fich, "Chi square: %g\n", LMM->chisq);
	fprintf(fich, "BSR\tPosition\tEfficiency\tExpansionRatio\tRealEffectiveSection\tFittedEffectiveSection\n");
	for (unsigned int i = 0; i < BSR.size(); i++) {
		fprintf(fich, "%g\t" , BSR[i]);
		fprintf(fich, "%g\t", Pos[i]);
		fprintf(fich, "%g\t", EFF[i]);
		fprintf(fich, "%g\t", ER[i]);
		fprintf(fich, "%g\t", AEFF[i]);
		fprintf(fich, "%g\n", Fun_Aeff->val(BSR[i], ER[i], EFF[i], Pos[i]));
	}
	fclose(fich);
}

void TAeffFunction_v2::faeff(dVector& x, dVector &a, double &y, dVector &dyda) {

	double AN = BladeMech->ToAN(x[1]);
	double b = a[1] * x[1] + a[2];
	double c = a[3] * x[1] + a[4];
	double d = a[5] * x[1] + a[6];
	double C7 = (x[3] - 1);
	double invRE_R = 1 / (1 + d * C7);
	double C1 = (1 + (b + pow2(x[0]) * diam) / rendMEAN);
	if(C1 < 0)
		C1 = 1e-10;
	else
		C1 = sqrt(C1);

	// double C1 = sqrt(1 + (b + pow2(blade) * diam) / rendMEAN);
	double C2 = (1 - x[2] * (1 - pow(invRE_R, g1)));
	double C4 = AR * a[0];

	double C6 = sqrt(pow2((AR * invRE_R * c) / (AN * C2)) + 1);

	y = C4 * C1 / C6;

	dyda[0] = y / a[0];

	dyda[1] = C4 / (2 * rendMEAN * C6 * C1) * x[1];

	dyda[2] = C4 / (2 * rendMEAN * C6 * C1);

	dyda[3] = -C4 * C1 / pow3(C6) * pow2(AR * invRE_R / C2 / AN) * c * x[1];

	dyda[4] = -C4 * C1 / pow3(C6) * pow2(AR * invRE_R / C2 / AN) * c;

	double C3 = 2 - C2;
	double C5 = pow2(c * AR / AN);

	dyda[5] = C4 * C1 * C5 * C7 * (pow3(invRE_R) / pow2(C3) - x[2] * g1 * pow(invRE_R,
								   -1 / g) * pow4(invRE_R) / pow3(C3)) / pow150(C5 * pow2(invRE_R) / pow2(C3) + 1) * x[1];

	dyda[6] = C4 * C1 * C5 * C7 * (pow3(invRE_R) / pow2(C3) - x[2] * g1 * pow(invRE_R,
								   -1 / g) * pow4(invRE_R) / pow3(C3)) / pow150(C5 * pow2(invRE_R) / pow2(C3) + 1);
}
#pragma package(smart_init)
