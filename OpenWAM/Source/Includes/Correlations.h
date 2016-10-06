// ---------------------------------------------------------------------------

#ifndef CorrelationsH
#define CorrelationsH

#include "Globales.h"

struct stBladeMechanism {

	double sVGTtoAlpha1;
	double sVGTtoAlpha2;
	double sR2geom;
	double sLTE;
	double sZ0;
	double sHBlade;
	double sGamma;

	stBladeMechanism() {
	}
	;

	double virtual Toz3(double Position) = 0;

	double virtual ToL2(double Position) = 0;

	double virtual ToAN(double Position) = 0;

	double virtual ToAlpha(double Position) = 0;

};

struct stBMVaned: stBladeMechanism {

	stBMVaned(double VGTtoAlpha1, double VGTtoAlpha2, double R2geom, double LTE, double Z0, double HBlade);

	double Toz3(double Position);

	double ToL2(double Position);

	double ToAN(double Position);

	double ToAlpha(double Position);

};

struct stBMVaneless: stBladeMechanism {

	double sAN;

	stBMVaneless(double R2geom, double HBlade);

	double Toz3(double Position);

	double ToL2(double Position);

	double ToAN(double Position);

	double ToAlpha(double Position);

};

struct stFun_z {

	double val(double bsr, double n, double vgt, dVector& z_) {
		return -(z_[0] * n + z_[1]) * bsr + (z_[2] * n + z_[3] * pow2(vgt) + z_[4] * vgt + z_[5]);
	}
};

struct stFun_Eff {
	dVector p;
	dVector a;
	stBladeMechanism *BladeMech;
	double Cp;
	double K3;
	double C;
	stFun_z Fun_z;

	stFun_Eff(double K1, double A0, double beta, double A3, double g, stBladeMechanism* BM, double C1, double C2,
			  double Dwheel) {
		p.resize(8);
		p[0] = K1;
		p[1] = A0;
		p[2] = beta;
		p[3] = A3;
		p[4] = g;
		p[5] = C1;
		p[6] = C2;
		p[7] = Dwheel;
		BladeMech = BM;
		Cp = __Gamma::Cp;
		BladeMech = BM;
	}
	;

	void Asign_a(dVector _a) {
		a = _a;
	}

	double val(double bsr, double aeff, double n, double vgt) {
		double y;
		K3 = pow2(n * __cons::Pi * p[7]) / 2 / Cp;
		C = p[5] * vgt + p[6];
		if(bsr <= 0.)
			y = 1e-6;
		else {

			double zeta = Fun_z.val(bsr, n * 60., vgt, a);
			double K2 = 2 * aeff / p[1] * (zeta * p[1] / p[3] * BladeMech->Toz3(vgt) * C * sin(BladeMech->ToAlpha(vgt)) + sqrt(
											   p[0] / 2) * tan(p[2]));
			double KK3 = 1 - (K3 / pow2(bsr));
			if(KK3 < 0) {
				y = -p[0] * pow2(bsr);
			} else {
				y = -p[0] * pow2(bsr) + K2 * (pow(KK3, 1 / (p[4] - 1))) * bsr;
			}
		}
		return y;
	}
};

struct stFun_Aeff {

	dVector p;
	dVector a;
	stBladeMechanism *BladeMech;
	double C1;
	double C2;
	double C22;
	double C31;
	double C33;
	double C3;
	double C4;
	double C5;
	double C6;
	double invPIR;
	double AN;
	double b;
	double c;
	double d;
	double g1;

	stFun_Aeff(double AR, double diam, double g, double rendMEAN, stBladeMechanism* BM) {
		p.resize(4);
		p[0] = AR;
		p[1] = diam;
		p[2] = g;
		p[3] = rendMEAN;
		BladeMech = BM;
		g1 = (g - 1) / g;
	}
	;

	void Asign_a(dVector _a) {
		a = _a;
	}
	;

	double val(double bsr, double ER, double EFF, double vgt) {
		AN = BladeMech->ToAN(vgt);
		b = a[1] * vgt + a[2];
		c = a[3] * vgt + a[4];
		d = a[5] * vgt + a[6];

		invPIR = 1 / (1 + d * (ER - 1));

		double C0 = (b + pow2(bsr) * p[1]) / p[3];
		if (C0 < -0.999)
			C0 = 0.001;
		C1 = sqrt(1 + C0);
		C22 = pow2(c * p[0] / AN);
		C2 = C22 * pow2(invPIR);
		C33 = EFF * (pow(invPIR, g1) - 1) + 1;
		C3 = pow2(C33);
		C4 = p[0] * a[0];
		C5 = C2 / C3 + 1;
		C6 = sqrt(C5);

		return C4 * C1 / C6;
	}
	;

};

struct stCorrelation {

	dVector p;

	stCorrelation() {

	}
	;

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {
	}
	;

	void virtual Fun(double & x, dVector & a, double&y, dVector & dyda) {
	}
	;

};

struct stCorrEffT: stCorrelation {

	stCorrEffT(double K1, double A0, double beta, double K3, double g) {
		p.resize(5);
		p[0] = K1;
		p[1] = A0;
		p[2] = beta;
		p[3] = K3;
		p[4] = g;
	}
	;

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {
		double C1 = a[1] - x[0] * a[0];
		double C2 = 1 / (p[4] - 1);
		double C3 = pow(1 - p[3] / pow2(x[0]), C2);
		double C4 = 2 * x[1] / p[1];

		y = -p[0] * pow2(x[0]) + ((C4 * (a[2] * sin(C1) + sqrt(p[0] / 2) * tan(p[2]))) * C3 * x[0]);

		dyda[0] = -C4 * pow2(x[0]) * a[2] * cos(C1) * C3;

		dyda[1] = C4 * x[0] * a[2] * cos(C1) * C3;

		dyda[2] = C4 * x[0] * sin(C1) * C3;

	}
	;

};

struct stCorrEffT2: stCorrelation {

	stFun_Eff *Fun_Eff;

	stCorrEffT2(stFun_Eff* _fun_eff) {
		Fun_Eff = _fun_eff;
		p = Fun_Eff->p;
	}
	;

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {

		Fun_Eff->Asign_a(a);
		y = Fun_Eff->val(x[0], x[1], x[2], x[3]);

		double KK3 = 1 - (Fun_Eff->K3 / pow2(x[0]));
		if(KK3 < 0.)
			KK3 = 0.;
		double Const = 2 * x[1] / p[1] * (pow(KK3,
											  1 / (p[4] - 1))) * x[0] * p[1] / p[3] * Fun_Eff->BladeMech->Toz3(x[3]) * Fun_Eff->C * sin(Fun_Eff->BladeMech->ToAlpha(
													  x[3]));

		dyda[0] = -Const * x[2] * x[0];

		dyda[1] = -Const * x[0];

		dyda[2] = Const * x[2];

		dyda[3] = Const * pow2(x[3]);

		dyda[4] = Const * x[3];

		dyda[5] = Const;

	}
	;

};

struct stCorrAeff: stCorrelation {

	stCorrAeff(double AR, double AN, double diam, double g, double rendMEAN) {
		p.resize(5);
		p[0] = AR;
		p[1] = AN;
		p[2] = diam;
		p[3] = g;
		p[4] = rendMEAN;
	}
	;

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {

		double C1 = sqrt(1 + (a[1] + pow2(x[0]) * p[2]) / p[4]);
		double C2 = pow2(p[0] * x[1] * a[2] * a[3]);
		double C31 = x[1] * a[3];
		double C32 = (p[3] - 1) / p[3];
		double C33 = x[2] * (pow(C31, C32) - 1) + 1;
		double C3 = pow2(p[1] * C33);
		double C4 = p[0] * a[0];
		double C5 = C2 / C3 + 1;
		double C6 = sqrt(C5);
		double C7 = pow3(C6);

		y = C4 * C1 / C6;

		dyda[0] = (p[0] * C1) / C6;

		dyda[1] = C4 / (2 * p[4] * C6 * C1);

		dyda[2] = -(C4 * C2 / a[2] * C1) / (C7 * C3);

		dyda[3] = -(C4 * (2 * C2 / a[3] / C3 - (2 * x[1] * C2 * x[2] * pow(C31,
												C32 - 1) * (p[3] - 1)) / (p[3] * C3 * C33)) * C1) / (2 * C7);

	}
	;

};

struct stCorrAeff2: stCorrelation {

	stFun_Aeff *Fun_Aeff;

	stCorrAeff2(stFun_Aeff *_fun_aeff) {
		Fun_Aeff = _fun_aeff;
		p = Fun_Aeff->p;
	}
	;

	void virtual Fun(dVector & x, dVector & a, double&y, dVector & dyda) {

		Fun_Aeff->Asign_a(a);
		y = Fun_Aeff->val(x[0], x[3], x[2], x[1]);

		double C7 = pow3(Fun_Aeff->C6);

		//y = Fun_Aeff->C4 * Fun_Aeff->C1 / Fun_Aeff->C6;

		dyda[0] = (p[0] * Fun_Aeff->C1) / Fun_Aeff->C6;

		dyda[2] = Fun_Aeff->C4 / (2 * p[3] * Fun_Aeff->C6 * Fun_Aeff->C1);

		dyda[1] = dyda[2] * x[1];

		dyda[4] = -(Fun_Aeff->C4 * Fun_Aeff->C1 * (Fun_Aeff->C5 - 1) / Fun_Aeff->c) / (C7 * Fun_Aeff->C3);

		dyda[3] = dyda[4] * x[1];

		dyda[6] = Fun_Aeff->C4 * Fun_Aeff->C1 * Fun_Aeff->C22 * (x[3] - 1) / Fun_Aeff->C3
				  * (pow2(Fun_Aeff->invPIR) - x[2] * Fun_Aeff->g1 * pow(Fun_Aeff->invPIR,
						  -(4 * p[2] + 1) / p[2]) / Fun_Aeff->C33) / pow150(Fun_Aeff->C22 * pow2(Fun_Aeff->invPIR) / Fun_Aeff->C3 + 1);

		dyda[5] = dyda[6] * x[1];

	}
	;

};

// ---------------------------------------------------------------------------
#endif
