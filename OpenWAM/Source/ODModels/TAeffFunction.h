// ---------------------------------------------------------------------------

#ifndef TAeffFunctionH
#define TAeffFunctionH
// ---------------------------------------------------------------------------

#include "Fit_MRQ.h"

class TAeffFunction {
  private:

	double a;
	double b;
	double c;
	double d;

	double AR;
	double AN;
	double diam;
	double g;
	double R;
	double rendMEAN;
	double Dwheel;

	stCorrelation *CorrAeff;

	FitmrqMult *LMM;

  public:

	TAeffFunction();

	TAeffFunction(double AR, double AN, double diam, double g, double R, double rendMEAN, double Dwheel);

	~TAeffFunction();

	double get_a() {
		return a;
	}
	;

	double get_b() {
		return b;
	}
	;

	double get_c() {
		return c;
	}
	;

	double get_d() {
		return d;
	}
	;

	double val(double blade, double REmod, double RendTS);

	void faeff(dVector& x, dVector &a, double &y, dVector &dyda);

	void fit(dVector BSR, dVector ER, dVector EFF, dVector MR);

};

class TAeffFunction_v2 {
  private:

	double a;
	double b1, b2;
	double c1, c2;
	double d1, d2;

	double AR;
	double diam;
	double g;
	double g1;
	double rendMEAN;

	stCorrelation *CorrAeff;
	stBladeMechanism *BladeMech;
	stFun_Aeff *Fun_Aeff;

	FitmrqMult *LMM;

  public:

	TAeffFunction_v2();

	TAeffFunction_v2(double AR, double diam, double g1, double rendMEAN, stBladeMechanism *BM);

	~TAeffFunction_v2();

	double get_a() {
		return a;
	}

	double get_b1() {
		return b1;
	}

	double get_b2() {
		return b2;
	}

	double get_c1() {
		return c1;
	}

	double get_c2() {
		return c2;
	}

	double get_d1() {
		return d1;
	}

	double get_d2() {
		return d2;
	}

	double val(double blade, double RE, double RendTS, double AN);

	void faeff(dVector& x, dVector &a, double &y, dVector &dyda);

	void fit(dVector BSR, dVector Pos, dVector EFF, dVector ER, dVector AEFF);

};
#endif
