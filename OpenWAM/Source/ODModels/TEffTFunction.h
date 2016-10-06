// ---------------------------------------------------------------------------

#ifndef TEffTFunctionH
#define TEffTFunctionH
// ---------------------------------------------------------------------------

#include "Fit_MRQ.h"
#include "TAeffFunction.h"

void fefft2(dVector & x, dVector & a, double&y, dVector & dyda);

class TEffTFunction {
  private:

	dVector z;
	dMatrix x;
	dVector y;
	dVector sig;

	double K1;
	double A0;
	double beta;
	double K3;
	double g;
	double Dwheel;

	FitmrqMult *LMM;

	stCorrelation *CorrEffT;

  public:

	TEffTFunction();

	TEffTFunction(double _K1, double _A0, double _beta, double _K3, double _g, double _Dwheel);

	~TEffTFunction();

	double get_z(int i) {
		return z[i];
	}
	;

	void put_z(int i, double val) {
		z[i] = val;
	}
	;

	double val(double blade, double Aeff, double ER);

	void fefft(dVector& x, dVector &a, double &y, dVector &dyda);

	void fit(dVector BSR, dVector ER, dVector EFF, dVector MR, TAeffFunction *AeffFun);

	void fit2();

	void FixParameter(int i, double val);

};

class TEffTFunction_v2 {
  private:
	dMatrix x;

	dVector z;
	double K1;
	double g;
	double Cp;
	double C1;
	double C2;
	double A0;
	double A3;
	double beta;
	double Dwheel;

	FitmrqMult *LMM;

	stCorrelation *CorrEffT;
	stBladeMechanism *BladeMech;
	stFun_Eff *Fun_Eff;

	double Fun_z(double bsr, double n, double vgt, dVector& z_);

  public:
	TEffTFunction_v2();

	TEffTFunction_v2(double _K1, double _A0, double _A3, double _beta, stBladeMechanism *BM, double C_1, double C_2,
					 double _g, double _Dwheel);

	~TEffTFunction_v2();

	double val(double bsr, double Aeff, double Nred, double VGT);

	void fefft(dVector& x, dVector &a, double &y, dVector &dyda);

	void fit(dVector BSR, dVector N, dVector EFF, dVector VGT, dVector Aeff);

	void fit2();

	double get_z(int i) {
		return z[i];
	}
	;

};
#endif
