/*
 * TPipeExtHeatW.cpp
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#include "TPipeExtHeatW.h"

TPipeExtHeatW::TPipeExtHeatW(int ncells, double vel, double mult, ArrayXd dext, double cellsize,
	double emis, double Tw) {

	Cooler = new TFluidWater();
	Text.setConstant(ncells, Tw);
	Visc.setZero(ncells);
	Cond.setZero(ncells);
	Rho.setZero(ncells);
	Pr.setZero(ncells);
	Vext.setConstant(ncells, vel);
	Dext = dext;
	h.setZero(ncells);
	Aext = Dext * __cons::Pi * cellsize;
	Emissivity = emis;
	HeatMultiplier = mult;
	Qtot.setZero(ncells);

}

TPipeExtHeatW::~TPipeExtHeatW() {
	// TODO Auto-generated destructor stub
}

Eigen::VectorXd TPipeExtHeatW::Heat(Eigen::VectorXd Twall) {

	if (HeatMultiplier > 0){
		double n1 = 0., n2 = 0.;

		Tavg = 0.5 * (Twall.array() + Text);
		DeltaT = Text - Twall.array();

		for (int i = 0; i < Tavg.size(); i++) {
			Rho(i) = Cooler->FunRHO(Tavg(i));
			if (Tavg(i) < 373) {
				Visc(i) = Cooler->FunVisc(Tavg(i));
				Cond(i) = Cooler->Funk(Tavg(i));
				Pr(i) = Cooler->FunPr(Tavg(i));
			}
			else {
				Visc(i) = 0.000282;
				Pr(i) = 1.76;
				Cond(i) = 0.6775;
			}
		}

		Re = Rho * Vext * Dext / Visc;

		// Termino de conveccion de Churchill Bernstein

		for (int i = 0; i < Re.size(); i++) {
			if ((2e4 < Re(i)) && (Re(i) < 4e5)) {
				n1 = 0.5;
				n2 = 1.;
			}
			else {
				n1 = 0.625;
				n2 = 0.8;
			}
			h(i) = 0.3 + 0.62 * sqrt(Re(i)) * cbrt(Pr(i)) / pow025(1 + pow(0.4 / Pr(i), 0.666666)) * pow(1 + pow(Re(i) / 282000,
				n1), n2) * Cond(i) / Dext(i);
		}

		Qconv = h * Aext * DeltaT;

		// Termino de radiación

		Qrad = __cons::Sigma * Emissivity * (Text.pow(4) - Twall.array().pow(4)) * Aext;

		Qtot = (Qconv + Qrad) * HeatMultiplier;
	}

	return Qtot.matrix();
}

