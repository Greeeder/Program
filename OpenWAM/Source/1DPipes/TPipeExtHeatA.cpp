/*
 * TPipeExtHeatA.cpp
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#include "TPipeExtHeatA.h"

TPipeExtHeatA::TPipeExtHeatA(int ncells, double vel, double mult, ArrayXd dext, double cellsize,
	double emis) {

	Cooler = new TFluidAir();
	Text.setConstant(ncells, __Ambient::T_K);
	Pext.setConstant(ncells, __Ambient::p_Pa);
	Visc.setZero(ncells);
	Cond.setZero(ncells);
	Vext.setConstant(ncells, vel);
	Dext = dext;
	h.setZero(ncells);
	Aext = Dext * __cons::Pi * cellsize;
	Emissivity = emis;
	HeatMultiplier = mult;

}

TPipeExtHeatA::~TPipeExtHeatA() {
	// TODO Auto-generated destructor stub
}

VectorXd TPipeExtHeatA::Heat(Eigen::VectorXd Twall) {

	double n1 = 0., n2 = 0.;

	Tavg = 0.5 * (Twall.array() + Text);
	DeltaT = Text - Twall.array();

	Rho = Pext / (Cooler->FunR() * Tavg);

	for(int i = 0; i < Tavg.size(); i++) {
		Visc(i) = Cooler->FunVisc(Tavg(i));
		Cond(i) = Cooler->Funk(Tavg(i));
	}

	double _Pr = 0.7;

	Re = Rho * Vext * Dext / Visc;

	// Termino de conveccion de Churchill Bernstein
	
	// TODO: Add natural convection

	for(int i = 0; i < Re.size(); i++) {
		if((2e4 < Re(i)) && (Re(i) < 4e5)) {
			n1 = 0.5;
			n2 = 1.;
		} else {
			n1 = 0.625;
			n2 = 0.8;
		}
		h(i) = 0.3 + 0.62 * sqrt(Re(i)) * cbrt(_Pr) / pow025(1 + pow(0.4 / _Pr, 0.666666)) * pow(1 + pow(Re(i) / 282000, n1),
				n2) * Cond(i) / Dext(i);
	}

	Qconv = h * Aext * DeltaT;

	// Termino de radiación

	Qrad = __cons::Sigma * Emissivity * Aext * (Text.pow(4) - Twall.array().pow(4));

	Qtot = (Qconv + Qrad) * HeatMultiplier;

	return Qtot.matrix();
}

