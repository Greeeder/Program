/*
 * TPortExtHeat.cpp
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#include "TPortExtHeat.h"

TPortExtHeat::TPortExtHeat() {

	TComponentArray_ptr Water;
	Water->resize(1);
	Water->at(0) = make_shared<TChH2Ol>("Water", __R::Universal, __PM::Air);

	Cooler = new TFluid(Water);

	Rho = Eigen::VectorXd::Ones(ncells) * 980; // Engine coolant density;
	Qrad = Eigen::VectorXd::Zero(ncells);

}

TPortExtHeat::~TPortExtHeat() {
	// TODO Auto-generated destructor stub
}

Eigen::VectorXd TPortExtHeat::Heat(Eigen::VectorXd Twall) {

	double n1 = 0., n2 = 0.;

	Tavg = Text;
	DeltaT = Text - Twall.array();

	double Te = Text.mean();  // Usually the external temperature is the same for all cells

	if(Te < 373) {
		Visc = Eigen::VectorXd::Ones(ncells) * Cooler->FunVisc(Te);
		Cond = Eigen::VectorXd::Ones(ncells) * Cooler->Funk(Te);
		Pr = Eigen::VectorXd::Ones(ncells) * Cooler->Funk(Te);
	} else {
		Visc = Eigen::VectorXd::Ones(ncells) * 0.000282;
		Pr = Eigen::VectorXd::Ones(ncells) * 1.76;
		Cond = Eigen::VectorXd::Ones(ncells) * 0.6775;
	}

	for(int i = 0; i < Twall.size(); i++) {
		if(Twall(i) < 373)
			ViscWall(i) = Cooler->FunVisc(Twall(i));
		else
			ViscWall(i) = 0.000282;
	}

	Re = Rho * Vext * Lchar / Visc;

	// Termino de conveccion de Churchill Bernstein

	for(int i = 0; i < Re.size(); i++) {

		h(i) = 0.027 * (1 + 24.2 / pow(2.3, 0.7) / pow025(Re(i))) * pow(Re(i), 0.8) * cbrt(Pr(i)) * pow(Visc(i) / ViscWall(i),
				0.14) * Cond(i) / (Lchar / 2.3);

	}

	Qconv = h * Aext * DeltaT;

	// Termino de radiación

	Qtot = Qconv + Qrad;

	return Qtot.matrix();
}

void TPortExtHeat::setVext(double n) {

	Vext = ArrayXd::Ones(ncells) * 5.64268e-7 * __units::RPMToRPS(n) * TorqueMaxPower / Ncil / pow2(Lchar);

}

