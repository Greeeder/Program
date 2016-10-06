/* --------------------------------------------------------------------------------*\
 ==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
  \\ |  X  | //  W ave     |
   \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
    \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*--------------------------------------------------------------------------------*/

/**
 * @file PipeMethod.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 *
 * @section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenWAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 * This file defines a pipe computation method.
 */

#include "PipeMethod.hpp"
#include "BoundaryCondition.hpp"

double TPipeMethod::ComputeCharacteristic(double& velocidadp, double& asonidop, int ind, double dist, int signo,
		double entropia, double DeltaTiempo) {
#ifdef usetry
	try {
#endif

		double characteristic = 0.;
		int ind1 = ind + signo;

		double w0 = FPipe->FU0(0, ind);
		double w1 = FPipe->FU0(1, ind);
		double w2 = FPipe->FU0(2, ind);
		double gammap = FPipe->FGamma(ind);
		double Rmezclap = FPipe->FR(ind);
		double diamep = FPipe->FD(ind);

// 		double tptubop = FTPTubo[0][ind];
// 		double hip = Fhi[ind];
// 		double rhop = Frho[ind];
// 		double Rep = FRe[ind];

		if(dist < 1. && dist > 0) {

			w0 = Interpola(FPipe->FU0(0, ind), FPipe->FU0(0, ind1), 1., dist);
			w1 = Interpola(FPipe->FU0(1, ind), FPipe->FU0(1, ind1), 1., dist);
			w2 = Interpola(FPipe->FU0(2, ind), FPipe->FU0(2, ind1), 1., dist);
			gammap = Interpola(FPipe->FGamma(ind), FPipe->FGamma(ind1), 1., dist);
			Rmezclap = Interpola(FPipe->FR(ind), FPipe->FR(ind1), 1., dist);
			diamep = Interpola(FPipe->FD(ind), FPipe->FD(ind1), 1., dist);
// 			if (FCoefAjusTC != 0 || FCoefAjusFric != 0) {
// 				tptubop = Interpola(FTPTubo[0][ind], FTPTubo[0][ind1], 1., dist);
// 				hip = Interpola(Fhi[ind], Fhi[ind1], 1., dist);
// 				rhop = Interpola(Frho[ind], Frho[ind1], 1., dist);
// 				Rep = Interpola(FRe[ind], FRe[ind1], 1., dist);
// 			}

		}
		double gamma1p = __Gamma::G1(gammap);
		double gamma3p = __Gamma::G3(gammap);
		double gamma5p = __Gamma::G5(gammap);
		velocidadp = w1 / w0;
		asonidop = Sqrt(gammap * gamma1p * (w2 / w0 - pow2(velocidadp) / 2));
		characteristic = asonidop - signo * gamma3p * velocidadp;

// 		// Las siguientes expresiones se pueden encontrar en la Tesis de Corberan
// 		// Pagina 22
// 		/* variacion debida a la transmision del calor */
// 		/* ------------------------------------------ */
// 		if (FCoefAjusTC != 0) {
//
// 			double q = 0;
//
// 			double tgasp = pow2(asonidop) / (gammap * Rmezclap);
//
// 			TransmisionCalor(tgasp, diamep, q, hip, rhop, tptubop);
//
// 			double dacal = gamma3p * gamma1p * DeltaTiempo * q * FCoefAjusTC / asonidop;
//
// 			caracteristica += dacal;
//
// 		}
		/* variacion debida a la variacion entropia */
		/* ---------------------------------------- */

		double presionp = __units::PaToBar((w2 - pow2(w1) / 2. / w0) * gamma1p / __geom::Circle_area(diamep));
		double entropiap = asonidop / pow(presionp, gamma5p);
		double increentropia = entropia * __cons::ARef - entropiap;
		double daen = asonidop * increentropia / entropiap;
		characteristic += daen;

		/* variacion debida al cambio de seccion */
		/* ------------------------------------- */
		if(FPipe->FD(ind1) != FPipe->FD(ind)) {
			double daar = 0.;
			if(signo == 1) {
				daar = -gamma3p * asonidop * velocidadp * 2 * (FPipe->FD(ind1) - FPipe->FD(ind)) * DeltaTiempo /
					   (diamep * FPipe->FXref);
			} else if(signo == -1) {
				daar = -gamma3p * asonidop * velocidadp * 2 * (FPipe->FD(ind) - FPipe->FD(ind1)) * DeltaTiempo /
					   (diamep * FPipe->FXref);
			}
			characteristic += daar;
		}

// 		/* variacion debida al termino de friccion */
// 		/* --------------------------------------- */
//
// 		if (velocidadp != 0. && FCoefAjusFric != 0.) {
//
// 			double velabs = fabs(velocidadp);
// 			double f = 0;
//
// 			Colebrook(FFriccion, diamep, f, Rep);
//
// 			double dafric = signo * gamma1p * (1. + signo * gamma1p * velocidadp / asonidop) * f * FCoefAjusFric
// 				* pow3(velocidadp) * DeltaTiempo / (diamep * velabs);
//
// 			caracteristica += dafric;
// 		}
		return characteristic;

#ifdef usetry
	} catch(exception & N) {
		std::cout << "ERROR: TPipeMethod::ComputeCharacteristic " << FPipe->FPipeNumber << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
#endif
}

double TPipeMethod::ComputeEntropy(double& velocidadp, int ind, double dist, int signo, double DeltaTiempo,
								   int indiceCC) {
#ifdef usetry
	try {
#endif

		double entropy = 0.;
		int ind1 = ind + signo;

		double w0 = FPipe->FU0(0, ind);
		double w1 = FPipe->FU0(1, ind);
		double w2 = FPipe->FU0(2, ind);
		double gammap = FPipe->FGamma(ind);
		double Rmezclap = FPipe->FR(ind);
		double diamep = FPipe->FD(ind);

// 		double tptubop = FTPTubo[0][ind];
// 		double hip = Fhi[ind];
		double rhop = FPipe->Frho[ind];
// 		double Rep = FRe[ind];

// 		for (int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
// 			FFraccionMasicaCC[indiceCC][j] = FFraccionMasicaEspecie[ind][j];
// 		}

		if(dist > 0. || dist < 1.) {
			w0 = Interpola(FPipe->FU0(0, ind), FPipe->FU0(0, ind1), 1., dist);
			w1 = Interpola(FPipe->FU0(1, ind), FPipe->FU0(1, ind1), 1., dist);
			w2 = Interpola(FPipe->FU0(2, ind), FPipe->FU0(2, ind1), 1., dist);
			gammap = Interpola(FPipe->FGamma(ind), FPipe->FGamma(ind1), 1., dist);
			Rmezclap = Interpola(FPipe->FR(ind), FPipe->FR(ind1), 1., dist);
			diamep = Interpola(FPipe->FD(ind), FPipe->FD(ind1), 1., dist);
// 			if (FCoefAjusTC != 0 || FCoefAjusFric != 0) {
// 				tptubop = Interpola(FTPTubo[0][ind], FTPTubo[0][ind1], 1., dist);
// 				hip = Interpola(Fhi[ind], Fhi[ind1], 1., dist);
// 				Rep = Interpola(FRe[ind], FRe[ind1], 1., dist);
// 			}

// 			for (int j = 0; j < FNumeroEspecies - FIntEGR; j++) {
// 				FFraccionMasicaCC[indiceCC][j] = Interpola(FFraccionMasicaEspecie[ind][j],
// 					FFraccionMasicaEspecie[ind1][j], 1., dist);
// 			}
		}

		double gamma1p = __Gamma::G1(gammap);
		double gamma3p = __Gamma::G3(gammap);
		double gamma5p = __Gamma::G5(gammap);
		velocidadp = w1 / w0;
		rhop = w0 / __geom::Circle_area(diamep);
		double asonidop = Sqrt(gammap * gamma1p * (w2 / w0 - pow2(velocidadp) / 2));
		double presionp = __units::PaToBar((w2 - pow2(w1) / 2. / w0) * gamma1p / __geom::Circle_area(diamep));
		double entropiap = asonidop / pow(presionp, gamma5p);
		entropy = entropiap;

// 		FAreaCC[indiceCC] = FArea[ind];
// 		FDensidadCC[indiceCC] = rhop;
// 		FVelocidadCC[indiceCC] = velocidadp; // Se que no hace falta hacerlo asi,pero me queda más ordenado.

// 		/* variacion de la entropia debida a la transmision del calor */
// 		/* ------------------------------------------ */
// 		if (!DoubEqZero(FCoefAjusTC)) {
//
// 			double q = 0;
// 			double tgasp = pow2(asonidop) / (gammap * Rmezclap);
//
// 			TransmisionCalor(tgasp, diamep, q, hip, rhop, tptubop);
//
// 			// Las siguientes expresiones estan en la Tesis de Corberan. Pagina 23
// 			double dacal = gamma3p * entropiap * q * FCoefAjusTC * DeltaTiempo / pow2(asonidop);
//
// 			entropia += dacal;
// 		}

// 		/* variacion de la entropia debida al termino de friccion */
// 		/* --------------------------------------- */
// 		if (!DoubEqZero(velocidadp) && !DoubEqZero(FCoefAjusFric)) {
//
// 			double f = 0.;
// 			double velabs = fabs(velocidadp);
// 			Colebrook(FFriccion, diamep, f, Rep);
// 			double dafric = gamma1p * FCoefAjusFric * f * entropiap * pow3(velabs) * DeltaTiempo
// 				/ (diamep * pow2(asonidop));
// 			entropia += dafric;
// 		}
		return entropy;
#ifdef usetry
	}

	catch(exception & N) {
		std::cout << "ERROR: TPipeMethod::ComputeEntropy " << FPipe->FPipeNumber << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
#endif
}

void TPipeMethod::ComputePipeEndCharacteristics(double dt) {
#ifdef usetry
	try {
#endif
		if(FPipe->getSpeed((unsigned int) 0) <= 0.) {
			FPipe->FLeftBC->setEntropy(InterpolateEntropy(nmLeft, dt));
		}
		FPipe->FLeftBC->setBeta(InterpolateCharacteristic(FPipe->FLeftBC->getEntropy(), 1, 0, dt));

		if(FPipe->getSpeed().tail(1)(0) >= 0) {
			FPipe->FRightBC->setEntropy(InterpolateEntropy(nmRight, dt));
		}
		FPipe->FRightBC->setLambda(InterpolateCharacteristic(FPipe->FRightBC->getEntropy(), -1, FPipe->getSpeed().cols() - 1,
								   dt));

#ifdef usetry
	} catch(exception & N) {
		std::cout << "ERROR: TPipeMethod::ComputePipeEndCharacteristics "
				  "in pipe number: " << FPipe->FPipeNumber << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
#endif
}

void TPipeMethod::Connect(const Pipe_ptr & pipe) {
	FPipe = pipe.get();
}


double TPipeMethod::getCourant() const {
	return FMaxCourant;
}

double TPipeMethod::getMaxTimeStep() {
// 	if(FPipe->FIsIntegrated == false) {
		ComputeMaxTimeStep();
		FPipe->FMaxTimeStep = FMaxTimeStep;
// 	}
	return FMaxTimeStep;
}

std::string TPipeMethod::getName() const {
	return FName;
}

void TPipeMethod::InterpolateBeta(nmPipeEnd pipe_end, double dt) {

}

double TPipeMethod::InterpolateCharacteristic(double entropia, int signo, int extremo, double DeltaTiempo) {
#ifdef usetry
	try {
#endif

		double dtdx = DeltaTiempo / FPipe->FXref;
		int ind = extremo;
		double characteristic = 0.;
		double velocidadp = 0.;
		double asonidop = 0.;

		if(DeltaTiempo < 1e-15) {
			characteristic = ComputeCharacteristic(velocidadp, asonidop, ind, 0., signo, entropia, DeltaTiempo);
		} else {

			dtdx = DeltaTiempo / FPipe->FXref;

			int ind1 = ind + signo;

			stCharOrigin CharOrigin(FPipe->FU0(0, ind), FPipe->FU0(1, ind), FPipe->FU0(2, ind), FPipe->FU0(0, ind1), FPipe->FU0(1,
									ind1), FPipe->FU0(2, ind1), FPipe->FGamma[ind], FPipe->FGamma(ind1),
									dtdx, signo);

			double dist = zbrent(CharOrigin, 0., 1., 1e-5);

			characteristic = ComputeCharacteristic(velocidadp, asonidop, ind, dist, signo, entropia, DeltaTiempo);
		}
		return characteristic / __cons::ARef;
#ifdef usetry
	} catch(exception & N) {
		std::cout << "ERROR: TPipeMethod::InterpolateCharacteristic " << FPipe->FPipeNumber << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
#endif
}

double TPipeMethod::InterpolateEntropy(nmPipeEnd pipe_end, double dt) {
#ifdef usetry
	try {
#endif

		int signo = 1;
		int extremo = 0;
		int indiceCC = 0;
		double entropy = 0.;

		if(pipe_end == nmRight) {   // PipeEnd Derecho
			signo = -1;
			extremo = FPipe->getArea().cols() - 1;
			indiceCC = 1;
		}
		double dtdx = dt / FPipe->FXref;
		int ind = extremo;
		double velocidadp = 0.;
		double asonidop = 0.;

		if(dt < 1e-15 || DoubEqZero(FPipe->getSpeed((unsigned int) extremo))) {

			entropy = ComputeEntropy(velocidadp, extremo, 0., signo, dt, indiceCC);

		} else {
			int ind1 = ind + signo;

			stPathOrigin PathOrigin(FPipe->FU0(0, ind), FPipe->FU0(1, ind), FPipe->FU0(0, ind1), FPipe->FU0(1, ind1), dtdx, signo);

			double dist = zbrent(PathOrigin, 0., 1., 1e-5);

			entropy = ComputeEntropy(velocidadp, ind, dist, signo, dt, indiceCC);
		}

		return entropy / __cons::ARef;
#ifdef usetry
	} catch(exception & N) {
		std::cout << "ERROR: TPipeMethod::InterpolateEntropy: " << std::endl;
		std::cout << "Tipo de error: " << N.what() << std::endl;
		throw Exception(N.what());
	}
#endif
}

void TPipeMethod::InterpolateLambda(nmPipeEnd pipe_end, double dt) {

}

void TBasicPipeMethod::setCourant(double C) {
	FMaxCourant = C;
}

void TPipeMethod::setPtTtU(double p_t, double T_t, double u) {
	auto T = T_t - u * u / 2. / __Gamma::Cp;
	auto p = p_t * pow(T / T_t, __Gamma::G9(__Gamma::G));
	setPTU(p, T, u);
}

void TPipeMethod::setPtTtU(const RowVector& p_t, const RowVector& T_t,
	const RowVector& u) {
	auto T = T_t - u * u / 2. / __Gamma::Cp;
	auto p = p_t * (T / T_t).pow(__Gamma::G9(__Gamma::G));
	setPTU(p, T, u);
}
