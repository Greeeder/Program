/**
 * @file TTCHTM2.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 22 de feb. de 2016
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
 * Turbocharger heat transfer model.
 */
#include "TTCHTM2.h"
#include <fstream>

TTC_HTM2::TTC_HTM2(TFluid *Oil) {
	// TODO Auto-generated constructor stub

	nnodes = 5;
	nboundy = 5;
	nparts = 3;

	NodeLabel.resize(nnodes);
	NodeLabel[0] = "T";
	NodeLabel[1] = "H1";
	NodeLabel[2] = "H2";
	NodeLabel[3] = "H3";
	NodeLabel[4] = "C";

	for(int i = 0; i < nnodes; i++) {
		MetalNode[NodeLabel[i]] = i;
	}

	BoundaryLabel.resize(nboundy);
	BoundaryLabel[0] = "GAS";
	BoundaryLabel[1] = "AIR";
	BoundaryLabel[2] = "OIL";
	BoundaryLabel[3] = "WAT";
	BoundaryLabel[4] = "AMB";

	for(int i = 0; i < nboundy; i++) {
		Boundary[BoundaryLabel[i]] = i;
	}

	PartLabel.resize(nparts);
	PartLabel[0] = "TURBINE";
	PartLabel[1] = "HOUSING";
	PartLabel[2] = "COMPRESSOR";

	for(int i = 0; i < nparts; i++) {
		Part[PartLabel[i]] = i;
	}

	dVector Y(1, 1);

	Fluid["GAS"] = new TFluidAir();
	Fluid["AIR"] = new TFluidAir();
	Fluid["AMB"] = new TFluidAir();
	Fluid["WAT"] = new TFluidWater();
	Fluid["OIL"] = Oil;

	D_in.setZero(nboundy);
	D_out.setZero(nboundy);

	ExternalDiameter.setZero(Part.size());
	ExternalLength.setZero(Part.size());
	ExternalArea.setZero(Part.size());

	T_oil = ArrayXd::Zero(3);
	T_cool = ArrayXd::Zero(2);

	FAlpha.resize(3);
	FAlpha << 0.58, 0.21, 0.21;

	VF_nodes.setZero(nnodes, nnodes);

	T1 = VectorXd::Zero(nnodes);
	T0 = VectorXd::Zero(nnodes);

	Qext_J = MatrixXd::Zero(nnodes, nboundy);
	Qext_W = MatrixXd::Zero(nnodes, nboundy);

	Qrad_J = MatrixXd::Zero(nnodes, nboundy);
	Qrad_W = MatrixXd::Zero(nnodes, nboundy);

	Kcond = MatrixXd::Zero(nnodes, nnodes);
	Krad = MatrixXd::Zero(nnodes, nnodes);

	Kconv = MatrixXd::Zero(nnodes, nboundy);
	Kradext = MatrixXd::Zero(nnodes, nboundy);

	ConvectiveConn = MatrixXi::Zero(nnodes, nboundy);
	ConvectiveConn(MetalNode["T"], Boundary["GAS"]) = 1;
	ConvectiveConn(MetalNode["C"], Boundary["AIR"]) = 1;
	ConvectiveConn(MetalNode["H2"], Boundary["OIL"]) = 1;
	ConvectiveConn(MetalNode["H2"], Boundary["WAT"]) = 1;
	ConvectiveConn.col(Boundary["AMB"]).setOnes();

	RadiativeConn = MatrixXi::Zero(nnodes, nboundy);
	RadiativeConn.col(Boundary["AMB"]).setOnes();

	Kfull1 = MatrixXd::Identity(nnodes + nboundy, nnodes + nboundy);
	Kfull0 = MatrixXd::Identity(nnodes + nboundy, nnodes + nboundy);

	Tb = VectorXd::Zero(nnodes + nboundy);

	FFitTurbExtHeat = 1.4;

	FInitialTemp = false;

	FImposedExtHeatCoefficient = false;
}

TTC_HTM2::~TTC_HTM2() {
	// TODO Auto-generated destructor stub
}

void TTC_HTM2::SolveNodeTemperatures(double dt, int mode) {

	VectorXd b;

	VectorXd Qe = (Qext_J) * VectorXd::Ones(nboundy);
	VectorXd Qini = Kcond * T1;

	if(mode == 0) {
		b = (Qext_J + Qrad_J) * VectorXd::Ones(nboundy);

		T1 = Kcond.ldlt().solve(Qe);
		T1 = (Kcond + Krad).ldlt().solve(b);
	} else {

		b = ((Kcond + Krad) / 2. - C / dt) * T0 + (Qext_J + Qrad_J) * VectorXd::Ones(nboundy);

		T1 = ((Kcond + Krad) / 2. + C / dt).ldlt().solve(b);
	}

	Qext_J.setZero();
	Qrad_J.setZero();

}

void TTC_HTM2::SolveExplicit(double dt, int mode) {

	MatrixXd K = Kcond;
	MatrixXd Q = Qext_J;
	T0 = T1;
	if(FExternalRadiation) {
		K += - Krad + Krad.rowwise().sum().asDiagonal().toDenseMatrix();
		Q += Qrad_J;
	}

	T1 = T0 + Cinv * (dt * (-K) * T0 + (Q) * VectorXd::Ones(nboundy));

	Qext_J.setZero();
	Qrad_J.setZero();

}

void TTC_HTM2::SolveImplicit(double dt, int mode) {

	Kfull1.block(nboundy, 0, nnodes, nboundy) = -(Kconv + Kradext);

	if (!Kfull1.block(nboundy, 0, nnodes, nboundy).isZero(1e-9)){

		Kfull1.block(nboundy, nnodes, nnodes, nboundy) = -Krad;
		Kfull1.block(nboundy, nnodes, nnodes, nboundy) -= (Kfull1.block(nboundy, 0, nnodes,
			nboundy + nnodes)).rowwise().sum().asDiagonal().toDenseMatrix();
		Kfull1.block(nboundy, nnodes, nnodes, nboundy) += Kcond;
		VectorXd T;
		if (mode == 0) {
			T = Kfull1.lu().solve(Tb);
		}
		else {
			Kfull1.block(nboundy, nnodes, nnodes, nboundy) /= 2.;
			Kfull0.block(nboundy, 0, nnodes, nboundy + nnodes) = -Kfull1.block(nboundy, 0, nnodes, nboundy + nnodes);
			Kfull1.block(nboundy, 0, nnodes, nboundy).setZero();
			Kfull1.block(nboundy, nnodes, nnodes, nboundy) += C / dt;
			Kfull0.block(nboundy, nnodes, nnodes, nboundy) += C / dt;

			Tb.tail(nnodes) = T0;

			T = Kfull1.lu().solve(Kfull0 * Tb);
		}

		T1 = T.tail(nnodes);
	}
	Qext_J.setZero();
	Qrad_J.setZero();

}

void TTC_HTM2::Read_HTMXML(xml_node node_ht) {

	double tmp = 0.;
	dVector vtmp;
	int np = 0;
	int water = 0;

	xml_node node_set = GetNodeChild(node_ht, "Htr_Settings");

	FIsWaterCooled = GetAttributeAsBool(node_set, "WaterCooled");

	//if (FIsWaterCooled) {
	//	FMatrix_Conn[FNode.H2][FNode.Water] = true;
	//}
	//else {
	//	FMatrix_Conn[FNode.H3][FNode.AirOut] = true;
	//}

	xml_node node_cond = GetNodeChild(node_ht, "Htr_Conduction");

	MetalConductances = MatrixXd::Zero(nnodes, nnodes);

	MetalConductances(MetalNode["T"], MetalNode["H1"]) = GetAttributeAsDouble(node_cond, "T_H1");
	MetalConductances(MetalNode["H1"], MetalNode["H2"]) = GetAttributeAsDouble(node_cond, "H1_H2");
	MetalConductances(MetalNode["H2"], MetalNode["H3"]) = GetAttributeAsDouble(node_cond, "H2_H3");
	MetalConductances(MetalNode["H3"], MetalNode["C"]) = GetAttributeAsDouble(node_cond, "H3_C");

	MatrixXd M_t = MetalConductances.transpose();

	MetalConductances += M_t;

	xml_node node_cap = GetNodeChild(node_ht, "Htr_Capacity");

	MetalCapacitances = MatrixXd::Zero(nnodes, nnodes);

	MetalCapacitances(MetalNode["T"], MetalNode["T"]) = GetAttributeAsDouble(node_cap, "T");
	MetalCapacitances(MetalNode["H1"], MetalNode["H1"]) = GetAttributeAsDouble(node_cap, "H1");
	MetalCapacitances(MetalNode["H2"], MetalNode["H2"]) = GetAttributeAsDouble(node_cap, "H2");
	MetalCapacitances(MetalNode["H3"], MetalNode["H3"]) = GetAttributeAsDouble(node_cap, "H3");
	MetalCapacitances(MetalNode["C"], MetalNode["C"]) = GetAttributeAsDouble(node_cap, "C");

	xml_node node_iconv = GetNodeChild(node_ht, "Htr_Convection");

	/** Correlation GAS - T * */

	ConvectiveNusselt["T"]["GAS"] = new TCConvCorr();

	xml_node node_gast = GetNodeChild(node_iconv, "Cnv_Gas_T");

	ConvectiveNusselt["T"]["GAS"]->K.push_back(GetAttributeAsDouble(node_gast, "k1"));
	ConvectiveNusselt["T"]["GAS"]->K.push_back(GetAttributeAsDouble(node_gast, "k2"));
	ConvectiveNusselt["T"]["GAS"]->K.push_back(GetAttributeAsDouble(node_gast, "k3"));
	ConvectiveNusselt["T"]["GAS"]->K.push_back(GetAttributeAsDouble(node_gast, "k4"));
	ConvectiveNusselt["T"]["GAS"]->K.push_back(1.467);

	ConvectiveNusselt["GAS"]["T"] = ConvectiveNusselt["T"]["GAS"];

	if(FIsWaterCooled) {
		/** Correlation Water - H2 * */

		ConvectiveNusselt["H2"]["WAT"] = new TCConvCorr();

		xml_node node_wath2 = GetNodeChild(node_iconv, "Cnv_Wat_H2");

		ConvectiveNusselt["H2"]["WAT"]->K.push_back(GetAttributeAsDouble(node_wath2, "k1"));
		ConvectiveNusselt["H2"]["WAT"]->K.push_back(GetAttributeAsDouble(node_wath2, "k2"));
		ConvectiveNusselt["H2"]["WAT"]->K.push_back(GetAttributeAsDouble(node_wath2, "k3"));

		ConvectiveNusselt["WAT"]["H2"] = ConvectiveNusselt["H2"]["WAT"];

	}

	/** Correlation Oil - H1 * */

	ConvectiveNusselt["H1"]["OIL"] = new TCConvCorr();

	xml_node node_oilh1 = GetNodeChild(node_iconv, "Cnv_Oil_H1");

	ConvectiveNusselt["H1"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh1, "k1"));
	ConvectiveNusselt["H1"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh1, "k2"));
	ConvectiveNusselt["H1"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh1, "k3"));
	ConvectiveNusselt["H1"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh1, "k4"));
	ConvectiveNusselt["H1"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh1, "k5"));


	///** Corelation Oil - H2 * */

	ConvectiveNusselt["OIL"]["H2"] = new TCConvCorr();

	xml_node node_oilh2c = GetNodeChild(node_iconv, "Cnv_Oil_H2_C");

	ConvectiveNusselt["OIL"]["H2"]->K.push_back(GetAttributeAsDouble(node_oilh2c, "k1"));
	ConvectiveNusselt["OIL"]["H2"]->K.push_back(GetAttributeAsDouble(node_oilh2c, "k2"));
	ConvectiveNusselt["OIL"]["H2"]->K.push_back(GetAttributeAsDouble(node_oilh2c, "k3"));
	ConvectiveNusselt["OIL"]["H2"]->K.push_back(GetAttributeAsDouble(node_oilh2c, "k4"));

	ConvectiveNusselt["H2"]["OIL"] = new TCConvCorr();

	xml_node node_oilh2h = GetNodeChild(node_iconv, "Cnv_Oil_H2_H");

	ConvectiveNusselt["H2"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh2h, "k1"));
	ConvectiveNusselt["H2"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh2h, "k2"));
	ConvectiveNusselt["H2"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh2h, "k3"));
	ConvectiveNusselt["H2"]["OIL"]->K.push_back(GetAttributeAsDouble(node_oilh2h, "k4"));

	///** Corelation Air - C * */

	ConvectiveNusselt["AIR"]["C"] = new TCConvCorr();

	xml_node node_aircc = GetNodeChild(node_iconv, "Cnv_Air_C_C");

	ConvectiveNusselt["AIR"]["C"]->K.push_back(GetAttributeAsDouble(node_aircc, "k1"));
	ConvectiveNusselt["AIR"]["C"]->K.push_back(GetAttributeAsDouble(node_aircc, "k2"));
	ConvectiveNusselt["AIR"]["C"]->K.push_back(GetAttributeAsDouble(node_aircc, "k3"));

	///** Corelation Air - C * */

	ConvectiveNusselt["C"]["AIR"] = new TCConvCorr();

	xml_node node_airch = GetNodeChild(node_iconv, "Cnv_Air_C_H");

	ConvectiveNusselt["C"]["AIR"]->K.push_back(GetAttributeAsDouble(node_airch, "k1"));
	ConvectiveNusselt["C"]["AIR"]->K.push_back(GetAttributeAsDouble(node_airch, "k2"));
	ConvectiveNusselt["C"]["AIR"]->K.push_back(GetAttributeAsDouble(node_airch, "k3"));

	xml_node node_rad = GetNodeChild(node_ht, "Htr_Radiation");

	FExternalRadiation = GetAttributeAsBool(node_rad, "Radiation");

	if(FExternalRadiation) {
		FFitRadiation = GetAttributeAsDouble(node_rad, "FitCoefficient");

		Emissivity["TUR"] = GetAttributeAsDouble(node_rad, "TurbineEmissivity");
		Emissivity["HOU"] = GetAttributeAsDouble(node_rad, "HousingEmissivity");
		Emissivity["COM"] = GetAttributeAsDouble(node_rad, "CompressorEmissivity");

		FTurbineShielded = GetAttributeAsBool(node_rad, "Shield");
		if(FTurbineShielded)
			Emissivity["SHI"] = GetAttributeAsDouble(node_rad, "ShieldEmissivity");
	}

	xml_node node_econv = GetNodeChild(node_ht, "Htr_ConvectionExt");

	FExternalConvection = GetAttributeAsBool(node_econv, "ExtConvection");

	if (FExternalConvection){
		FFitConvection = GetAttributeAsDouble(node_econv, "FitCoefficient");

		xml_node node_ehtc = GetNodeChild(node_econv, "Ecv_HeatCoefficients");
		if (node_ehtc){
			FImposedExtHeatCoefficient = true;
			FTExtHeatCoefficient = GetAttributeAsDouble(node_ehtc, "Turbine");
			FHExtHeatCoefficient = GetAttributeAsDouble(node_ehtc, "Housing");
			FCExtHeatCoefficient = GetAttributeAsDouble(node_ehtc, "Compressor");
		}
		else{
			if (node_econv.attribute("ExtVelocity")){
				FExtVelocity = GetAttributeAsDouble(node_econv, "ExtVelocity");
				FTurbExtVelocity = FExtVelocity;
				FHousExtVelocity = FExtVelocity;
				FCompExtVelocity = FExtVelocity;
			}
			else{
				xml_node node_extvel = GetNodeChild(node_econv, "Ecv_Velocities");
				FTurbExtVelocity = GetAttributeAsDouble(node_extvel, "Turbine");
				FHousExtVelocity = GetAttributeAsDouble(node_extvel, "Housing");
				FCompExtVelocity = GetAttributeAsDouble(node_extvel, "Compressor");
			}
		}
	}
	
	if(node_econv.attribute("TurbineHeatMultiplier")) {
		FFitTurbExtHeat = GetAttributeAsDouble(node_econv, "TurbineHeatMultiplier");
	}

	xml_node node_temp = GetNodeChild(node_ht, "Htr_Temperatures");

	if (node_temp){
		Tinit.resize(nnodes);
		FInitialTemp = true;
		Tinit(0) = GetAttributeAsDouble(node_temp, "T");
		Tinit(1) = GetAttributeAsDouble(node_temp, "H1");
		Tinit(2) = GetAttributeAsDouble(node_temp, "H2");
		Tinit(3) = GetAttributeAsDouble(node_temp, "H3");
		Tinit(3) = GetAttributeAsDouble(node_temp, "C");
	}

	setKCondMatrix();
	setCapMatrix();

}

void TTC_HTM2::setHeatFromGas(double dt, double T, double m) {

	int NodeG = Boundary["GAS"];
	int NodeT = MetalNode["T"];

	setKFromGas(T, m);

	Qext_W(NodeT, NodeG) = Kconv(NodeT, NodeG) * (T - T1(NodeT));

	Qext_J(NodeT, NodeG) = Qext_J(NodeT, NodeG) + Qext_W(NodeT, NodeG) * dt;

}

void TTC_HTM2::setKFromGas(double T, double m) {

	int NodeG = Boundary["GAS"];
	int NodeT = MetalNode["T"];

	Tb(NodeG) = T;

	dVector Input;

	TFluid* Gas = Fluid["GAS"];

	double mu = Gas->FunVisc(T);
	double mu_w = Gas->FunVisc(T1(NodeT));
	double k = Gas->Funk(T);
	// Re - Input(1)
	Input.push_back(fabs(m) / (__cons::Pi_4 * mu * D_in(NodeG)));
	// Pr - Input(2)
	Input.push_back(Gas->FunPr(T));
	// mu / mu_w - Input(3)
	Input.push_back(mu / mu_w);
	// Eff_max - Input(4)
	Input.push_back(TurbMaxEfficiency);

	Kconv(NodeT, NodeG) = ConvectiveNusselt["T"]["GAS"]->Nusselt(Input) * k *
						  __geom::Circle_area(ExternalDiameter(Part["TURBINE"])) / D_in[NodeG];

}

void TTC_HTM2::setHeatFromAir(double dt, double T, double m) {

	int NodeB = Boundary["AIR"];
	int NodeM = MetalNode["C"];

	setKFromAir(T, m);

	Qext_W(NodeM, NodeB) = Kconv(NodeM, NodeB) * (T - T1(NodeM));

	Qext_J(NodeM, NodeB) = Qext_J(NodeM, NodeB) + Qext_W(NodeM, NodeB) * dt;
}

void TTC_HTM2::setKFromAir(double T, double m) {

	int NodeB = Boundary["AIR"];
	int NodeM = MetalNode["C"];

	dVector Input;

	TFluid* Air = Fluid["AIR"];

	Tb(NodeB) = T;

	double mu = Air->FunVisc(T);
	double k = Air->Funk(T);
	// Re - Input(1)
	Input.push_back(fabs(m) / (__cons::Pi_4 * mu * D_in(NodeB)));
	// Pr - Input(2)
	Input.push_back(Air->FunPr(T));

	if(T > T1(NodeM)) {
		Kconv(NodeM, NodeB) = ConvectiveNusselt["AIR"]["C"]->Nusselt(Input) * k *
							  __cons::Pi * ExternalDiameter(Part["COMPRESSOR"]);
	} else {
		Kconv(NodeM, NodeB) = ConvectiveNusselt["C"]["AIR"]->Nusselt(Input) * k *
							  __cons::Pi * ExternalDiameter(Part["COMPRESSOR"]);
	}
}

void TTC_HTM2::setHeatFromOil(double dt, double T, double m, double mech_pow) {

	int NodeB = Boundary["OIL"];
	int NodeM = MetalNode["H2"];

	setKFromOil(T, m, mech_pow);

	Qext_W(NodeM, NodeB) = Kconv(NodeM, NodeB) * (T_oilH2 - T1(NodeM));

	Qext_J(NodeM, NodeB) = Qext_J(NodeM, NodeB) + Qext_W(NodeM, NodeB) * dt;
}

void TTC_HTM2::setKFromOil(double T, double m, double mech_pow) {

	int NodeB = Boundary["OIL"];

	dVector Input;

	TFluid* Oil = Fluid["OIL"];

	// OIL from/to H2

	int NodeM = MetalNode["H2"];

	T_oil(0) = T;

	T_oil(1) = T_oil(0) + mech_pow / m / Oil->FunCp(T_oilMech);

	T_oilMech = (T_oil(0) + T_oil(1)) / 2.;

	T_oil(2) = T_oil(1) - Qext_W(NodeM, NodeB) / m / Oil->FunCp(T_oilH2);
	if(Qext_W(NodeM, NodeB) < 0 && T_oil(2) > T1(NodeM))
		T_oil(2) = T1(NodeM);
	if(Qext_W(NodeM, NodeB) > 0 && T_oil(2) < T1(NodeM))
		T_oil(2) = T1(NodeM);

	T_oilH2 = (T_oil(1) + T_oil(2)) / 2;

	Tb(NodeB) = T_oilH2;

	double mu = Oil->FunVisc(T_oilH2);
	double mu_w = Oil->FunVisc(T1(NodeM));
	double k = Oil->Funk(T_oilH2);
	// Re - Input(1)
	Input.push_back(fabs(m) / (__cons::Pi_4 * mu * D_in(NodeB)));
	// Pr - Input(2)
	Input.push_back(Oil->FunPr(T_oilH2));
	// mu / mu_w - Input(3)
	Input.push_back(mu / mu_w);

	if(T > T1(NodeM)) {
		Kconv(NodeM, NodeB) = ConvectiveNusselt["OIL"]["H2"]->Nusselt(Input) * k *
							  __cons::Pi * ExternalLength(Part["HOUSING"]);
	} else {
		Kconv(NodeM, NodeB) = ConvectiveNusselt["H2"]["OIL"]->Nusselt(Input) * k *
							  __cons::Pi * ExternalLength(Part["HOUSING"]);
	}
}

void TTC_HTM2::setHeatFromCoolant(double dt, double T, double m) {

	if(FIsWaterCooled) {

		int NodeB = Boundary["WAT"];
		int NodeM = MetalNode["H2"];

		setKFromCoolant(T, m);

		Qext_W(NodeM, NodeB) = Kconv(NodeM, NodeB) * (T - T1(NodeM));

		T_cool(0) = T;
		T_cool(1) = T - Qext_W(NodeM, NodeB) / Fluid["WAT"]->FunCp((T_cool(0) + T_cool(1)) / 2) / m;

		Qext_J(NodeM, NodeB) = Qext_J(NodeM, NodeB) + Qext_W(NodeM, NodeB) * dt;
	}

}

void TTC_HTM2::setKFromCoolant(double T, double m) {

	if(FIsWaterCooled) {
		int NodeB = Boundary["WAT"];
		int NodeM = MetalNode["H2"];

		Tb(NodeB) = T;

		dVector Input;

		TFluid* Coolant = Fluid["WAT"];

		double mu = Coolant->FunVisc(T);
		double k = Coolant->Funk(T);
		// Re - Input(1)
		Input.push_back(fabs(m) / (__cons::Pi_4 * mu * D_in(NodeB)));
		// Pr - Input(2)
		Input.push_back(Coolant->FunPr(T));

		Kconv(NodeM, NodeB) = ConvectiveNusselt["WAT"]["H2"]->Nusselt(Input) * k *
							  __cons::Pi * ExternalLength(Part["HOUSING"]);
	}
}

void TTC_HTM2::setHeatFromAmbient(double dt, double T, double p) {

	int NodeB = Boundary["AMB"];

	if (FExternalConvection) {

		setKFromAmbient(T, p);

		for(int NodeM = 0; NodeM < nnodes; NodeM++) {
			Qext_W(NodeM, NodeB) = Kconv(NodeM, NodeB) * (T - T1(NodeM));

			Qext_J(NodeM, NodeB) = Qext_J(NodeM, NodeB) + Qext_W(NodeM, NodeB) * dt;
		}

	}
	else{
		Qext_W.col(NodeB).setZero();
	}

}

void TTC_HTM2::setKFromAmbient(double T, double p) {

	int NodeB = Boundary["AMB"];
	Tb(NodeB) = T;

	if(FExternalConvection) {
		if (FImposedExtHeatCoefficient){
			int NodeM = MetalNode["T"];
			int TC_part = Part["TURBINA"];

			Kconv(NodeM, NodeB) = ExternalArea(TC_part) * FTExtHeatCoefficient * FFitConvection;

			NodeM = MetalNode["H1"];
			TC_part = Part["HOUSING"];

			Kconv(NodeM, NodeB) = FAlpha(0) * ExternalArea(TC_part) * FHExtHeatCoefficient * FFitConvection;

			NodeM = MetalNode["H2"];

			Kconv(NodeM, NodeB) = FAlpha(1) * ExternalArea(TC_part) * FHExtHeatCoefficient * FFitConvection;

			NodeM = MetalNode["H3"];

			Kconv(NodeM, NodeB) = FAlpha(2) * ExternalArea(TC_part) * FHExtHeatCoefficient * FFitConvection;

			NodeM = MetalNode["C"];
			TC_part = Part["COMPRESSOR"];

			Kconv(NodeM, NodeB) = FAlpha(0) * ExternalArea(TC_part) * FHExtHeatCoefficient * FFitConvection;
		}
		else{

			double Tprop = 0., Beta = 0., rho = 0., mu = 0., visc = 0., Gr = 0., Re = 0., Pr = 0., Nu = 0.;
			double g = 9.8; // [m/s**2]

			TFluid* Ambient = Fluid["AMB"];

			// Turbine to ambient

			int NodeM = MetalNode["T"];
			int TC_part = Part["TURBINA"];

			Tprop = (T1(NodeM) + T) / 2;

			Beta = 1 / Tprop;
			rho = p / Ambient->FunR() / Tprop;
			mu = Ambient->FunVisc(Tprop);
			visc = mu / rho;
			Pr = Ambient->FunPr(Tprop);

			Gr = g * Beta * fabs(T1(NodeM) - T) * pow3(ExternalDiameter(TC_part) * FFitTurbExtHeat) / pow2(visc);

			Re = rho * FTurbExtVelocity * FFitTurbExtHeat * ExternalDiameter(TC_part) / mu;

			if (Gr > 10 * pow2(Re)) {
				Nu = NusseltFreeConv(Gr, Pr);
			}
			else if (Gr < 0.1 * pow2(Re)) {
				Nu = NusseltForcConv(Re, Pr);
			}
			else {
				Nu = cbrt(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)));
			}
			Kconv(NodeM, NodeB) = __cons::Pi * ExternalLength(TC_part) * Nu * Ambient->Funk(Tprop) * FFitConvection;

			// H1 to ambient

			NodeM = MetalNode["H1"];
			TC_part = Part["HOUSING"];

			Tprop = (T1(NodeM) + T) / 2;

			Beta = 1 / Tprop;
			rho = p / Ambient->FunR() / Tprop;
			mu = Ambient->FunVisc(Tprop);
			visc = mu / rho;
			Pr = Ambient->FunPr(Tprop);

			Gr = g * Beta * fabs(T1(NodeM) - T) * pow3(ExternalDiameter(TC_part)) / pow2(visc);

			Re = rho * FHousExtVelocity * ExternalDiameter(TC_part) / mu;

			if (Gr > 10 * pow2(Re)) {
				Nu = NusseltFreeConv(Gr, Pr);
			}
			else if (Gr < 0.1 * pow2(Re)) {
				Nu = NusseltForcConv(Re, Pr);
			}
			else {
				Nu = cbrt(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)));
			}
			Kconv(NodeM, NodeB) = __cons::Pi * ExternalLength(TC_part) * FAlpha(0) * Nu * Ambient->Funk(Tprop) * FFitConvection;

			// H2 to ambient

			NodeM = MetalNode["H2"];
			TC_part = Part["HOUSING"];

			Tprop = (T1(NodeM) + T) / 2;

			Beta = 1 / Tprop;
			rho = p / Ambient->FunR() / Tprop;
			mu = Ambient->FunVisc(Tprop);
			visc = mu / rho;
			Pr = Ambient->FunPr(Tprop);

			Gr = g * Beta * fabs(T1(NodeM) - T) * pow3(ExternalDiameter(TC_part)) / pow2(visc);

			Re = rho * FHousExtVelocity * ExternalDiameter(TC_part) / mu;

			if (Gr > 10 * pow2(Re)) {
				Nu = NusseltFreeConv(Gr, Pr);
			}
			else if (Gr < 0.1 * pow2(Re)) {
				Nu = NusseltForcConv(Re, Pr);
			}
			else {
				Nu = cbrt(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)));
			}
			Kconv(NodeM, NodeB) = __cons::Pi * ExternalLength(TC_part) * FAlpha(1) * Nu * Ambient->Funk(Tprop) * FFitConvection;

			// H3 to ambient

			NodeM = MetalNode["H3"];
			TC_part = Part["HOUSING"];

			Tprop = (T1(NodeM) + T) / 2;

			Beta = 1 / Tprop;
			rho = p / Ambient->FunR() / Tprop;
			mu = Ambient->FunVisc(Tprop);
			visc = mu / rho;
			Pr = Ambient->FunPr(Tprop);

			Gr = g * Beta * fabs(T1(NodeM) - T) * pow3(ExternalDiameter(TC_part)) / pow2(visc);

			Re = rho * FHousExtVelocity * ExternalDiameter(TC_part) / mu;

			if (Gr > 10 * pow2(Re)) {
				Nu = NusseltFreeConv(Gr, Pr);
			}
			else if (Gr < 0.1 * pow2(Re)) {
				Nu = NusseltForcConv(Re, Pr);
			}
			else {
				Nu = cbrt(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)));
			}
			Kconv(NodeM, NodeB) = __cons::Pi * ExternalLength(TC_part) * FAlpha(2) * Nu * Ambient->Funk(Tprop) * FFitConvection;

			// C to ambient

			NodeM = MetalNode["C"];
			TC_part = Part["COMPRESSOR"];

			Tprop = (T1(NodeM) + T) / 2;

			Beta = 1 / Tprop;
			rho = p / Ambient->FunR() / Tprop;
			mu = Ambient->FunVisc(Tprop);
			visc = mu / rho;
			Pr = Ambient->FunPr(Tprop);

			Gr = g * Beta * fabs(T1(NodeM) - T) * pow3(ExternalDiameter(TC_part)) / pow2(visc);

			Re = rho * FCompExtVelocity * ExternalDiameter(TC_part) / mu;

			if (Gr > 10 * pow2(Re)) {
				Nu = NusseltFreeConv(Gr, Pr);
			}
			else if (Gr < 0.1 * pow2(Re)) {
				Nu = NusseltForcConv(Re, Pr);
			}
			else {
				Nu = cbrt(pow3(NusseltFreeConv(Gr, Pr)) + pow3(NusseltForcConv(Re, Pr)));
			}
			Kconv(NodeM, NodeB) = __cons::Pi * ExternalLength(TC_part) * Nu * Ambient->Funk(Tprop) * FFitConvection;
		}
	}

}

double TTC_HTM2::NusseltFreeConv(double Gr, double Pr) {

	double Ra = Gr * Pr;

	return pow2(0.6 + 0.387 * cbrt(cbrt(Ra)) / pow(1 + 0.559 / pow(Pr, 9. / 16.), 8. / 27.));

}

double TTC_HTM2::NusseltForcConv(double Re, double Pr) {

	return 0.3 + 0.62 * sqrt(Re) * cbrt(Pr) / sqrt(sqrt(1 + 0.4 / pow2(cbrt(Pr)))) * pow(1 + pow(Re / 282000, 5. / 8.),
			4 / 5);

}

void TTC_HTM2::setViewFactorBetweenNodes(double Dt, double Dc, double Dh, double L) {

	double A1 = FAlpha(0);
	double A2 = FAlpha(1);
	double A3 = FAlpha(2);

	double Dt2 = pow2(Dt);
	double Dc2 = pow2(Dc);
	double Dh2 = pow2(Dh);

	double DtDh = Dt2 - Dh2;
	double DcDh = Dc2 - Dh2;
	double LDh4 = L * Dh * 4;

	int _c = MetalNode["C"];
	int _t = MetalNode["T"];
	int _h1 = MetalNode["H1"];
	int _h2 = MetalNode["H2"];
	int _h3 = MetalNode["H3"];

	VF_nodes(_c, _t) = ViewFactorDiskDisk(Dc / 2, Dt / 2, Dh / 2, L);
	VF_nodes(_t, _c) = VF_nodes(_c, _t) * DcDh / DtDh;
	VF_nodes(_h1, _t) = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L * A1);
	VF_nodes(_t, _h1) = VF_nodes(_h1, _t) * LDh4 * A1 / DtDh;
	VF_nodes(_t, _h2) = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L * (A1 + A2)) * (A1 + A2) * LDh4 / DtDh - VF_nodes(_t, _h1);
	VF_nodes(_h2, _t) = VF_nodes(_t, _h2) * DtDh / A2 / LDh4;
	VF_nodes(_t, _h3) = ViewFactorDiskCylinder(Dh / 2, Dt / 2, L) * LDh4 / DtDh - VF_nodes(_t, _h1) - VF_nodes(_t, _h2);
	VF_nodes(_h3, _t) = VF_nodes(_t, _h3) * DtDh / LDh4 / A3;
	VF_nodes(_h3, _c) = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L * A3);
	VF_nodes(_c, _h3) = VF_nodes(_h3, _c) * LDh4 * A3 / DcDh;
	VF_nodes(_c, _h2) = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L * (A3 + A2)) * (A3 + A2) * LDh4 / DcDh - VF_nodes(_c, _h3);
	VF_nodes(_h2, _c) = VF_nodes(_c, _h2) * DcDh / LDh4 / A2;
	VF_nodes(_c, _h1) = ViewFactorDiskCylinder(Dh / 2, Dc / 2, L) * LDh4 / DcDh - VF_nodes(_c, _h3) - VF_nodes(_c, _h2);
	VF_nodes(_h1, _c) = VF_nodes(_c, _h1) * DcDh / LDh4 / A1;

	if(FTurbineShielded) {
		VF_nodes(_t, _c) = VF_nodes(_t, _c) * Emissivity["SHI"] /
						   (2 * (Emissivity["SHI"] + VF_nodes(_t, _c) * (1 - Emissivity["SHI"])));
		VF_nodes(_c, _t) = VF_nodes(_t, _c) * DtDh / DcDh;

		VF_nodes(_t, _h1) = VF_nodes(_t, _h1) * Emissivity["SHI"] /
							(2 * (Emissivity["SHI"] + VF_nodes(_t, _h1) * (1 - Emissivity["SHI"])));
		VF_nodes(_h1, _t) = VF_nodes(_t, _h1) / (LDh4 * A1 / DtDh);

		VF_nodes(_t, _h2) = VF_nodes(_t, _h2) * Emissivity["SHI"] /
							(2 * (Emissivity["SHI"] + VF_nodes(_t, _h2) * (1 - Emissivity["SHI"])));
		VF_nodes(_h2, _t) = VF_nodes(_t, _h2) / (LDh4 * A2 / DtDh);

		VF_nodes(_t, _h3) = VF_nodes(_t, _h3) * Emissivity["SHI"] /
							(2 * (Emissivity["SHI"] + VF_nodes(_t, _h3) * (1 - Emissivity["SHI"])));
		VF_nodes(_h3, _t) = VF_nodes(_t, _h3) / (LDh4 * A3 / DtDh);

	}

}

void TTC_HTM2::setViewFactorToBoundaries(double Dt, double Dc, double Dh, double L) {

	int _a = Boundary["AMB"];
	int _c = MetalNode["C"];
	int _t = MetalNode["T"];
	int _h1 = MetalNode["H1"];
	int _h2 = MetalNode["H2"];
	int _h3 = MetalNode["H3"];

	VF_toBo = MatrixXd::Zero(nnodes, nboundy);

	VF_toBo.col(_a) = VectorXd::Ones(nnodes) - VF_nodes.rowwise().sum();

	if(FTurbineShielded) {
		VF_toBo(_t, _a) = VF_nodes(_t, _a) * Emissivity["SHI"] /
						  (2 * (Emissivity["SHI"] + VF_nodes(_t, _a) * (1 - Emissivity["SHI"])));
	}

}

double TTC_HTM2::ViewFactorDiskDisk(double r1, double r2, double rc, double h) {

	double R1 = r1 / h;
	double R2 = r2 / h;
	double Rc = rc / h;
	double A = pow2(R1) - pow2(Rc);
	double B = pow2(R2) - pow2(Rc);
	double C = R1 + R2;
	double D = R2 - R1;
	double Y = sqrt(A) + sqrt(B);

	double VF1 = A / 2 * acos(Rc / R2) + B / 2 * acos(Rc / R1) + 2 * Rc * (atan(Y) - atan(sqrt(A)) - atan(sqrt(B)));
	double VF2 = sqrt((1 + pow2(C)) * (1 + pow2(D))) * atan(sqrt((1 + pow2(C)) * (pow2(Y) - pow2(D)) / (1 + pow2(D)) /
				 (pow2(C) - pow2(Y))));
	double VF3 = sqrt((1 + pow2(R1 + Rc)) * (1 + pow2(R1 - Rc))) * atan(sqrt((1 + pow2(R1 + Rc)) * (R1 - Rc) / (1 + pow2(
					 R1 - Rc)) / (R1 + Rc)));
	double VF4 = sqrt((1 + pow2(R2 + Rc)) * (1 + pow2(R2 - Rc))) * atan(sqrt((1 + pow2(R2 + Rc)) * (R2 - Rc) / (1 + pow2(
					 R2 - Rc)) / (R2 + Rc)));

	return (VF1 - VF2 + VF3 + VF4) / A / __cons::Pi;

}

double TTC_HTM2::ViewFactorDiskCylinder(double r1, double r2, double h) {

	double R = r1 / r2;
	double H = h / r2;
	double A = pow2(H) + pow2(R) - 1;
	double B = pow2(H) - pow2(R) + 1;

	return B / 8 / R / H + 1 / (__cons::Pi_x_2) * (acos(A / B) - 1 / (2 * H) * sqrt(pow2(A + 2) / pow2(R) - 4) * acos(
				A * R / B) - A / 2 / R / H * asin(R));

}

void TTC_HTM2::setRadiatioFromAmbient(double dt, double T) {

	if(FExternalRadiation) {
		int b = Boundary["AMB"];

		setKRadiatioFromAmbient(T);

		for(int m = 0; m < nnodes; m++) {
			Qrad_W(m, b) = Kradext(m, b) * (T - T1(m));
			Qrad_J(m, b) = Qrad_J(m, b) + Qrad_W(m, b) * dt;
		}

	}
}

void TTC_HTM2::setKRadiatioFromAmbient(double T) {

	if(FExternalRadiation) {

		int b = Boundary["AMB"];

		int m = MetalNode["T"];
		int p = Part["TURBINE"];
		Kradext(m, b) = KRadiationExternal(ExternalLength(p), ExternalDiameter(p) * FFitTurbExtHeat, Emissivity["TUR"],
			VF_nodes(m, b), T1(m), T) * FFitRadiation;

		m = MetalNode["C"];
		p = Part["COMPRESSOR"];
		Kradext(m, b) = KRadiationExternal(ExternalLength(p), ExternalDiameter(p) * FFitTurbExtHeat, Emissivity["COM"],
			VF_nodes(m, b), T1(m), T) * FFitRadiation;

		m = MetalNode["H1"];
		p = Part["HOUSING"];
		Kradext(m, b) = KRadiationExternalHous(ExternalArea(p) * FAlpha(0), Emissivity["HOU"], VF_nodes(m, b),
			T1(m), T) * FFitRadiation;

		m = MetalNode["H2"];
		p = Part["HOUSING"];
		Kradext(m, b) = KRadiationExternalHous(ExternalArea(p) * FAlpha(1), Emissivity["HOU"], VF_nodes(m, b),
			T1(m), T) * FFitRadiation;

		m = MetalNode["H3"];
		p = Part["HOUSING"];
		Kradext(m, b) = KRadiationExternalHous(ExternalArea(p) * FAlpha(2), Emissivity["HOU"], VF_nodes(m, b),
			T1(m), T) * FFitRadiation;

	}
	else{
		Kradext.setZero();
	}
}

double TTC_HTM2::KRadiationExternal(double L, double D, double e, double F, double T, double Tamb) {

	return (D * (1 + F / (e * (1 - F) + F)) + 4 * L) * __cons::Pi_4 * D * e * __cons::Sigma * (pow2(T) + pow2(Tamb)) *
		   (T + Tamb);

}

double TTC_HTM2::KRadiationExternalHous(double A, double e, double F, double T, double Tamb) {

	return e * F / (e * (1 - F) + F) * A * __cons::Sigma * (pow2(T) + pow2(Tamb)) * (T + Tamb);

}

void TTC_HTM2::setKCondMatrix() {

	Kcond = -MetalConductances + (MetalConductances.rowwise().sum()).asDiagonal().toDenseMatrix();

}

void TTC_HTM2::setCapMatrix() {

	C = MetalCapacitances;
	Cinv = C.inverse();
}

double TTC_HTM2::KRadiationSurfaces(double A1, double A2, double e1, double e2, double F12, double T1, double T2) {

	double R1 = (1 - e1) / A1 / e1;
	double R2 = (1 - e2) / A2 / e2;
	double R12 = 1 / A1 / F12;

	double K = __cons::Sigma * (pow2(T1) + pow2(T2)) * (T1 + T2) / (R1 + R12 + R2);

	return K;

}

void TTC_HTM2::setKRadMatrix() {

	if(FExternalRadiation) {

		int _c = MetalNode["C"];
		int _t = MetalNode["T"];
		int _h1 = MetalNode["H1"];
		int _h2 = MetalNode["H2"];
		int _h3 = MetalNode["H3"];
		int pt = Part["TURBINE"];
		int ph = Part["HOUSING"];
		int pc = Part["COMPRESSOR"];

		Krad(_t, _h1) = KRadiationSurfaces(ExternalArea(pt), ExternalArea(ph) * FAlpha(0), Emissivity["TUR"], Emissivity["HOU"],
			VF_nodes(_t, _h1), T1(_t), T1(_h1)) * FFitRadiation;
		Krad(_t, _h2) = KRadiationSurfaces(ExternalArea(pt), ExternalArea(ph) * FAlpha(1), Emissivity["TUR"], Emissivity["HOU"],
			VF_nodes(_t, _h2), T1(_t), T1(_h2)) * FFitRadiation;
		Krad(_t, _h3) = KRadiationSurfaces(ExternalArea(pt), ExternalArea(ph) * FAlpha(2), Emissivity["TUR"], Emissivity["HOU"],
			VF_nodes(_t, _h3), T1(_t), T1(_h3)) * FFitRadiation;
		Krad(_t, _c) = KRadiationSurfaces(ExternalArea(pt), ExternalArea(pc), Emissivity["TUR"], Emissivity["COM"],
			VF_nodes(_t, _c), T1(_t), T1(_c)) * FFitRadiation;

		Krad(_c, _t) = KRadiationSurfaces(ExternalArea(pc), ExternalArea(pt), Emissivity["COM"], Emissivity["TUR"],
			VF_nodes(_c, _t), T1(_c), T1(_t)) * FFitRadiation;
		Krad(_c, _h1) = KRadiationSurfaces(ExternalArea(pc), ExternalArea(ph) * FAlpha(0), Emissivity["COM"], Emissivity["HOU"],
			VF_nodes(_c, _h1), T1(_c), T1(_h1)) * FFitRadiation;
		Krad(_c, _h2) = KRadiationSurfaces(ExternalArea(pc), ExternalArea(ph) * FAlpha(1), Emissivity["COM"], Emissivity["HOU"],
			VF_nodes(_c, _h2), T1(_c), T1(_h2)) * FFitRadiation;
		Krad(_c, _h3) = KRadiationSurfaces(ExternalArea(pc), ExternalArea(ph) * FAlpha(2), Emissivity["COM"], Emissivity["HOU"],
			VF_nodes(_c, _h3), T1(_c), T1(_h3)) * FFitRadiation;

		Krad(_h1, _t) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pt) * FAlpha(0), Emissivity["HOU"], Emissivity["TUR"],
			VF_nodes(_h1, _t), T1(_h1), T1(_t)) * FFitRadiation;
		Krad(_h1, _c) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pc) * FAlpha(0), Emissivity["COM"], Emissivity["TUR"],
			VF_nodes(_h1, _c), T1(_h1), T1(_c)) * FFitRadiation;

		Krad(_h2, _t) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pt) * FAlpha(1), Emissivity["HOU"], Emissivity["TUR"],
			VF_nodes(_h2, _t), T1(_h2), T1(_t)) * FFitRadiation;
		Krad(_h2, _c) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pc) * FAlpha(1), Emissivity["COM"], Emissivity["TUR"],
			VF_nodes(_h2, _c), T1(_h2), T1(_c)) * FFitRadiation;

		Krad(_h3, _t) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pt) * FAlpha(2), Emissivity["HOU"], Emissivity["TUR"],
			VF_nodes(_h3, _t), T1(_h3), T1(_t)) * FFitRadiation;
		Krad(_h3, _c) = KRadiationSurfaces(ExternalArea(ph), ExternalArea(pc) * FAlpha(2), Emissivity["COM"], Emissivity["TUR"],
			VF_nodes(_h3, _c), T1(_h3), T1(_c)) * FFitRadiation;
	}
	else{
		Krad.setZero();
	}

}

double TTC_HTM2::AmbientHeatFlowConv() {

	return Qext_W.col(Boundary["AMB"]).sum();
}

double TTC_HTM2::AmbientToTurbineHeatFlowConv() {

	return Qext_W(Boundary["AMB"], MetalNode["T"]);
}

double TTC_HTM2::AmbientToHousingHeatFlowConv() {

	return Qext_W(Boundary["AMB"], MetalNode["H1"]) + Qext_W(Boundary["AMB"], MetalNode["H2"]) + Qext_W(Boundary["AMB"], MetalNode["H3"]);
}

double TTC_HTM2::AmbientToCompressorHeatFlowConv() {

	return Qext_W(Boundary["AMB"], MetalNode["C"]);
}

double TTC_HTM2::AmbientHeatFlowRad() {

	return Qrad_W.col(Boundary["AMB"]).sum();
}

double TTC_HTM2::AmbientToTurbineHeatFlowRad() {

	return Qrad_W(Boundary["AMB"], MetalNode["T"]);
}

double TTC_HTM2::AmbientToHousingHeatFlowRad() {

	return Qrad_W(Boundary["AMB"], MetalNode["H1"]) + Qrad_W(Boundary["AMB"], MetalNode["H2"]) + Qrad_W(Boundary["AMB"], MetalNode["H3"]);
}

double TTC_HTM2::AmbientToCompressorHeatFlowRad() {

	return Qrad_W(Boundary["AMB"], MetalNode["C"]);
}

double TTC_HTM2::AmbientHeatFlowTotal() {

	return AmbientHeatFlowConv() + AmbientHeatFlowRad();
}

double TTC_HTM2::WaterHeatFlow() {

	return Qext_W.col(Boundary["WAT"]).sum();
}

double TTC_HTM2::TurbineHeatFlow() {

	return Qext_W.col(Boundary["GAS"]).sum();
}

double TTC_HTM2::CompressorHeatFlow() {

	return Qext_W.col(Boundary["AIR"]).sum();
}

double TTC_HTM2::OilHeatFlow() {

	return Qext_W.col(Boundary["OIL"]).sum();
}

double TTC_HTM2::getAdiabaticEfficiency(string Case, TurboMachine T, TurboMachine C) {

	C.Cp = Fluid["AIR"]->FunCp(C.it);
	C.g = Fluid["AIR"]->FunGamma(C.it);
	T.Cp = Fluid["GAS"]->FunCp(T.it);
	T.g = Fluid["GAS"]->FunGamma(T.it);
	double O_Cp = Fluid["OIL"]->FunCp(T_oil(0));
	double rtc;
	double M_eff;
	bool error = true;
	double err = 1.;
	double ip;
	double DeltaTMech;
	double DT_max = 0.05;
	int itera = 0;

	T.sin = __geom::Circle_area(D_in(Boundary["GAS"]));
	C.sin = __geom::Circle_area(D_in(Boundary["AIR"]));
	T.sout = __geom::Circle_area(D_out(Boundary["GAS"]));
	C.sout = __geom::Circle_area(D_out(Boundary["AIR"]));

	T.sig = -1;
	C.sig = 1;

	T.op = __Ambient::p_Pa;
	T.ip0 = T.op * T.pr;
	C.ip0 = __Ambient::p_Pa;
	C.op0 = C.ip0 * C.pr;

	if(Case == "Turbine") {
		if(C.m == 0)
			C.m = T.m;
		if(C.eff == 0)
			C.eff = 0.6;
		T.dit = T.it;
		T.ip = T.ip0;
		while(err > 0.0001) {
			T.it0 = T.get_it0();
			ip = T.get_p_static(T.ip0, T.it, T.it0, T.g);
			err = fabs(ip - T.ip);
			T.ip = ip;
		}
		T.ots = T.get_ts();
		T.pows = T.get_pows();
		T.powr = T.get_powr();
		C.pr = pow(T.powr * C.eff / (C.m * C.Cp * C.it) + 1, C.g / (C.g - 1));
		double Lim1 = C.pr - 0.1;
		if(Lim1 < 1)
			Lim1 = 1;
		double Lim2 = C.pr + 0.1;
		stComPower CompressorPower(T.powr, C.m, C.eff, C.it, C.Cp, (C.g - 1) / C.g);

		if(zbrac(CompressorPower, Lim1, Lim2)) {
			C.pr = FindRoot(CompressorPower, Lim1, Lim2);
		}
		rtc = T.rtc;
		C.op0 = C.ip0 * C.pr;
		err = 1;
		C.ip = C.ip0;
		while(err > 0.0001) {
			C.it0 = C.get_it0();
			ip = C.get_p_static(C.ip0, C.it, C.it0, C.g);
			err = fabs(ip - C.ip);
			C.ip = ip;
		}
		C.ots = C.get_ts();
		C.pows = C.get_pows();

	} else if(Case == "Compressor") {
		if (T.m == 0){
			T.m = C.m * C.pr * sqrt(C.it / T.it);

			double umax = 100;
			double T0max = T.it + pow2(umax) / 2 / T.Cp;
			double pmax = (C.pr + 0.5) * pow(T0max / T.it, T.g / (1 - T.g));
			double rhomax = __units::BarToPa(pmax) / Fluid["GAS"]->getR() / T.it;
			double mmax = rhomax * umax * T.sin;
			if (T.m > mmax)
				T.m = mmax;
		}

		if(T.eff == 0)
			T.eff = 0.6;
		C.ip = C.ip0;
		err = 1.;
		while(err > 0.0001) {
			C.it0 = C.get_it0();
			ip = C.get_p_static(C.ip0, C.it, C.it0, C.g);
			err = fabs(ip - C.ip);
			C.ip = ip;
		}
		//C.it0 = C.get_it0();
		C.ots = C.get_ts();
		C.pows = C.get_pows();
		C.powr = C.get_powr();
		T.pr = pow(1 - C.powr / (T.eff * T.m * T.Cp * T.it), -T.g / (T.g - 1));
		double Lim1 = T.pr - 0.1;
		if(Lim1 < 1)
			Lim1 = 1;
		double Lim2 = T.pr + 0.1;
		stTurPower TurbinePower(C.powr, T.m, T.eff, T.it, T.Cp, (T.g - 1) / T.g);
		if(zbrac(TurbinePower, Lim1, Lim2)) {
			T.pr = FindRoot(TurbinePower, Lim1, Lim2);
		}
		rtc = C.rtc;
		T.ip0 = T.pr * T.op;
		err = 1;
		T.ip = T.ip0;
		while(err > 0.0001) {
			T.it0 = T.get_it0();
			ip = T.get_p_static(T.ip0, T.it, T.it0, T.g);
			err = fabs(ip - T.ip);
			T.ip = ip;
		}
		T.dit = T.it;
	} else {

	}

	double O_pow = FMechLosses->P_oil(T_oil(0), __units::RPMToRad_s(rtc), 1.0, C.pr, T.pr, 1.0, m_oil);
	if(Case == "Turbine")
		M_eff = 1 - O_pow / T.powr;
	else if(Case == "Compressor")
		M_eff = C.powr / (C.powr + O_pow);

	C.ot0 = C.get_ot0();
	C.ot = C.get_t_static(C.ot0, C.op0, C.m, C.sout, C.Cp, C.g);
	C.op = C.get_p_static(C.op0, C.ot, C.ot0, C.g);
	C.aot = C.ot + 10;

	if(C.ot < T_oil(0)) {
		C.aot = C.ot - 10;
		if(C.aot < C.ots)
			C.aot = C.ots;
	}

	T.ot = T.get_ts();


	//O_pow = T.pows * (1 - M_eff);
	T_oil(2) = T_oil(0) + O_pow / m_oil / O_Cp;

	double DOT0 = C.aot;
	double TOT0 = T.ot;
	double OOT0 = T_oil(2);

	T1(0) = 0.80 * (T.it - C.ot) + C.ot;
	T1(1) = 0.50 * (T.it - C.ot) + C.ot;
	T1(2) = 0.15 * (T.it - C.ot) + C.ot;
	T1(3) = 0.10 * (T.it - C.ot) + C.ot;
	T1(4) = 0.05 * (T.it - C.ot) + C.ot;

	//T1 = VectorXd::Constant(nnodes, T_oil(0));

	T0 = T1;

	Tb(Boundary["GAS"]) = T.dit;
	Tb(Boundary["AIR"]) = C.ot;
	Tb(Boundary["OIL"]) = T_oil(0);
	Tb(Boundary["WAT"]) = T_cool(0);
	Tb(Boundary["AMB"]) = C.it;

	while(error) {
		Qext_J = MatrixXd::Zero(nnodes, nboundy);
		Qrad_J = MatrixXd::Zero(nnodes, nboundy);

		setHeatFromGas(1, T.dit, T.m);
		setHeatFromAir(1, (C.ot + C.aot) / 2., C.m);
		setHeatFromOil(1, T_oil(0), m_oil, O_pow);
		setHeatFromCoolant(1, T_cool(0), m_cool);
		setHeatFromAmbient(1, __Ambient::T_K, __Ambient::p_Pa);

		setRadiatioFromAmbient(1., __Ambient::T_K);

		setKRadMatrix();

		Tb(Boundary["OIL"]) = T_oilH2;

		SolveImplicit(1., 0);

		T.it = T.dit - TurbineHeatFlow() / T.m / T.Cp;

		C.aot = C.ot + CompressorHeatFlow() / C.m / C.Cp;
		if(C.aot < C.ots)
			C.aot = C.ots;

		C.ot0 = C.aot + pow2(C.m * 287 * C.aot / (C.op * C.sout)) / (2 * C.Cp);
		C.it0 = C.it + pow2(C.m * 287 * C.it / (C.ip * C.sin)) / (2 * C.Cp);
		C.powr = C.m * C.Cp * (C.ot0 - C.it0);

		T.powr = C.powr + O_pow;

		T.it0 = T.it + pow2(T.m * 287 * T.it / (T.ip * T.sin)) / (2 * T.Cp);
		T.ot = T.it0 - T.powr / T.m / T.Cp;

		if((fabs(C.aot - DOT0) > DT_max || fabs(T.ot - TOT0) > DT_max || fabs(T_oil(2) - OOT0) > DT_max) && itera < 1000) {
			T0 = T1;
			DOT0 = C.aot;
			TOT0 = T.ot;
			OOT0 = T_oil(2);
			itera += 1;
		} else {
			error = false;
			if(itera == 1000) {
				printf("Process does not converge. Error: Turbine %lf K, Compressor %lf K, Oil %lf K\n", T.ot - TOT0, C.aot - DOT0,
					   T_oil(2) - OOT0);
			}
		}

	}

	if(Case == "Compressor") {
		C.aeff = C.pows / C.powr;

		PrintCompressorAdiMap(T, C);

		return C.aeff;
	} else if(Case == "Turbine") {
		T.ots = T.get_ts();
		T.pows = T.get_pows();
		T.aeff = T.powr / T.pows;

		PrintTurbineAdiMap(T, C);

		return T.aeff;
	}
	return 0.;
}

void TTC_HTM2::PrintTurbineAdiMap(TurboMachine T, TurboMachine C) {

	FILE *fres = fopen("TurbRes.dat", "a");

	fprintf(fres, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", T.rtc, T.m, T.pr, T.aeff, T.eff,
			TurbineHeatFlow(), CompressorHeatFlow(), OilHeatFlow(), WaterHeatFlow(), AmbientHeatFlowTotal());

	for(int i = 0; i < nnodes; i++) {
		fprintf(fres, "\t%lf", T1(i));
	}
	for(int i = 0; i < nboundy; i++) {
		fprintf(fres, "\t%lf", Tb(i));
	}

	fprintf(fres, "\t%lf", T_oil(2));
	fprintf(fres, "\t%lf", C.ot);
	fprintf(fres, "\t%lf", T_cool(1));

	fprintf(fres, "\t%lf\t%lf\t%lf\t%lf", C.powr, T.powr, C.powr / T.powr, T.m * T.Cp * T.dit);

	fprintf(fres, "\n");
	fclose(fres);
}

void TTC_HTM2::PrintCompressorAdiMap(TurboMachine T, TurboMachine C) {

	FILE *fres = fopen("CompRes.dat", "a");

	fprintf(fres, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", C.rtc, C.m, C.pr, C.aeff, C.eff, TurbineHeatFlow(),
			CompressorHeatFlow(), OilHeatFlow(), WaterHeatFlow(), AmbientHeatFlowTotal());

	for(int i = 0; i < nnodes; i++) {
		fprintf(fres, "\t%lf", T1(i));
	}
	for(int i = 0; i < nboundy; i++) {
		fprintf(fres, "\t%lf", Tb(i));
	}

	fprintf(fres, "\t%lf", T_oil(2));
	fprintf(fres, "\t%lf", T.ot);
	fprintf(fres, "\t%lf", T_cool(1));

	fprintf(fres, "\t%lf\t%lf\t%lf\t%lf", C.powr, T.powr, C.powr / T.powr, T.m * T.Cp * T.din);

	fprintf(fres, "\n");
	fclose(fres);

}

void TTC_HTM2::TurbochargerData(double Doil, double DWater, double DT, double LT, double DC,
								double LC, double DH, double LH) {

	D_in(Boundary["OIL"]) = Doil;
	D_in(Boundary["WAT"]) = DWater;

	ExternalDiameter(Part["TURBINE"]) = DT;
	ExternalLength(Part["TURBINE"]) = LT;
	ExternalArea(Part["TURBINE"]) = __geom::Ring_area(DH, DT);
	ExternalDiameter(Part["HOUSING"]) = DH;
	ExternalLength(Part["HOUSING"]) = LH;
	ExternalArea(Part["HOUSING"]) = __cons::Pi * DH * LH;
	ExternalDiameter(Part["COMPRESSOR"]) = DC;
	ExternalLength(Part["COMPRESSOR"]) = LC;
	ExternalArea(Part["COMPRESSOR"]) = __geom::Ring_area(DH, DC);

	if(FExternalRadiation) {
		setViewFactorBetweenNodes(DT * FFitTurbExtHeat, DC, DH, LH);

		setViewFactorToBoundaries(DT * FFitTurbExtHeat, DC, DH, LH);
	}

}

void TTC_HTM2::setOilConditions(double m, double T) {
	m_oil = m;
	T_oil = ArrayXd::Constant(3, T);
	T_oilMech = T;
	T_oilH2 = T;
}

void TTC_HTM2::setCoolantConditions(double m, double T) {
	m_cool = m;
	T_cool(0) = T;
	T_cool(1) = T;
}

void TTC_HTM2::HeaderInsTemperatures(stringstream & insoutput, int num) {

	std::string Label;

	for(int i = 0; i < nnodes; i++) {
		Label = "\t" + PutLabel(5007) + "/" + std::to_string(num) + "/" + PutLabel(4005) + PutLabel(910) + "/" + NodeLabel[i];
		insoutput << Label.c_str();
	}

}

void TTC_HTM2::PrintInsTemperatures(stringstream & insoutput) {

	for(int i = 0; i < nnodes; i++) {
		insoutput << "\t" << __units::KTodegC(T1(i));
	}
}

void TTC_HTM2::PrintInsHeatFlow(stringstream & insoutput) {

	for(int i = 0; i < nnodes; i++) {
		for(int j = 0; j < nboundy; j++) {
			if(ConvectiveConn(i, j) == 1)
				insoutput << "\t" << Qext_W(i, j);
			if(RadiativeConn(i, j) == 1)
				insoutput << "\t" << Qrad_W(i, j);
		}
	}

}

void TTC_HTM2::HeaderInsHeatFlow(stringstream & insoutput, int num) {

	std::string Label;

	for(int i = 0; i < nnodes; i++) {
		for(int j = 0; j < nboundy; j++) {
			if(ConvectiveConn(i, j) == 1) {
				Label = "\t" + PutLabel(5007) + "/" + std::to_string(num) + "/" + PutLabel(4010) + PutLabel(
							903) + "/" + "Convected_From_" + BoundaryLabel[j] + "_To_" + NodeLabel[i];
				insoutput << Label.c_str();
			}
			if(RadiativeConn(i, j) == 1) {
				Label = "\t" + PutLabel(5007) + "/" + std::to_string(num) + "/" + PutLabel(4010) + PutLabel(
							903) + "/" + "Radiated_From_" + BoundaryLabel[j] + "_To_" + NodeLabel[i];
				insoutput << Label.c_str();
			}
		}
	}

}

void TTC_HTM2::InitializeTemp(double T3, double T2, double Toil, double Twater){

	if (FInitialTemp){
		T1 = Tinit;
	}
	else{
		T1(0) = 0.80 * (T3 - T2) + T2;
		T1(1) = 0.50 * (T3 - T2) + T2;
		T1(2) = 0.15 * (T3 - T2) + T2;
		T1(3) = 0.10 * (T3 - T2) + T2;
		T1(4) = 0.05 * (T3 - T2) + T2;
	}

	T0 = T1;

	T_cool.setConstant(Twater);

	T_oil.setConstant(Toil);
	T_oilH2 = Toil;
	T_oilMech = Toil;

}

