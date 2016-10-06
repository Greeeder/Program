/* --------------------------------------------------------------------------------*\
|===========================|
 | \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 |  \\ |  X  | //  W ave     |
 |   \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |    \\/   \//    M odel    |
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


 \*-------------------------------------------------------------------------------- */

/**
 * @file TPipeHeatT.cpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
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
 * This file defines the heat transfer in the pipes and the calculation of the
 * wall temperatures.
 */

// ---------------------------------------------------------------------------
#include "TPipeHeatT.h"

TPipeHeatT::TPipeHeatT() {

}

TPipeHeatT::~TPipeHeatT() {

}

void TPipeHeatT::ReadHeatTransferData(xml_node node_pipe, int nc, double cs, std::map<string, TSolid*> mdb) {

	string MatName;
	
	xml_node node_ht = GetNodeChild(node_pipe, "Pip_HeatTransfer");

	nlayers = CountNodes(node_ht, "Pht_Layer");
	string WallCalc = node_ht.attribute("WallCalculation").as_string();
	if (WallCalc != "Constant" && nlayers == 0){
		cout << "ERROR: One or more layer must be defined to calculate pipe wall temperature" << endl;
	}
	Material.resize(nlayers);

	ncells = nc;
	CellSize = cs;
	nnodes = ncells * nlayers;

	Rint = VectorXd::Zero(nnodes);
	Rext = VectorXd::Zero(nnodes);
	Conductance = VectorXd::Zero(nlayers);
	double T = GetAttributeAsDouble(node_ht, "WallTemperature");

	int i = 0;
	for(xml_node node_layer = GetNodeChild(node_ht, "Pht_Layer"); node_layer;
		node_layer = node_layer.next_sibling("Pht_Layer")) {

		MatName = node_layer.attribute("Material").as_string();
		Material[i] = new TSolid(mdb[MatName]);
		Material[i]->setTemperature(T);
		Thickness.push_back(GetAttributeAsDouble(node_layer, "Thickness"));
		IsFluid.push_back(GetAttributeAsBool(node_layer, "Fluid"));
		i++;

	}

	Initiate(T);
}

void TPipeHeatT::BuildMatrix(Eigen::VectorXd D) {

	VectorXd ri = D / 2;
	VectorXd Ones = VectorXd::Ones(ncells);
	Rint.segment(0, ncells) = ri;
	Rext.segment(0, ncells) = ri + Ones * Thickness[0];

	// Vector for internal and external radiud.
	//
	for(int j = 1; j < nlayers; j++) {
		Rint.segment(j * ncells, ncells) = Rext.segment((j - 1) * ncells, ncells);
		Rext.segment(j * ncells, ncells) = Rint.segment(j * ncells, ncells) + Ones * Thickness[j];
	}

	// Build capacitances matrix.

	MatrixXd c = (Rext.segment(0, ncells).array().pow(2) - Rint.segment(0,
						  ncells).array().pow(2)).matrix() * Material[0]->getDensity() * Material[0]->getHeatCapacity() * __cons::Pi * CellSize;
	C.block(0, 0, ncells, ncells) = c.asDiagonal();

	if(nlayers > 1) {
		for(int j = 1; j < nlayers; j++) {
			c = (Rext.segment(j * ncells, ncells).array().pow(2) - Rint.segment(j * ncells,
						 ncells).array().pow(2)).matrix() * Material[j]->getDensity() * Material[j]->getHeatCapacity() * __cons::Pi * CellSize;
			C.block(j * ncells, j * ncells, ncells, ncells) = c.asDiagonal().toDenseMatrix();
		}
	}

	Cinv = C.inverse();

	// Build conductances matrix
	
	VectorXd k;
	VectorXd con1;

	for(int j = 0; j < nlayers; j++) {
		K_layer.segment(j * ncells, ncells) = (2 * __cons::Pi * Material[j]->getConductivity() * CellSize / (Rext.segment(j * ncells,
											   ncells).array() / Rint.segment(j * ncells, ncells).array()).log()).matrix();
	}

	for(int j = 0; j < nlayers; j++) {
		k = __cons::Pi_2 * ((Rext.segment(j * ncells, ncells).array().pow(2) - Rint.segment(j * ncells, ncells).array().pow(2)).head(ncells - 1)
			   + (Rext.segment(j * ncells, ncells).array().pow(2) - Rint.segment(j * ncells,
					   ncells).array().pow(2)).tail(ncells - 1)).matrix() * Material[j]->getConductivity() / CellSize;

		K1.block(j * ncells, j * ncells + 1, ncells - 1, ncells - 1) = k.asDiagonal();
		if(j < nlayers - 1) {
			con1 = (K_layer.segment(j * ncells, ncells) + K_layer.segment((j + 1) * ncells, ncells)) * 0.5;
			K1.block(j * ncells, (j + 1) * ncells, ncells, ncells) = con1.matrix().asDiagonal();
		}
	}

	MatrixXd K_tras = K1.transpose();
	K1 = K1 + K_tras;

	K1 = K1.rowwise().sum().asDiagonal().toDenseMatrix() - K1;

}

void TPipeHeatT::Initiate(double T) {

	Twall0.setConstant(nnodes, T);
	Twall1.setConstant(nnodes, T);

	Twallint.setConstant(ncells, T);
	Twallext.setConstant(ncells, T);

	Qint.setZero(ncells);

	K0.setZero(nnodes, nnodes);
	K1.setZero(nnodes, nnodes);

	K_layer.setZero(nnodes);

	C.setZero(nnodes, nnodes);


	Rint.setZero(nnodes);
	Rext.setZero(nnodes);

}

void TPipeHeatT::Solve() {

	VectorXd b;
	VectorXd Q = VectorXd::Zero(nnodes);
	Twall0 = Twall1;
	Q.head(ncells) = Qint;
	Q.tail(ncells) = Q.tail(ncells) + ExtHeat->Heat(Twallext);

	b = Q + (K1 - C / TimeStep) * Twall0;

	Twall1 = (K1 + C / TimeStep).ldlt().solve(b);

	Twallint = Twall1.head(ncells) + (Qint.array() / (0.5 * K_layer.head(ncells).array())).matrix();
	Twallext = Twall1.tail(ncells) + (Q.tail(ncells).array() / (0.5 * K_layer.tail(ncells).array())).matrix();

	Qint.setZero();
}

void TPipeHeatT::SolveExplicit(double dt) {

	VectorXd Q = VectorXd::Zero(nnodes);
	VectorXd Qext = ExtHeat->Heat(Twallext) * dt;
	Q.head(ncells) = Qint;
	Q.tail(ncells) = Q.tail(ncells) + Qext;
	Twall0 = Twall1;
	Twall1 = Twall0 + Cinv * (dt * (-K1) * Twall0 + (Q));

	Twallint = Twallint + 4 * Cinv.topLeftCorner(ncells, ncells) * (((Twall0.head(ncells) - Twallint).array()
		* 2 * K_layer.head(ncells).array() * dt + Qint.array())).matrix();
	Twallext = Twallext + 4 * Cinv.bottomRightCorner(ncells, ncells) * (((Twall0.tail(ncells) - Twallext).array()
		* 2 * K_layer.tail(ncells).array() * dt + Qext.array())).matrix();

	Qint.setZero();

}

void TPipeHeatT::Solve(double dt) {

}

void TPipeHeatT::AddInternalHeat(int i, int q) {
	Qint(i) = Qint(i) + q;
}

void TPipeHeatT::AddInternalHeat(Eigen::VectorXd q) {
	Qint = Qint + q;
}

