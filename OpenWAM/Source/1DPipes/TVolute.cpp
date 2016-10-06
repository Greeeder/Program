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

/*!
 * \file TVolute.cpp
 * \author Francisco Jose Arnau <farnau@mot.upv.es>
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 *
 * \section LICENSE
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
 * \section DESCRIPTION
 * This file defines a quasi-2D volute.
 */

#include "TVolute.hpp"

double TVolute::getLateralEnthalpyFlow(double x) const {
	return linear_interp(Fx, FLateralEnthalpyFlow, x);
}

double TVolute::getLateralEnthalpyFlow(Uint i) const {
	if(i < FLateralEnthalpyFlow.size()) {
		return FLateralEnthalpyFlow(i);
	} else {
		return FLateralEnthalpyFlow.tail(1)(0);
	}
}

RowVector TVolute::getLateralEnthalpyFlow() const {
	return FLateralEnthalpyFlow;
}

double TVolute::getLateralMassFlow(double x) const {
	return linear_interp(Fx, FLateralMassFlow, x);
}

double TVolute::getLateralMassFlow(Uint i) const {
	if(i < FLateralMassFlow.size()) {
		return FLateralMassFlow(i);
	} else {
		return FLateralMassFlow.tail(1)(0);
	}
}

RowVector TVolute::getLateralMassFlow() const {
	return FLateralMassFlow;
}

RowVector TVolute::getLateralNozzleArea() const {
	return dynamic_cast<TLateralNozzles*>(FSource.get())->getArea();
}

double TVolute::getLateralNozzleArea(Uint i) const
{
	return dynamic_cast<TLateralNozzles*>(FSource.get())->getArea(i);
}

double TVolute::getWindowArea(double x) const {
	return linear_interp(Fx, FWindowArea, x);
}

double TVolute::getWindowArea(Uint i) const {
	if(i < FWindowArea.size()) {
		return FWindowArea(i);
	} else {
		return FWindowArea.tail(1)(0);
	}
}

RowVector TVolute::getWindowArea() const {
	return FWindowArea;
}

void TVolute::setGeometry(const RowVector& x, double dx, const RowVector& D) {
	double n_nodes = round((x.tail(1)(0) - x(0)) / dx) + 1;
	dx = (x.tail(1)(0) - x(0)) / n_nodes;
	Fdx.setConstant(n_nodes - 1, dx);
	FXref = Fdx(0);
	Fx.setLinSpaced(n_nodes, x(0), x.tail(1)(0));
	FMx = (Fx.head(n_nodes - 1) + Fx.tail(n_nodes - 1)) / 2.;
	FArea = linear_interp(x, D.unaryExpr(&__geom::Circle_area), Fx);
	FD = FArea.unaryExpr(&__geom::Circle_diameter);
	FMArea = (FArea.head(n_nodes - 1) + FArea.tail(n_nodes - 1)) / 2.;
	FDcell = FMArea.unaryExpr(&__geom::Circle_diameter);
	FVolume = dx * (FArea.head(n_nodes - 1) + FArea.tail(n_nodes - 1)) / 2.;
	FDerLinArea = (FArea.tail(n_nodes - 1) - FArea.head(n_nodes - 1)) / Fdx;
	FU0.setZero(3, n_nodes);
	FNin = n_nodes;
}

void TVolute::setLateralNozzleArea(RowVector A) {
	dynamic_cast<TLateralNozzles*>(FSource.get())->setArea(A);
}

void TVolute::setLateralNozzleArea(double A) {
	dynamic_cast<TLateralNozzles*>(FSource.get())->setArea(A);
}

void TVolute::setLateralNozzleArea(double A, Uint i) {
	dynamic_cast<TLateralNozzles*>(FSource.get())->setArea(A, i);
}

void TVolute::UpdateStateVector() {
	TPipe::UpdateStateVector();
	FLateralMassFlow = -FSource->getSource().row(0);
	FLateralEnthalpyFlow = -FSource->getSource().row(2);
}

Volute_ptr create_volute(const xml_node& node,
	const TComponentArray_ptr& fluid,
	const std::map<string, TSolid*> &MDB)
{
	auto pipe = make_shared<TVolute>();
	pipe->setName(node.attribute("Name").as_string());
	pipe->FNodeLeft = GetAttributeAsInt(node, "NodeL_ID");
	pipe->FNodeRight = GetAttributeAsInt(node, "NodeR_ID");

	pipe->FNIdenticalPipes = GetAttributeAsInt(node, "ParallelPipes");

	auto node_geometry = GetNodeChild(node, "Pip_Geometry");
	auto node_num_method = GetNodeChild(node, "Pip_NumericalMethod");
	std::string numerical_method = node_num_method.attribute("Scheme").as_string();
	
	pipe->FMeshSize = GetAttributeAsDouble(node_num_method, "MeshSize");

	pipe->FNumStretch = CountNodes(node_geometry, "Geo_Stretch");

	pipe->FDStretch.resize(pipe->FNumStretch + 1);
	pipe->FLStretch.resize(pipe->FNumStretch + 1);

	pipe->FDStretch[0] = GetAttributeAsDouble(node_geometry, "Diameter");
	pipe->FLStretch[0] = 0.;

	auto id = 0;
	for(auto node_stretch = GetNodeChild(node_geometry, "Geo_Stretch");
		node_stretch; node_stretch = node_stretch.next_sibling("Geo_Stretch")) {
		id = GetAttributeAsInt(node_stretch, "Stretch_ID");
		if(id > pipe->FNumStretch) {
			std::cout << "ERROR: The stretches of the pipe are not correctly ordered in pipe: "
				<< pipe->FPipeNumber << std::endl;
			throw Exception("");
		}
		pipe->FLStretch[id] = GetAttributeAsDouble(node_stretch, "Length");
		pipe->FDStretch[id] = GetAttributeAsDouble(node_stretch, "Diameter");
	}

	ArrayXd LStr(pipe->FLStretch.size());
	LStr(0) = pipe->FLStretch(0);
	for(auto i = 1; i < pipe->FLStretch.size(); ++i) {
		LStr(i) = pipe->FLStretch(i) + LStr(i - 1);
	}

	pipe->setGeometry(LStr, pipe->FMeshSize, pipe->FDStretch);

	auto node_prop = GetNodeChild(node, "GasProperties");

	double Pini, Tini, Vmean;
	RowVector Y;
	Y.setZero(fluid->size());
	ReadGasProperties(node_prop, Pini, Tini, Vmean, Y, fluid);

	pipe->WorkingFluid.resize(pipe->FD.size() - 1);
	for (auto i = 0; i < pipe->WorkingFluid.size(); i++) {
		pipe->WorkingFluid[i] = make_shared<TFluid>(fluid);
		pipe->WorkingFluid[i]->SetComposition(Y);
	}

	set_numerical_method(pipe, numerical_method);
	set_null_source(pipe);

	pipe->setPTU(Pini, Tini, Vmean);
	pipe->FLateralMassFlow.setZero(pipe->getNin());
	pipe->FLateralEnthalpyFlow.setZero(pipe->getNin());

	close_pipe_end(pipe, nmRight);
	pipe->setCourant(GetAttributeAsDouble(node_num_method, "Courant"));
	pipe->ReadInsResults(node);

	// HEAT TRANSFER -------
	xml_node node_ht = GetNodeChild(node, "Pip_HeatTransfer");
	double IntHeatMult = GetAttributeAsDouble(node_ht, "IntMultiplier");
	string Type = node_ht.attribute("HT_Type").as_string();

	if (Type == "IntakePipe"){
		pipe->FHeatTransfer = make_unique<TPipeIntHeatIntake>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "ExhaustPipe"){
		pipe->FHeatTransfer = make_unique<TPipeIntHeatExhaust>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "IntakePort"){
		pipe->FHeatTransfer = make_unique<TPortIntHeatIntake>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else if (Type == "ExhaustPort"){
		pipe->FHeatTransfer = make_unique<TPortIntHeatExhaust>(pipe->FNin, IntHeatMult, pipe->FDcell, pipe->FXref);
	}
	else{
		cout << "Internal heat transfer in pipe not correctly defined" << endl;
	}

	string WallCalc = node_ht.attribute("WallCalculation").as_string();
	pipe->FWallHT = make_shared<TPipeHeatT>();
	pipe->FWallHT->ReadHeatTransferData(node, pipe->FNin, pipe->FXref, MDB);
	

	if (WallCalc != "Constant"){

		pipe->FWallHT->BuildMatrix(pipe->FDcell);
		xml_node node_eheat = GetNodeChild(node_ht, "Pht_External");
		double velocity = GetAttributeAsDouble(node_eheat, "Velocity");
		double extMult = GetAttributeAsDouble(node_eheat, "ExtMultiplier");
		double emis = GetAttributeAsDouble(node_eheat, "Emissivity");
		string Type2 = node_eheat.attribute("Type").as_string();
		TPipeExtHeat *extheatobj;
		if (Type2 == "AirCooled"){
			extheatobj = new TPipeExtHeatA(pipe->FNin, velocity, extMult, pipe->FWallHT->getExtDiameter().array(),
				pipe->FXref, emis);
		}
		else if (Type2 == "WaterCooled"){
			double Tw = GetAttributeAsDouble(node_eheat, "WaterTemp");
			extheatobj = new TPipeExtHeatW(pipe->FNin, velocity, extMult, pipe->FWallHT->getExtDiameter().array(),
				pipe->FXref, emis, Tw);
		}
		else if (Type2 == "Port"){

		}
		pipe->FWallHT->AddExtHeatObject(extheatobj);
	}

	xml_node node_friction = GetNodeChild(node, "Pip_Friction");
	string FricType = node_friction.attribute("Type").as_string();
	if (FricType == "Colebrook"){
		pipe->FPipeFriction = make_unique<TPipeFrictionColebrook>();
	}
	else{
		pipe->FPipeFriction = make_unique<TPipeFriction>();
	}
	pipe->FPipeFriction->ReadFrictionData(node);
	pipe->FPipeFriction->setRelativeRugosity(pipe->FDcell);
	pipe->FFric.setZero(pipe->FNin);

	return pipe;
}
