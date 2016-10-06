#include "Globales.h"
#include "TPipeHeatT.h"
#include "TPipeIntHeatIntake.h"
#include "TPipeIntHeatExhaust.h"
#include "TPortIntHeatIntake.h"
#include "TPortIntHeatExhaust.h"
#include "TPipeExtHeatA.h"
#include "TFluidAir.h"
#include "TSldAl.h"
#include "TMultiWiebe.h"
#include "THRL.h"
#include <chrono>
#include <fstream>

namespace __Ambient {
	double p_Pa = 100000;	   //!< The ambient pressure [Pa]
	double p_bar = 1;		   //!< The ambient pressure [bar]
	double T_K = 300;		   //!< The ambient temperature [K]
	double T_degC = 27;		   //!< The ambient temperature [degC]
	double HR = 50;			   //!< The humidity [%]
}
;

void TestHeatTransferInPipes(xml_node node_openwam);

void BuildMaterialsDB(xml_node node_mat, std::map<string, TSolid*> &MaterialsDB);

void TestCombustion(xml_node node_openwam);

int main(){

	xml_document File;
	if (!File.load_file("D:/OpenWAM/TCMVSProj/bin/debug/Input.xml"))
		cout << "The input file does not exist" << endl;
	xml_node node_openwam = File.child("OpenWAM");

	//TestHeatTransferInPipes(node_openwam);
	
	TestCombustion(node_openwam);

	return EXIT_SUCCESS;
}

void TestCombustion(xml_node node_openwam){

	TCombustion* Comb;
	TCombustion* Comb2;

	xml_node node_eng = GetNodeChild(node_openwam, "EngineBlock");
	xml_node node_cmb = GetNodeChild(node_eng, "Combustion");
	//xml_node node_wb = GetNodeChild(node_cmb, "Wiebe");
	xml_node node_mwb = GetNodeChild(node_cmb, "MultiWiebe");
	xml_node node_hrl = GetNodeChild(node_cmb, "HRL");
	//Comb = new TMultiWiebe();
	
	Comb = new THRL();
	Comb->ReadCombustionData(node_hrl);

	Comb2 = new THRL(dynamic_cast<THRL*>(Comb));

	for (int i = -10; i < 10; i++){
		cout << (double)i << " " << Comb->getHRL((double)i) << endl;
	}
	double hrl = Comb->getHRL(5);
	//double hrl2 = Comb2->getHRL(5);
}

void TestHeatTransferInPipes(xml_node node_openwam){

	int ncells = 5;
	double cellsize = 0.02;
	fstream fout;
	std::map<string, TSolid*> MaterialsDB;

	xml_node node_GD = GetNodeChild(node_openwam, "GeneralData");
	xml_node node_mat = GetNodeChild(node_GD, "Materials");

	BuildMaterialsDB(node_mat, MaterialsDB);

	xml_node node_pipeblock = GetNodeChild(node_openwam, "BlockOfPipes");
	xml_node np = GetNodeChild(node_pipeblock, "Bop_Pipe");

	VectorXd D = VectorXd::Ones(5) * 0.05;
	VectorXd Tgas = VectorXd::Constant(5, 530);
	VectorXd Re = VectorXd::Constant(5, 12000);
	std::vector<TFluid_ptr> Fluid;
	Fluid.resize(ncells);
	for (int i = 0; i < ncells; i++){
		Fluid[i] = make_shared<TFluidAir>();
		Fluid[i]->Funk(Tgas(i));
		Fluid[i]->FunVisc(Tgas(i));
	}


	xml_node node_ht = GetNodeChild(np, "Pip_HeatTransfer");
	double IntHeatMult = GetAttributeAsDouble(node_ht, "IntMultiplier");
	string Type = node_ht.attribute("HT_Type").as_string();


	TPipeIntHeat *intheatobj;
	TPipeHeatT *htobj;
	TPipeExtHeat *extheatobj;

	if (Type == "IntakePipe"){
		intheatobj = new TPipeIntHeatIntake(ncells, IntHeatMult, D, cellsize);
	}
	else if (Type == "ExhaustPipe"){
		intheatobj = new TPipeIntHeatExhaust(ncells, IntHeatMult, D, cellsize);
	}
	else if (Type == "IntakePort"){
		intheatobj = new TPortIntHeatIntake(ncells, IntHeatMult, D, cellsize);
	}
	else if (Type == "ExhaustPort"){
		intheatobj = new TPortIntHeatExhaust(ncells, IntHeatMult, D, cellsize);
	}
	else{
		cout << "Internal heat transfer in pipe not correctly defined" << endl;
	}

	string WallCalc = node_ht.attribute("WallCalculation").as_string();
	if (WallCalc != "Constant"){
		htobj = new TPipeHeatT();
		htobj->ReadHeatTransferData(np, ncells, 0.02, MaterialsDB);
		htobj->BuildMatrix(D);

		xml_node node_eheat = GetNodeChild(node_ht, "Pht_External");
		double velocity = GetAttributeAsDouble(node_eheat, "Velocity");
		double extMult = GetAttributeAsDouble(node_eheat, "ExtMultiplier");
		double emis = GetAttributeAsDouble(node_eheat, "Emissivity");
		string Type2 = node_eheat.attribute("Type").as_string();
		if (Type2 == "AirCooled"){
			extheatobj = new TPipeExtHeatA(ncells, velocity, extMult, htobj->getExtDiameter().array(),
				cellsize, emis);
		}
		else if (Type2 == "WaterCooled"){
			double Tw = GetAttributeAsDouble(node_eheat, "WaterTemp");
			extheatobj = new TPipeExtHeatW(ncells, velocity, extMult, htobj->getExtDiameter().array(),
				cellsize, emis, Tw);
		}
		else if (Type2 == "Port"){

		}
		htobj->AddExtHeatObject(extheatobj);
	}

	fout.open("Output.dat", fstream::out);
	VectorXd Qi;
	double dt = 0.1;
	auto start = chrono::steady_clock::now();
	for (int i = 0; i < 3000; i++){
		Qi = intheatobj->Heat(Tgas.array(), htobj->getTwallint().array(), Re, Fluid);


		Qi *= dt;

		htobj->AddInternalHeat(Qi);

		htobj->SolveExplicit(dt);

		fout << htobj->getTwallint()(0) << "\t";
		fout << htobj->getTnodes()(0) << "\t";
		fout << htobj->getTwallext()(0);
		fout << endl;
	}
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << chrono::duration <double, milli>(diff).count() << " ms" <<
		endl;

	fout.close();

	cout << endl << htobj->getTnodes() << endl;

	cout << htobj->getTwallint() << endl;
	cout << htobj->getTwallext();
}

void BuildMaterialsDB(xml_node node_mat, std::map<string, TSolid*> &MaterialsDB){

	for (xml_node node_m = GetNodeChild(node_mat, "Material"); node_m; 
		node_m = node_m.next_sibling("Material")){
		string matname = node_m.attribute("Name").as_string();
		if (matname == "Aluminium"){
			MaterialsDB[matname] = new TSldAl();
		}
		else{
			if (MaterialsDB.count(matname) == 0){
				MaterialsDB[matname] = new TSolid(matname);
				xml_node node_prop;
				int ncoef;
				int i;
				ArrayXd Coefs;

				node_prop = GetNodeChild(node_m, "Conductivity");
				ncoef = CountNodes(node_prop, "ConductivityCoef");
				Coefs.setZero(ncoef);
				i = 0;
				for (xml_node node_c = GetNodeChild(node_prop, "ConductivityCoef"); node_c;
					node_c = node_c.next_sibling("ConductivityCoef")){
					Coefs(i) = GetAttributeAsDouble(node_c, "Value");
					i++;
				}
				MaterialsDB[matname]->setCoefCond(Coefs);

				node_prop = GetNodeChild(node_m, "Density");
				ncoef = CountNodes(node_prop, "DensityCoef");
				Coefs.setZero(ncoef);
				i = 0;
				for (xml_node node_c = GetNodeChild(node_prop, "DensityCoef"); node_c;
					node_c = node_c.next_sibling("DensityCoef")){
					Coefs(i) = GetAttributeAsDouble(node_c, "Value");
					i++;
				}
				MaterialsDB[matname]->setCoefDens(Coefs);

				node_prop = GetNodeChild(node_m, "HeatCapacity");
				ncoef = CountNodes(node_prop, "HeatCapacityCoef");
				Coefs.setZero(ncoef);
				i = 0;
				for (xml_node node_c = GetNodeChild(node_prop, "HeatCapacityCoef"); node_c;
					node_c = node_c.next_sibling("HeatCapacityCoef")){
					Coefs(i) = GetAttributeAsDouble(node_c, "Value");
					i++;
				}
				MaterialsDB[matname]->setCoefHeCap(Coefs);
			}
			else{
				cout << "ERROR: Material " << matname << " is defined twice" << endl;
			}
		}
	}
}

