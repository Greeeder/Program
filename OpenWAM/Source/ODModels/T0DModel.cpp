/**
* @file T0DModel.cpp
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @date 30 de mar. de 2016
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
* This file include the different methods to solve the mass and energy balance in OD elements.
*/

#include "T0DModel.h"

T0DModel::T0DModel() {
	
	dFQL = 0.;
	EnthalpyIn = 0.;
	EnthalpyOut = 0.;
	Rcil = __R::Air;
	dUsen = 0.;
	dWi = 0.;
	Cv = __Gamma::Cv;
	dEbb = 0.;
	dEfinj = 0.;
	dEAinj = 0.;
	dQconv = 0.;
	dQrad = 0.;

}

T0DModel::~T0DModel() {
	// TODO Auto-generated destructor stub
}

void T0DModel::AppendMass(double m){
	PMass += m;
}

void T0DModel::Initialize(double Pressure, double Temperature, double Mass, double Volume){
	PPressure = Pressure;
	PPressure0 = Pressure;
	PTemperature = Temperature;
	PTemperature0 = Temperature;
	PMass = Mass;
	PMass0 = Mass;
	PVolume = Volume;
	PVolume0 = Volume;
}

void T0DModel::Initialize(double Pressure, double Temperature, double Mass, double Volume, TFluid *FluidCyl, BasicHeatTransfer *HeatTransfer)
{
	PPressure = Pressure;
	PPressure0 = Pressure;
	PTemperature = Temperature;
	PTemperature0 = Temperature;
	PMass = Mass;
	PMass0 = Mass;
	PVolume = Volume;
	PVolume0 = Volume;
	PHeatTransfer = HeatTransfer;
}

//Calculate Mass Balance
void T0DModel::SolveMassBalance(std::vector<double> *MassIn){

	//New mass is calculated
	for (int i = 0; i < MassIn->size(); i++){ 
		PMass += MassIn->at(i);
	}

}



/// devuelve la presión y la temperatura en un instante partiendo de los valores correspondientes en el instante anterior
///CICLO ABIERTO (CALMEC) - CICLO ABIERTO + CERRADO (SICILO)
void T0DModel::SolveNewPressure(double Volume, double dFQLp, std::vector<double> *MassIn, std::vector<double> *EnthIn, std::vector<TFluid*>*FluidIn, TFluid *FluidCyl){

	TFluid_ptr FlCyl = make_shared<TFluid>(*FluidCyl);
	//TFluidArray_ptr FluIn = make_shared<TFluidArray>(*FluidIn);
	TFluidArray Flarr;
	for (int i = 0; i < FluidIn->size(); i++){

		Flarr.push_back(make_shared<TFluid>(FluidIn->at(i)));
	}
	TFluidArray_ptr FluIn = make_shared<TFluidArray>(Flarr);

	SolveNewPressure(Volume, dFQLp, MassIn, EnthIn, FluIn, FlCyl);

	UpdateStateVector();
}



void T0DModel::SolveNewPressure(double Volume, double dFQLp, std::vector<double> *MassIn, std::vector<double> *EnthIn, TFluidArray_ptr FluidIn, TFluid_ptr FluidCyl){
	//The method calculates the T and p and the next step

	double Check = 100;// Variable to check the precision of the calculation at a determined iteration
	double Precision = 0.0001; //Limit to exit from the iterative loop 1 por mil de d_U
	double EnergyError = 0; //Unbalance of the energy equation in [J]
	double Tmed;//Mean temperatrure [K]
	double delta;
	double PTemperature_aux;

	dFQL = dFQLp;

	PVolume = Volume; //Chamber volume in [m3] at the current step
	Rcil = FluidCyl->FunR();
	//Rcil = 286.98;

	double Mmed = (PMass0 + PMass) / 2.;//Average mass in [K]

	do{

		Tmed = (PTemperature0 + PTemperature) / 2.;
		PPressureStep = (PPressure0 + PPressure) / 2.;

		delta = 0.;
		EnthalpyIn = 0.; //Term with the difference (he+ve2/2-ue) 
		EnthalpyOut = 0.;//Term due to the outlet flows
		//Flow terms are added. Inlet mass is positive, outlet negative.
		for (int i = 0; i < MassIn->size(); i++){ //En Calmec: [1] MBB, [2] GastoAdm*deltat, [3] GastoEsc*deltat
			if (MassIn->at(i)> 0){//Flujo entrante 
				delta = getDeltaFlujoEntrante(EnthIn->at(i), *FluidIn->at(i), MassIn->at(i), Tmed);
				EnthalpyIn += delta; //Terms addition
			}
			else{//Flujo saliente
				//delta=-1*Rcil * Tmed * MassIn->at(i);//Energy term due to the gas leaving  the chamber is calculated at chamber temperature
				delta = getDeltaFlujoSaliente(Rcil, MassIn->at(i), Tmed);
				EnthalpyOut += delta;//Terms addition
			}
		}

		//PHeatTransfer->Calcula_Calor_Woschni(Tmed, PPressureStep);
		//dQconv = PHeatTransfer->getQW(); //Convective HT [J]
		//dQrad = PHeatTransfer->getQRAD();//Recupera valor del objeto TransCalor // Radiative HT  [J] 

		dWi = PPressureStep * (PVolume - PVolume0); //Work [J]

		Cv = FluidCyl->FunCv(Tmed);//Specific heat 

		dUsen = Mmed * Cv *(PTemperature - PTemperature0);// Variation of sensible specific internal energy [J] (first time =0) //REVISAR

		EnergyError = dUsen + dWi + dQconv + dQrad - dFQL + EnthalpyIn + EnthalpyOut;

		PTemperature_aux = PTemperature - EnergyError / Cv / Mmed;//New temperature

		PTemperature = PTemperature_aux;
		PPressure = Rcil * PMass * PTemperature / PVolume;//New pressure
		
		if (dUsen == 0)
			Check = 1;
		else Check = fabs(EnergyError / dUsen);//Error a comprobar

	} while (Check > Precision);
}

void T0DModel::SolveNewPressure(std::vector<double> *MassIn, std::vector<double> *EnthIn, TFluidArray_ptr FluidIn, TFluid_ptr FluidCyl){
	//The method calculates the T and p and the next step

	double Check = 100;// Variable to check the precision of the calculation at a determined iteration
	double Precision = 0.0001; //Limit to exit from the iterative loop 1 por mil de d_U
	double EnergyError = 0; //Unbalance of the energy equation in [J]
	double Tmed;//Mean temperatrure [K]
	double delta;
	double PTemperature_aux;


	Rcil = FluidCyl->FunR();
	//Rcil = 286.98;

	double Mmed = (PMass0 + PMass) / 2.;//Average mass in [K]

	do{

		Tmed = (PTemperature0 + PTemperature) / 2.;

		EnthalpyIn = 0.; //Term with the difference (he+ve2/2-ue) 
		EnthalpyOut = 0.;//Term due to the outlet flows
		//Flow terms are added. Inlet mass is positive, outlet negative.
		for (int i = 0; i < MassIn->size(); i++){ //En Calmec: [1] MBB, [2] GastoAdm*deltat, [3] GastoEsc*deltat
			if (MassIn->at(i)> 0){//Flujo entrante 
				delta = getDeltaFlujoEntrante(EnthIn->at(i), *FluidIn->at(i), MassIn->at(i), Tmed);
				EnthalpyIn += delta; //Terms addition
			}
			else{//Flujo saliente
				//delta=-1*Rcil * Tmed * MassIn->at(i);//Energy term due to the gas leaving  the chamber is calculated at chamber temperature
				delta = getDeltaFlujoSaliente(Rcil, MassIn->at(i), Tmed);
				EnthalpyOut += delta;//Terms addition
			}
		}

		Cv = FluidCyl->FunCv(Tmed);//Specific heat 

		dUsen = Mmed * Cv *(PTemperature - PTemperature0);// Variation of sensible specific internal energy [J] (first time =0) //REVISAR

		EnergyError = dUsen + EnthalpyIn + EnthalpyOut;

		PTemperature -= EnergyError / Cv / Mmed;//New temperature

		if (dUsen == 0)
			Check = fabs(EnergyError);
		else Check = fabs(EnergyError / dUsen);//Error a comprobar

	} while (Check > Precision);

	PPressure = Rcil * PMass * PTemperature / PVolume;//New pressure
	PPressureStep = (PPressure0 + PPressure) / 2.;

}


///Devuelve la FQL para ciclo cerrado
//Si está puesta la opción de ciclo completo también calculará la FQL en ciclo abierto
void T0DModel::SolveNewFQL(double Volume, double Pressure, std::vector<double> *MassIn, std::vector<double> *EnthIn, std::vector<TFluid*>*FluidIn, TFluid *FluidCyl){

	double Pressure0 = PPressure;//Pressure in [Pa] at the previous step
	double Temperature0 = PTemperature;//Chamber temperature in [K] at the previous step
	double Volume0 = PVolume;//Chamber volume in [m3] at the previous step
	//double dQconv;//Convective HT in [J]
	//double dQrad; //RadiativHT in [J]
	double Tmed;//Mean temperatrure [K]
	double Pmed;//Mean pressure [Pa]
	double Mmed;//Mean mass [kg]
	double delta; //Auxiliar variable

	PPressure = Pressure;
	PVolume = Volume; //Chamber volume in [m3] at the current step
	Rcil = FluidCyl->FunR();//R of the cylinder gas



	PTemperature = PPressure * PVolume / PMass / Rcil;//New temperature is calculated

	if (PTemperature < 0)
	{
		PTemperature = PTemperature;
	}

	Pmed = (Pressure0 + PPressure) / 2.;
	Tmed = (Temperature0 + PTemperature) / 2.;
	Mmed = (PMass0 + PMass) / 2.;


	delta = 0.;
	EnthalpyIn = 0.; //Term with the difference (he+ve2/2-ue)
	EnthalpyOut = 0.;//Term due to the outlet flows
	//Flow terms are added. Inlet mass is positive, outlet negative.
	for (int i = 0; i < MassIn->size(); i++){// En Calmec: [1] MBB, [2] MFEVAP, [3] MAINY
		if (MassIn->at(i) > 0){//Flujo entrante 
			//delta=-1*(EnthIn->at(i) - FluidIn->at(i)->FunU(Tmed)) * MassIn->at(i);//Energy term due to the gas entering to the chamber. Inlet enthalpy is an input but U is calculated here
			delta = getDeltaFlujoEntrante(EnthIn->at(i), *FluidIn->at(i), MassIn->at(i), Tmed);
			EnthalpyIn += delta; //Terms addition
		}
		else{//Flujo saliente
			//delta=-1*Rcil * Tmed * MassIn->at(i);//Energy term due to the gas leaving  the chamber is calculated at chamber temperature
			delta = getDeltaFlujoSaliente( Rcil, MassIn->at(i), Tmed);
			EnthalpyOut += delta;//Terms addition
		}
	}

	PHeatTransfer->Heat(Tmed, Pmed);
	dQconv = PHeatTransfer->ConvectiveHeat(); //Convective HT [J]
	dQrad = PHeatTransfer->RadiativeHeat();//Recupera valor del objeto TransCalor // Radiative HT  [J] JUANMA

	dWi = Pmed * (PVolume - Volume0); //Work [J]

	Cv = FluidCyl->FunCv(Tmed);//Specific heat 

	dUsen = Mmed * Cv *(PTemperature - Temperature0);// Variation of sensible specific internal energy [J]

	dFQL = dUsen + dWi + dQconv + dQrad + EnthalpyIn + EnthalpyOut; //Heat release  [J] 

}

double T0DModel::getDeltaFlujoEntrante(double EnthIn, TFluid FluidIn, double MassIn, double Tmed)
{
	return - (EnthIn - FluidIn.FunCv(Tmed) * Tmed) * MassIn;
}


double T0DModel::getDeltaFlujoSaliente(double R, double MassIn, double Tmed)
{
	return -R * Tmed * MassIn;
}


