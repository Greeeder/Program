//#include "stdafx.h"
#include "TChemicalComponent.h"

 TChemicalComponent::TChemicalComponent() {
 	// TODO Auto-generated constructor stub
 	
 }

TChemicalComponent::~TChemicalComponent() {
	// TODO Auto-generated destructor stub
}

void TChemicalComponent::setName(std::string nombre) {
	_nombre = nombre;
}

std::string TChemicalComponent::getName() {
	return _nombre;
}

void TChemicalComponent::setCoeficients(string Parameter, dVector vcoef) {

	if(Parameter == "Viscosity") {
		//ViscCoef.resize(vcoef.size());
		ViscCoef = vcoef;
	}
}

