/**
 * @file TChGasolina.h
 * @author F.J. Arnau <farnau@mot.upv.es>
 * @version 0.1.0
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
 * The TChGasolina class represents the chemical component gasoline.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the Gasoline.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TCHGASOLINA_H
#define __TCHGASOLINA_H

/*!
 * \class	TChGasolina
 *
 * \brief	Chemical component: Gasoline.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	05/12/2015
 */

class TChGasolina: public TChemicalComponent {
  private:

  public:

	/*!
	 * \fn	TChGasolina(std::string nombre, double CAHFL, double CBHFL, double RU, double PM,
	 * 		double HC, double OC, double HV, double DENSF, double PMA, double XO2);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	CAHFL 	The coefficient A for the calculation of the enthalpy [J / kg].
	 * \param	CBHFL 	The coefficient B for the calculation of the enthalpy [J / (kg K)].
	 * \param	RU	  	The universal gas constant [J / (mol K)].
	 * \param	PM	  	The molecular weigth [gr / mol].
	 * \param	HC	  	The hc.
	 * \param	OC	  	The oc.
	 * \param	HV	  	The hv.
	 * \param	DENSF 	The density [kg / m^3].
	 * \param	PMA   	The molecular weight of the air [gr / mol].
	 * \param	XO2   	The molar fraction of the oxigen in the air [-].
	 */

	TChGasolina(std::string nombre, double CAHFL, double CBHFL, double RU, double PM, double HC, double OC, double HV,
				double DENSF, double PMA, double XO2);

	/*!
	 * \fn	~TChGasolina();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	~TChGasolina();

	/*!
	 * \fn	double FunCp(double T);
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat C_p [J / (kg K)].
	 */

	virtual double FunCp(double T);

	/*!
	 * \fn	double FunCv(double T);
	 *
	 * \brief	Function for C_v.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat C_v [J / (kg K)].
	 */

	virtual double FunCv(double T);

	/*!
	 * \fn	double FunU(double T);
	 *
	 * \brief	Function for internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	virtual double FunU(double T);

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T);

};

#endif
