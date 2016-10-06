/**
 * @file TChBlend.h
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
 * The TChGasolina class represents the chemical component fuel blend.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of a fuel blend.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"
#include "TChGasoil.h"
#include "TChGasolina.h"
#include <string>

#ifndef __TCHBLEND_H
#define __TCHBLEND_H

/*!
 * \class	TChBlend
 *
 * \brief	Chemical component: Fuel blend.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	05/12/2015
 */

class TChBlend: public TChemicalComponent {

  private:
	TComponentArray_obj _Fuel;	//!< Vector that contain the different fuels of the blend.
	std::vector<double> _YFuel;//!< Mass fraction for the different fuels [-]

  public:

	/*!
	 * \fn	TChBlend(std::string nombre, std::vector<TChemicalComponent *> Fuel,
	 * 		std::vector<double> Fraction, double RU);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	nombre			The name of the fuel blend.
	 * \param [in,out]	Fuel	If non-null, vector with the different fuels.
	 * \param	Fraction		Mass fraction of the different fuels [-].
	 * \param	RU				The universal gas constant [J / (mol K)].
	 */

	TChBlend(std::string nombre, TComponentArray_obj Fuel, std::vector<double> Fraction, double RU);

	/*!
	 * \fn	~TChBlend();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

// 	~TChBlend();

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
	 * \return	Sepecific heat C_p [J / (kg K)].
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

// 	/*!
// 	 * \fn	double FunR();
// 	 *
// 	 * \brief	Gets the gas constant.
// 	 *
// 	 * \author	F.J. Arnau (farnau@mot.upv.es)
// 	 * \date	05/12/2015
// 	 *
// 	 * \return	Gas constant [J / (kg K)].
// 	 */
//
// 	virtual double FunR();

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

// 	/*!
// 	 * \fn	double FunPM();
// 	 *
// 	 * \brief	Gets the molecula weight.
// 	 *
// 	 * \author	F.J. Arnau (farnau@mot.upv.es)
// 	 * \date	05/12/2015
// 	 *
// 	 * \return	Molecular weight [gr / mol].
// 	 */
//
// 	virtual double FunPM();
//
// 	/*!
// 	 * \fn	double FunRelacionHC();
// 	 *
// 	 * \brief	Gets the relation H-C.
// 	 *
// 	 * \author	F.J. Arnau (farnau@mot.upv.es)
// 	 * \date	05/12/2015
// 	 *
// 	 * \return	Hidrogen to Carbon ratio [-].
// 	 */
//
// 	virtual double FunRelacionHC();
//
// 	/*!
// 	 * \fn	double FunRelacionOC();
// 	 *
// 	 * \brief	Gets the relation O-C
// 	 *
// 	 * \author	F.J. Arnau (farnau@mot.upv.es)
// 	 * \date	05/12/2015
// 	 *
// 	 * \return	Oxigen to Carbon ratio [-].
// 	 */
//
// 	virtual double FunRelacionOC();
};

#endif
