/**
 * @file TChN2.h
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
 * The TChN2 class represents the chemical component N2.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the N2.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TCHN2_H
#define __TCHN2_H

/*!
 * \class	TChN2
 *
 * \brief	Chemical component: N_2
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/12/2015
 */

class TChN2: public TChemicalComponent {

  private:

  public:

	/*!
	 * \fn	TChN2(std::string nombre, double RU, double PMA);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	RU	  	Universal gas constant [J / (mol K)].
	 * \param	PMA   	Molecular weight [gr / mol].
	 */

	TChN2(std::string nombre, double RU, double PMA);

	/*!
	 * \fn	~TChN2();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 */

	~TChN2();

	/*!
	 * \fn	double FunCp(double T);
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
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
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat C_v [J / (kg K)].
	 */

	virtual double FunCv(double T);

	/*!
	 * \fn	double FunU(double T);
	 *
	 * \brief	Function for the internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	virtual double FunU(double T);

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T);

	/*!
	* \brief	Function for the viscosity.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	03/02/2016
	*
	* \param	T	Temperature [K].
	*
	* \return	Viscosity [Pa s].
	*/

	virtual double FunVisc(double T);

};

#endif
