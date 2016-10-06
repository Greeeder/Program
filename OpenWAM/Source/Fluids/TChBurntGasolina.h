/**
 * @file TChBurntGasolina.h
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
 * The TChBurntGasolina class represents the chemical component burnt gas produced by a
 * combustion with gasoline.
 *
 * This class include different methods to determine the chemical and
 * thermodynamics properties of a burnt gas.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TCHBURNTGASOLINA_H
#define __TCHBURNTGASOLINA_H

/*!
 * \class	TChBurntGasolina
 *
 * \brief	Chemical component: Burnt gas by gasoline.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	07/12/2015
 */

class TChBurntGasolina: public TChemicalComponent {
  private:

  public:

	/*!
	 * \fn	TChBurntGasolina(std::string nombre, double RU, double PM);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	RU	  	The universal gas constant [J / (mol K)].
	 * \param	PM	  	The pmolecular weight [gr / mol].
	 */

	TChBurntGasolina(std::string nombre, double RU, double PM);

	/*!
	 * \fn	~TChBurntGasolina();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 */

	~TChBurntGasolina();

	/*!
	 * \fn	double FunCp(double T);
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
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
	 * \date	07/12/2015
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
	 * \date	07/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / (kg K)].
	 */

	virtual double FunU(double T);

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	T	The double to process.
	 *
	 * \return	Enthalpy [J /(kg K)].
	 */

	virtual double FunH(double T);

};

#endif
