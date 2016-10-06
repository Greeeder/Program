/**
 * @file TChIdealAir.h
 * @author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
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
 * The TChIdealAir class represents ideal air.
 *
 * This class declares different methods to determine the chemical and
 * thermodynamics properties of ideal air.
 *
 */

#include "Constantes.h"
#include "TChemicalComponent.h"

#ifndef __TCHIDEALAIR_H
#define __TCHIDEALAIR_H

/*!
 * \class	TChIdealAir
 *
 * \brief	Chemical component: Ideal air.
 *
 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
 * \date	05/12/2015
 */
class TChIdealAir: public TChemicalComponent {

  private:

  public:

	/*!
	 *
	 * \brief	Constructor.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	RU	  	Universal gas constant [J / (mol K)].
	 * \param	PMA   	The molecula weight [gr / mol].
	 */
	TChIdealAir(std::string nombre, double RU, double PMA);

	/*!
	 *
	 * \brief	Destructor.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 */
	virtual ~TChIdealAir();

	/*!
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat capacity at constant pressure. [J / (kg * K)]
	 */
	virtual double FunCp(double T);

	/*!
	 *
	 * \brief	Function for C_v.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat capacity at constant volume. [J / (kg * K)]
	 */
	virtual double FunCv(double T);

	/*!
	 *
	 * \brief	Function for the internal energy.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */
	virtual double FunU(double T);

	/*!
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */
	virtual double FunH(double T);

	/*!
	 *
	 * \brief	Function for the thermal conductivity.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2016
	 *
	 * \param	T	Temperature [T].
	 *
	 * \return	Thermal conductivity [W / (m K)].
	 */
	virtual double Funk(double T);

	/*!
	 *
	 * \brief	Function for the dynamic viscosity.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Dynamic viscosity [Pa s].
	 */
	virtual double FunVisc(double T);
};

#endif
