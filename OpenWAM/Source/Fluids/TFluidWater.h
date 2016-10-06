/**
 * @file TFluidWater.h
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
 * The TFluidWater class is used for a fluid composed only for liquid water.
 *
 * It represent a fluid composed by water.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"
#include "TFluid.h"
#include <vector>

#ifndef __TFLUIDWATER_H
#define __TFLUIDWATER_H

/*!
 * \class	TFluidWater
 *
 * \brief	Fluid composed by water.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	05/12/2015
 */

class TFluidWater : public TFluid {

  private:

	TComponentArray_obj comp;	//!< The chemical components that compose the fluid (only water)

  public:

	/*!
	 * \fn	TFluidWater();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	TFluidWater();

	/*!
	 * \fn	~TFluidWater();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	~TFluidWater();

	/*!
	 * \fn	double FunCp(double T);
	 *
	 * \brief	Function for the specific heat at constant pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / (kg K)].
	 */

	virtual double FunCp(double T);

	/*!
	 * \fn	double FunU(double T);
	 *
	 * \brief	Function for the internal energy.
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
	 * \fn	double FunRHO(double T);
	 *
	 * \brief	Function for the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Density [kg / m ** 3].
	 */

	virtual double FunRHO(double T);

	/*!
	 * \fn	double FunCv(double T);
	 *
	 * \brief	Function for the specific heat at constant volume.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / kg].
	 */

	virtual double FunCv(double T);

	/*!
	 * \fn	double FunR();
	 *
	 * \brief	Gets the gas constant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \return	Gas constant [J / kg].
	 */

	virtual double FunR();

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T);

	/*!
	 * \fn	double Funk(double T);
	 *
	 * \brief	Function for the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Conductivity [W / (m K)].
	 */

	virtual double Funk(double T);

	/*!
	 * \fn	double FunVisc(double T);
	 *
	 * \brief	Function for the viscosity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Viscosity [Pa s].
	 */

	virtual double FunVisc(double T);

	/*!
	 * \fn	double FunPr(double T);
	 *
	 * \brief	Function for the Prandtl number
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Prandtl number [-].
	 */

	virtual double FunPr(double T);


};



#endif
