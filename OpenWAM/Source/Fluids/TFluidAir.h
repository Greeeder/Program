/**
 * @file TFluidAir.h
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
 * The TFluidAir class is used for a fluid composed only by air.
 *
 * It represent a fluid composed by air.
 *
 */

//#include "stdafx.h"
#include "TChemicalComponent.h"
#include "TFluid.h"
#include "TChAir.h"
#include <vector>

#ifndef __TFLUIDAIR_H
#define __TFLUIDAIR_H

/*!
 * \class	TFluidAir
 *
 * \brief	Fluid composed by air.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	05/12/2015
 */

class TFluidAir : public TFluid {

  private:

	TComponentArray_obj comp;	//!< The chemila components array (only air)

  public:

	/*!
	 * \fn	TFluidAir();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	TFluidAir();

	/*!
	 * \fn	~TFluidAir();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	virtual ~TFluidAir();

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

	double FunU(double T);

	/*!
	 * \fn	double FunRHO(double T);
	 *
	 * \brief	Function for density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Density [kg / m ** 3].
	 */

	virtual double FunRHO(double T);

	/*!
	 * \fn	double FunCv(double T);
	 *
	 * \brief	Function for the specific heat a constant volume.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / (kg K)].
	 */

	virtual double FunCv(double T);

	/*!
	 * \fn	double FunR();
	 *
	 * \brief	Function for the gas constant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Gas constant [J / (kg K)].
	 */

	virtual double FunR();

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
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
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Conductivity [W / (m K)].
	 */

	virtual double Funk(double T);

	/*!
	 * \fn	double FunVisc(double T);
	 *
	 * \brief	Function for the dynamic viscosity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Viscosity [Pa s].
	 */

	virtual double FunVisc(double T);

	/*!
	 * \fn	double FunPr(double T);
	 *
	 * \brief	Function for the Prandtl number.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Prandtl number [-].
	 */

	virtual double FunPr(double T);


};



#endif
