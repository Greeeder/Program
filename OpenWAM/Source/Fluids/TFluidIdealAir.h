/**
 * @file TFluidIdealAir.h
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
 * The TFluidIdealAir class is used for a fluid composed only by ideal air.
 *
 * It represents a fluid composed by ideal air. This file includes the
 * declaration of such class.
 *
 */

#include "TChemicalComponent.h"
#include "TFluid.h"
#include "TChIdealAir.h"
#include <vector>

#ifndef __TFLUIDIDEALAIR_H
#define __TFLUIDIDEALAIR_H

/*!
 * \class	TFluidIdealAir
 *
 * \brief	Fluid composed by ideal air.
 *
 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
 * \date	30/05/2015
 */
class TFluidIdealAir : public TFluid {

  private:

	TComponentArray_obj comp;	//!< The chemila components array (only air)

  public:

	/*!
	 *
	 * \brief	Default constructor.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2015
	 */
	TFluidIdealAir();

	/*!
	 *
	 * \brief	Destructor.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2015
	 */
	virtual ~TFluidIdealAir();

	/*!
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat capacity at constant pressure [J / (kg K)].
	 */
	virtual double FunCp(double T);

	/*!
	 *
	 * \brief	Function for the internal energy.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	30/05/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */
	double FunU(double T);

	/*!
	 *
	 * \brief	Function for density.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Density [kg / m ** 3].
	 */
	virtual double FunRHO(double T);

	/*!
	 *
	 * \brief	Function for the specific heat a constant volume.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / (kg K)].
	 */
	virtual double FunCv(double T);

	/*!
	 *
	 * \brief	Function for the gas constant.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Gas constant [J / (kg K)].
	 */
	virtual double FunR();

	/*!
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
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
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Thermal conductivity [W / (m K)].
	 */
	virtual double Funk(double T);

	/*!
	 *
	 * \brief	Function for the dynamic viscosity.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Viscosity [Pa s].
	 */
	virtual double FunVisc(double T);

	/*!
	 *
	 * \brief	Function for the Prandtl number.
	 *
	 * \author	L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Prandtl number [-].
	 */
	virtual double FunPr(double T);
};

#endif
