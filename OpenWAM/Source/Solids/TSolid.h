/**
 * @file TSolid.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 9 de mar. de 2016
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
 * Class for the calculation of any generic solid material.
 */
#ifndef SOURCE_SOLIDS_TSOLID_H_
#define SOURCE_SOLIDS_TSOLID_H_


#include "Math_wam.h"
#include <string>

/*!
 * \class	TSolid
 *
 * \brief	Material: Generic material user defined.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	11/03/2016
 */

class TSolid {

protected:
	double Conductivity;			   //!< The conductivity [W / (m K)]
	double HeatCapacity;			   //!< The heat capacity [J / (kg K)]
	double Density;					   //!< The density [kg / m ** 3]
	double Elasticity;				   //!< The elasticity [N / m]
	ArrayXd CoefCond;				   //!< The coefficients for the conductivity correlation
	ArrayXd ExpCond;				   //!< The temperature exponents for the conductivity correlation
	ArrayXd CoefHeCap;				   //!< The coefficients for the heat capacity correlation
	ArrayXd ExpHeCap;				   //!< The temperature exponents for the heat capacity correlation
	ArrayXd CoefDens;				   //!< The coefficients for the density correlation
	ArrayXd ExpDens;				   //!< The temperature exponents for the density correlation
	ArrayXd CoefElas;				   //!< The coefficients for the elasticity correlation
	ArrayXd ExpElas;				   //!< The temperature exponents for the elasticity correlation
	string Name;					   //!< The name of the material
public:

	/*!
	 * \fn	TSolid();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 */

	TSolid();

	/*!
	 * \fn	TSolid(const TSolid* obj);
	 *
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	obj	The original object.
	 */

	TSolid(const TSolid* obj);

	/*!
	 * \fn	TSolid(string nme);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	nme	The name of the material.
	 */

	TSolid(string nme);

	/*!
	 * \fn	virtual ~TSolid();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 */

	virtual ~TSolid();

	/*!
	 * \fn	double getConductivity() const
	 *
	 * \brief	Gets the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \return	The conductivity [W / (m K)].
	 */

	double getConductivity() const {
		return Conductivity;
	}

	/*!
	 * \fn	double getDensity() const
	 *
	 * \brief	Gets the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \return	The density [kg / m ** 3].
	 */

	double getDensity() const {
		return Density;
	}

	/*!
	 * \fn	double getElasticity() const
	 *
	 * \brief	Gets the elasticity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \return	The elasticity [N / m].
	 */

	double getElasticity() const {
		return Elasticity;
	}

	/*!
	 * \fn	double getHeatCapacity() const
	 *
	 * \brief	Gets the heat capacity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \return	The heat capacity [J / (kg K)].
	 */

	double getHeatCapacity() const {
		return HeatCapacity;
	}

	/*!
	 * \fn	virtual double getConductivity(double T);
	 *
	 * \brief	Gets the conductivity as a function of the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	T	The temperature [K].
	 *
	 * \return	The conductivity [W / (m K)].
	 */

	virtual double getConductivity(double T);

	/*!
	 * \fn	virtual double getDensity(double T);
	 *
	 * \brief	Gets the density as a function of the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	T	The temperature [K].
	 *
	 * \return	The density [kg / m ** 3].
	 */

	virtual double getDensity(double T);

	/*!
	 * \fn	virtual double getElasticity(double T);
	 *
	 * \brief	Gets the elasticity as a function of the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	T	The temperature [K].
	 *
	 * \return	The elasticity [N / m].
	 */

	virtual double getElasticity(double T);

	/*!
	 * \fn	virtual double getHeatCapacity(double T);
	 *
	 * \brief	Gets the heat capacity as a function of the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	T	The temperature [K].
	 *
	 * \return	The heat capacity [J / (kg K)].
	 */

	virtual double getHeatCapacity(double T);

	/*!
	 * \fn	void setCoefCond(const ArrayXd& coefCond)
	 *
	 * \brief	Sets the coefficients for the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	coefCond	The array of coefficients.
	 */

	void setCoefCond(const ArrayXd& coefCond) {
		CoefCond = coefCond;
		ExpCond.setLinSpaced(CoefCond.size(), 0, CoefCond.size() - 1);
	}

	/*!
	 * \fn	void setCoefDens(const ArrayXd& coefDens)
	 *
	 * \brief	Sets the coefficients for the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	coefDens	The array of coefficients.
	 */

	void setCoefDens(const ArrayXd& coefDens) {
		CoefDens = coefDens;
		ExpDens.setLinSpaced(CoefDens.size(), 0, CoefDens.size() - 1);
	}

	/*!
	 * \fn	void setCoefElas(const ArrayXd& coefElas)
	 *
	 * \brief	Sets the coefficients for the elasticity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	coefElas	The array of coefficients.
	 */

	void setCoefElas(const ArrayXd& coefElas) {
		CoefElas = coefElas;
		ExpElas.setLinSpaced(CoefElas.size(), 0, CoefElas.size() - 1);
	}

	/*!
	 * \fn	void setCoefHeCap(const ArrayXd& coefHeCap)
	 *
	 * \brief	Sets the coefficients for the heat capacity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	coefHeCap	The array of coefficients.
	 */

	void setCoefHeCap(const ArrayXd& coefHeCap) {
		CoefHeCap = coefHeCap;
		ExpHeCap.setLinSpaced(CoefHeCap.size(), 0, CoefHeCap.size() - 1);
	}

	/*!
	 * \fn	void setTemperature(double T);
	 *
	 * \brief	Sets the perties for a temperature value.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 *
	 * \param	T	The temperature [K].
	 */

	void setTemperature(double T);
};

#endif /* SOURCE_SOLIDS_TSOLID_H_ */
