/**
 * @file TPipeIntHeat.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 feb. 2016
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
 * This class calculates the heat transfer between the gas and the wall.
 */

#ifndef TPIPEINTHEAT_H_
#define TPIPEINTHEAT_H_

#include <memory>
#include "Globales.h"
#include "TFluid.h"

class TPipeIntHeat {

  protected:

	ArrayXd h;						   //!< Heat transfer coefficient [W / (m^2 K)]
	ArrayXd Cond;					   //!< Internal fluid conductivity [W / (m K)]
	ArrayXd Dint;					   //!< Internal diameter [m]
	ArrayXd Aint;					   //!< Internal area [m^2]
	ArrayXd Visc;					   //!< Gas viscosity [Pa s]
	ArrayXd ViscWall;				   //!< Gas viscosity as wall temperature [Pa s]

	ArrayXd Qint;					   //!< Internal heat transfer [W]

	double HeatMultiplier;			   //!< Internal heat transfer multiplier [-]

	int ncells;						   //!< Number of cells of the pipe
	double MeshSize;				   //!< Cell length [m]

	ArrayXd Turbulence;				   //!< Turbulence coefficient [-]

  public:

	/*!
	 * \fn	TPipeIntHeat();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 * 			
	 * \param	HeatMult		Internal heat multiplier;
	 * \param	D				Internal diameter of each cell [m]
	 * \param	cellsize		Length of the cells [m]
	 */

	TPipeIntHeat(double HeatMult, VectorXd D, double cellsize);

	/*!
	 * \fn	virtual ~TPipeIntHeat();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	virtual ~TPipeIntHeat();

	/*!
	 * \fn	virtual ArrayXd Heat(ArrayXd Tgas, ArrayXd Twall, ArrayXd Re, std::vector<TFluid *> Gas)
	 *
	 * \brief	Heats.
	 * 			
	 *	Heat is positive when the gas losses energy.
	 *	
	 * \f[ \dot{Q} = h \cdot A \cdot \left( T_{gas} - T_{wall} \right) \f]
	 * 
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 *
	 * \param	Tgas	   	The gas temperature [K].
	 * \param	Twall	   	The wall temperature [K].
	 * \param	Re		   	The Reynolds number [-].
	 * \param [in,out]	Gas	The object that includes the gas properties.
	 *
	 * \return	.
	 */

	virtual ArrayXd Heat(const ArrayXd &Tgas, const ArrayXd &Twall, const ArrayXd &Re, const std::vector<TFluid_ptr>& Gas) {

		return ArrayXd::Zero(ncells);

	}
};

typedef unique_ptr<TPipeIntHeat> TPipeIntHeat_ptr;

#endif /* TPIPEINTHEAT_H_ */
