/**
 * @file TPortIntHeatIntake.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 de feb. de 2016
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
 * Object to calculate internal heat transfer in intake pipes.
 */
#ifndef SOURCE_1DPIPES_TPortIntHeatIntake_H_
#define SOURCE_1DPIPES_TPortIntHeatIntake_H_

#include "TPipeIntHeat.h"

class TPortIntHeatIntake: public TPipeIntHeat {

  public:

	/*!
	 * \fn	TPortIntHeatIntake(double HeatMult);
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 * 
	 * \param	ncells		Number of cells of the pipe
	 * \param	HeatMult	Internal heat multiplier
	 * \param	D			Internal diameter of each cell [m]
	 * \param	cellsize	Length of the cells [m]
	 */

	TPortIntHeatIntake(int ncells, double HeatMult, VectorXd D, double cellsize);

	/*!
	 * \fn	virtual ~TPortIntHeatIntake();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	virtual ~TPortIntHeatIntake();

	/*!
	 * \fn	ArrayXd Heat(ArrayXd Tgas, ArrayXd Twall, ArrayXd Re, std::vector<TFluid *> Gas);
	 *
	 * \brief	Calculate internal heat transfer in ports.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 *
	 * \param	Tgas	   	The gas temperature [K].
	 * \param	Twall	   	The wall temperature [K].
	 * \param	Re		   	The Reynolds number [-].
	 * \param [in,out]	Gas	The object that contains the properties of the fluid.
	 *
	 * \return	.
	 */

	ArrayXd Heat(const ArrayXd &Tgas, const ArrayXd &Twall, const ArrayXd &Re, const std::vector<TFluid_ptr>& Gas);

};

#endif /* SOURCE_1DPIPES_TINTAKEPIPEINTHEAT_H_ */
