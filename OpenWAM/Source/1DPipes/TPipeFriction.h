/**
 * @file TPipeFriction.h
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
 * This class calculates the friction between the gas and the wall.
 */

#ifndef TPipeFriction_H_
#define TPipeFriction_H_

#include <memory>
#include "Globales.h"
#include "TFluid.h"

class TPipeFriction {

  protected:

	double FrictionMultiplier;		   //!< Friction multiplier [-]
	double Rugosity;				   //!< The rugosity of the inner surface of the pipe [m]
	ArrayXd RelativeRugosity;		   //!< The relative rugosity (r / 3.7 / D) [-]

	ArrayXd f;
	ArrayXd D;

	int ncells;						   //!< Number of cells of the pipe
	double MeshSize;				   //!< Cell length [m]

  public:

	/*!
	 * \fn	TPipeFriction();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	TPipeFriction();

	/*!
	 * \fn	virtual ~TPipeFriction();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	virtual ~TPipeFriction();

	/*!
	 *
	 * \brief	Heats.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 *
	 * \param	Re		   	The Reynolds number [-].
	 *
	 * \return	Friction coefficient [-].
	 */
	virtual ArrayXd FrictionCoefficient(const ArrayXd & Re);

	double getMultiplier(){
		return FrictionMultiplier;
	}

	void ReadFrictionData(xml_node node_pipe);

	virtual void setRelativeRugosity(const RowVector& D){};
};

typedef unique_ptr<TPipeFriction> TPipeFriction_ptr;

#endif /* TPipeFriction_H_ */
