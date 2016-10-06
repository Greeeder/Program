/**
 * @file TPipeFrictionColebrook.h
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

#ifndef TPipeFrictionColebrook_H_
#define TPipeFrictionColebrook_H_

#include "Globales.h"
#include "TFluid.h"
#include "TPipeFriction.h"

class TPipeFrictionColebrook: public TPipeFriction {

  private:


  public:

	/*!
	 * \fn	TPipeFrictionColebrook();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	TPipeFrictionColebrook();

	/*!
	 * \fn	virtual ~TPipeFrictionColebrook();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	08/02/2016
	 */

	virtual ~TPipeFrictionColebrook();

	/*!
	 * \fn	virtual ArrayXd FrictionCoefficient(ArrayXd Re)
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

	virtual ArrayXd FrictionCoefficient(const ArrayXd& Re);

	virtual void setRelativeRugosity(const RowVector& D);
};

#endif /* TPipeFrictionColebrook_H_ */
