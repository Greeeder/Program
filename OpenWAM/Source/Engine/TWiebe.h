/**
 * @file TWiebe.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 11 de mar. de 2016
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
 * Calculates the heat released by means of a Wiebe function.
 */
#ifndef SOURCE_ENGINE_TWIEBE_H_
#define SOURCE_ENGINE_TWIEBE_H_

#include "TCombustion.h"

/*!
 * \brief	Calculates the heat released by means of a Wiebe function.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	14/03/2016
 */

class TWiebe: public TCombustion {

private:

	double Duration;				   //!< Combustion duration
	double SOC;						   //!< Start of combustion
	double m;						   //!< Form parameter
	double C;						   //!< Completeness parameter

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	TWiebe();

	/*!
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param [in,out]	obj	If non-null, the object.
	 */

	TWiebe(TWiebe* obj);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	virtual ~TWiebe();

	/*!
	 * \brief	Gets the heat released.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param	angle	The angle [deg].
	 *
	 * \return	The heat released [-].
	 */

	virtual double getHRL(double angle);

	/*!
	 * \brief	Reads the combustion data.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param	node_comb	The combustion xml node.
	 */

	virtual void ReadCombustionData(xml_node node_comb);
};

#endif /* SOURCE_ENGINE_TWIEBE_H_ */
