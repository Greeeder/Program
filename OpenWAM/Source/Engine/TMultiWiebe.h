/**
 * @file TMultiWiebe.h
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
 * Calculates heat released by means of a multi-Wiebe function.
 */
#ifndef SOURCE_ENGINE_TMULTIWIEBE_H_
#define SOURCE_ENGINE_TMULTIWIEBE_H_

#include "TCombustion.h"
#include "TWiebe.h"

/*!
 * \class	TMultiWiebe
 *
 * \brief	Calculates heat released by means of a multi-Wiebe function.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	14/03/2016
 */

class TMultiWiebe: public TCombustion {

private:

	std::vector<TWiebe *> Wiebe;	   //!< The array of Wiebe function
	std::vector<double> Beta;		   //!< The weight of the Wiebes

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	TMultiWiebe();

	/*!
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param [in,out]	obj	If non-null, the object.
	 */

	TMultiWiebe(TMultiWiebe * obj);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	virtual ~TMultiWiebe();

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
	 * \param	node_comb	The node comb.
	 */

	virtual void ReadCombustionData(xml_node node_comb);
};

#endif /* SOURCE_ENGINE_TMULTIWIEBE_H_ */
