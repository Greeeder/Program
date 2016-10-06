/**
 * @file THRL.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 14 de mar. de 2016
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
 * Impose heat realeas by means of a 1D look-up table (Angle-HRL).
 */
#ifndef SOURCE_ENGINE_THRL_H_
#define SOURCE_ENGINE_THRL_H_

#include "TCombustion.h"

/*!
 * \brief	Imposes heat realeas by means of a 1D look-up table (Angle-HRL).
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	14/03/2016
 */

class THRL: public TCombustion {

private:

	dVector Angle;					   //!< The angle array [deg]
	dVector HeatRL;					   //!< The heat released arry

	double SOC;						   //!< The start of combustion

	Hermite_interp* interpHRL;		   //!< The interpolator object

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	THRL();

	/*!
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param [in,out]	obj	If non-null, the object.
	 */

	THRL(THRL* obj);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 */

	virtual ~THRL();

	/*!
	 * \brief	Gets a hrl.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param	angle	The angle [deg].
	 *
	 * \return	The current heat release.
	 */

	virtual double getHRL(double angle);

	/*!
	 * \brief	Reads combustion data.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	14/03/2016
	 *
	 * \param	node_comb	The combustion xml node.
	 */

	virtual void ReadCombustionData(xml_node node_comb);

};

#endif /* SOURCE_ENGINE_THRL_H_ */
