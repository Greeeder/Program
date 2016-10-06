/**
 * @file TSldAl.h
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
 * Class for the calculation of the Aluminium properties.
 */
#ifndef SOURCE_SOLIDS_TSLDAL_H_
#define SOURCE_SOLIDS_TSLDAL_H_

#include "TSolid.h"

/*!
 * \class	TSldAl
 *
 * \brief	Material: Aluminium.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	11/03/2016
 */

class TSldAl: public TSolid {

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 */

	TSldAl();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	11/03/2016
	 */

	virtual ~TSldAl();
};

#endif /* SOURCE_SOLIDS_TSLDAL_H_ */
