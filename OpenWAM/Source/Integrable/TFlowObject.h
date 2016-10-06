/**
 * @file T0DObject.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 11 de abr. de 2016
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
 * Virtual flow object.
 */
#ifndef SOURCE_ODMODELS_TFLOWOBJECT_H_
#define SOURCE_ODMODELS_TFLOWOBJECT_H_

#include "TFluid.h"
#include <memory>

/*!
 * \brief	Virtual flow object.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	10/05/2016
 */
class TFlowObject {
public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 */
	TFlowObject();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 */
	virtual ~TFlowObject();

	/*!
	* \brief Gets the mesh size.
	*
	* \return The mesh size. [m]
	*/
	virtual double getDeltaX() const = 0;

	/*!
	* \brief Gets the current time.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*
	* \return Current time. [s]
	*/
	virtual double getCurrentTime() const = 0;

	/*!
	 * \brief	Gets the pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \param	i	Index of the cell.
	 *
	 * \return	The pressure [Pa].
	 */
	virtual double getPressure(Uint i) const = 0;

	/*!
	 * \brief	Gets the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \param	i	Index of the cell.
	 *
	 * \return	The temperature [K].
	 */
	virtual double getTemperature(Uint i) const = 0;

	/*!
	 * \brief	Gets the speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \param	i	Index of the cell.
	 *
	 * \return	The speed [m/s].
	 */
	virtual double getSpeed(Uint i) const = 0;

	/*!
	 * \brief	Gets the area.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \return	The area [m ** 2].
	 */
	virtual RowVector getArea() const = 0;

	/*!
	 * \brief	Gets the number of cells.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \return	The number of cells.
	 */
	virtual Uint getNin() const = 0;

	/*!
	 * \brief	Gets the fluid object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	10/05/2016
	 *
	 * \param	i	Index of the cell.
	 *
	 * \return	The fluid.
	 */
	virtual TFluid_ptr getFluid(Uint i) const = 0;

	/*!
	 * \brief Returns the fluid components array.
	 * 
	 * \author L.M. García-Cuevas (luiga12@mot.upv.es)
	 * \date 24/05/2016
	 * 
	 * This function returns an array of pointers to the chemical components
	 * of the working fluid. It doesn't return the mass or molar fractions.
	 * 
	 * \return The fluid components array.
	 */
	virtual TComponentArray_ptr getFluidComponents() const = 0;
};

typedef std::shared_ptr<TFlowObject> TFlowObject_ptr;

#endif /* SOURCE_ODMODELS_T0DOBJECT_H_ */
