// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

/*!
 * \file T1DAcousticTurbine.h
 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
 * \author Francisco Jose Arnau Martinez <farnau@mot.upv.es>
 * \version 0.1
 *
 * \section LICENSE
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
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \section DESCRIPTION
 * The T1DAcousticTurbine class represents the acoustics of a classic
 * one-dimensional turbine.
 *
 * This file declares a TQ2DAcousticTurbine class.
 */

#ifndef T1DAcousticTurbineH
#define T1DAcousticTurbineH

#include <memory>
#include "TAcousticTurbine.h"
#include "TPipe.hpp"
#include "TVolute.hpp"
#include "ConstantConditionsBC.hpp"
#include "TQ2DAcousticTurbine.h"


/*!
 * \brief A 1D turbine, acoustic part.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 */
class T1DAcousticTurbine: public TQ2DAcousticTurbine {

	friend class T1DTurbine;

  public:

	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	T1DAcousticTurbine();

	/*!
	 * \brief TQ2DAcousticTurbine constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * It sets the inlet pipe, the volute and the outlet pipe.
	 * 
	 * \param inlet Inlet pipe.
	 * \param volute Volute.
	 * \param outlet Outlet pipe.
	 */
	T1DAcousticTurbine(Pipe_ptr inlet, Volute_ptr volute, Pipe_ptr outlet);

	/*!
	 * \brief Sets the effective area for one volute lateral nozzle.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param A Effective area. [m ** 2]
	 * \param i Cell number.
	 */
	virtual void setVoluteNozzleEffectiveArea(double A, unsigned int i);

	/*!
	 * \brief Sets the volute outlet conditions.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param p Volute outlet pressure. [Pa]
	 * \param T Volute outlet temperature. [K]
	 */
	virtual void setVoluteOutletConditions(double p, double T);

	/*!
	 * \brief Volute outlet mass flow.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet mass flow. [kg / s]
	 */
	virtual RowVector VoluteOutletMassFlow() const;

	/*!
	 * \brief Volute outlet pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet pressure. [Pa]
	 */
	virtual RowVector VoluteOutletp() const;

	/*!
	 * \brief Volute outlet total pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet total pressure. [Pa]
	 */
	virtual RowVector VoluteOutletp0() const;

	/*!
	 * \brief Volute outlet temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet temperature. [K]
	 */
	virtual RowVector VoluteOutletT() const;

	/*!
	 * \brief Volute outlet total temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Volute outlet total temperature. [K]
	 */
	virtual RowVector VoluteOutletT0() const;
};

/*!
 * \brief Shared pointer to a TQ2DAcousticTurbine object.
 */
typedef std::shared_ptr<T1DAcousticTurbine> OneDAcousticTurbine_ptr;

#endif
