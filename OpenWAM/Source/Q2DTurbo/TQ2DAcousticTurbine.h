// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

/*!
 * \file TQ2DAcousticTurbine.h
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
 * The TQ2DAcousticTurbine class represents the acoustics of a
 * quasi-two-dimensional turbine.
 *
 * This file declares a TQ2DAcousticTurbine class.
 */

#ifndef TQ2DAcousticTurbineH
#define TQ2DAcousticTurbineH

#include <memory>
#include "TAcousticTurbine.h"
#include "TPipe.hpp"
#include "TVolute.hpp"
#include "ConstantConditionsBC.hpp"


/*!
 * \brief A Q2D turbine, acoustic part.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 */
class TQ2DAcousticTurbine: public TAcousticTurbine, public TIntegrableGroup {

	friend class TQ2DTurbine;

  protected:

	Pipe_ptr FInlet; ///< Inlet pipe.
	Volute_ptr FVolute; ///< Volute.
	Pipe_ptr FOutlet; ///< Outlet pipe.

  public:

	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	TQ2DAcousticTurbine();

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
	TQ2DAcousticTurbine(Pipe_ptr inlet, Volute_ptr volute, Pipe_ptr outlet);

	/*!
	 * \brief Turbine inlet diameter.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 *
	 * \param i Inlet pipe number.
	 * \return The turbine inlet diameter for a given inlet. [m]
	 */
	virtual double DIn(int i) const;

	/*!
	 * \brief Turbine inlet diameter.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 *
	 * \return The turbine inlet diameter. [m]
	 */
	virtual double DInTot() const;

	/*!
	 * \brief Turbine outlet diameter.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return The turbine outlet diameter. [m]
	 */
	virtual double DOut() const;

	/*!
	 * \brief Return the turbine expansion ratio.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine total to static expansion ratio.
	 */
	virtual double ExpRatio() const;

	/*!
	 * \brief Gets the turbine inlet mass flow for one of its inlet ducts.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param i Inlet duct index.
	 * 
	 * \return The turbine inlet mass flow. [kg / s]
	 */
	virtual double MassIn(int i) const;

	/*!
	 * \brief Gets the turbine inlet mass flow.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return The turbine inlet mass flow. [kg / s]
	 */
	virtual double MassIn() const;

	/*!
	 * \brief Gets the turbine outlet mass flow.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return The turbine outlet mass flow. [kg / s]
	 */
	virtual double MassOut() const;

	/*!
	 * \brief Outlet pipe pressure, measured at the connection with the outlet nozzle.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Outlet pipe pressure. [Pa]
	 */
	virtual double OutletNozzlep() const;

	/*!
	 * \brief Returns the turbine inlet pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet pressure. [bar]
	 */
	virtual double P3() const;

	/*!
	 * \brief Returns the turbine inlet total pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet total pressure. [bar]
	 */
	virtual double P30() const;

	/*!
	 * \brief Returns the turbine inlet total pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param i Turbine inlet duct.
	 * \return Turbine inlet total pressure. [bar]
	 */
	virtual double P30(int i) const;

	/*!
	 * \brief Returns the turbine inlet pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet pressure. [bar]
	 */
	virtual double P4() const;

	/*!
	 * \brief Returns the turbine inlet total pressure.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet total pressure. [bar]
	 */
	virtual double P40() const;

	/*!
	 * \brief Gas constant at the turbine inlet.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Gas constant at the turbine inlet. [J / (kg * K)]
	 */
	double R() const;

	/*!
	 * \brief Returns the rotor inlet mass flow.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return Rotor inlet mass flow. [kg / s]
	 */
	double RotorInletMassFlow() const;

	/*!
	 * \brief Sets the rotor outlet nozzle effective area.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param A Effective area. [m ** 2]
	 */
	void setRotorOutletEffectiveArea(double A);

	/*!
	 * \brief Sets the rotor nozzle inlet conditions.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param p Rotor inlet pressure. [Pa]
	 * \param T Rotor inlet temperature. [K]
	 */
	void setRotorInletConditions(double p, double T);

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
	 * \brief Turbine inlet area.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return The turbine inlet area. [m ** 2]
	 */
	virtual double SIn() const;

	/*!
	 * \brief Turbine outlet area.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return The turbine outlet area. [m ** 2]
	 */
	virtual double SOut() const;

	/*!
	 * \brief Returns the turbine inlet temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet temperature. [K]
	 */
	virtual double T3() const;

	/*!
	 * \brief Returns the turbine inlet total temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine inlet total temperature. [K]
	 */
	virtual double T30() const;

	/*!
	 * \brief Returns the turbine inlet total temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param i Turbine inlet duct.
	 * \return Turbine inlet total temperature. [K]
	 */
	virtual double T30(int i) const;

	/*!
	 * \brief Returns the turbine outlet temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine outlet temperature. [K]
	 */
	virtual double T4() const;

	/*!
	 * \brief Returns the turbine outlet total temperature.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Turbine outlet total temperature. [K]
	 */
	virtual double T40() const;

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
typedef std::shared_ptr<TQ2DAcousticTurbine> Q2DAcousticTurbine_ptr;

#endif
