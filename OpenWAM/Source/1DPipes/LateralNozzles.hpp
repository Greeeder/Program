/* --------------------------------------------------------------------------------*\
 ==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
  \\ |  X  | //  W ave     |
   \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
    \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*--------------------------------------------------------------------------------*/

/*!
 * \file LateralNozzles.hpp
 * \author Francisco Jose Arnau <farnau@mot.upv.es>
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
 * along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \section DESCRIPTION
 * This file declares a lateral nozzles source term for a TPipe object.
 * It also includes helper functions.
 */

#ifndef TLateralNozzles_hpp
#define TLateralNozzles_hpp

#include <memory>
#include "TGenericSource.hpp"

/*!
 * \brief A lateral nozzles source. It produces flow through a lateral window.
 */
class TLateralNozzles: public TGenericSource {
  protected:
	RowArray FArea; //!< Nozzle effective area.
	double FOutletPressure; //!< Outlet pressure. [Pa]
	double FOutletTemperature; //!< Outlet temperature. [K]

  public:
	/*!
	 * \brief Creates a lateral nozzles source object.
	 * 
	 * \param pipe The pipe affected by this source.
	 */
	TLateralNozzles(const Pipe_ptr & pipe);

	/*!
	 * \brief Computes the source terms.
	 * 
	 * \param t Current time. [s]
	 * \param dt Time-step. [s]
	 * \return The source terms. 
	 */
	virtual RowArray ComputeSource(double t, double dt);

	/*!
	 * \brief Gets the nozzle effective areas.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Effective areas. [m ** 2]
	 */
	RowVector getArea() const;

	/*!
	 * \brief Gets the effective area for a nozzle.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param i Nozzle number.
	 * 
	 * \return Effective area. [m ** 2]
	 */
	double getArea(unsigned int i) const;

	/*!
	 * \brief Sets the nozzle effective areas.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param area Effective areas. [m ** 2]
	 */
	void setArea(const RowVector & area);

	/*!
	 * \brief Sets the effective area for a nozzle.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param area Effective area. [m ** 2]
	 * \param i Nozzle number.
	 */
	void setArea(double area, unsigned int i);

	/*!
	 * \brief Sets the same effective area for all the nozzles.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param area Effective area. [m ** 2]
	 */
	void setArea(double area);

	/*!
	 * \brief Sets the nozzle outlet plenum pressure and temperature.
	 * 
	 * \param p Outlet plenum pressure. [Pa]
	 * \param T Outlet plenum temperature. [K]
	 */
	void setOutletConditions(double p, double T);
};

/*!
 * \brief Computes the mass flow rate through a nozzle.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * \date 2016/05/04
 * 
 * \param p_1 Nozzle inlet total pressure. [Pa]
 * \param p_2 Nozzle outlet pressure. [Pa]
 * \param T_1 Nozzle inlet total temperature. [K]
 * \param T_2 Nozzle outlet temperature. [K]
 * \param R Gas constant. [J / (kg * K)]
 * \param cp Specific heat capacity at constant pressure. [J / (kg / K)]
 * \param g Specific heat capacities ratio.
 * \param A Nozzle area. [m ** 2]
 */
double nozzle_mass_flow(double p_1, double p_2, double T_1, double T_2,
	double R, double cp, double g, double A);

/*!
 * \brief Sets a TLateralNozzles object to a pipe.
 * 
 * \param pipe Pipe.
 */
void set_lateral_nozzles(const Pipe_ptr & pipe);

#endif
