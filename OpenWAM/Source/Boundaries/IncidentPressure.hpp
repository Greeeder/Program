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
 * \file IncidentPressure.hpp
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
 * This file declares a boundary condition that imposes an incident pressure
 * wave.
 */

#ifndef IncidentPressureBC_hpp
#define IncidentPressureBC_hpp

#include "ConstantConditionsBC.hpp"
#include "CSVLoader.hpp"

/*!
 * \brief A BC that imposes an incident pressure wave.
 */
class TIncidentPressureBC: public TConstantConditionsBC {

	friend void attach_to_incident_pressure_BC(const Pipe_ptr& pipe,
		nmPipeEnd pipe_end, const RowVector& t, const RowVector& p,
		const RowVector& A_A, TFluid_ptr fluid);

  protected:
	/*!
	 * \brief Default constructor.
	 */
	TIncidentPressureBC();

	/*!
	 * \brief Constructor, only used from friend functions.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param pipe Pipe attached to this BC.
	 * \param pipe_end Pipe end connected to this BC.
	 * \param pipe_cc Constant conditions pipe attached to this BC.
	 * \param virtual_pipe Virtual pipe used by this BC.
	 */
	TIncidentPressureBC(const Pipe_ptr & pipe, nmPipeEnd pipe_end,
		const ConstantConditionsPipe_ptr & pipe_cc,
		const VirtualPipe_ptr & virtual_pipe);

	RowVector FTimeVector; //!< Time vector. [s]
	RowVector FIncidentPressure; //!< Incident pressure vector. [Pa]
	RowVector FA_A; //!< Entropy level vector. [-]

  public:
	/*!
	 * \brief Computes the boundary condition flux vector.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/09
	 * 
	 * \param t Time for which the flux is going to be computed. [s]
	 * \param dt Time-step. [s]
	 * \return The BC flux vector.
	 */
	virtual ColVector Flux(double t, double dt);

	/*!
	 * \brief Sets the pressure and entropy level vectors.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param t Time vector. [s]
	 * \param p Incident pressure. [Pa]
	 * \param A_A Entropy level. [-]
	 */
	void setPA_A(const RowVector& t, const RowVector& p, const RowVector& A_A);
};

/*!
 * \brief Unique pointer to a TIncidentPressureBC object.
 */
typedef unique_ptr<TIncidentPressureBC> IncidentPressureBC_ptr;

/*!
 * \brief Joins a pipe to an incident pressure BC.
 *
 * \author Luis Miguel Garcia-Cuevas González <luiga12@mot.upv.es>
 * \date 2016/06/02
 * 
 * \param pipe Pipe connected to this boundary condition.
 * \param pipe_end Pipe end connected to the boundary condition.
 * \param t Time vector. [s]
 * \param p Incident pressure at the BC. [Pa]
 * \param A_A Entropy level. [K]
 * \param fluid Pointer to the fluid
 */
void attach_to_incident_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
	const RowVector& t, const RowVector& p, const RowVector& A_A, TFluid_ptr fluid);

/*!
 * \brief Loads incident pressure data from a CSV-like file.
 * 
 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
 * \date 2016/06/03
 * 
 * Loads data from a CSV-like file. The different data collumns are
 * separator by the given column separator. The labels for the time,
 * pressure and entropy level columns are also needed. The units
 * for the time and pressure columns are also used.
 * 
 * \param file_name File name.
 * \param separator Column separator.
 * \param t_label Time label.
 * \param p_label Pressure label.
 * \param A_A_label Entropy level lavel.
 * \param time_unit Time unit.
 * \param p_unit Pressure unit.
 * \param t Time vector, output. [s]
 * \param p Pressure vector, output. [Pa]
 * \param A_A Entropy level vector, output. [-]
 */
void load_incident_pressure_data(const std::string& file_name, char separator,
	const std::string& t_label, const std::string& p_label,
	const std::string& A_A_label, const std::string& t_unit,
	const std::string& p_unit, RowVector* t, RowVector* p, RowVector* A_A);

#endif
