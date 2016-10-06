/*--------------------------------------------------------------------------------*\
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

/**
 * @file TQ2DTurbine.hpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
 * This file declares a quasi-two-dimensional turbine class.
 */

#ifndef TQ2DTurbineH
#define TQ2DTurbineH

#include <iostream>
#ifdef __BORLANDC__
#include <vcl.h>
#endif
#include "Constantes.h"
#include "TIntegrable.hpp"
#include "TTurbina.h"
#include "TQ2DAcousticTurbine.h"
#include "T1DAcousticTurbine.h"

/**
 * @brief Q2D turbine.
 */
class TQ2DTurbine: public TTurbina, public TIntegrableGroup {
  protected:
	Q2DAcousticTurbine_ptr FQ2DAcTurb; ///< Q2D acoustic turbine objetc.
	RowVector FBSR; ///< Blade speed ratio.
	RowVector Fcp_med; ///< Specific heat capacity at constant pressure. [J / (kg * K)]
	RowVector FExpansionRatio; ///< Total to static pressure ratio.
	RowVector FInletPressure; ///< Stator inlet pressure. [Pa]
	RowVector FInletTemperature; ///< Stator inlet temperature. [K]
	RowVector FInletTotalPressure; ///< Stator inlet total pressure. [Pa]
	RowVector FInletTotalTemperature; ///< Stator inlet total temperature. [K]
	RowVector FInletMassFlow; ///< Stator inlet mass flow. [kg / s]
	RowVector FOutletIsentropicTemperature; ///< Rotor outlet isentropic temperature. [K]
	double FOutletMassFlow; ///< Rotor outlet mass flow. [kg / s]
	double FOutletPressure; ///< Rotor outlet pressure.
	RowVector FReducedSpeed; ///< Reduced speed. [rpm / sqrt(K)]
	RowVector FReducedMassFlow; ///< Reduced flow rate. [kg * sqrt(K) / (s * Pa)]
	RowVector FTurbineEfficiency; ///< Isentropic efficiency.
	RowVector FIsentropicWork; ///< Turbine isentropic work. [J]

	/**
	 * @brief Updates the plenum properties.
	 * 
	 * @param t New time. [s]
	 */
	void UpdatePlenum(double t);

	/**
	 * @brief Updates the state vector with the already-computed results.
	 *
	 * After solving the system and computing the state vector update, this
	 * function is called to update it.
	 */
	virtual void UpdateStateVector();

	/**
	 * @brief Updates the turbine working points.
	 * 
	 * @param t New time. [s]
	 */
	void UpdateWorkingPoint(double t);

  public:
	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/17
	 */
    TQ2DTurbine();

	virtual void AcumulaMedias(double Tiempo);

	virtual void AsignaEntradaSalidaCC();

	virtual void CabeceraResultadosMedTurb(stringstream & medoutput);

	virtual void CabeceraResultadosInstantTurb(stringstream & insoutput);

	virtual void CalculaCondicionTurbina(double TimeCalculo);

	virtual void CalculaResultadosMediosTurb();

	/*!
	 * \brief Gets the turbine efficiency, non-const version.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return Turbine adiabatic efficiency. [-]
	 */
	virtual double GetEfficiency();

	/*!
	 * \brief Gets the turbine efficiency.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return Turbine adiabatic efficiency. [-]
	 */
	virtual double GetEfficiency() const;

	/*!
	 * \brief Gets a vector with all the inlet pipes.
	 * 
	 * In the particular case of a single-entry turbine, returns a vector with
	 * only one component.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/04
	 * 
	 * \return A vector with pointers to all the inlet pipes.
	 */
	std::vector<Pipe_ptr> getInletPipes() const;

	/*!
	 * \brief Gets a pointer to the outlet pipe.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/04
	 * 
	 * \return A pointer to the outlet pipe.
	 */
	Pipe_ptr getOutletPipe() const;

	/*!
	 * \brief Gets a vector with all the volutes.
	 * 
	 * In the particular case of a single-entry turbine, returns a vector with
	 * only one component.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/03
	 * 
	 * \return A vector with pointers to all the volutes.
	 */
	std::vector<Volute_ptr> getVolutes() const;

	virtual void ImprimeResultadosInstantTurb(stringstream & insoutput);

	virtual void ImprimeResultadosMediosPantalla();

	virtual void ImprimeResultadosMedTurb(stringstream & medoutput);

	virtual void IniciaMedias();

	virtual void LeeResultadosInstantTurb(const char *FileWAM, fpos_t &filepos);

	virtual void LeeResultadosInstantTurbXML(xml_node node_turb);

	/*!
	* \brief Initialises a TQ2DTurbine object from an XML file.
	* 
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	* \date 2016/03/17
	* 
	* \param xml_node XML node that contains the object data.
	* \param fluid Working fluid.
	*/
	void ReadXML(const xml_node & node,
		const TComponentArray_ptr& fluid, const std::map<string, TSolid*> &MDB);

	virtual void ReadAverageResultsTurb(const char *FileWAM, fpos_t &filepos);

	virtual void ReadAverageResultsTurbXML(xml_node node_turb);

	/*!
	 * \brief Reads the instantaneous results setup from an XML node.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param node XML node with the object data.
	 */
	virtual void ReadInsResults(const pugi::xml_node& node);

	virtual void ResultadosInstantTurb();

	/*!
	 * \brief Sets the turbine inlet BC.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/21
	 * 
	 * \param bc Inlet boundary condition.
	 */
	void setInletBC(BoundaryCondition_ptr bc);

	/*!
	 * \brief Sets the turbine outlet BC.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/21
	 * 
	 * \param bc Outlet boundary condition.
	 */
	void setOutletBC(BoundaryCondition_ptr bc);

	/**
	 * @brief Computes the turbine, using the already-set time-step.
	 * 
	 * Integrates the flow evolution inside the turbine, updating the state vector
	 * at the end. Uses an already-set time-step.
	 */
	virtual void Solve();

	/**
	 * @brief Computes the turbine until a given time is reached.
	 *
	 * Integrates the flow evolution inside the turbine, updating the state vector
	 * at the end. The computation is performed until the time passed as an argument
	 * is reached.
	 *
	 * @param t New time. [s]
	 */
	virtual void Solve(double t);

	/*!
	 * \brief Writes the instantaneous results header.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsHeader(std::stringstream & output) const;

	/*!
	 * \brief Writes the instantaneous results.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsResults(std::stringstream & output) const;
};

/**
 * @brief Shared pointer to a TQ2DTurbine.
 */
typedef shared_ptr<TQ2DTurbine> Q2DTurbine_ptr;

#endif
