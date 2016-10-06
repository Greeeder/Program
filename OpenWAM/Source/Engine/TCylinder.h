/**
 * @file TCylinder.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 30 de mar. de 2016
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
 * This file include the different methods to solve the thermodynamics in the cylinder.
 */
#ifndef SOURCE_ENGINE_TCYLINDER_H_
#define SOURCE_ENGINE_TCYLINDER_H_

#include <memory>

#include "T0DModel.h"
#include "TInjector.h"
#include "AllCombustionType.h"
#include "TCrankMechanism.h"
#include "TBasicPlenum.h"
#include "HeatTransfer.h"

#include "BoundaryCondition.hpp"

struct CylinderOutput {
	stOutput Angle;				//!< Angle output
	stOutput Pressure;			//!< Pressure output
	stOutput Temperature;		//!< Temperature output
	stOutput Mass;				//!< Mass output
	stOutput Volume;			//!< Volume output
	stOutput IntakeMass;		//!< Intake mass flow output
	stOutput ExhaustMass;		//!< Exhaust mass flow output
	stOutput FuelMass;			//!< Fuel mass rate output
	stOutput BlowBy;			//!< Blow by output
	stOutput HeatPower;			//!< Heat power output
	stOutput HeatReleased;		//!< Heat release rate output
};

class TCylinder : public TBasicPlenum {

	/*!
	 * \brief	Creates a cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	node		 	The node.
	 * \param	node_cylinder	The node cylinder.
	 * \param	WorkFluid	 	The work fluid.
	 *
	 * \return	The new cylinder.
	 */

	friend shared_ptr<TCylinder> create_cylinder(xml_node node, xml_node node_cylinder, TComponentArray_ptr WorkFluid);

private:

	int FCylinder_ID;				   //!< Identifier for the cylinder

	int FCycle;

	bool FInj;
	bool FComb;
	bool FFirstCycle;

	bool FClosedLoop;

	TInjector_ptr FInjector;		   //!< The injector object
	TCombustion_ptr FCombustion;	   //!< The combustion object
	TCrankMechanism_ptr FCrankMech;	   //!< The engine mechanism object
	//HeatTransfer *FHeatTransferModel;   //!< The heat transfer model

	double FAngleStep;				   //!< The angle step [rad]
	double FCombAngle;				   //!< The current angle used for combustion [rad]. Range 4S [-2pi 2pi], 2S [-pi pi]
	double FIVCAngle;
	double FEVOAngle;

	double FHRL0;					   //!< The current heat relesed [-]
	double FHRL1;					   //!< The next heat relesed [-]

	double FPhaseDiference;			   //!< Phase difference between this cyclinder and the cylinder 1 [rad]
	double FAngleCycle;
	double FOrder;

	TChemicalComponent_ptr FFuel;	   //!< The fuel component

	vector<TBoundaryCondition*> IntakeValve;   //!< The array of intake valves boundaries
	vector<TBoundaryCondition*> ExhaustValve;   //!< The array of exhaust valves boundaries
	
	stResMediosCilindro FAVGOutput;	   //!< The average results output
	stResInstantCilindro FINSOutput;	//!< The instantaneous results output

	double FNetTorque;					   //!< The torque [Nm]
	double FWork;					   //!< The work [J]

	bool FPrintAVG;					   //!< true to print average results
	bool FPrintINS;					   //!< true to print instantaneous results

	mutable CylinderOutput FCylinderOutput;		//!< Cylinder output manager

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 */

	TCylinder();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 */

	virtual ~TCylinder();

	/*!
	 * \brief	Integrates one time-step without updating the state vector.
	 *
	 * \author	L.M. García-Cuevas (luiga12@mot.upv.es)
	 * \date	27/05/2016
	 */

	virtual void IntegrateWithoutUpdating();

	TCrankMechanism* getCrankMechanism() {
		return FCrankMech.get();
	}

	TInjector_ptr getInjector(){
		return FInjector;
	}

	/*!
	 * \brief	Gets the torque.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \return	The torque [Nm].
	 */

	double getNetTorque(){
		return FNetTorque;
	}

	/*!
	 * \brief	Gets displaced volume.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \return	The volume displaced [m ** 3].
	 */

	double getVolumeDisplaced(){
		return FCrankMech->getVolumeDisplaced();
	}

	/*!
	 * \brief	Gets the work.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \return	The work [J].
	 */

	double getWork(){
		return FWork;
	}

	/*!
	 * \brief	Sets initial conditions.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	p		  	The pressure [Pa].
	 * \param	T		  	The temperature [T].
	 * \param	EngineType	Type of the engine.
	 */

	void setInitialConditions(double p, double T, string EngineType);

	/*!
	 * \brief	Sets average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 */

	void SetAverageOutput();

	/*!
	 * \brief	Sets a fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	com	The components array.
	 * \param	Y  	The mass fraction array [-].
	 */

	void setFluid(const TComponentArray_ptr& com, RowVector Y);

	/*!
	* \brief	Append a pointer to an exhaust valve of the cylinder.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	27/05/2016
	*
	* \param	ech_valv	Pointer to the exhaust valve.
	*/

	void appendExhaustValve(TBoundaryCondition* ech_valv) {
		IntakeValve.push_back(ech_valv);
	}

	/*!
	* \brief	Append a pointer to an intake valve of the cylinder.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	27/05/2016
	*
	* \param	int_valv	Pointer to the intake valve.
	*/

	void appendIntakeValve() {
		IntakeValve.push_back(Boundaries.back().get());
	}

	/*!
	 * \brief	Sets instantaneous output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 */

	void SetInstantaneousOutput();

	/*!
	 * \brief	Header average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param [in,out]	medoutput	The average results container.
	 */

	void HeaderAverageOutput(stringstream& medoutput);

	/*!
	 * \brief	Header instantaneous output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param [in,out]	insoutput	The instantaneous results container.
	 */

	void HeaderInstantaneousOutput(stringstream& insoutput);

	/*!
	 * \brief	Integrate average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	TActual	The current time.
	 */

	void IntegrateAverageOutput(double TActual);

	/*!
	 * \brief	Print average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param [in,out]	medoutput	The average results container.
	 */

	void PrintAverageOutput(stringstream& medoutput);

	/*!
	 * \brief	Print instantaneuos output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param [in,out]	insoutput	The instantaneous results container.
	 */

	void PrintInstantaneuosOutput(stringstream& insoutput);

	/*!
	 * \brief	Reads average output (XML format).
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	node_cyl	The node cylinder.
	 */

	void ReadAverageOutputXML(xml_node node_cyl);

	/*!
	 * \brief	Reads instantaneous output (XML format).
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	node_cyl	The node cylinder.
	 */

	void ReadInstantaneousOutputXML(xml_node node_cyl);

	void setAngle();

	/*!
	 * \brief	Sets an angle.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	angle	The angle [deg].
	 */

	void setAngle(double angle);

	/*!
	 * \brief	Sets a combustion object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	Comb	The combustion object.
	 */

	void setCombustion(TCombustion_ptr Comb);

	/*!
	 * \brief	Sets engine mechanism.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	CrMch	The engine mechanism.
	 */

	void setCrankMechanism(TCrankMechanism_ptr CrMch);

	/*!
	* \brief	Sets the angle step.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	27/05/2016
	*
	* \param	da	The angle step [rad].
	*/

	void setDeltaAngle(double da);

	/*!
	 * \brief	Sets engine speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	n	The speed [rpm].
	 */

	void setEngineSpeed(double n);

	/*!
	 * \brief	Sets the time-step for the integration of this pipe.
	 *
	 * \author	L.M. García-Cuevas (luiga12@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	dt	Time-step. [s].
	 */

	virtual void setTimeStep(double dt);

	/*!
	 * \brief	Integrate the flow, updating the state vector.
	 *
	 * \author	L.M. García-Cuevas (luiga12@mot.upv.es)
	 * 			
	 * 			Integrates the flow evolution inside the duct, updating the state vector at the end.
	 * 			Uses an already-set time-step.
	 * \date	27/05/2016
	 */

	virtual void Solve();

	/*!
	 * \brief	Integrate the flow, updating the state vector.
	 *
	 * \author	L.M. García-Cuevas (luiga12@mot.upv.es)
	 * 			
	 * 			Integrates the flow evolution inside the duct, updating the state vector at the end.
	 * 			Computes until the time passed is reached.
	 * \date	27/05/2016
	 *
	 * \param	t	Time at the end of the integration. [s].
	 */

	virtual void Solve(double t);

	/*!
	 * \brief	Updates the state vector with the already-computed results.
	 *
	 * \author	L.M. García-Cuevas (luiga12@mot.upv.es)
	 * 			
	 * 			After solving the system and computing the state vector update, this function is
	 * 			called to update it.
	 * \date	27/05/2016
	 */

	virtual void UpdateStateVector();

	/*!
	* \brief Build the header for the cycle average results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Array of string where the headers are stored
	*/
	virtual void HeaderForCycleOutput(std::vector<std::string> & output) const;

	/*!
	* \brief Build the header for the instantaneous plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Array of string where the headers are stored
	*/
	virtual void HeaderForPlots(std::vector<std::string> & output) const;

	/*!
	* \brief Integrate the last cycle instantaneous results to obtain the average.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param Output	Matrix that contains the instantaneous results.
	* \param Cicle	Array where the cycle average results are stored
	*
	*/
	virtual void IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const;

	/*!
	* \brief Store the instantaneous results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Matrix where the instantenous results are stored
	*/
	virtual void StoreOutput(std::vector<std::vector<float>>& output) const;

	/*!
	* \brief Reads the output results setup from an XML node.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*
	* \param node XML node with the object data.
	*/
	virtual void ReadOutput(const xml_node& node);

};

typedef shared_ptr<TCylinder> TCylinder_ptr;

TCylinder_ptr create_cylinder(xml_node node, xml_node node_cylinder, TComponentArray_ptr WorkFluid);

#endif /* SOURCE_ENGINE_TCYLINDER_H_ */
