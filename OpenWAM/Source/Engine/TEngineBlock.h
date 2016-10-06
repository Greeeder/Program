/**
 * @file TEngineBlock.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 5 de abr. de 2016
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
 * Includes all properties and methos to solve the physics of the engine block and its components.
 */
#ifndef SOURCE_ENGINE_TENGINEBLOCK_H_
#define SOURCE_ENGINE_TENGINEBLOCK_H_

#include "Globales.h"
#include "TCylinder.h"
#include "TInjectionSystem.h"

struct EngineBlockOuput {
	EngineBlockOuput();
	stOutput VolumetricEfficiency;
};

class TEngineBlock : public TIntegrableGroup, public TGraphicable{
private:

	double FAngle;					   //!< The angle [deg]
	double FAngleStep;				   //!< The angle step [deg]
	double FSpeed;					   //!< The speed [rpm]

	double FTotalVolume;			   //!< The total volume [m ** 3]

	string FEngineType;				   //!< Type of the engine

	double FAngleCycle;				   //!< The total angle of the cycle [rad]

	std::vector<TCylinder_ptr> Cylinder;	//!< The cylinder object array

	shared_ptr<TController> FRPMController; //!< The engine speed controller

	InjectionSystem_ptr FInjectionSystem;

	int FRPMControllerID;			   //!< Identifier for the engine speed controller
	bool FRPMControlled;			   //!< true if the engine speed is controlled

	stResMediosMotor FAverageOutput;	//!< The average output

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 */

	TEngineBlock();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 */

	virtual ~TEngineBlock();

	/*!
	 * \brief	Reads engine block XML.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param	node_openwam	The node openwam.
	 * \param	WorkFluid   	The working fluid components array.
	 */

	void ReadEngineBlockXML(xml_node node_openwam, TComponentArray_ptr WorkFluid);

	/*!
	 * \brief	Gets a cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param	i	Zero-based index of the cylinder.
	 *
	 * \return	The cylinder object.
	 */

	shared_ptr<TCylinder> getCylinder(int i){
		return Cylinder[i];
	}

	std::vector<TCylinder_ptr> getCylinders() {
		return Cylinder;
	}

	double* getSpeed_ptr() {
		return &FSpeed;
	}

	/*!
	 * \brief	Reads the average output XML.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param	node_engine	The node engine.
	 */

	void ReadAVGOutputXML(xml_node node_engine);

	/*!
	 * \brief	Header for the average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param [in,out]	medoutput	The average output container.
	 */

	void HeaderAVGOutput(std::stringstream& medoutput);

	/*!
	 * \brief	Integrate a vg output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param	TActual		  	The current time.
	 * \param	CilindroActual	The current cylindrer.
	 */

	void IntegrateAVGOutput(double TActual, int CilindroActual);

	void AppendInjectionSystem(InjectionSystem_ptr injsys) {
		FInjectionSystem = move(injsys);
	}

	/*!
	 * \brief	Calculates average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 */

	void ComputeAVGOutput();

	/*!
	 * \brief	Print average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param [in,out]	medoutput	The average output container.
	 */

	void PrintAVGOutput(std::stringstream& medoutput);

	/*!
	 * \brief	Sets the engine speed controller controller.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	01/06/2016
	 *
	 * \param	Controller	The controller.
	 */

	void SetRPMController(vector<shared_ptr<TController>> Controller);

	void Solve();

	void UpdateWorkingPoint();

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
	* \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	*
	* \param node XML node with the object data.
	*/
	virtual void ReadOutput(const xml_node& node);

};

typedef shared_ptr<TEngineBlock> EngineBlock_ptr;

#endif /* SOURCE_ENGINE_TENGINEBLOCK_H_ */
