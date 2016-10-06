/**
 * @file TBasicPlenum.h
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
 * This file include the different methods to solve a basic 0D element.
 */
#ifndef SOURCE_ENGINE_TBASICPLENUM_H_
#define SOURCE_ENGINE_TBASICPLENUM_H_

#include <memory>

#include "TIntegrable.hpp"
#include "T0DModel.h"
#include "TFlowObject.h"
#include "TGraphicable.h"

#include "BoundaryCondition.hpp"

struct stPlenumEntry{
	double p;
	double T;
	double p0;
	double T0;
	double u;
	double m;
};

struct PlenumOutput {
	PlenumOutput();				//!< Constructor
	stOutput Mass;			//!< Mass output
	stOutput Pressure;		//!< Pressure output
	stOutput Temperature;	//!< Temperature output
	stOutput Volume;		//!< Volume output
};


class TBasicPlenum: public TIntegrable, public TFlowObject, public TGraphicable {

	/*!
	 * \brief	Creates a basicplenum.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	node 	The node.
	 * \param	fluid	The fluid.
	 *
	 * \return	The new basicplenum.
	 */

	friend shared_ptr<TBasicPlenum> create_basicplenum(const xml_node& node,
	const TComponentArray_ptr& fluid);


protected:

	int FBasicPlenum_ID;			   //!< Identifier for the basic plenum

	T0DModel_ptr FZeroD;			   //!< The 0-dimensional thermodynamic system

	TFluid_ptr FWorkingFluid;		   //!< The working fluid

	double FPreviousTime;			   //!< Time of the previous [s]
	
	double FVolume;					   //!< The volume [m ** 3]

	vector<BoundaryCondition_ptr> Boundaries;   //!< The array of boundaries
	
	bool FPrintAVG;					   //!< true to print average results
	bool FPrintINS;					   //!< true to print instantaneous results

	dVector FMass_in;
	dVector FEnth_in;
	TFluidArray FFluid_in;
	std::vector<stPlenumEntry> PlenumEntry;

	mutable PlenumOutput FPlenumOutput;			//!< Plenum output results

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 */

	TBasicPlenum();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 */

	virtual ~TBasicPlenum();

	/*!
	 * \brief	Integrates one time-step without updating the state vector.
	 *
	 * \author	L.M. Garc�a-Cuevas &lt;luiga12@mot.upv.es&gt;
	 * \date	25/05/2016
	 */

	virtual void IntegrateWithoutUpdating();

	void AppendNewEntry();

	/*!
	 * \brief	Gets the area.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The area [m ** 2].
	 */

	virtual RowVector getArea() const {
		//! \todo return the area of the connection i;
		return RowVector::Constant(1, 10.);
	}

	/*!
	 * \brief	Gets the current time.
	 *
	 * \author	L.M. Garc�a-Cuevas &lt;luiga12@mot.upv.es&gt;
	 * \date	25/05/2016
	 *
	 * \return	Current time. [s].
	 */

	virtual double getCurrentTime() const {
		return FCurrentTime;
	}

	/*!
	 * \brief	Gets the fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The fluid object.
	 */

	TFluid_ptr getFluid(){
		return FWorkingFluid;
	}

	/*!
	 * \brief	Gets a fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	i	Index of the cell (not used).
	 *
	 * \return	The fluid object.
	 */

	virtual TFluid_ptr getFluid(unsigned int i) const {
		return FWorkingFluid;
	}

	/*!
	* \brief Returns the fluid components array.
	*
	* \author L.M. Garc�a-Cuevas (luiga12@mot.upv.es)
	* \date 24/05/2016
	*
	* This function returns an array of pointers to the chemical components
	* of the working fluid. It doesn't return the mass or molar fractions.
	*
	* \return The fluid components array.
	*/
	virtual TComponentArray_ptr getFluidComponents() const;

	virtual int getID() const { return FBasicPlenum_ID; }

	/*!
	 * \brief	Gets the mass.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The mass [kg].
	 */

	double getMass(){
		return FZeroD->getMass();
	}

	/*!
	 * \brief	Gets maximum time step.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The maximum time step [s].
	 */

	virtual double getMaxTimeStep();

	/*!
	 * \brief	Gets delta x coordinate.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The delta x coordinate [m].
	 */

	virtual double getDeltaX() const {
		return 1.;
	}

	/*!
	 * \brief	Gets the number of cells.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The number of cells.
	 */

	virtual unsigned int getNin() const {
		return PlenumEntry.size();
	}

	/*!
	 * \brief	Gets the pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The pressure [Pa].
	 */

	double getPressure(){
		return FZeroD->getPressure();
	}

	/*!
	 * \brief	Gets a pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	i	Index of the cell (not used).
	 *
	 * \return	The pressure [Pa].
	 */

	virtual double getPressure(unsigned int i) const {
		return PlenumEntry[i].p;
	}

	/*!
	 * \brief	Gets the temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \return	The temperature [K].
	 */

	double getTemperature(){
		return FZeroD->getTemperature();
	}

	/*!
	 * \brief	Gets a temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	i	Index of the cell (not used).
	 *
	 * \return	The temperature [K].
	 */

	virtual double getTemperature(unsigned int i) const {
		return PlenumEntry[i].T;
	}

	/*!
	 * \brief	Gets a speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	i	Index of the cell (not used).
	 *
	 * \return	The speed [m / s].
	 */

	virtual double getSpeed(unsigned int i) const {
		return PlenumEntry[i].u;
	}

	void setBoundary(BoundaryCondition_ptr Boundary){
		Boundaries.push_back(move(Boundary));
	}

	/*!
	 * \brief	Sets initial conditions.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	p	The pressure [Pa].
	 * \param	T	The temperature [K].
	 */

	void setInitialConditions(double p, double T);

	/*!
	 * \brief	Sets total flow conditions at the entry i.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	i	Index of the entry.
	 * \param	p	Total pressure [Pa].
	 * \param	T	Total temperature [K].
	 * \param	u	Speed [u].
	 */

	void setPtTtU(int i, double p, double T, double u);

	/*!
	 * \brief	Sets total flow conditions at the entry i.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	27/05/2016
	 *
	 * \param	i	Index of the entry.
	 * \param	p	Static pressure [Pa].
	 * \param	T	Static temperature [K].
	 * \param	u	Speed [u].
	 */

	void setPTU(int i, double p, double T, double u);

	/*!
	 * \brief	Sets average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 */

	void SetAverageOutput();

	/*!
	 * \brief	Sets a fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	com	The components array.
	 * \param	Y  	The mass fraction array [-].
	 */

	void setFluid(const TComponentArray_ptr& com, RowVector Y);

	/*!
	 * \brief	Sets instantaneous output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 */

	void SetInstantaneousOutput();

	/*!
	 * \brief	Header average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param [in,out]	medoutput	The medoutput.
	 */

	void HeaderAverageOutput(stringstream& medoutput);

	/*!
	 * \brief	Header instantaneous output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param [in,out]	insoutput	The insoutput.
	 */

	void HeaderInstantaneousOutput(stringstream& insoutput);

	/*!
	 * \brief	Integrate average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	TActual	The current time [s].
	 */

	void IntegrateAverageOutput(double TActual);

	/*!
	 * \brief	Print average output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param [in,out]	medoutput	The medoutput.
	 */

	void PrintAverageOutput(stringstream& medoutput);

	/*!
	 * \brief	Print instantaneuos output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param [in,out]	insoutput	The insoutput.
	 */

	void PrintInstantaneuosOutput(stringstream& insoutput);

	/*!
	 * \brief	Reads average output XML.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	node_cyl	The node cylinder.
	 */

	void ReadAverageOutputXML(xml_node node_cyl);

	/*!
	 * \brief	Reads instantaneous output XML.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	node_cyl	The node cylinder.
	 */

	void ReadInstantaneousOutputXML(xml_node node_cyl);

	/*!
	 * \brief	Sets 0-d thermodynamic model.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	25/05/2016
	 *
	 * \param	zerod	The 0-d thermodynamic model.
	 */

	void set0DModel(T0DModel_ptr zerod);

	/*!
	 * \brief	Sets the time-step for the integration of this plenum.
	 *
	 * \author	L.M. Garc�a-Cuevas (luiga12@mot.upv.es);
	 * \date	25/05/2016
	 *
	 * \param	dt	Time-step. [s].
	 */

//	virtual void setTimeStep(double dt);

	/*!
	 * \brief	Integrate the flow, updating the state vector.
	 *
	 * \author	L.M. Garc�a-Cuevas (luiga12@mot.upv.es);
	 * 			
	 * 			Integrates the flow evolution inside the duct, updating the state vector at the end.
	 * 			Uses an already-set time-step.
	 * \date	25/05/2016
	 */

	virtual void Solve();

	/*!
	 * \brief	Integrate the flow, updating the state vector.
	 *
	 * \author	L.M. Garc�a-Cuevas (luiga12@mot.upv.es);
	 * 			
	 * 			Integrates the flow evolution inside the duct, updating the state vector at the end.
	 * 			Computes until the time passed is reached.
	 * \date	25/05/2016
	 *
	 * \param	t	Time at the end of the integration. [s].
	 */

	virtual void Solve(double t);

	/*!
	 * \brief	Updates the state vector with the already-computed results.
	 *
	 * \author	L.M. Garc�a-Cuevas &lt;luiga12@mot.upv.es&gt;
	 * 			
	 * 			After solving the system and computing the state vector update, this function is
	 * 			called to update it.
	 * \date	25/05/2016
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
	* \author L.M. Garc�a-Cuevas <luiga12@mot.upv.es>
	*
	* \param node XML node with the object data.
	*/
	virtual void ReadOutput(const xml_node& node);
};

/*!
 * \brief	A shared pointer to a TBasisPlenum object.
 */

typedef shared_ptr<TBasicPlenum> BasicPlenum_ptr;

/*!
 * \brief	Creates a basic plenum.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	25/05/2016
 *
 * \param	node 	The node.
 * \param	fluid	The fluid components array.
 *
 * \return	The new basic plenum.
 */

BasicPlenum_ptr create_basicplenum(const xml_node& node, const TComponentArray_ptr& fluid);

#endif /* SOURCE_ENGINE_TCYLINDER_H_ */
