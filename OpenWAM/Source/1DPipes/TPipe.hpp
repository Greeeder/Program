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
 * \file TPipe.hpp
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
 * This file declares a basic one-dimensional pipe.
 */

#ifndef TPipe_hpp
#define TPipe_hpp

#include "TCondicionContorno.h"
#include "Math_wam.h"
#include "BasicPipeMethod.hpp"
#include "TGenericSource.hpp"
#include "TIntegrable.hpp"
#include "TFlowObject.h"
#include "TGraphicable.h"
#include "TPipeIntHeat.h"
#include "TPipeFriction.h"
#include "TPipeFrictionColebrook.h"
#include "TPipeIntHeatIntake.h"
#include "TPipeIntHeatExhaust.h"
#include "TPortIntHeatIntake.h"
#include "TPortIntHeatExhaust.h"
#include "TPipeHeatT.h"
#include <memory>
#include "AllValves.hpp"

class TLaxWendroff;
class TBoundaryCondition;
typedef unique_ptr<TBoundaryCondition> BoundaryCondition_ptr;
typedef unique_ptr<TBasicPipeMethod> PipeMethod_ptr;

/*!
 * \brief A structure with instantaneous results at one point inside a pipe.
 */
struct PipeInsResults {
	PipeInsResults(); //!< Constructor.
	double Distance; //!< Distance from the pipe inlet. [m]
	bool Cp; //!< Heat capacity ratio at constant pressure status.
	bool Cv; //!< Heat capacity ratio at constant volume status.
	bool Gamma; //!< Heat capacities ratio status.
	bool LeftPressure; //!< Left-travelling pressure wave status.
	bool MassFlow; //!< Mass flow status.
	bool Pressure; //!< Pressure measurement status.
	bool R; //!< Gas constant status.
	bool RightPressure; //!< Right-travelling pressure wave status.
	bool Temperature; //!< Static temperature status.
	bool TotalPressure; //!< Total pressure measurement status.
	bool TotalTemperature; //!< Total temperature status.
	bool Velocity; //!< Velocity measurement status.
};

/*!
* \brief A structure with output results at one point inside a pipe.
*/
struct PipeOutput {
	PipeOutput();				//!< Constructor
	double Location;			//!< Location for the results [m]
	stOutput LeftPressure;      //!< Left-travelling pressure wave output.
	stOutput MassFlow;			//!< Mass flow output
	stOutput Pressure;			//!< Static pressure output
	stOutput RightPressure;     //!< Right-travelling pressure wave output.
	stOutput Temperature;		//!< Static temperature output
	stOutput TotalPressure;		//!< Total pressure output
	stOutput TotalTemperature;	//!< Total temperature output
	stOutput Velocity;			//!< Gas velocity output
};
/*!
 * \brief A one-dimensional pipe.
 *
 * This class implements a pipe, including its geometry, state vector and
 * variables such as its pressure or entropy level.  It also owns a flow
 * integration method object of class TBasicPipeMethod, as well as boundary
 * condition objects of class TCondicionContorno.
 */
class TPipe: public TIntegrable, public TFlowObject, public TGraphicable {
	friend class TBasicPipeMethod;
	friend class TPipeMethod;
	friend class TLaxWendroff;
	friend class TGodunov;
	friend class TMUSCL;

	friend Pipe_ptr create_pipe(const xml_node & node,
		const TComponentArray_ptr& components,
		const std::map<string, TSolid*> &MDB);

	friend void attach_to_pressure_BC(const Pipe_ptr& pipe, nmPipeEnd pipe_end,
		double p, double T, TFluid_ptr fluid);

	friend void attach_to_incident_pressure_BC(const Pipe_ptr& pipe,
		nmPipeEnd pipe_end, const RowVector& t, const RowVector& p,
		const RowVector& A_A, TFluid_ptr fluid);

	friend void attach_to_0Dto1Dconnection(const int id, const TFlowObject_ptr& pipe, nmPipeEnd pipe_end,
		const TFlowObject_ptr& zeroD, TTipoValvula_ptr valve);

  protected:
	RowArray FU0; ///< State vector at current time.
	double FXref; ///< Cell size. [m]
	RowVector Fa; ///< Speed of sound. [m / s]
	RowVector FArea; ///< Node or cell interface area. [m ** 2]
	RowVector FA_A; ///< Non-dimensional entropy level. [-]
	RowVector Fbeta; ///< Left-travelling non-dimensional characteristic. [-]
	RowVector FD; ///< Node or cell interface diameter. [m]
	RowVector FDcell; ///< Cell diameter. [m]
	RowVector Fhi; ///< Interior heat transfer coefficient. [W / (m ** 2 * K)]
	RowVector Frho; ///< Density. [kg / m ** 3]
	RowVector FRe; //!< Reynolds number.
	RowVector FDerLinArea; ///< First derivative of the area. [m]
	RowArray FTWPipe; ///< Pipe wall temperature. [K]
	RowVector FGamma; ///< Specific heat capacities ratio.
	RowVector Flambda; ///< Right-travelling non-dimensional characteristic. [-]
	RowVector FR; ///< Gas constant. [J / (kg * K)]
	RowVector Fcv; ///< Specific heat capacity at constant volume. [J / (kg * K)]
	RowVector Fcp; ///< Specific heat capacity at constant pressure. [J / (kg * K)]
	RowVector FViscosity; ///< Viscosity of the gas. [kg / (m * s)]
	RowVector Fdx; ///< Distance between nodes or cell length. [m]
	RowVector FGamma1; ///< Gamma - 1.
	RowVector FMArea; ///< Cell center or in-between node area. [m ** 2]
	RowVector FMx; ///< Cell center or in-between node position. [m]
	double FNin; ///< Number of cells or nodes.
	RowArray FU1; ///< State vector at next time-step.
	RowVector Fpressure; ///< Pressure vector. [Pa]
	RowVector Fspeed; ///< Flow speed vector. [m / s]
	RowVector Ftemperature; ///< Temperature vector. [K]
	RowVector FTotalPressure; ///< Total pressure vector. [Pa]
	RowVector FTotalTemperature; ///< Total temperature vector. [K]
	RowVector FVolume; ///< Cell volume. [m ** 2]
	RowVector Fx; ///< Node or cell interface position. [m]
	RowVector Fx_sv; //!< Points where the state vector is known. [m]
	double FCoefAjusFric; ///< Friction fitting coefficient.
	double FCoefAjusTC; ///< Heat transfer fitting coefficient.
	double FFriction; ///< Pipe friction.
	bool FIsIntegrated; ///< Whether or not the pipe is already integrated.
	BoundaryCondition_ptr FLeftBC; ///< Left boundary condition.
	BoundaryCondition_ptr FRightBC; ///< Right boundary condition.
	Uint FPipeNumber; ///< Pipe number.
	PipeMethod_ptr FMethod; ///< Pipe computation method.
	RowVector FQint; ///< Internal heat transfer in the pipe [W]
	TPipeIntHeat_ptr FHeatTransfer; ///< Internal heat transfer object in the pipe
	RowVector FFric; ///< Friction coefficient [-]
	TPipeFriction_ptr FPipeFriction; ///< Friction between the gas and the walls
	Source_ptr FSource; ///< Extra source terms object.

	TPipeHeatT_ptr FWallHT; ///< Heat trasnfer object in the wall

	int FNodeLeft;					   //!< The node left ID
	int FNodeRight;					   //!< The node right ID

	int FNIdenticalPipes;			   //!< Number of identical pipes
	int FNumStretch;				   //!< Number of stretches of the pipe
	double FMeshSize;				   //!< User defined mesh size;
	RowVector FDStretch;			   //!< The stretches diameter
	RowVector FLStretch;			   //!< The stretches length

	std::vector<PipeInsResults> InsResults; //!< Istantaneous results storage.
	mutable std::vector<PipeOutput> FPipeOutput; //!< Array with points for output results
	std::vector<TFluid_ptr> WorkingFluid; //!< The working fluid for each cell

  public:
	/*!
	 * \brief Creates a pipe object, with some default parameters.
	 */
	TPipe();

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TPipe();
	
	/*!
	 * \fn	void ReadDataXML(xml_node node_pipe);
	 *
	 * \brief	Read the input data from a xml node.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/02/2016
	 *
	 * \param	node_pipe	An xml node containing the pipe data.
	 */

	/*!
	 * \brief Gets the node area or cell interface area.
	 *
	 * \return Node area or cell interface area. [m ** 2]
	 */
	virtual RowVector getArea() const;

	/*!
	 * \brief Gets the node area or cell interface area for a given node or interface.
	 *
	 * \param i Node id or cell interface id.
	 * 
	 * \return Node area or cell interface area. [m ** 2]
	 */
	double getArea(Uint i) const;

	/*!
	 * \brief Gets the pipe entropy level vector.
	 *
	 * Gets the entropy level inside the pipe:
	 *
	 * @f[
	 * A_A = \cfrac{a}{a_0} / \left( \cfrac{p}{p_0} \right)
	 * ^ \frac{\gamma - 1}{2 \cdot \gamma}
	 * @f]
	 *
	 * where @f$ A_A @f$ is the entropy level, @f$ a @f$ is the speed of
	 * sound, @f$ a_0 @f$ is the reference speed of sound #ARef,
	 * @f$ p @f$ is the flow pressure, @f$ p_0 @f$ is the reference pressure
	 * #PRef and @f$ \gamma @f$ is the specific heat capacities ratio.
	 *
	 * Its value is returned for each node/cell of the pipe.
	 *
	 * \return The pipe entropy level vector. [-]
	 */
	RowVector getA_A() const;

	/*!
	 * \brief Gets the pipe entropy level at a given cell.
	 *
	 * Gets the entropy level inside the pipe:
	 *
	 * @f[
	 * A_A = \cfrac{a}{a_0} / \left( \cfrac{p}{p_0} \right)
	 * ^ \frac{\gamma - 1}{2 \cdot \gamma}
	 * @f]
	 *
	 * where @f$ A_A @f$ is the entropy level, @f$ a @f$ is the speed of
	 * sound, @f$ a_0 @f$ is the reference speed of sound #ARef,
	 * @f$ p @f$ is the flow pressure, @f$ p_0 @f$ is the reference pressure
	 * #PRef and @f$ \gamma @f$ is the specific heat capacities ratio.
	 *
	 * Its value is returned for ith node/cell of the pipe.
	 *
	 * \param i Cell number.
	 * \return The pipe entropy level at a given cell. [-]
	 */
	double getA_A(Uint i) const;

	/*!
	 * \brief Gets the pipe entropy level at a given distance from the inlet.
	 *
	 * Gets the entropy level inside the pipe:
	 *
	 * @f[
	 * A_A = \cfrac{a}{a_0} / \left( \cfrac{p}{p_0} \right)
	 * ^ \frac{\gamma - 1}{2 \cdot \gamma}
	 * @f]
	 *
	 * where @f$ A_A @f$ is the entropy level, @f$ a @f$ is the speed of
	 * sound, @f$ a_0 @f$ is the reference speed of sound #ARef,
	 * @f$ p @f$ is the flow pressure, @f$ p_0 @f$ is the reference pressure
	 * #PRef and @f$ \gamma @f$ is the specific heat capacities ratio.
	 *
	 * Its value is returned for a node/cell that is at a distance @f$ x @f$
	 * of the pipe inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe entropy level at a given point. [-]
	 */
	double getA_A(double x) const;

	/*!
	 * \brief Gets the pipe beta vector.
	 *
	 * Gets the left-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \beta = \cfrac{a}{a_0} - \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \beta @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for each node/cell of the pipe.
	 *
	 * \return The pipe beta vector. [-]
	 */
	RowVector getBeta() const;

	/*!
	 * \brief Gets the pipe beta at a given cell.
	 *
	 * Gets the left-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \beta = \cfrac{a}{a_0} - \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \beta @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for ith node/cell of the pipe.
	 *
	 * \param i Cell number.
	 * \return The pipe beta at a given cell. [-]
	 */
	double getBeta(Uint i) const;

	/*!
	 * \brief Gets the pipe beta at a given distance from the inlet.
	 *
	 * Gets the left-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \beta = \cfrac{a}{a_0} - \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \beta @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for a node/cell that is at a distance @f$ x @f$
	 * of the pipe inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe beta at a given point. [-]
	 */
	double getBeta(double x) const;

	/*!
	 * \brief Gets the pipe specific heat capacity at constant pressure vector.
	 *
	 * \return The pipe specific heat capacity at constant pressure vector. [J / (kg * K)]
	 */
	RowVector getcp() const;

	/*!
	 * \brief Gets the pipe specific heat capacity at constant pressure at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe specific heat capacity at constant pressure at a given cell. [J / (kg * K)]
	 */
	double getcp(Uint i) const;

	/*!
	 * \brief Gets the pipe specific heat capacity at constant pressure at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe specific heat capacity at constant pressure at a given point. [J / (kg * K)]
	 */
	double getcp(double x) const;

	/*!
	* \brief Gets the current time.
	*
	* \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	*
	* \return Current time. [s]
	*/

	virtual double getCurrentTime() const {
		return FCurrentTime;
	}

	/*!
	 * \brief Gets the node or inter-cell diameter.
	 * 
	 * \return Diameter. [m]
	 */
	RowVector getD() const;

	/*!
	 * \brief Gets the node or inter-cell diameter for a given node or interface.
	 * 
	 * \param i Node number or cell interface number.
	 * 
	 * \return Diameter. [m]
	 */
	double getD(Uint i) const;

	/*!
	 * \brief Gets the mesh size.
	 *
	 * \return The mesh size. [m]
	 */
	double getDeltaX() const;

	/*!
	 * \brief Gets the pipe density vector.
	 *
	 * \return The pipe density vector. [kg / m ** 3]
	 */
	RowVector getDensity() const;

	/*!
	 * \brief Gets the pipe density at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe density at a given cell. [kg / m ** 3]
	 */
	double getDensity(Uint i) const;

	/*!
	 * \brief Gets the pipe density at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe density at a given point. [kg / m ** 3]
	 */
	double getDensity(double x) const;

	/*!
	 * \brief Gets the extra source terms object.
	 * 
	 * Note for users of this function: the value returned by this function
	 * shall not be used to construct a new managed pointer!
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return The extra source terms object.
	 */
	TGenericSource * getExtraSources() const;

	/*!
	 * \brief Gets the working fluid at a cell.
	 * 
	 * \param i Cell number.
	 * \return Working fluid.
	 */
	virtual TFluid_ptr getFluid(Uint i) const;

	/*!
	 * \brief Returns the fluid components array.
	 * 
	 * \author L.M. Garcia-Cuevas (luiga12@mot.upv.es)
	 * \date 24/05/2016
	 * 
	 * This function returns an array of pointers to the chemical components
	 * of the working fluid. It doesn't return the mass or molar fractions.
	 * 
	 * \return The fluid components array.
	 */
	virtual TComponentArray_ptr getFluidComponents() const;

	/*!
	 * \brief Gets the pipe specific heat capacities ratio vector.
	 * 
	 * \return The pipe specific heat capacities ratio vector. [-]
	 */
	RowVector getGamma() const;

	/*!
	 * \brief Gets the pipe specific heat capacities ratio at a given cell.
	 * 
	 * \param i Cell number.
	 * \return The pipe specific heat capacities ratio at a given cell. [-]
	 */
	double getGamma(Uint i) const;

	/*!
	 * \brief Gets the pipe specific heat capacities ratio at a given distance from the inlet.
	 * 
	 * \param x Distance from the inlet. [m]
	 * \return The pipe specific heat capacities ratio at a given point. [-]
	 */
	double getGamma(double x) const;

	/*!
	 * \brief Gets the pipe lambda vector.
	 *
	 * Gets the right-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \lambda = \cfrac{a}{a_0} + \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \lambda @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for each node/cell of the pipe.
	 *
	 * \return The pipe lambda vector. [-]
	 */
	RowVector getLambda() const;

	/*!
	 * \brief Gets the pipe lambda at a given cell.
	 *
	 * Gets the right-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \lambda = \cfrac{a}{a_0} + \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \lambda @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for ith node/cell of the pipe.
	 *
	 * \param i Cell number.
	 * \return The pipe lambda at a given cell. [-]
	 */
	double getLambda(Uint i) const;

	/*!
	 * \brief Gets the pipe lambda at a given distance from the inlet.
	 *
	 * Gets the right-travelling non-dimensional characteristic:
	 *
	 * @f[
	 * \lambda = \cfrac{a}{a_0} + \cfrac{\gamma - 1}{2} \cdot \cfrac{u}{a_0}
	 * @f]
	 *
	 * where @f$ \lambda @f$ is the non-dimensional characteristic,
	 * @f$ a @f$ is the speed of sound, @f$ a_0 @f$ is the reference speed of
	 * sound #ARef, @f$ \gamma @f$ is the specific heat capacities ratio and
	 * @f$ u @f$ is the flow speed.
	 *
	 * Its value is returned for a node/cell that is at a distance @f$ x @f$
	 * of the pipe inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe lambda at a given point. [-]
	 */
	double getLambda(double x) const;

	/*!
	 * \brief Gets the left BC.
	 * 
	 * Note for users of this function: the value returned by this function
	 * shall not be used to construct a new managed pointer!
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return The left BC.
	 */
	TBoundaryCondition * getLeftBC() const;

	/*!
	 * \brief Returns the left node id.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 * 
	 * \return Left node id.
	 */
	int getLeftNode() const;

	/*!
	 * \brief Gets the left-travelling pressure wave at a point.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The left-travelling pipe pressure at a given point. [Pa]
	 */
	double getLeftPressure(double x) const;

	/*!
	 * \brief Gets the inter-node area or cell center area.
	 *
	 * \return Inter-node area or cell center area. [m ** 2]
	 */
	RowVector getMArea() const;

	/*!
	 * \brief Gets the mass flow rate at a given distance from the inlet.
	 * 
	 * \author L.M Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param x Distance from the inlet. [m]
	 * \return The mass flow rate. [kg / s]
	 */
	double getMassFlow(double x) const;

	/*!
	 * \brief Returns the maximum allowable time-step.
	 *
	 * Returns the maximum allowable time-step due to stability criteria.
	 *
	 * \return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep();

	/*!
	 * \brief Gets the number of cells or nodes.
	 * 
	 * \author Luis Miguel Garcia-Cuevas González <luiga12@mot.upv.es>
	 * \date 15/03/2016
	 * 
	 * \return Number of cells or nodes.
	 */
	Uint getNin() const;

	/*!
	 * \brief Gets the pipe number.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/21
	 * 
	 * \return Pipe number.
	 */
	Uint getPipeNumber() const;

	/*!
	 * \brief Gets the pipe pressure vector.
	 *
	 * \return The pipe pressure vector. [Pa]
	 */
	RowVector getPressure() const;

	/*!
	 * \brief Gets the pipe pressure at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe pressure at a given cell. [Pa]
	 */
	virtual double getPressure(Uint i) const;

	/*!
	 * \brief Gets the pipe pressure at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe pressure at a given point. [Pa]
	 */
	double getPressure(double x) const;

	/*!
	 * \brief Gets the pipe gas constant vector.
	 *
	 * \return The pipe gas constant vector. [J /(kg * K)]
	 */
	RowVector getR() const;

	/*!
	 * \brief Gets the pipe gas constant at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe gas constant at a given cell. [J /(kg * K)]
	 */
	double getR(Uint i) const;

	/*!
	 * \brief Gets the pipe gas constant at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe gas constant at a given point. [J /(kg * K)]
	 */
	double getR(double x) const;

	/*!
	 * \brief Gets the right-travelling pressure wave at a point.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The right-travelling pipe pressure at a given point. [Pa]
	 */
	double getRightPressure(double x) const;

	/*!
	 * \brief Gets the right BC.
	 * 
	 * Note for users of this function: the value returned by this function
	 * shall not be used to construct a new managed pointer!
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return The right BC.
	 */
	TBoundaryCondition * getRightBC() const;

	/*!
	 * \brief Returns the right node id.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 * 
	 * \return Right node id.
	 */
	int getRightNode() const;

	/*!
	 * \brief Gets the pipe speed vector.
	 *
	 * \return The pipe speed vector. [m / s]
	 */
	RowVector getSpeed() const;

	/*!
	 * \brief Gets the pipe speed at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe speed at a given cell. [m / s]
	 */
	virtual double getSpeed(Uint i) const;

	/*!
	 * \brief Gets the pipe speed at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe speed at a given point. [m / s]
	 */
	double getSpeed(double x) const;

	/*!
	 * \brief Gets the pipe speed vector, non-blocking version.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 * 
	 * Returns the speed vector, without blocking the mutex.
	 * 
	 * \return The pipe speed vector. [m / s]
	 */
	RowVector getSpeedNB() const;

	/*!
	 * \brief Gets the pipe speed of sound vector.
	 *
	 * \return The pipe speed vector. [m / s]
	 */
	RowVector getSpeedOfSound() const;

	/*!
	 * \brief Gets the pipe speed of sound at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe speed of sound at a given cell. [m / s]
	 */
	double getSpeedOfSound(Uint i) const;

	/*!
	 * \brief Gets the pipe speed of sound at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe speed of sound at a given point. [m / s]
	 */
	double getSpeedOfSound(double x) const;

	/*!
	 * \brief Gets the pipe speed of sound vector, non-blocking version.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 * 
	 * Returns the speed of sound vector, without blocking the mutex.
	 */
	RowVector getSpeedOfSoundNB() const;

	/*!
	 * \brief Gets the state vector.
	 * 
	 * \return The state vector.
	 */
	RowArray getStateVector() const;

	/*!
	 * \brief Gets the state vector at a given cell.
	 * 
	 * \param i Cell number.
	 * \return The state vector.
	 */
	ColArray getStateVector(Uint i) const;

	/*!
	 * \brief Gets the pipe temperature vector.
	 *
	 * \return The pipe temperature vector. [K]
	 */
	RowVector getTemperature() const;

	/*!
	 * \brief Gets the pipe temperature at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe temperature at a given cell. [K]
	 */
	virtual double getTemperature(Uint i) const;

	/*!
	 * \brief Gets the pipe temperature at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe temperature at a given point. [K]
	 */
	double getTemperature(double x) const;

	/*!
	 * \brief Gets the pipe total pressure vector.
	 *
	 * \return The pipe total pressure vector. [Pa]
	 */
	RowVector getTotalPressure() const;

	/*!
	 * \brief Gets the pipe total pressure at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe total pressure at a given cell. [Pa]
	 */
	double getTotalPressure(Uint i) const;

	/*!
	 * \brief Gets the pipe total pressure at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe total pressure at a given point. [Pa]
	 */
	double getTotalPressure(double x) const;

	/*!
	 * \brief Gets the pipe total temperature vector.
	 *
	 * \return The pipe total temperature vector. [K]
	 */
	RowVector getTotalTemperature() const;

	/*!
	 * \brief Gets the pipe total temperature at a given cell.
	 *
	 * \param i Cell number.
	 * \return The pipe total temperature at a given cell. [K]
	 */
	double getTotalTemperature(Uint i) const;

	/*!
	 * \brief Gets the pipe total temperature at a given distance from the inlet.
	 *
	 * \param x Distance from the inlet. [m]
	 * \return The pipe total temperature at a given point. [K]
	 */
	double getTotalTemperature(double x) const;

	/*!
	 * \brief Gets the cell volume vector.
	 *
	 * \return Cell volume vector. [m ** 3]
	 */
	RowVector getVolume() const;

	/*!
	 * \brief Gets the cell volume for a given node or interface.
	 *
	 * \param i Cell id.
	 * 
	 * \return Cell volume. [m ** 3]
	 */
	double getVolume(Uint i) const;

	/*!
	 * \brief	Gets working fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	18/05/2016
	 *
	 * \return	The working fluid pointer.
	 */
	std::vector<TFluid_ptr> getWorkingFluid() const{
		return WorkingFluid;
	}


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
	 * \brief Integrates one time-step without updating the state vector.
	 */
	virtual void IntegrateWithoutUpdating();

	/*!
	 * \brief Gets the node or cell interface positions.
	 *
	 * \return Node or cell interface positions vector. [m]
	 */
	RowVector getX() const;

	/*!
	 * \brief Reads the instantaneous results setup from an XML node.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \param node XML node with the object data.
	 */
	virtual void ReadInsResults(const pugi::xml_node& node);

	/*!
	* \brief Reads the output results setup from an XML node.
	*
	* \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	*
	* \param node XML node with the object data.
	*/
	virtual void ReadOutput(const xml_node& node);

	/*!
	 * \brief Sets the boundary conditions of the pipe.
	 *
	 * To compute the flow inside the pipe, the initial conditions are needed,
	 * but also a pair of boundary conditions: one for the left end and another
	 * one for the right end.
	 *
	 * Objects of type ::TBoundaryCondition may have state depending on the
	 * duct that they affect. To ensure that they have a single owner (i.e.,
	 * they are assigned to a single pipe), the pointer passed as an argument
	 * to this function has to be of type ::BoundaryCondition_ptr. Even a
	 * junction boundary condition connecting two ::TPipe objects is also
	 * unique to each pipe: some information is shared between both boundary
	 * conditions, but they are not the same.
	 *
	 * \param leftBC Left end boundary condition.
	 * \param rightBC Right end boundary condition.
	 */
	void setBCs(BoundaryCondition_ptr leftBC, BoundaryCondition_ptr rightBC);

	/*!
	 * \brief Sets the Courant number.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * \date 2016-05-04
	 *
	 * Sets the Courant number. The maximum allowable time-step is limited
	 * due to the Courant-Friedrichs-Lewy condition, so the eigenvector of the
	 * problem travelling at maximum speed (having the maximum eigenvalue)
	 * doesn't travel more than a given fraction of the mesh size, being this
	 * fraction the Courant number.
	 *
	 * \param C Courant number.
	 */
	void setCourant(double C);

	/*!
	 * \brief Sets the object generating extra source terms.
	 *
	 * Sets an object generating extra source terms. Useful for extra heat
	 * flow or for imposing a lateral window in the pipe.
	 *
	 * \param source Extra source terms object.
	 */
	void setExtraSource(Source_ptr source);

	/*!
	 * \brief Sets the pipe geometry.
	 *
	 * Sets the pipe geometry, using several sections with linear variations
	 * of area. Tries to keep the objective distance between nodes, or cell length
	 * (when using a finite-volume method), keeping a minimum of 3 nodes (or 2
	 * cells). The final node distance will be rounded to keep the pipe length.
	 * The diameter keeps a linear relationship with the axial coordinate at each
	 * stretch.
	 *
	 * \param x Known node (or cell interface) positions. [m]
	 * \param dx Node distance (or cell length) objective. [m]
	 * \param D Known node (or cell interface) diameter. [m]
	 */
	virtual void setGeometry(const RowVector& x, double dx, const RowVector& D);

	/*!
	 * \brief Sets the left boundary condition.
	 *
	 * To compute the flow inside the pipe, the initial conditions are needed,
	 * but also a pair of boundary conditions: one for the left end and another
	 * one for the right end.
	 *
	 * Objects of type ::TBoundaryCondition may have state depending on the
	 * duct that they affect. To ensure that they have a single owner (i.e.,
	 * they are assigned to a single pipe), the pointer passed as an argument
	 * to this function has to be of type ::BoundaryCondition_ptr. Even a
	 * junction boundary condition connecting two ::TPipe objects is also
	 * unique to each pipe: some information is shared between both boundary
	 * conditions, but they are not the same.
	 *
	 * \param leftBC Left end boundary condition.
	 */
	void setLeftBC(BoundaryCondition_ptr leftBC);

	/*!
	 * \brief Sets the integration method.
	 *
	 * In order to integrate the flow inside the duct, several integration
	 * methods can be used. This functions sets the integration method for
	 * this pipe.
	 *
	 * A pipe integration method has some internal state that
	 * depends on some attributes of the pipe, so it can be used only with
	 * one pipe. To ensure that the pipe integration method has a single
	 * pipe as an owner (i.e., is assigned to a single pipe), the pointer
	 * passed as an argument to this functions has to be of type
	 * ::PipeMethod_ptr.
	 *
	 * \param method Pipe integration method.
	 */
	void setMethod(PipeMethod_ptr method);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given total pressure, total temperature
	 * and flow speed.
	 *
	 * \param p_t Total pressure. [Pa]
	 * \param T_t Total temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	void setPtTtU(double p_t, double T_t, double u);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given total pressure, total temperature
	 * and flow speed, one set of values for each node/cell.
	 *
	 * \param p_t Total pressure. [Pa]
	 * \param T_t Total temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	void setPtTtU(const RowVector& p_t, const RowVector& T_t, const RowVector& u);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given pressure, temperature and flow speed.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	void setPTU(double p, double T, double u);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given pressure, temperature and flow speed,
	 * one set of values for each node/cell.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 */
	void setPTU(const RowVector& p, const RowVector& T, const RowVector& u);

	/*!
	 * \brief Sets the state vector.
	 *
	 * Sets the state vector with a given pressure, temperature and flow speed.
	 * These values are known at different points of the duct, so the final
	 * state vector is interpolated and set in its node/cells.
	 *
	 * \param p Pressure. [Pa]
	 * \param T Temperature. [K]
	 * \param u Flow speed. [m / s]
	 * \param x Points where the state vector is known. [m]
	 */
	void setPTU(const RowVector& p, const RowVector& T, const RowVector& u,
			const RowVector& x);

	/*!
	 * \brief Sets the right boundary condition.
	 *
	 * To compute the flow inside the pipe, the initial conditions are needed,
	 * but also a pair of boundary conditions: one for the left end and another
	 * one for the right end.
	 *
	 * Objects of type ::TBoundaryCondition may have state depending on the
	 * duct that they affect. To ensure that they have a single owner (i.e.,
	 * they are assigned to a single pipe), the pointer passed as an argument
	 * to this function has to be of type ::BoundaryCondition_ptr. Even a
	 * junction boundary condition connecting two ::TPipe objects is also
	 * unique to each pipe: some information is shared between both boundary
	 * conditions, but they are not the same.
	 *
	 * \param rightBC Right end boundary condition.
	 */
	void setRightBC(BoundaryCondition_ptr rightBC);

	/*!
	 * \brief Saves the object state into an XML node.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 *
	 * \param base Base node where the element state will be saved.
	 */
	virtual void saveState(xml_node* base);

	/*!
	 * \brief Integrate the flow.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct.
	 */
	virtual void Solve();

	/*!
	 * \brief Integrate the flow, updating the state vector.
	 *
	 * Integrates the flow evolution inside the duct, updating the state
	 * vector at the end. Computes until the time passed is reached.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param t Time at the end of the integration. [s]
	 */
	virtual void Solve(double t);

	/*!
	* \brief Store the instantaneous results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Matrix where the instantenous results are stored
	*/
	virtual void StoreOutput(std::vector<std::vector<float>>& Output) const;

	/*!
	 * \brief Updates the state vector with the already-computed results.
	 *
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 * 
	 * After solving the system and computing the state vector update, this
	 * function is called to update it.
	 */
	virtual void UpdateStateVector();

	/*!
	 * \brief Writes the instantaneous results header.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsHeader(std::stringstream & output) const;

	/*!
	* \brief Writes the instantaneous results header.
	*
	* \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	*/
	virtual void WriteInsHeader(std::vector<std::string> & output) const;

	/*!
	 * \brief Writes the instantaneous results.
	 * 
	 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsResults(std::stringstream & output) const;

};

/*!
 * \brief A shared pointer to a TPipe object.
 */
typedef shared_ptr<TPipe> Pipe_ptr;

/*!
 * \brief Creates a new pipe.
 * 
 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
 * \date 2016/03/21
 * 
 * \param node xml_node with the pipe data.
 * \param fluid Pointer to the species vector.
 */
Pipe_ptr create_pipe(const xml_node & node,
	const TComponentArray_ptr& components,
	const std::map<string, TSolid*> &MDB);

#endif
