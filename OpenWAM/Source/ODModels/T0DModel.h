/**
* @file T0DModel.h
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
* This file include the different methods to solve the mass and energy balance in OD elements.
*/

#include "TFluid.h"
#include "BasicHeatTransfer.h"
#include <vector>
#include <memory>

#ifndef __T0DMODEL_H
#define __T0DMODEL_H

class T0DModel {

	private:
		//Calculate and saves previous instant
		double PPressure;			   //!< The pressure in the next step [Pa]
		double PPressure0;			   //!< The pressure in the current step [Pa]
		double PPressureStep;		   //!< The average pressure during the last integration step [Pa]
		double PTemperature;		   //!< The temperature in next step [K]
		double PTemperature0;		   //!< The temperature in current step [K]
		double PMass0;				   //!< The mass in [kg] at the current step
		double PMass;				   //!< The mass in [kg] at the next step
		double PVolume0;			   //!< The volume in the current step [m ** 3]
		double PVolume;				   //!< The volume in the next step [m ** 3]
		double dFQL;				   //!< The rate of heat release in [J]
		double EnthalpyIn;			   //!< Addition on inlet flow terms in [J]
		double EnthalpyOut;			   //!< Addition on outlet flow terms in [J]
		double Rcil;				   //!< Specific gas constant in [J/kg K]
		double dUsen;				   //!< Variation of sensible internal energy in [J]
		double dWi;					   //!< Work in [J]
		double Cv;					   //!< Specific heat at constant volume in [J/kg K]
		double dEbb;				   //!< Energy term due to blow-by in [J]
		double dEfinj;				   //!< Energy term due to fuel injection in  [J]
		double dEAinj;				   //!< Energy term due to air injection in [J]
		double dQconv;				   //!< Convective HT in [J]
		double dQrad;				   //!< RadiativHT in [J]
		
		BasicHeatTransfer *PHeatTransfer;	//!< The heat transfer object
		
		shared_ptr<TFluid> Fluid;				   //!< The fluid object

		/*!
		 * \brief	Gets the solve enthalpy.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	.
		 */

		double SolveEnthalpy();

	public:

		/*!
		 * \brief	Default constructor.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 */

		T0DModel();

		/*!
		 * \brief	Destructor.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 */

		~T0DModel();

		/*!
		 * \brief	Appends a mass quantity.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	m	The mass [kg].
		 */

		void AppendMass(double m);

		void Initialize(double Pressure, double Temperature, double Mass, double Volume);

		/*!
		 * \brief	Initializes this object.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	Pressure				The pressure [Pa].
		 * \param	Temperature				The temperature [K].
		 * \param	Mass					The mass [kg].
		 * \param	Volume					The volume [m ** 3].
		 * \param [in,out]	FluidCyl		If non-null, the fluid in the cylinder.
		 * \param [in,out]	HeatTransfer	If non-null, the heat transfer object.
		 */

		void Initialize(double Pressure, double Temperature, double Mass, double Volume, TFluid *FluidCyl, BasicHeatTransfer *HeatTransfer);

		/*!
		 * \brief	Gets the pressure.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The pressure [Pa].
		 */

		double getPressure() const { return PPressure0; }

		double getPressureStep() const { return PPressureStep; }

		/*!
		 * \brief	Gets the temperature.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The temperature [K].
		 */

		double getTemperature() const { return PTemperature0; }

		/*!
		 * \brief	Gets the mass.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The mass [kg].
		 */

		double getMass() const { return PMass0; } 

		/*!
		 * \brief	Gets the getd rate of heat release.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	Rate of heat release [-].
		 */

		double getdFQL() const { return dFQL; } 

		/*!
		 * \brief	Gets the inlet enthalpy.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The enthalpy [J].
		 */

		double getEnthalpyIn() const { return EnthalpyIn; }

		/*!
		 * \brief	Gets the outlet enthalpy.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The enthalpy [J].
		 */

		double getEnthalpyOut() const { return EnthalpyOut; }

		shared_ptr<TFluid> getFluid() { return Fluid; }

		/*!
		 * \brief	Gets the gas constant.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The gas constant [J / (kg K)].
		 */

		double getR() const { return Rcil; }

		/*!
		 * \brief	Gets the sensible heat.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	Sensible heat [J].
		 */

		double getdUsen() const { return dUsen; }

		/*!
		 * \brief	Gets the work.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The work [J].
		 */

		double getdWi() const { return dWi; } 

		/*!
		 * \brief	Gets the specific heat at constant volume.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	The specific heat at constant volume [J / (kg K)].
		 */

		double getCv() const { return Cv; }

		/*!
		 * \brief	Gets the convective heat transfer.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	Heat [J].
		 */

		double getdQconv() const { return dQconv; }

		/*!
		 * \brief	Gets the radiative heat transfer
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \return	Heat [J].
		 */

		double getdQrad() const { return dQrad; }

		/*!
		 * \brief	Gets energy variation due to outlet flow.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	R			Gas constant of the fluid in the cylinder [J / (kg K)].
		 * \param	MassIn  	The mass that has gone out [kg].
		 * \param	Tmed		The average temperature in the cylinder [K].
		 *
		 * \return	Energy [J].
		 */

		double getDeltaFlujoSaliente(double R, double MassIn, double Tmed);

		/*!
		 * \brief	Gets energy variation due to the inlet flow.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	EnthIn 	The enthalpy of the inlet flow [J].
		 * \param	FluidIn	The fluid that has gone in.
		 * \param	MassIn 	The mass that has gone in [kg].
		 * \param	Tmed   	The tmed.
		 *
		 * \return	The delta flujo entrante.
		 */

		double getDeltaFlujoEntrante(double EnthIn, TFluid FluidIn, double MassIn, double Tmed);

		BasicHeatTransfer* getHeatTransfer() { return PHeatTransfer; };

		void setComposition(RowVector Y){
			Fluid->SetComposition(Y);
		}

		void setFluid(TFluid_ptr fluid){
			Fluid = fluid;
		}

		void setHeatTransfer(BasicHeatTransfer* HT) {
			PHeatTransfer = HT;
		}

		/*!
		 * \brief	Solve the mass balance.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param [in,out]	MassIn	Array with all inlet/outlet mass flows [kg].
		 */

		void SolveMassBalance(std::vector<double> *MassIn);

		/*!
		 * \brief	Solve the energy balance imposing heat released.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	Volume				The volume [m ** 3].
		 * \param	dFQLp				The rate of heat release [J].
		 * \param [in,out]	MassIn  	If non-null, the inlet/outlet mass array [kg].
		 * \param [in,out]	EnthIn  	If non-null, the inlet enthalpy array [J / kg].
		 * \param [in,out]	FluidIn 	If non-null, the inlet fluid array.
		 * \param [in,out]	FluidCyl	If non-null, the fluid in the cylinder.
		 */

		void SolveNewPressure(double Volume, double dFQLp, std::vector<double> *MassIn, std::vector<double> *EnthIn, std::vector<TFluid*>*FluidIn, TFluid *FluidCyl);

		/*!
		* \brief	Solve the energy balance imposing heat released.
		*
		* \author	F.J. Arnau (farnau@mot.upv.es)
		* \date	30/03/2016
		*
		* \param [in,out]	MassIn  	If non-null, the inlet/outlet mass array [kg].
		* \param [in,out]	EnthIn  	If non-null, the inlet enthalpy array [J / kg].
		* \param [in,out]	FluidIn 	If non-null, the inlet fluid array.
		* \param [in,out]	FluidCyl	If non-null, the fluid in the cylinder.
		*/

		void SolveNewPressure(std::vector<double> *MassIn, std::vector<double> *EnthIn, TFluidArray_ptr FluidIn, TFluid_ptr FluidCyl);

		/*!
		* \brief	Solve the energy balance imposing heat released.
		*
		* \author	F.J. Arnau (farnau@mot.upv.es)
		* \date	30/03/2016
		*
		* \param	Volume				The volume [m ** 3].
		* \param	dFQLp				The rate of heat release [J / kg].
		* \param [in,out]	MassIn  	If non-null, the inlet/outlet mass array [kg].
		* \param [in,out]	EnthIn  	If non-null, the inlet enthalpy array [J].
		* \param [in,out]	FluidIn 	If non-null, the pointer to the inlet fluid array.
		* \param [in,out]	FluidCyl	If non-null, the pointer to the fluid in the cylinder.
		*/

		void SolveNewPressure(double Volume, double dFQLp, std::vector<double> *MassIn, std::vector<double> *EnthIn, TFluidArray_ptr FluidIn, TFluid_ptr FluidCyl);

		/*!
		 * \brief	Solve the energy balance imposing pressure.
		 *
		 * \author	F.J. Arnau (farnau@mot.upv.es)
		 * \date	30/03/2016
		 *
		 * \param	Volume				The volume [m ** 3].
		 * \param	Pressure			The pressure [Pa].
		 * \param [in,out]	MassIn  	If non-null, the inlet/outlet mass array [kg].
		 * \param [in,out]	EnthIn  	If non-null, the inlet enthalpy [J / kg].
		 * \param [in,out]	FluidIn 	If non-null, the inlet fluid array.
		 * \param [in,out]	FluidCyl	If non-null, the fluid in the cylinder.
		 */

		void SolveNewFQL(double Volume, double Pressure, std::vector<double> *MassIn, std::vector<double> *EnthIn, std::vector<TFluid*>*FluidIn, TFluid *FluidCyl);

		void UpdateStateVector(){
			PTemperature0 = PTemperature;
			PPressure0 = PPressure;
			PMass0 = PMass;
			PVolume0 = PVolume;
		}
};

typedef unique_ptr<T0DModel> T0DModel_ptr;

#endif

