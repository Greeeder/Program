/**
 * @file TFluid.h
 * @author F.J. Arnau <farnau@mot.upv.es>
 * @version 0.1.0
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
 * The TFluid class is used for the different fluid to be used in OpenWAM.
 *
 * It represent a fluid composed by several chemical components and calculates
 * the thermodynamic properties depending on the composition.
 *
 */

//#include "stdafx.h"

#include "TChemicalComponent.h"
#include "Globales.h"
#include <vector>
#include <memory>

#ifndef __TFLUID_H
#define __TFLUID_H

/*!
 * \class	TFluid
 *
 * \brief	Fluid composed by different chemical components.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	05/12/2015
 */

class TFluid {

  protected:

	RowVector Y;					   //!< Mass fraction of the different chemical components that composes the fluid [-].

//	std::vector<TChemicalComponent_ptr> *_Component;//!< Vector with the different chemical components that composes the fluid.
	TComponentArray_ptr _Component;//!< Vector with the different chemical components that composes the fluid.


	double R;//!< The gas constant of the fluid [J / (kg K)]

	double Viscosity;//!< The viscosity [Pa s]
	double Conductivity;//!< The conductivity [W / (m K)]
	double Cp;//!< The specific heat [J / (kg K)]

  public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	TFluid();

	/*!
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/04/2016
	 *
	 * \param [in,out]	fluid	If non-null, the fluid.
	 */

	TFluid(TFluid* fluid);

	/*!
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	componentes	Array with the chemical components that composes the fluid.
	 */

	TFluid(TComponentArray_obj *componentes);

	/*!
	* \brief	Constructor.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	05/12/2015
	*
	* \param [in,out]	componentes	Pointer to an array with the chemical components that composes the fluid.
	*/

	TFluid(TComponentArray_ptr componentes);

	/*!
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	name	The name of the fluid.
	 */

	TFluid(string name);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 */

	virtual ~TFluid();

	/*!
	 * \brief	Get the mass fraction of the chemical component called "componente".
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	componente	The name of the chemical component.
	 *
	 * \return	Mass fraction [-].
	 */

	double GetY(std::string componente);

	/*!
	 * \brief	Get mass fration for the chemical component i.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	i	Index of the chemical component.
	 *
	 * \return	Mass fraction [-].
	 */

	double GetY(int i);

	/*!
	 * \brief	Sets the mass fraction of the chemical component called "nombre".
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	valor 	The mass fraction [-].
	 */

	void SetY(std::string nombre, double valor);

	/*!
	 * \brief	Join two fluid and determines the new composition.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	fluid1	The fluid 2.
	 * \param	m1			  	Mass of the fluid 2 [kg].
	 * \param	m0			  	Mass of the original fluid [kg].
	 */

	void AppendFluid(shared_ptr<TFluid> fluid1, double m1, double m0);

	/*!
	 * \brief	Function for C_v.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / (kg K)].
	 */

	virtual double FunCv(double T);

	/*!
	 * \brief	Function for C_v as a function of the composicion.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	comp	Vector with the mass fraction of the components [-].
	 * \param	T				Temperature [K].
	 *
	 * \return	Specific heat C_v [J / (kg K)].
	 */

	double FunCvComposicion(std::vector<double>* comp, double T);

	/*!
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat [J / (kg K)].
	 */

	virtual double FunCp(double T);

	/*!
	 * \brief	Function for the internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	virtual double FunU(double T);

	/*!
	* \brief	Function for the sensible internal energy.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	05/12/2015
	*
	* \param	T	Temperature [K].
	*
	* \return	Internal energy [J / kg].
	*/

	virtual double FunUsens(double T);

	/*!
	 * \brief	Function for the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Density [kg / m ** 3].
	 */

	virtual double FunRHO(double T);

	/*!
	 * \brief	Function internal energy as a fucntion of the input composition.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	comp	Vector with the mass fraction of the different components [-].
	 * \param	T				Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	double FunUComposicion(std::vector<double>* comp, double T);

	/*!
	 * \brief	Function for the enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T);

	/*!
	 * \brief	Gets the gas constant R.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \return	Gas constant [J / (kg K)].
	 */

	virtual double FunR();

	/*!
	 * \brief	Get gas constant as a function of the input composition.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	comp	Mass fraction of the components.
	 *
	 * \return	Gas constant [J / (kg K)].
	 */

	double FunRComposicion(std::vector<double>* comp);

	/*!
	 * \brief	Function for the specific heat ratio.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat ratio [-].
	 */

	virtual double FunGamma(double T);

	/*!
	 * \brief	Function for the specific heat ratio as a function of the input composition.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	comp	Mass fraction of the components.
	 * \param	T				Temperature [K].
	 *
	 * \return	Specific heat ratio [-].
	 */

	double FunGammaComposicion(std::vector<double>* comp, double T);

	/*!
	 * \brief	Gets the chemical component of the fluid as a function of its name.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	nombre	Name of the chemical component.
	 *
	 * \return	null if it fails, else the chemical component.
	 */

	TChemicalComponent_ptr getChemicalComponent(std::string nombre);

	/*!
	 * \brief	Sets the components of the fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param [in,out]	componentes	If non-null, the components.
	 */

	void SetComponentes(TComponentArray_obj *componentes);

	/*!
	* \brief	Sets the components of the fluid.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	05/12/2015
	*
	* \param [in,out]	componentes	If non-null, the components pointer.
	*/

	void SetComponentes(TComponentArray_ptr componentes);

	/*!
	 * \brief	Gets a component as a function of its name.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/12/2015
	 *
	 * \param	i	The index of the chemical component.
	 *
	 * \return	null if it fails, else the component.
	 */

	TChemicalComponent_ptr GetComponent(int i);

	/*!
	 * \brief Gets an array with pointers to the chemical components.
	 * 
	 * \author L.M. García-Cuevas (luiga12@mot.upv.es)
	 * \date 24/05/2016
	 * 
	 * \return An array with pointers to the chemical components.
	 */
	TComponentArray_ptr GetComponents() const;

	/*!
	 * \brief	Gets the temperature as a function of the internal energy [J / kg].
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param	u	The internal energy [J / kg]
	 *
	 * \return	The temperature [K].
	 */

	double getTemperature(double u);

	/*!
	 * \brief	Gets temperature from the Equation of State for ideal gases.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param	rho	The density [kg / s].
	 * \param	p  	The pressure [Pa].
	 *
	 * \return	The temperature [K].
	 */

	double getTfromEqOfState(double rho, double p);

	/*!
	 * \brief	Gets the pressure from the Equation of State for idel gases.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param	rho	The density [kg / s].
	 * \param	T  	The temperature [K].
	 *
	 * \return	The pressure [Pa].
	 */

	double getPfromEqOfState(double rho, double T);

	/*!
	 * \brief	Gets the density from the Equation of State for idel gases.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param	p	The pressure [Pa]
	 * \param	T	The temperature [K].
	 *
	 * \return	The density [kg / s].
	 */

	double getRhofromEqOfState(double p, double T);

	/*!
	 * \brief	Gets the gas constant of the fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \return	The gas constant [J / (kg K)].
	 */

	double getR() {
		return R;
	}

	/*!
	 * \brief	Update the fluid composition and its gas constant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param	Ynew	The new compposition [-].
	 */

	void SetComposition(RowVector Ynew);

	/*!
	 * \brief	Gets the current composition of the fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \return	The composition [-].
	 */

	RowVector GetComposition() {
		return Y;
	}

	/*!
	 * \brief	Sets the gas constant as a function of the current composition.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 */

	void SetR() {
		R = FunR();
	}

	/*!
	 * \brief	Function for the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature [T].
	 *
	 * \return	Conductivity [W / (m K)].
	 */

	virtual double Funk(double T);

	/*!
	 * \brief	Function for the viscosity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Viscosity [Pa s].
	 */

	virtual double FunVisc(double T);

	/*!
	 * \brief	Function for the Prandtl.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Prandtl [-].
	 */

	virtual double FunPr(double T);

	/*!
	 * \brief	Gets the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \return	The conductivity [W / (m K)].
	 */

	double getConductivity() {
		return Conductivity;
	}

	/*!
	 * \brief	Gets the cp.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \return	The specific heat [J / (kg K)].
	 */

	double getCp() {
		return Cp;
	}

	/*!
	 * \brief	Gets molar fraction by chemical component name.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	02/06/2016
	 *
	 * \param	nombre	The name of the chemical component.
	 *
	 * \return	The molar fraction [-].
	 */

	double getMolarFraction(std::string nombre);

	/*!
	 * \brief	Gets molecular weight of the fluid.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	02/06/2016
	 *
	 * \return	The molecular weight [g / mol].
	 */

	double getMolecularWeight();

	/*!
	 * \brief	Gets the viscosity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \return	The viscosity [Pa s].
	 */

	double getViscosity() {
		return Viscosity;
	}

	/*!
	 * \brief	Sets a conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	k	Conductivity [W / (m K)].
	 */

	void setConductivity(double k) {
		Conductivity = k;
	}

	/*!
	 * \brief	Sets a cp.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	c_p	Specific heat [J / (kg K)].
	 */

	void setCp(double c_p) {
		Cp = c_p;
	}

	/*!
	 * \brief	Sets a viscosity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	mu	Viscosity [Pa s].
	 */

	void setViscosity(double mu) {
		Viscosity = mu;
	}

	/*!
	 * \brief	Gets the size.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/04/2016
	 *
	 * \return	Number of components of the fluid.
	 */

	unsigned int size(){ return _Component->size(); }
};

/*!
 * \struct	stTforRealGas
 *
 * \brief	Error function between the internal energy provided and the internal energy obtained as a
 * 			function of the temperature.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	17/12/2015
 */

struct stTforRealGas {
	TFluid *Fluid;					   //!< The fluid
	double u;//!< The internal energy [J / kg]

	/*!
	 * \fn	stTforRealGas(TFluid* _Fluid, double _u)
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \param [in,out]	_Fluid	If non-null, the fluid.
	 * \param	_u			  	The internal energy [J / kg].
	 */

	stTforRealGas(TFluid* _Fluid, double _u) {
		u = _u;
		Fluid = _Fluid;
	}

	/*!
	 * \fn	double operator()(const double x)
	 *
	 * \brief	 Return the error in internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	17/12/2015
	 *
	 * \return	The error.
	 */

	double operator()(const double x) {
		return u - Fluid->FunUsens(x);
	}

};

typedef shared_ptr<TFluid> TFluid_ptr;
typedef vector<TFluid_ptr> TFluidArray;
typedef shared_ptr<TFluidArray> TFluidArray_ptr;

#endif
