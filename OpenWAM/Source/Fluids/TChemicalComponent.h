//#include "stdafx.h"
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "Math_wam.h"

#ifndef __TCHEMICALCOMPONENT_H
#define __TCHEMICALCOMPONENT_H

/*!
 * \class	TChemicalComponent
 *
 * \brief	Chemical component.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	04/12/2015
 */

class TChemicalComponent {

  private:
	std::string _nombre;			   //!< The name of the chemical component

  protected:
	double _CAHFL;//!< Coefficient A for enthalpy calculation [J / kg]
	double _CBHFL;//!< Coefficient B for enthalpy calculation [J / (kg K)]

	double _R;//!< Gas constant [J / (kg K)]
	double _PM;//!< Molecular Weight [gr / mol]
	double _HC;//!< The hidrogen/carbono ratio.
	double _OC;//!< The oxigen/carbono ratio
	double _HV;//!< The heating value [J / kg]
	double _DENSF;//!< Density [kg / m^3]
	double _FE;//!< Stoichiometric fuel to air ratio [-]

	double _RU;//!< Universal constant [J / (mol K)]

	dVector ViscCoef;				   //!< Coefficients for viscosity correlation

  public:

 	/*!
 	 * \brief	Default constructor.
 	 *
 	 * \author	F.J. Arnau (farnau@mot.upv.es)
 	 * \date	04/12/2015
 	 */

 	TChemicalComponent();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 */

	virtual ~TChemicalComponent();

	/*!
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat Cp [J / (kg K)].
	 */

	virtual double FunCp(double T) = 0;

	/*!
	 * \brief	Function for C_v.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat Cv [J / (kg K)].
	 */

	virtual double FunCv(double T) = 0;

	/*!
	 * \brief	Gets the gas constant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	Gas constant [J / (kg K)].
	 */

	virtual double FunR() {
		return _R;
	}

	/*!
	 * \brief	Function for the internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	virtual double FunU(double T) = 0;

	/*!
	 * \brief	Function for the Enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	T	Temperature [T].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T) = 0;

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

	virtual double Funk(double T) {

		std::cout << "The conductivity of this chemical component is not defined" << std::endl;
		return 0.;
	};

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

	virtual double FunVisc(double T) {

		std::cout << "The viscosity of this chemical component is not defined" << std::endl;
		return 0.;
	};

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

	virtual double FunPr(double T) {

		return FunCp(T) * FunVisc(T) / Funk(T);
	}

	/*!
	 * \brief	Gets the fun hv.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	.
	 */

	virtual double FunHV() {
		return _HV;
	}

	/*!
	 * \brief	Gets the molecula weight.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	Molecular weight [gr / mol].
	 */

	virtual double FunPM() {
		return _PM;
	}

	/*!
	 * \brief	Gets the fun relacion hc.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	.
	 */

	virtual double FunRelacionHC() {
		return _HC;
	}

	/*!
	 * \brief	Gets the fun relacion oc.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	.
	 */

	virtual double FunRelacionOC() {
		return _OC;
	}

	/*!
	 * \brief	Gets the fe.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	The fe.
	 */

	double getFE() {
		return _FE;
	}

	/*!
	 * \brief	Gets the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	Density [kg / m ** 3].
	 */

	virtual double FunRHO() {
		return _DENSF;
	}

	/*!
	 * \brief	Gets the density.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Density [kg / m ** 3].
	 */

	virtual double FunRHO(double T) {
		return _DENSF;
	};

	/*!
	 * \brief	Gets the fun cahfl.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	.
	 */

	virtual double Fun_CAHFL() {
		return _CAHFL;
	}

	/*!
	 * \brief	Gets the fun cbhfl.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	.
	 */

	virtual double Fun_CBHFL() {
		return _CBHFL;
	}

	/*!
	 * \brief	Sets a name of the chemical component.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \param	nombre	The name of the chemical species.
	 */

	void setName(std::string nombre);

	/*!
	 * \brief	Gets the name of the chemical component.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/12/2015
	 *
	 * \return	The name of the chemical component.
	 */

	std::string getName();

	/*!
	 * \brief	Sets the coeficients of the correlation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	Parameter	The parameter of the correlation.
	 * \param	vcoef	 	The list of coefficients.
	 */

	void setCoeficients(string Parameter, dVector vcoef);

};

typedef shared_ptr<TChemicalComponent> TChemicalComponent_ptr;
typedef vector<TChemicalComponent_ptr> TComponentArray_obj;
typedef shared_ptr<TComponentArray_obj> TComponentArray_ptr;

#endif
