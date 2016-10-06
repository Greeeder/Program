//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TChH2Ol_H
#define __TChH2Ol_H

/*!
 * \class	TChH2Ol
 *
 * \brief	Chemical component: H2O (liquid).
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	29/01/2016
 */

class TChH2Ol: public TChemicalComponent {

  private:

  public:

	/*!
	 * \fn	TChH2Ol(std::string nombre, double RU, double PMA);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	RU	  	The universal gas constan [J / (mol K)].
	 * \param	PMA   	The molecular weigth [gr / mol].
	 */

	TChH2Ol(std::string nombre, double RU, double PMA);

	/*!
	 * \fn	~TChH2Ol();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 */

	~TChH2Ol();

	/*!
	 * \fn	double FunCp(double T);
	 *
	 * \brief	Function for C_p.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat C_p [J / (kg K)].
	 */

	virtual double FunCp(double T);

	/*!
	 * \fn	double FunCv(double T);
	 *
	 * \brief	Function for C_v.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Specific heat C_v [J / (kg K)].
	 */

	virtual double FunCv(double T);

	/*!
	 * \fn	double FunU(double T);
	 *
	 * \brief	Function for the internal energy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Internal energy [J / kg].
	 */

	virtual double FunU(double T);

	/*!
	 * \fn	double FunH(double T);
	 *
	 * \brief	Function for the enthalpy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Enthalpy [J / kg].
	 */

	virtual double FunH(double T);

	/*!
	 * \fn	double Funk(double T);
	 *
	 * \brief	Function for the conductivity.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature [K].
	 *
	 * \return	Conductivity [W / (m K)].
	 */

	virtual double Funk(double T);

	/*!
	 * \fn	double Funk(double T);
	 *
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
	 * \fn	double Funk(double T);
	 *
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
	* \fn	double FunRHO(double T);
	*
	* \brief	Function for the density.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	03/02/2016
	*
	* \param	T	Temperature [K].
	*
	* \return	Density [kg / m ** 3].
	*/
	virtual double FunRHO(double T);

};

#endif
