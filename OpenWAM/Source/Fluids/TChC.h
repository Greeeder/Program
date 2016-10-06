//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TCHC_H
#define __TCHC_H

/*!
 * \class	TChC
 *
 * \brief	Chemical component: C (Soot)
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	29/01/2016
 */

class TChC: public TChemicalComponent {

  private:

  public:

	/*!
	 * \fn	TChC(std::string nombre, double RU, double PMA);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	nombre	The name of the chemical component.
	 * \param	RU	  	The universal gas constant [J / (mol K)].
	 * \param	PMA   	The molecular weight [gr / mol].
	 */

	TChC(std::string nombre, double RU, double PMA);

	/*!
	 * \fn	~TChCO();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 */

	~TChC();

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

};

#endif
