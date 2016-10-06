//#include "stdafx.h"
#include "TChemicalComponent.h"

#ifndef __TCHNO_H
#define __TCHNO_H

/*!
 * \class	TChNO
 *
 * \brief	Chemical Component: NO
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	29/01/2016
 */

class TChNO: public TChemicalComponent {

  private:

  public:

	/*!
	 * \fn	TChNO(std::string nombre, double RU, double PMA);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 *
	 * \param	nombre	The name of the chemical NOmponent.
	 * \param	RU	  	The universal gas NOnstant [J / (mol K)].
	 * \param	PMA   	The molecular weight [gr / mol].
	 */

	TChNO(std::string nombre, double RU, double PMA);

	/*!
	 * \fn	~TChNO();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/01/2016
	 */

	~TChNO();

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
