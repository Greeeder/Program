/*
 * TPipeExtHeat.h
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#ifndef SOURCE_1DPIPES_TPIPEEXTHEAT_H_
#define SOURCE_1DPIPES_TPIPEEXTHEAT_H_

#include "Globales.h"
#include "TFluid.h"

/*!
 * \class	TPipeExtHeat
 *
 * \brief	Pipe external heat calculation object.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	03/02/2016
 */

class TPipeExtHeat {

  protected:

	int ncells;				   //!< Number of cells in the pipe

	TFluid *Cooler;			   //!< The cooler fluid object

	ArrayXd Dext;			   //!< External pipe diameter array [m]
	ArrayXd Aext;			   //!< External pipe area array [m^2]
	ArrayXd Tavg;			   //!< Average temperature between wall and infinity [K]
	ArrayXd Pext;			   //!< External pressure [Pa]
	ArrayXd Rho;				   //!< Density at average temperature [kg/m^3]
	ArrayXd DeltaT;			   //!< Temperature difference between cooler and wall temperature [K]

	ArrayXd Text;			   //!< Cooler tempeerature [K]
	ArrayXd Visc;			   //!< Viscosity at average temperature [Pa s]
	ArrayXd Cond;			   //!< Conductivity at average temperature [W / (m K)]

	ArrayXd Vext;			   //!< External velocity [m / s]
	ArrayXd Re;				   //!< Reynolds number [-]
	ArrayXd Pr;				   //!< Prandtl number [-]

	ArrayXd h;				   //!< Heat transfer coeficient [W / (m^2 K)]

	ArrayXd Qconv;			   //!< Heat due to convection [W]
	ArrayXd Qrad;			   //!< Heat due to radiation [W]
	ArrayXd Qtot;			   //!< Total heat (convection + radiation) [W]

	double Emissivity;		   //!< External emissivity [-]
	double HeatMultiplier;	   //!< Multiplier for the external heat [-]	

  public:

	/*!
	 * \fn	TPipeExtHeat();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 */

	TPipeExtHeat();

	/*!
	 * \fn	virtual ~TPipeExtHeat();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 */

	virtual ~TPipeExtHeat();

	/*!
	 * \fn	virtual VectorXd Heat(VectorXd Twall);
	 *
	 * \brief	Heats the given twall.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	Twall	Wall temperature [K].
	 *
	 * \return	Total heat exchanged [W].
	 */

	virtual VectorXd Heat(VectorXd Twall){
		return VectorXd::Zero(Twall.size());
	}

	/*!
	 * \fn	void setText(double T);
	 *
	 * \brief	Sets a constant temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature [K].
	 */

	void setText(double T);

	/*!
	 * \fn	void setText(ArrayXd T);
	 *
	 * \brief	Sets temperature array for each cell.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	T	Temperature array [K].
	 */

	void setText(ArrayXd T);
};

#endif /* SOURCE_1DPIPES_TPIPEEXTHEAT_H_ */
