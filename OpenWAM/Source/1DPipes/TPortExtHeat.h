/*
 * TPortExtHeat.h
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#ifndef SOURCE_1DPIPES_TPortExtHeat_H_
#define SOURCE_1DPIPES_TPortExtHeat_H_

#include "TPipeExtHeat.h"
#include "TChH2Ol.h"

class TPortExtHeat: public TPipeExtHeat {

  private:

	double TorqueMaxPower;			   //!< The torque at maximum power [Nm]
	double EngineSpeed;				   //!< The engine speed [rpm]
	double Lchar;					   //!< The characteristic length [m]
	double Ncil;					   //!< The number of cylinders

	ArrayXd ViscWall;		   //!< The coolant viscosity at wall temperature conditions [Pa s]

  public:

	/*!
	 * \fn	TPortExtHeat();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 */

	TPortExtHeat();

	/*!
	 * \fn	virtual ~TPortExtHeat();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 */

	virtual ~TPortExtHeat();

	/*!
	 * \fn	virtual VectorXd Heat(VectorXd Twall);
	 *
	 * \brief	Heats as a funtion of a given wall temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	Twall	Wall temperature [K].
	 *
	 * \return	Heat transfer [W].
	 */

	virtual VectorXd Heat(VectorXd Twall);

	/*!
	 * \fn	void setTorqueMaxPower(double val)
	 *
	 * \brief	Sets torque at maximum power.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	val	The torque [N m].
	 */

	void setTorqueMaxPower(double val) {
		TorqueMaxPower = val;
	}

	/*!
	 * \fn	void setLchar(double val)
	 *
	 * \brief	Sets the characteristic length as a function of the cylinder diameter.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	val	Engine cylinder diameter [m].
	 */

	void setLchar(double val) {
		Lchar = val * 1.5;
	}

	/*!
	 * \fn	void setVext(double val);
	 *
	 * \brief	Sets the coolant velocity as a function of the engine speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/02/2016
	 *
	 * \param	val	Engine speed [rpm].
	 */

	void setVext(double val);
};

#endif /* SOURCE_1DPIPES_TPortExtHeat_H_ */
