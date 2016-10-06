/*
 * TPipeExtHeatW.h
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#ifndef SOURCE_1DPIPES_TPIPEEXTHEATW_H_
#define SOURCE_1DPIPES_TPIPEEXTHEATW_H_

#include "TPipeExtHeat.h"
#include "TChH2Ol.h"

/*!
 * \class	TPipeExtHeatW
 *
 * \brief	External heat transfered from the water to the pipe wall.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	09/03/2016
 */

class TPipeExtHeatW: public TPipeExtHeat {

  public:

	  /*!
	   * \fn	TPipeExtHeatW(int ncells, double vel, double mult, ArrayXd dext, double cellsize,
	   * 		double emis, double Tw);
	   *
	   * \brief	Constructor.
	   *
	   * \author	F.J. Arnau (farnau@mot.upv.es)
	   * \date	09/03/2016
	   *
	   * \param	ncells  	The number of the cells of the pipe.
	   * \param	vel			The external fluid velocity [m/s].
	   * \param	mult		The external heat multiplier.
	   * \param	dext		The external diameter [m].
	   * \param	cellsize	The length of the cells [m].
	   * \param	emis		The external emissivity.
	   * \param	Tw			The external fluid temperture [K].
	   */

	  TPipeExtHeatW(int ncells, double vel, double mult, ArrayXd dext, double cellsize,
		  double emis, double Tw);

	/*!
	 * \fn	virtual ~TPipeExtHeatW();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	09/03/2016
	 */

	virtual ~TPipeExtHeatW();

	/*!
	 * \fn	virtual VectorXd Heat(VectorXd Twall);
	 *
	 * \brief	Heats extracted from the water.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	09/03/2016
	 *
	 * \param	Twall	External wall temperature [K].
	 *
	 * \return	Heat [W].
	 */

	virtual VectorXd Heat(VectorXd Twall);
};

#endif /* SOURCE_1DPIPES_TPIPEEXTHEATW_H_ */
