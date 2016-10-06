/*
 * TPipeExtHeatA.h
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#ifndef SOURCE_1DPIPES_TPIPEEXTHEATA_H_
#define SOURCE_1DPIPES_TPIPEEXTHEATA_H_

#include "TPipeExtHeat.h"
#include "TFluidAir.h"

/*!
 * \class	TPipeExtHeatA
 *
 * \brief	Pipe external heat calculation object for air cooled pipes.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	03/02/2016
 */

class TPipeExtHeatA: public TPipeExtHeat {

  public:

	/*!
	 * \fn	TPipeExtHeatA();
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 * 
	 * \param	ncells		Number of cells		
	 * \param	vel			External air velocity [m/s]
	 * \param	mult		External heat multiplier [-]	
	 * \param	dext		External diameter of the pipe [m]
	 * \param	cellsize	Length of the cells [m]	
	 * \param	emis		Emissivity [-]					
	 */

	TPipeExtHeatA(int ncells, double vel, double mult, ArrayXd dext, double cellsize, double emis);

	/*!
	 * \fn	virtual ~TPipeExtHeatA();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 */

	virtual ~TPipeExtHeatA();

	/*!
	 * \fn	virtual VectorXd Heat(VectorXd Twall);
	 *
	 * \brief	Heats as a function of twall.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/02/2016
	 *
	 * \param	Twall	Wall temperature [K].
	 *
	 * \return	Heat exchanged [K].
	 */

	virtual VectorXd Heat(VectorXd Twall);
};

#endif /* SOURCE_1DPIPES_TPIPEEXTHEATA_H_ */
