/*!
* \file TGraphicable.h
* \author Francisco Jose Arnau <farnau@mot.upv.es>
* \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
*
* \section LICENSE
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
* \section DESCRIPTION
* This file declares basic function for output results.
*/

// ---------------------------------------------------------------------------

#ifndef TGraphicableH
#define TGraphicableH

#include "Globales.h"
#include <sstream>


// ---------------------------------------------------------------------------

/*!
* \brief	Base class with the functions to manage the output results
*/

class TGraphicable {
  protected:

	  bool Is_plot;			//!< true if the object will provide results for plots
	  bool Is_cycleout;		//!< true if the object will provide cycle average results

  public:

	/*!
	* \brief Default constructor.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*/
	TGraphicable();

	/*!
	* \brief Default destructor.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*/
	virtual ~TGraphicable();

	/*!
	* \brief Build the header for the cycle average results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Array of string where the headers are stored
	*/
	virtual void HeaderForCycleOutput(std::vector<std::string> & output) const {};

	/*!
	* \brief Build the header for the instantaneous plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Array of string where the headers are stored
	*/
	virtual void HeaderForPlots(std::vector<std::string> & output) const {};

	/*!
	* \brief Integrate the last cycle instantaneous results to obtain the average.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param Output	Matrix that contains the instantaneous results.
	* \param Cicle	Array where the cycle average results are stored
	*
	*/
	virtual void IntegrateOutput(std::vector<std::vector<float>>& Output, std::vector<float>& Cicle) const {};

	/*!
	* \brief Store the instantaneous results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param output	Matrix where the instantenous results are stored
	*/
	virtual void StoreOutput(std::vector<std::vector<float>>& output) const {};

	/*!
	* \brief true if the object provides results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \return	Boolean to indicated if has plots
	*/
	virtual bool is_plot(){
		return Is_plot;
	}

	/*!
	* \brief true if the object provides cycle average results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \return	Boolean to indicated if has cycle average results
	*/
	virtual bool is_cycleout(){
		return Is_cycleout;
	}

};

typedef shared_ptr<TGraphicable> Graphicable_ptr;
#endif
