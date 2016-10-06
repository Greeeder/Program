/* --------------------------------------------------------------------------------*\
 ==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
  \\ |  X  | //  W ave     |
   \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
    \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*--------------------------------------------------------------------------------*/

/**
 * @file TVirtualPipe.hpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
 * This file declares a pipe connection and helper functions.
 */

#ifndef TVirtualPipe_hpp
#define TVirtualPipe_hpp

#include "TPipe.hpp"
#include "TFlowObject.h"
#include "Godunov.hpp"
#include <utility>
#include "AllValves.hpp"

class TVirtualPipe: public TPipe {
  protected:
	vector<ColVector> FFluxSign; ///< Flux sign correction.
	vector<unsigned int> FCells; ///< Cells connected.
	vector<TFlowObject *> FFlowObj; ///< Flow elements connected
	mutable std::mutex vp_mtx; //!< Virtual pipe mutex.

	TTipoValvula_ptr FValve;
	double FDischCoef;
	double FGeometricalArea;

	/**
	 * @brief Updates the state vector.
	 *
	 * @param time New current time for this virtual pipe. [s]
	 */
	void UpdateStateVector(double time);

	void UpdateStateVector();

  public:

	unique_ptr<TGodunov> FMethod; ///< Virtual pipe computation method.

	/**
	 * @brief Default constructor.
	 */
	TVirtualPipe();

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TVirtualPipe();

	/**
	 * @brief Creates a virtual pipe for a connection between two flow elements.
	 *
	 * @param fobj_0 Left flow element.
	 * @param pipe_end_0 Left pipe end connected using this virtual pipe.
	 * @param fobj_1 Right flow element.
	 * @param pipe_end_1 Right pipe end connected using this virtual pipe.
	 */

	TVirtualPipe(const TFlowObject_ptr & fobj_0, nmPipeEnd pipe_end_0, const TFlowObject_ptr & fobj_1, nmPipeEnd pipe_end_1);

	/**
	 * @brief Returns the flow produced at the middle of this virtual pipe.
	 *
	 * @param t Current time. [s]
	 * @param pipe_number Pipe for which the flow is to be returned.
	 * @return The flow at the middle of this virtual pipe.
	 */
	ColVector Flow(double t, unsigned int pipe_number);

	ColVector Flow();

	/*!
	 * \brief Gets the virtual pipe area.
	 * 
	 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
	 * \date 15/03/2016
	 * 
	 * \return Virtual pipe area. [m ** 2]
	 */
	double getVirtualArea() const;

	double CalculateDCin(double t);

	double CalculateDCout(double t);

	/*!
	 * \brief Returns whether the virtual pipe needs an update or not.
	 * 
	 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
	 * 
	 * \param t Time. [s]
	 * \return True if needs an update.
	 */
	bool needsUpdate(double t) const;

	/*!
	 * \brief Sets the virtual pipe area.
	 * 
	 * \author Luis Miguel García-Cuevas González <luiga12@mot.upv.es>
	 * \date 15/03/2016
	 * 
	 * \param A New virtual pipe area. [m ** 2]
	 */
	void setArea(double A);

	void setArea();

	TFlowObject* getFlowObject(int i){
		return FFlowObj[i];
	}

	int getCell(int i){
		return FCells[i];
	}

	void setCell(int i, int c){
		FCells[i] = c;
	}

	void setValve(TTipoValvula_ptr valve){
		FValve = move(valve);
	}
};

/**
 * @brief A shared poitner to a TVirtualPipe object.
 */
typedef shared_ptr<TVirtualPipe> VirtualPipe_ptr;


#endif
