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
 * @file TGenericSource.hpp
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
 * This file declares a generic source term for a TPipe object.
 */

#ifndef TGenericSource_hpp
#define TGenericSource_hpp

#include "Math_wam.h"
#include <memory>

class TPipe;

/**
 * @brief A shared pointer to a TPipe object.
 */
typedef std::shared_ptr<TPipe> Pipe_ptr;

/**
 * @brief A generic source.
 * 
 * Computes extra source terms for a TPipe object.
 */
class TGenericSource {
  protected:
	TPipe * FPipe; ///< TPipe object affected by this object.
	RowArray FSource; ///< Source terms array, for memoisation purposes.
  public:
	/**
	 * @brief Creates a generic source object.
	 * 
	 * @param pipe The pipe affected by this source.
	 */
	TGenericSource(const Pipe_ptr & pipe);

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TGenericSource();

	/**
	 * @brief Returns the source terms.
	 */
	RowArray getSource() const;

	/**
	 * @brief Computes the source terms.
	 * 
	 * @param t Current time. [s]
	 * @param dt Time-step. [s]
	 * @return The source terms. 
	 */
	virtual RowArray ComputeSource(double t, double dt) = 0;
};

/**
 * @brief A shared pointer to a TGenericSource object.
 */
typedef std::unique_ptr<TGenericSource> Source_ptr;

#endif
