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
 * @file TNullSource.hpp
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
 * This file declares a null source term for a TPipe object.
 */

#ifndef TNullSource_hpp
#define TNullSource_hpp

#include <memory>
#include "TGenericSource.hpp"

/**
 * @brief A null source. It produces no mass, momentum or energy.
 */
class TNullSource: public TGenericSource {
  public:
	/**
	 * @brief Creates a null source object.
	 * 
	 * @param pipe The pipe affected by this source.
	 */
	TNullSource(const Pipe_ptr & pipe);

	/**
	 * @brief Computes the source terms.
	 * 
	 * @param t Current time. [s]
	 * @param dt Time-step. [s]
	 * @return The source terms. 
	 */
	virtual RowArray ComputeSource(double t, double dt);
};

/**
 * @brief Sets a TNullSource object to a pipe.
 * 
 * @param pipe Pipe.
 */
void set_null_source(const Pipe_ptr & pipe);

#endif
