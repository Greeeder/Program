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
 * @file AllBoundaryConditions.hpp
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
 * This file includes all the boundary condition interfaces and declares a BC
 * creator function.
 */

#ifndef AllBoundaryConditions_hpp
#define AllBoundaryConditions_hpp

#include "BoundaryCondition.hpp"
#include "ClosedEnd.hpp"
#include "ConstantConditionsBC.hpp"
#include "ExternalConnectionBC.hpp"
#include "IncidentPressure.hpp"
#include "PipeConnection.hpp"
#include "PressureBC.hpp"
#include "CheckXML.h"
#include "TBasicPlenum.h"
#include "TControlUnit.h"
#include "TEngineBlock.h"
#include <sstream>

/*!
 * \brief Creates BCs using XML data.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * \date 2016/04/10
 * 
 * \param block XML block with the BCs data.
 * \param pipes Vector of pipes.
 * \param fluid Pointer to the components array.
 */
void create_BCs(const xml_node & openwam, const std::vector<Pipe_ptr> & pipes, const std::vector<BasicPlenum_ptr>& Plenum, 
	const TComponentArray_ptr& fluid, const ControlUnit_ptr& ControlUnit, const EngineBlock_ptr& engine);

#endif
