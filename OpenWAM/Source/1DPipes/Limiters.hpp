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
 * @file Limiters.hpp
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
 * The functions in this file represent limiters that can be used
 * in MUSCL-like schemes.
 *
 * This file contains the interface of such functions.
 */

#ifndef Limiters_hpp
#define Limiters_hpp

#include "Math_wam.h"


/**
 * @brief Limiter function.
 *
 * Returns the value to limit the slopes in a high-resolution scheme,
 * which is a function of the ratio of consecutive gradients of the
 * state vector @f$ w @f$.
 *
 * The ratio of consecutive gradients is defined as:
 *
 * @f[
 * r_{i,r} = \cfrac{w_{i+1} - w_i}{w_i - w_{i-1}}
 * @f]
 *
 * for the right extrapolation and
 *
 * @f[
 * r_{i,l} = \cfrac{w_i - w_{i-1}}{w_{i+1} - w_i}
 * @f]
 *
 * for the left extrapolation.  This way, the right extrapolated value
 * for the boundary between the cell @f$ i @f$ and the cell
 * @f$ i + 1 @f$ becomes:
 *
 * @f[
 * w_{i,r} = w_i + \cfrac{\phi \left( r_{i,r} \right)}{2}
 * \cdot \left( w_i - w_{i - 1} \right)
 * @f]
 *
 * and the left extrapolated value:
 *
 * @f[
 * w_{i+1,l} = w_{i+1} - \cfrac{\phi \left( r_{i+1,l} \right)}{2}
 * \cdot \left( w_{i+2} - w_{i+1} \right)
 * @f]
 *
 * where @f$ \phi @f$ is the value returned by the limiter function.
 *
 * @param ratio Ratio of consecutive gradients.
 * @param phi Limiter value for the slopes extrapolation.
 */
typedef void (*Limiter_pt)(const RowArray &, RowArray *);


/**
 * @brief Monotonized central limiter.
 *
 * A monotonized central limiter function by van Leer.  Symmetric,
 * second-order TVD. See Limiter_pt().
 *
 * @param ratio Ratio of consecutive gradients.
 * @param phi Limiter value for the slopes extrapolation.
 */
void MC(const RowArray & ratio, RowArray * phi);


/**
 * @brief Minmod central limiter.
 *
 * The minmod limiter function by Roe.  Symmetric, second-order TVD.
 * See Limiter_pt().
 *
 * @param ratio Ratio of consecutive gradients.
 * @param phi Limiter value for the slopes extrapolation.
 */
void Minmod(const RowArray & ratio, RowArray * phi);


/**
 * @brief Superbee limiter.
 *
 * The superbee limiter function by Roe.  Symmetric, second-order TVD.
 * See Limiter_pt().
 *
 * @param ratio Ratio of consecutive gradients.
 * @param phi Limiter value for the slopes extrapolation.
 */
void Superbee(const RowArray & ratio, RowArray * phi);


/**
 * @brief Van Leer limiter.
 *
 * The van Leer limiter function by van Leer.  Symmetric, second-order TVD.
 * See Limiter_pt().
 *
 * @param ratio Ratio of consecutive gradients.
 * @param phi Limiter value for the slopes extrapolation.
 */
void VanLeer(const RowArray & ratio, RowArray * phi);

#endif
