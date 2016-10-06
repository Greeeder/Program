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

/*!
 * \file TVolute.hpp
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
 * This file declares a quasi-2D volute.
 */

#ifndef TVolute_hpp
#define TVolute_hpp

#include "TPipe.hpp"
#include "AllPipeMethods.hpp"
#include "AllExtraSourceTerms.hpp"
#include "AllBoundaryConditions.hpp"

class TVolute;

/*!
 * \brief Shared pointer to a TVolute object.
 */
typedef shared_ptr<TVolute> Volute_ptr;

/*!
 * \brief A quasi-2D volute.
 * 
 * This volute class is a basic TPipe with extra data about the flow through its
 * lateral window.
 */
class TVolute: public TPipe {
	friend class TBasicPipeMethod;
	friend class TPipeMethod;
	friend class TGodunov;
	friend class TMUSCL;

	friend Volute_ptr create_volute(const xml_node & node,
		const TComponentArray_ptr& components,
		const std::map<string, TSolid*> &MDB);

  protected:
	RowVector FLateralEnthalpyFlow; ///< Lateral enthalpy flow. [W]
	RowVector FLateralMassFlow; ///< Lateral mass flow. [kg / s]
	RowVector FWindowArea; ///< Lateral window area. [m ** 2]

  public:

	/*!
	 * \brief Gets the lateral enthalpy flow.
	 * 
	 * \return The lateral enthalpy flow. [W]
	 */
	RowVector getLateralEnthalpyFlow() const;

	/*!
	 * \brief Gets the lateral enthalpy flow at a given cell.
	 * 
	 * \param i Cell number.
	 * \return The lateral enthalpy flow. [W]
	 */
	double getLateralEnthalpyFlow(Uint i) const;

	/*!
	 * \brief Gets the lateral enthalpy flow at a given distance from the inlet.
	 * 
	 * \param x Distance from the inlet. [m]
	 * \return The lateral enthalpy flow. [W]
	 */
	double getLateralEnthalpyFlow(double x) const;

	/*!
	 * \brief Gets the lateral mass flow.
	 * 
	 * \return The lateral mass flow. [kg / s]
	 */
	RowVector getLateralMassFlow() const;

	/*!
	 * \brief Gets the lateral mass flow at a given cell.
	 * 
	 * \param i Cell number.
	 * \return The lateral mass flow. [kg / s]
	 */
	double getLateralMassFlow(Uint i) const;

	/*!
	 * \brief Gets the lateral mass flow at a given distance from the inlet.
	 * 
	 * \param x Distance from the inlet. [m]
	 * \return The lateral mass flow. [kg / s]
	 */
	double getLateralMassFlow(double x) const;

	/*!
	 * \brief Gets the effective area of all the lateral nozzles.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \return The lateral nozzle effective area. [m ** 2]
	 */
	RowVector getLateralNozzleArea() const;

	/*!
	 * \brief Gets the lateral nozzle area at a given cell.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param i Cell number.
	 * \return The lateral nozzle effective area. [m ** 2]
	 */
	double getLateralNozzleArea(Uint i) const;

	/*!
	 * \brief Gets the lateral window area.
	 * 
	 * \return The lateral window area. [m ** 2]
	 */
	RowVector getWindowArea() const;

	/*!
	 * \brief Gets the lateral window area at a given cell.
	 * 
	 * \param i Cell number.
	 * \return The lateral window area. [m ** 2]
	 */
	double getWindowArea(Uint i) const;

	/*!
	 * \brief Gets the lateral window area at a given distance from the inlet.
	 * 
	 * \param x Distance from the inlet. [m]
	 * \return The lateral window area. [m ** 2]
	 */
	double getWindowArea(double x) const;

	/*!
	 * \brief Sets the pipe geometry.
	 *
	 * Sets the pipe geometry, using several sections with linear variations
	 * of area. Tries to keep the objective distance between nodes, or cell length
	 * (when using a finite-volume method), keeping a minimum of 3 nodes (or 2
	 * cells). The final node distance will be rounded to keep the pipe length.
	 * The area changes has a linear relationship with the axial coordinate at each
	 * stretch.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/05/03
	 * 
	 * \param x Known node (or cell interface) positions. [m]
	 * \param dx Node distance (or cell length) objective. [m]
	 * \param D Known node (or cell interface) diameter. [m]
	 */
	virtual void setGeometry(const RowVector & x, double dx, const RowVector & D);

	/*!
	 * \brief Sets the lateral nozzle area for one cell.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param A Effective area. [m ** 2]
	 * \param i Cell number.
	 */
	void setLateralNozzleArea(double A, Uint i);

	/*!
	 * \brief Sets the same lateral nozzle area for all the cells.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param A Effective area. [m ** 2]
	 */
	void setLateralNozzleArea(double A);

	/*!
	 * \brief Sets the lateral nozzle area for all the cells.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 * 
	 * \param A Effective area. [m ** 2]
	 */
	void setLateralNozzleArea(RowVector A);

	/*!
	 * \brief Updates the state vector with the already-computed results.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * After solving the system and computing the state vector update, this
	 * function is called to update it.
	 */
	virtual void UpdateStateVector();
};

/*!
 * \brief Creates a new volute.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * \date 2016/03/21
 * 
 * \param node xml_node with the volute data.
 * \param fluid Pointer to the species vector.
 * \param MDB Materials database.
 */
Volute_ptr create_volute(const xml_node & node,
	const TComponentArray_ptr& components,
	const std::map<string, TSolid*> &MDB);

#endif
