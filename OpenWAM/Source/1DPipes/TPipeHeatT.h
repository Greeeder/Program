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
 * @file TTubo.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
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
 * This file declares a finite differences pipe.
 */

//---------------------------------------------------------------------------
#ifndef TPipeHeatTH
#define TPipeHeatTH

#include <memory>
#include "Math_wam.h"
#include "TPipeExtHeatA.h"
#include "TPipeExtHeatW.h"
#include "TSldAl.h"

class TPipeHeatT {
  private:
	int nlayers;					   //!< The nlayers
	int ncells;						   //!< The ncells
	int nnodes;						   //!< The nnodes

	double CellSize;				   //!< Size of the cell
	double TimeStep;				   //!< The time step
	double InternalHeatMult;		   //!< The internal heat multiplier

	bool FixedTimeStep;				   //!< true to fixed time step

	std::vector<TSolid *> Material;		//!< Material object of the layer
	std::vector<double> Thickness;	   //!< The thickness of the layer [m]
	std::vector<bool> IsFluid;		   //!< true if the layer is fluid

	VectorXd Rint;					   //!< The internal radious of the nodes [m]
	VectorXd Rext;					   //!< The external radious of the nodes [m]
	VectorXd Conductance;			   //!< The conductance of the nodes [W / K]

	VectorXd Twall1;				   //!< The temperature of the nodes at the current time step [K]
	VectorXd Twall0;				   //!< The temperature of the nodes at the previous time step [K]

	VectorXd Twallint;				   //!< The internal wall temperature of the pipe [K]
	VectorXd Twallext;				   //!< The external wall temperature of the pipe [K]

	VectorXd Qint;					   //!< The heat through the internal wall [J]
	MatrixXd K1;					   //!< The matrix conductances between wall nodes [W / K]
	MatrixXd C;						   //!< The matrix capacitances of the wall nodes [J / K]
	MatrixXd Cinv;					   //!< The inverse of the capacitances matrix
	MatrixXd K0;					   //!< The matrix conductances between wall nodes [W / K]
	VectorXd K_layer;				   //!< The layer

	TPipeExtHeat* ExtHeat;			   //!< The extent heat

  public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 */

	TPipeHeatT();

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 */

	~TPipeHeatT();

	/*!
	 * \brief	Reads heat transfer data.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	node_pipe	The node pipe.
	 * \param	ncells		Number of pipe cells.
	 * \param	cs			Cell size [m]
	 * \param	mdb			Data base with the material to be used in the model
	 */

	void ReadHeatTransferData(xml_node node_pipe, int ncells, double cs, std::map<string, TSolid*> mdb);

	/*!
	 * \brief	Initiates this object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 * 			
	 * \param	T	Initial temperature [K]
	 */

	void Initiate(double T);

	/*!
	 * \brief	Builds a matrix.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	D	The VectorXd to process.
	 */

	void BuildMatrix(VectorXd D);

	/*!
	 * \brief	Solves this object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 */

	void Solve();

	/*!
	 * \brief	Solves.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	dt	The dt.
	 */

	void Solve(double dt);

	/*!
	 * \brief	Solve explicit.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	dt	The dt.
	 */

	void SolveExplicit(double dt);

	/*!
	 * \brief	Adds an internal heat to 'q'.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	i	Zero-based index of the.
	 * \param	q	The int to process.
	 */

	void AddInternalHeat(int i, int q);

	/*!
	 * \brief	Adds an internal heat.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	05/03/2016
	 *
	 * \param	q	The VectorXd to process.
	 */

	void AddInternalHeat(VectorXd q);

	/*!
	 * \brief	Gets the internal wall temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	21/03/2016
	 *
	 * \return	The internal wall temperature [K].
	 */

	VectorXd getTwallint(){
		return Twallint;
	}

	/*!
	 * \brief	Gets the external wall temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	21/03/2016
	 *
	 * \return	The external wall temperature [K].
	 */

	VectorXd getTwallext(){
		return Twallext;
	}

	/*!
	 * \brief	Adds an extent heat object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	21/03/2016
	 *
	 * \param [in,out]	exHeat	If non-null, the external heat object.
	 */

	void AddExtHeatObject(TPipeExtHeat* exHeat){
		ExtHeat = exHeat;
	}

	/*!
	 * \brief	Gets external diameter of the pipe.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	21/03/2016
	 *
	 * \return	The external diameter [m].
	 */

	VectorXd getExtDiameter(){
		return Rext.tail(ncells);
	}

	/*!
	 * \brief	Gets the temperature of the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	21/03/2016
	 *
	 * \return	The temperature [K].
	 */

	MatrixXd getTnodes(){
		return Twall1;
	}
};

typedef std::shared_ptr<TPipeHeatT> TPipeHeatT_ptr;

#endif

