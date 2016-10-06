/**
 * @file TCrankMechanism.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 30 de mar. de 2016
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
 * This file includes the methods to solve the crank mechanism (Kinematics and dynamics) in the cylinder.
 */


#ifndef SOURCE_ENGINE_TCRANKMECHANISM_H_
#define SOURCE_ENGINE_TCRANKMECHANISM_H_

#include <memory>
#include "Math_wam.h"
#include "Globales.h"

class TCrankMechanism {

private:

	double FCurrentAngle;			   //!< The current angle [rad]
	double FCurrentRotSpeed;		   //!< The current rotational speed [rad / s]
	double FCurrentRotAccel;		   //!< The current rotational acceleration [rad / s^2]
	double FCurrentTime;			   //!< The current time [s]

	double FPreviousAngle;			   //!< The previous angle [rad]
	double FPreviousRotSpeed;		   //!< The previous rotational speed [rad / s]
	double FPreviousTime;			   //!< The previous time [s]

	double FTimeStep;				   //!< The time step [s]

	double FCycleDuration;			   //!< Duration of the cycle (2 or 4 strokes) [rad]

	double FCompressionRatio;		   //!< The compression ratio of the cylinder [-]
	double FVolumeCC;				   //!< The volume of the combustion chamber [m^3]
	double FVolumeDisplaced;		   //!< The total volume displaced [m^3]

	double FCrankRadius;			   //!< The crank radius [m]
	double FConnectingRodLength;	   //!< The connecting rod lenght [m]
	double FEccentricity;			   //!< The eccentricity [m]
	double FMaxDistance;			   //!< The maximum vertical distance between the bulon and the crankshaft center [m]
	double FStroke;					   //!< The crank diameter [m]
	double FAngleOffset;

	double FPistonDisplacement;		   //!< The piston displacement from the TDC [m]

	double FPressureCoefficient;	   //!< The pressure coefficient for deformations [m / Pa]
	double FInertiaCoefficient;		   //!< The inertia coefficient for deformations [s^2]

	double FDiameter;				   //!< The diameter of the cylinder [m]
	double FArea;					   //!< The transversal area of the cylinder [m^2]
	double FPerimeter;				   //!< The perimeter of the cylinder [m]

	double FVolume0;				   //!< The previous cylinder volume [m^3]
	double FVolume;					   //!< The current cylinder volume [m^3]

	/*!
	 * \brief	Gets the acceleration.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \param	angle   	The angle [rad].
	 * \param	rotspeed	The rotational speed [rad / s].
	 *
	 * \return	The acceleration [m / s^2].
	 */

	double getAcceleration(double angle, double rotspeed);

	/*!
	 * \brief	Gets the piston displacement.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \param	angle	The angle [rad].
	 *
	 * \return	The displacement [m].
	 */

	double getDisplacement(double angle);

	/*!
	 * \brief	Gets the deformation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \param	pressure		The pressure [Pa].
	 * \param	acceleration	The acceleration [m / s^2].
	 *
	 * \return	The deformation [m].
	 */

	double getDeformation(double pressure, double acceleration);

public:

	/*!
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 */

	TCrankMechanism();

	TCrankMechanism(double speed);

	/*!
	 * \brief	Copy constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \param [in,out]	CM	If non-null, the original object.
	 */

	TCrankMechanism(TCrankMechanism *CM);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 */

	virtual ~TCrankMechanism();

	/*!
	 * \brief	Gets wet cylinder area for heat transfer.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \return	The area [m^2].
	 */

	double getAreaCyl();

	double getCurrentAngle(){
		return FCurrentAngle;
	}

	double getPreviousAngle(){
		return FPreviousAngle;
	}

	void setCurrentAngle(double angle){
		FCurrentAngle = angle;
	}

	void setPreviousAngle(double angle){
		FPreviousAngle = angle;
	}

	void setTime(double dt){
		FTimeStep = dt;
		FPreviousTime = FCurrentTime;
		FCurrentTime += dt;
	}



	/*!
	 * \brief	Gets the diameter of the cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \return	The diameter [m].
	 */

	double getDiameter(){ return FDiameter; }

	/*!
	 * \brief	Gets the cylinder volume.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	30/03/2016
	 *
	 * \param	angle   	The angle [rad].
	 * \param	time		The time [s].
	 * \param	pressure	The pressure [Pa].
	 *
	 * \return	The volume [m^3].
	 */

	double getVolume(double angle, double pressure);

	double getVolume(double pressure);

	double getVolumeStatic(double angle);

	double getVolume(){ return FVolume; }

	void setVolume(double v){
		FVolume = v;
	}

	double getVolume0(){ return FVolume0; }

	void setVolume0(double v){
		FVolume0 = v;
	}

	double getVolumeDisplaced(){ return FVolumeDisplaced; }

	double getProjectedArea(double angel);

	double getRotationalSpeed(){
		return FCurrentRotSpeed;
	}

	void setRotationalSpeed(double value){
		FPreviousRotSpeed = FCurrentRotSpeed;
		FCurrentRotSpeed = value;
	}

	void ReadInputDataXML(xml_node node_geometry, string EngineType);

	void UpdateStateData();

};

typedef std::unique_ptr<TCrankMechanism> TCrankMechanism_ptr;

#endif /* SOURCE_ENGINE_TCRANKMECHANISM_H_ */
