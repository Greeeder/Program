/**
 * @file TCrankMechanism.cpp
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
#include "TCrankMechanism.h"

TCrankMechanism::TCrankMechanism() {


}

TCrankMechanism::TCrankMechanism(double speed){
	FCurrentRotSpeed = speed;
	FPreviousRotSpeed = speed;
	FPreviousAngle = 0.;
	FCurrentTime = 0.;
	FPreviousTime = 0.;
	FCurrentRotAccel = 0.;
	FTimeStep = 0.;


}

TCrankMechanism::TCrankMechanism(TCrankMechanism* CM) {

	FCycleDuration = CM->FCycleDuration;

	FCompressionRatio = CM->FCompressionRatio;
	FVolumeCC = CM->FVolumeCC;
	FVolumeDisplaced = CM->FVolumeDisplaced;

	FCrankRadius = CM->FCrankRadius;
	FConnectingRodLength = CM->FConnectingRodLength;
	FEccentricity = CM->FEccentricity;
	FMaxDistance = CM->FMaxDistance;
	FStroke = CM->FStroke;

	FPressureCoefficient = CM->FPressureCoefficient;
	FInertiaCoefficient = CM->FInertiaCoefficient;

	FDiameter = CM->FDiameter;
	FArea = CM->FArea;
	FPerimeter = CM->FPerimeter;

}


TCrankMechanism::~TCrankMechanism() {
	
}

double TCrankMechanism::getAcceleration(double angle, double rotspeed){
	//! \todo Define the acceleration of the piston
	return 0.;
}

double TCrankMechanism::getAreaCyl(){
	return FPerimeter * FPistonDisplacement;
}

double TCrankMechanism::getDeformation(double pressure, double acceleration){
	return pressure * FPressureCoefficient - acceleration * FInertiaCoefficient;
}

double TCrankMechanism::getDisplacement(double angle){
	return FMaxDistance - (FCrankRadius * cos(angle) + sqrt(pow2(FConnectingRodLength) - pow2(FCrankRadius *
		sin(angle) - FEccentricity)));
}

double TCrankMechanism::getProjectedArea(double angle){

	double b = asin((FEccentricity - FCrankRadius * sin(angle)) / FConnectingRodLength);
	return FArea * FCrankRadius * sin(angle - b) / cos(b);
}

double TCrankMechanism::getVolume(double angle, double pressure){

	FCurrentRotAccel = (FCurrentRotSpeed - FPreviousRotSpeed) / FTimeStep;

	double acceleration = getAcceleration(FCurrentAngle, FCurrentRotSpeed);

	FPistonDisplacement = getDisplacement(angle + FAngleOffset) + getDeformation(pressure, acceleration);

	FVolume = FVolumeCC + FArea * FPistonDisplacement;
	return FVolume;

}

double TCrankMechanism::getVolume(double pressure){

	FVolume = getVolume(FCurrentAngle, pressure);

	return FVolume;

}

double TCrankMechanism::getVolumeStatic(double angle){

	FPistonDisplacement = getDisplacement(angle + FAngleOffset);
	return FVolumeCC + FArea * FPistonDisplacement;


}

void TCrankMechanism::ReadInputDataXML(xml_node node_geom, string EngineType){

	FDiameter = GetAttributeAsDouble(node_geom, "Diameter");
	FPerimeter = __cons::Pi * FDiameter;
	FArea = __geom::Circle_area(FDiameter);
	FCrankRadius = GetAttributeAsDouble(node_geom, "CrankRadius");
	FConnectingRodLength = GetAttributeAsDouble(node_geom, "ConnectingRodLength");
	FEccentricity = GetAttributeAsDouble(node_geom, "Eccentricity");
	FStroke = sqrt(pow2(FConnectingRodLength + FCrankRadius) - pow2(FEccentricity)) - 
		sqrt(pow2(FConnectingRodLength - FCrankRadius) - pow2(FEccentricity));
	FVolumeDisplaced = __geom::Cylinder_volume(FDiameter, FStroke);
	FCompressionRatio = GetAttributeAsDouble(node_geom, "CompressionRatio");
	FVolumeCC = FVolumeDisplaced / (FCompressionRatio - 1);
	FPressureCoefficient = GetAttributeAsDouble(node_geom, "PressureCoefficient");
	FInertiaCoefficient = GetAttributeAsDouble(node_geom, "InertiaCoefficient");
	FAngleOffset = asin(FEccentricity / (FCrankRadius + FConnectingRodLength));
	FMaxDistance = sqrt(pow2(FCrankRadius + FConnectingRodLength) - pow2(FEccentricity));
	if (EngineType == "4Strokes"){
		FCycleDuration = 4 * __cons::Pi;
	}
	else{
		FCycleDuration = __cons::Pi_x_2;
	}

}

void TCrankMechanism::UpdateStateData()
{
	FPreviousAngle = FCurrentAngle;
	FPreviousTime = FCurrentTime;
	FPreviousRotSpeed = FCurrentRotSpeed;
	FVolume0 = FVolume;
}


