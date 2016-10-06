/*--------------------------------------------------------------------------------*\
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
 * @file Constantes.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 *
 * @section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU <<<<<<< HEAD
 General Public License as published by
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
 * This file defines several constants used in OpenWAM.
 */

//---------------------------------------------------------------------------
#ifndef ConstantesH
#define ConstantesH
//#include "TBloqueMotor.h"
#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include <limits>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Math_wam.h"

//---------------------------------------------------------------------------

/*!
 * \brief Gas constant
 */

namespace __R {
	const double Universal = 8314.4; //!< Universal gas constant [J / (kmol *K)]
	const double Air = 287.; //!< Gas constant for the air [J/(kg*K)]
	const double Fuel = 55.95; //!< Gas constant for the fuel [J/(kg*K)]
	const double Burnt = 285.4; //!< Gas constant for the burnt gas [J/(kg*K)]
	const double O2 = 259.825; //!< Gas constant for the O2 [J/(kg*K)]
	const double CO2 = 188.9207; //!< Gas constant for the CO2 [J/(kg*K)]
	const double H2O = 461.398; //!< Gas constant for the water [J/(kg*K)]
	const double N2 = 296.837; //!< Gas constant for the N2 [J/(kg*K)]
	const double Ar = 208.12; //!< Gas constant for the Ar [J/(kg*K)]
	const double Diesel = 55.95; //!< Gas constant for the diesel fuel [J/(kg*K)]
	const double Gasoline = 72.42; //!< Gas constant for the gasoline fuel [J/(kg*K)]
}


/*!
 * \brief Molecular weight
 */

namespace __PM {
	const double O2 = 32.; //!< Molar mass of the O2 [g/mol]
	const double CO2 = 44.01; //!< Molar mass of the CO2 [g/mol]
	const double H2O = 18.02; //!< Molar mass of the water [g/mol]
	const double N2 = 28.01; //!< Molar mass of the N2 [g/mol]
	const double Ar = 39.95; //!< Molar mass of the Ar [g/mol]
	const double C = 12.01; //!< Molar mass of the C [g/mol]
	const double NO2 = 46; //!< Molar mass of the NO2 [g/mol]
	const double NO = 30; //!< Molar mass of the NO [g/mol]
	const double CO = 28.01; //!< Molar mass of the CO [g/mol]
	const double UHC = 55.04; //!< Molar mass of the unburned fuel [g/mol]
	const double Diesel = 148.4; //!< Molar mass of the diesel fuel [g/mol]
	const double Gasoline = 114.8; //!< Molar mass of the gasoline fuel [g/mol]
	const double Air = 28.85; //!< Molar mass of the air (79% N2 + 21% O2) [g/mol]
	const double IdealAir = 28.97; //!< Molar mass of ideal air. [g/mol]
}

/*!
* \brief General constants
*/
namespace __cons {
	const double Pi = 3.14159265358979323846; //!< Pi number [-]
	const double Pi_2 = Pi / 2; //!< Pi / 2 [-]
	const double Pi_4 = Pi / 4; //!< Pi / 4 [-]
	const double Pi_x_2 = 2 * Pi; //!< 2 * Pi [-]
	const double _1_Pi = 1 / Pi; //!< 1 / Pi [-]
	const double _2_Pi = 2 / Pi; //!< 2 / Pi [-]
	const double _4_Pi = 4 / Pi; //!< 4 / Pi [-]
	const double SQR_4_Pi = sqrt(_4_Pi); //!< sqrt(4 / Pi) [-]
	const double Sigma = 5.670373e-8; //!< Stephan-Boltzman constant [W / (m^2 * K^4)]
	const double ARef = 343.11; //!< Reference speed of sound [m / s]
	const double ARef2 = ARef * ARef; //!< A_ref^2 [m^2 / s^2]
	const double TRef = 292.99271; //!< Reference temperature [K]
	const double PRef = 1.0; //!< Reference pressure [bar]
	const double Kb = 1.38054e-23; ///< Boltzmann constant. [J / K]
}

/*!
* \brief Units converters
*/
namespace __units {

	const double _DegToRad = 2 * __cons::Pi / 360;  //!< The degrees to radians
	const double _RadToDeg = 1 / _DegToRad; //!< The radians to degrees

	/*!
	 * \brief	Convert any parameter to kilo...
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Parameter.
	 *
	 * \return	Parameter in kilo...
	 */

	inline double To_kilo(double p) {
		return p * 0.001;
	}

	/*!
	 * \brief	Convert any parameter from kilo ...
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Parameter in kilo...
	 *
	 * \return	Parameter.
	 */

	inline double From_kilo(double p) {
		return p * 1000;
	}

	/*!
	 * \brief	Convert degree in radians.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Parameter [deg].
	 *
	 * \return	Parameter [rad].
	 */

	inline double DegToRad(double p) {
		return p * _DegToRad;
	}

	/*!
	* \brief	Convert radians in degrees.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	04/03/2016
	*
	* \param	p	Parameter [rad].
	*
	* \return	Parameter [deg].
	*/

	inline double RadToDeg(double p) {
		return p * _RadToDeg;
	}

	/*!
	 * \brief	Convert Bar to Pascals.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Pressure [bar].
	 *
	 * \return	Pressure [Pa].
	 */

	inline double BarToPa(double p) {
		return p * 1e5;
	}

	/*!
	 * \brief	Convert Pascals to bar.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Pressure [Pa].
	 *
	 * \return	Pressure [bar].
	 */

	inline double PaToBar(double p) {
		return p * 1e-5;
	}

	/*!
	 * \brief	Convert rotational speed from RPM to RPS.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Rotational speed [rpm].
	 *
	 * \return	Rotational speed [rps].
	 */

	inline double RPMToRPS(double p) {
		return p * 0.016666666666667;
	}

	/*!
	 * \brief	Convert rotational speed from RPM to rad/s.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Rotational speed [rpm].
	 *
	 * \return	Rotational speed [rad / s].
	 */

	inline double RPMToRad_s(double p) {
		return p * 0.104719755119660;
	}

	/*!
	 * \brief	Convert cotational speed from rad/s to rpm.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Rotational speed [rad / s].
	 *
	 * \return	Rotational speed [rpm].
	 */

	inline double Rad_sToRPM(double p) {
		return p * 9.549296585513721;
	}

	/*!
	 * \brief	Convert speed from m/s to km/h.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Speed [m / s].
	 *
	 * \return	Speed [kg / h].
	 */

	inline double m_sTokm_h(double p) {
		return p * 3.6;
	}

	/*!
	 * \brief	Convert temperature from Celsius to Kelvin.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Temperature [degC].
	 *
	 * \return	Temperature [K].
	 */

	inline double degCToK(double p) {
		return p + 273.15;
	}

	/*!
	 * \brief	Convert temperature form Kelvin to Celsius.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	p	Temperature [K].
	 *
	 * \return	Temperature [degC].
	 */

	inline double KTodegC(double p) {
		return p - 273.15;
	}
}


/*!
* \brief Enthalpy of formation for different chemical components
*/

namespace __HFormacion {
	const double CO2 = -393510; ///< Enthalpy of formation of CO2. [J / mol]
	const double CO = -110530; ///< Enthalpy of formation of CO. [J / mol]
	const double NO2 = 33100; ///< Enthalpy of formation of NO2. [J / mol]
	const double NO = 90290; ///< Enthalpy of formation of NO. [J / mol]
	const double H2O = -241830; ///< Enthalpy of formation of H2O. [J / mol]
}

/*!
* \brief Geometrical functios
*/

namespace __geom {

	/*!
	 * \brief	Calculate the area of a circle.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	d	Diameter [m].
	 *
	 * \return	The area [m ** 2].
	 */

	inline double Circle_area(double d) {
		return d * d * __cons::Pi_4;
	}

	/*!
	 * \brief	Calculate the diameter of a circle.
	 *
	 * \author	L.M. García-Cuevas (luiga12.upv.es)
	 * \date	15/03/2016
	 *
	 * \param	A Circle area [m ** 2].
	 *
	 * \return	Diameter [m].
	 */
	inline double Circle_diameter(double A) {
		return sqrt(A / __cons::Pi_4);
	}

	/*!
	 * \brief	Compute the volume of a cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	d	Diameter of a cylinder [m].
	 * \param	l	Length of the cylinder [m].
	 *
	 * \return	The volume [m ** 3].
	 */

	inline double Cylinder_volume(double d, double l) {
		return Circle_area(d) * l;
	}

	/*!
	 * \brief	Compute the area of a ring.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	din 	Inlet diameter [m].
	 * \param	dout	Outlet diameter [m].
	 *
	 * \return	The area [m ** 2].
	 */

	inline double Ring_area(double din, double dout) {
		return (dout * dout - din * din) * __cons::Pi_4;
	}
}

/*!
* \brief Ambient conditions
*/

namespace __Ambient {
	extern double p_Pa;				//!< The ambient pressure [Pa]
	extern double p_bar;			//!< The ambient pressure [bar]
	extern double T_K;				//!< The ambient temperature [K]
	extern double T_degC;			//!< The ambient temperature [degC]
	extern double HR;				//!< The humidity [%]
	extern RowVector Y_amb;			//!< The ambient mass fraction [-]
}

/*!
* \brief Specific heat ratio and derived coeficients
*/

namespace __Gamma {

// FOR AIR AT AMBIENT CONDITIONS
	const double Cp = 1004.5;		   //!< The specific heat \f$C_p\f$ [J / (kg K)]
	const double G = Cp / (Cp - __R::Air);  //!< The specific heat ratio \f$\gamma\f$ [-]
	const double G_1 = G - 1;		   //!< \f$\gamma - 1\f$
	const double G_2 = G + 1;		   //!< \f$\gamma + 1\f$
	const double G_3 = (G - 1) / 2;	   //!< \f$\frac{\gamma - 1}{2}\f$
	const double G_4 = 2 * G / (G - 1); //!< \f$\frac{2 \gamma}{\gamma - 1}\f$
	const double G_5 = (G - 1) / G / 2; //!< \f$\frac{\gamma - 1}{2 \gamma}\f$
	const double G_6 = 1 / (G - 1);	   //!< \f$\frac{1}{\gamma - 1}\f$
	const double G_7 = (3 - G) / (G + 1);   //!< \f$\frac{3 - \gamma}{\gamma + 1}\f$ 
	const double G_8 = (G - 1) / G;	   //!< \f$\frac{\gamma - 1}{\gamma}\f$ 
	const double G_9 = G / (G - 1);	   //!< \f$\frac{\gamma}{\gamma - 1}\f$ 
	const double Cp_x2 = 2 * Cp;	   //!< \f$2 C_p\f$
	const double Cv = Cp - __R::Air;	//!< \f$C_v = C_p - R\f$
	const double gxR = G * __R::Air;	//!< \f$\gamma R\f$

	/*!
	 * \fn	inline double GG(double Cp, double Cv)
	 *
	 * \brief	Return specific heat ratio.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	Cp	Specific heat at constant pressure [J / (kg K)].
	 * \param	Cv	Specific heat at constant volume [J / (kg K)].
	 *
	 * \return	Specific heat ratio [-].
	 */

	inline double GG(double Cp, double Cv) {
		return Cp / Cv;
	}

	/*!
	 * \fn	inline double G1(double g)
	 *
	 * \brief	Get \f$\gamma - 1\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\gamma - 1\f$.
	 */

	inline double G1(double g) {
		return g - 1;
	}

	/*!
	 * \fn	inline double G2(double g)
	 *
	 * \brief	Get \f$\gamma + 1\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\gamma + 1\f$.
	 */

	inline double G2(double g) {
		return g + 1;
	}

	/*!
	 * \fn	inline double G3(double g)
	 *
	 * \brief	Get \f$\frac{\gamma - 1}{2}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{\gamma - 1}{2}\f$.
	 */

	inline double G3(double g) {
		return (g - 1) * 0.5;
	}

	/*!
	 * \fn	inline double G4(double g)
	 *
	 * \brief	Get \f$\frac{2 \gamma}{\gamma - 1}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{2 \gamma}{\gamma - 1}\f$.
	 */

	inline double G4(double g) {
		return 2. * g / (g - 1);
	}

	/*!
	 * \fn	inline double G5(double g)
	 *
	 * \brief	Get \f$\frac{\gamma - 1}{2 \gamma}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{\gamma - 1}{2 \gamma}\f$.
	 */

	inline double G5(double g) {
		return (g - 1) / 2. / g;
	}

	/*!
	 * \fn	inline double G6(double g)
	 *
	 * \brief	Get \f$\frac{1}{\gamma - 1}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{1}{\gamma - 1}\f$.
	 */

	inline double G6(double g) {
		return 1 / (g - 1);
	}

	/*!
	 * \fn	inline double G7(double g)
	 *
	 * \brief	Get \f$\frac{3 - \gamma}{\gamma + 1}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{3 - \gamma}{\gamma + 1}\f$.
	 */

	inline double G7(double g) {
		return (3 - g) / (g + 1);
	}

	/*!
	 * \fn	inline double G8(double g)
	 *
	 * \brief	Get \f$\frac{\gamma - 1}{\gamma}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{\gamma - 1}{\gamma}\f$.
	 */

	inline double G8(double g) {
		return (g - 1) / g;
	}

	/*!
	 * \fn	inline double G9(double g)
	 *
	 * \brief	Get \f$\frac{\gamma}{\gamma - 1}\f$ as a function of \f$\gamma\f$.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	04/03/2016
	 *
	 * \param	g	Specific heat ratio [-].
	 *
	 * \return	\f$\frac{\gamma}{\gamma - 1}\f$.
	 */

	inline double G9(double g) {
		return g / (g - 1);
	}
}

#endif
