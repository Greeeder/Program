/**
 * @file TTCHTM2.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 22 de feb. de 2016
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
 * TODO Insert here the description.
 */
#ifndef SOURCE_TURBOCOMPRESSOR_TTCHTM2_H_
#define SOURCE_TURBOCOMPRESSOR_TTCHTM2_H_

#include "Globales.h"
#include "turbo_bearings.hpp"


struct TCConvCorr {

	dVector K;						   //!< Vector with the correlation constants

	/*!
	 * \fn	TCConvCorr()
	 *
	 * \brief	Default constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 */

	TCConvCorr() {

	}

	/*!
	 * \fn	double Nusselt(dVector Input)
	 *
	 * \brief	Nusselt for a given input.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	Input	The input parameters (Re, Pr, mu ...).
	 *
	 * \return	Nusselt number [-].
	 */

	double Nusselt(dVector Input) {
		double nu = K[0];
		for(int i = 1; i < K.size(); i++) {
			nu = nu * pow(Input[i - 1], K[i]);
		}
		return nu;
	}
};

struct TurboMachine {
	double m = 0;					   //!< The mass flow [kg/s]
	double pr = 0;					   //!< The pressure ratio [-]
	double ip = 0;					   //!< The inlet pressure [Pa]
	double ip0 = 0;					   //!< The total inlet pressure [Pa]
	double op = 0;					   //!< The outlet pressure [Pa]
	double op0 = 0;					   //!< The total outlet pressure [Pa]
	double eff = 0;					   //!< The efficiency [-]
	double aeff = 0;				   //!< The adiabatic efficiency [-]
	double rtc = 0;					   //!< The rotational speed [rpm]
	double powr = 0;				   //!< The real power [W]
	double pows = 0;				   //!< The isentropic power [W]
	double it = 0;					   //!< The inlet temperature [K]
	double it0 = 0;					   //!< The total inlet temperature [K]
	double ait0 = 0;				   //!< The total adiabatic inlet temperature [K]
	double ot = 0;					   //!< The outlet temperature [K]
	double ot0 = 0;					   //!< The total outlet temperature [K]
	double aot = 0;					   //!< The adiabatic outlet temperature [K]
	double dit = 0;					   //!< The diabatic inlet temperature [K]
	double ots = 0;					   //!< The isentropic outlet temperature [K]
	double din = 0;					   //!< The inlet diameter [m]
	double sin = 0;					   //!< The inlet section [m ** 2]
	double dout = 0;				   //!< The outlet diameter [m]
	double sout = 0;				   //!< The outlet secction [m ** 2]
	double g = 0;					   //!< The specific heat reatio [-]
	double Cp = 0;					   //!< The constant pressure specific heat [J / (kg K)]
	double sig = 0;					   //!< The process (Expansion = -1, Compression = 1)

	/*!
	 * \fn	double get_it0()
	 *
	 * \brief	Gets the total inlet temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	The total inlet temperature [K].
	 */

	double get_it0() {
		return it + pow2(m * 287 * it / ip / sin) / 2 / Cp;
	}

	/*!
	 * \fn	double get_ts()
	 *
	 * \brief	Gets the isentropic outlet temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	The isentropic outlet temperature [K].
	 */

	double get_ts() {
		return it0 * pow(pr, sig * (g - 1) / g);
	}

	/*!
	 * \fn	double get_pows()
	 *
	 * \brief	Gets the isentropic power.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	The isentropic power [W].
	 */

	double get_pows() {
		return sig * m * Cp * (ots - it0);
	}

	/*!
	 * \fn	double get_powr()
	 *
	 * \brief	Gets the real power.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	The real power [W].
	 */

	double get_powr() {
		return sig * m * Cp * (ots - it0) * eff;
	}

	/*!
	 * \fn	double get_ot0()
	 *
	 * \brief	Gets total outlet temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	The total outlet temperature [K].
	 */

	double get_ot0() {
		return it0 + sig * pows / m / Cp / eff;
	}

	/*!
	 * \fn	double get_t_static(double T0, double P0, double M, double A, double Cp, double g)
	 *
	 * \brief	Gets the static temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T0	The total temperature [K].
	 * \param	P0	The total pressure [Pa].
	 * \param	M 	The mass flow [kg/s].
	 * \param	A 	The pipe section [m ** 2].
	 * \param	Cp	The constant pressure specific heat [J / (kg K)].
	 * \param	g 	The specific heat ratio [-].
	 *
	 * \return	The static temperature [K].
	 */

	double get_t_static(double T0, double P0, double M, double A, double Cp, double g) {
		double T;
		double p = P0;
		double error = 1;
		int iter = 0;
		do {
			T = T0 - pow2(M * 287 * T0 / P0 / A) / 2 / Cp;
			error = get_p_static(P0, T, T0, g) - p;
			p = p + error;
			iter++;
		} while(iter < 100 && fabs(error) < 0.001);
		return T;
	}

	/*!
	 * \fn	double get_p_static(double P0, double T, double T0, double g)
	 *
	 * \brief	Gets the static pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	P0	The total pressure [Pa].
	 * \param	T 	The static temperature [K].
	 * \param	T0	The total temperature [K].
	 * \param	g 	The specific heat ratio [-].
	 *
	 * \return	The static pressure [Pa].
	 */

	double get_p_static(double P0, double T, double T0, double g) {
		return P0 * pow(T / T0, g / (g - 1));
	}


};

struct stTurPower {
	double Power;					   //!< The power [W]
	double m;						   //!< The mass flow [kg / s]
	double efft;					   //!< The efficiency [-]
	double Cp;						   //!< The specific heat C_p [J / (kg K)]
	double Tin;						   //!< The inlet temperature [K]
	double ga;						   //!< The specific heat ratio [-]

	/*!
	* \fn	stTurPower(double pow, double ma, double efftur, double Ti, double Cpgas, double gam)
	*
	* \brief	Constructor.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	07/12/2015
	*
	* \param	pow   	The power [W].
	* \param	ma	  	The mass flow [kg / s].
	* \param	efftur	The turbine efficiency [-].
	* \param	Ti	  	The inlet temperature [K].
	* \param	Cpgas 	The specific heat C_p [J / (kg K)].
	* \param	gam   	The specific heat ratio [-].
	*/

	stTurPower(double pow, double ma, double efftur, double Ti, double Cpgas, double gam) {
		Power = pow;
		m = ma;
		efft = efftur;
		Tin = Ti;
		Cp = Cpgas;
		ga = gam;

	}

	/*!
	* \fn	double operator()(const double er)
	*
	* \brief	 Calculate the difference between the turbine power and the calculated one as a function of the expansion ratio.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	07/12/2015
	*
	* \return	The error [W].
	*/

	double operator()(const double er) {
		return Power - efft * m * Cp * Tin * (1 - pow(er, -ga));
	}

};

/*!
* \struct	stCompressorPower
*
* \brief	Function for compressor power calculation.
*
* \author	F.J. Arnau (farnau@mot.upv.es)
* \date	07/12/2015
*/

struct stComPower {
	double Power;					   //!< The power [W]
	double m;						   //!< The mass flow [kg / s]
	double effc;					   //!< The compressor efficiency [-]
	double Cp;						   //!< The specific heat C_p [J / (kg K)]
	double Tin;						   //!< The inlet temperature [K]
	double ga;						   //!< The specific heat ratio [-]

	/*!
	* \fn	stComPower(double pow, double ma, double effcom, double Ti, double Cpgas,
	* 		double gam)
	*
	* \brief	Constructor.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	07/12/2015
	*
	* \param	pow   	The power [W].
	* \param	ma	  	The mass flow [kg / s].
	* \param	effcom	The compressor efficiency [-].
	* \param	Ti	  	The inlet temperature [K].
	* \param	Cpgas 	The specific heat C_p [J / (kg K)].
	* \param	gam   	The specific heat ratio [-].
	*/

	stComPower(double pow, double ma, double effcom, double Ti, double Cpgas, double gam) {
		Power = pow;
		m = ma;
		effc = effcom;
		Tin = Ti;
		Cp = Cpgas;
		ga = gam;
	}

	/*!
	* \fn	double operator()(const double cr)
	*
	* \brief	 Difference between the compressor power and the power calculated as a function of the compressor ratio.
	*
	* \author	F.J. Arnau (farnau@mot.upv.es)
	* \date	07/12/2015
	*
	* \return	The error [W].
	*/

	double operator()(const double cr) {
		return Power - m * Cp * Tin * (pow(cr, ga) - 1) / effc;
	}

};


class TTC_HTM2 {

  private:

	bool FIsWaterCooled;			   //!< true if the turbocharger is water cooled
	bool FInitialTemp;				   //!< true if initial metal temperatures are imposed by the use
	int nnodes;						   //!< Number of metal nodes
	int nboundy;					   //!< Number of boundaries
	int nparts;						   //!< Number of parts that compose the turbocharger

	double TurbMaxEfficiency;		   //!< The turbine maximum efficiency

	std::map<string, int> MetalNode;	//!< The metal nodes
	std::map<string, int> Boundary;	   //!< The boundaries
	std::map<string, int> Part;		   //!< The parts

	std::vector<string> NodeLabel;	   //!< The metal node labels
	std::vector<string> BoundaryLabel;  //!< The boundary labels
	std::vector<string> PartLabel;	   //!< The part labels


	std::map<string, TFluid*> Fluid;	//!< The array of fluids that exchange heat with the turbocharger

	std::map<string, std::map<string, TCConvCorr* >>
			ConvectiveNusselt; //!< The arrary of convective correlations to determine the Nusselt

	std::map<string, double> Emissivity;	//!< The emissivity [-]
	ArrayXd ExternalDiameter;		   //!< The external diameter [m]
	ArrayXd ExternalLength;			   //!< The external length [m]
	ArrayXd ExternalArea;			   //!< The external surface [m ** 2]

	double FFitTurbExtHeat;			   //!< Multiplier for external turbine diameter [-]
	bool FExternalRadiation;		   //!< true if external radiation is considered
	bool FTurbineShielded;			   //!< true if turbine is shielded
	double FFitRadiation;			   //!< Multiplier for the external radiation [-]

	bool FExternalConvection;		   //!< true if external convection is considered
	double FFitConvection;			   //!< Multiplier for the external convection
	double FExtVelocity;			   //!< External air velocity [m / s]
	double FTurbExtVelocity;		   //!< External air velocity in the turbine [m / s].
	double FHousExtVelocity;	       //!< External air velocity in the housing [m / s].
	double FCompExtVelocity;	       //!< External air velocity in the compressor [m / s].

	double FTExtHeatCoefficient;	   //!< The turbine external heat coefficient [W / (m ** 2 K)]
	double FHExtHeatCoefficient;	   //!< The housing external heat coefficient [W / (m ** 2 K)]
	double FCExtHeatCoefficient;	   //!< The compressor external heat coefficient [W / (m ** 2 K)]

	bool FImposedExtHeatCoefficient;	//!< true if external heat coeficient is imposed

	ArrayXd FAlpha;					   //!< Housing geometry distribution.

	MatrixXd Kcond;					   //!< Matrix with the conductive conductances between metal nodes [W / K]
	MatrixXd Krad;					   //!< Matrix with the radiative conductances between metal nodes [W / K]
	MatrixXd Kconv;					   //!< Matrix with the convective conductances between metal nodes and boundaries [W / K]
	MatrixXd Kradext;				   //!< Matrix with the radiative conductances between metal nodes and boundaries [W / K]
	MatrixXd Kfull1;				   //!< Full matrix with conductances [W / K]
	MatrixXd Kfull0;				   //!< Full matrix with condcutances for previous time step [W / K]
	MatrixXd VF_nodes;				   //!< Matrix with the view factors between nodes [-]
	MatrixXd VF_toBo;				   //!< Matrix with the view factors between nodes and boundaries [-]
	MatrixXd C;						   //!< Matrix with the capacitances of metal nodes [J / K]
	MatrixXd Cinv;					   //!< Inverse of the matrix with the capacitances of metal nodes [J / K]
	MatrixXd Qext_J;				   //!< Convective accumulated heat from the boundaries [J]
	MatrixXd Qrad_J;				   //!< Radiative accumulated heat from the boundaries [J]
	MatrixXd Qext_W;				   //!< Convective heat power from the boundaries [W]
	MatrixXd Qrad_W;				   //!< Radiative heat power from the boundaries [W]

	MatrixXd MetalConductances;		   //!< Matrix with the metal conductances between nodes [W / K]
	MatrixXd MetalCapacitances;		   //!< Matrix with the metal capacitances of the nodes [J / K]

	MatrixXi ConvectiveConn;		   //!< Matrix for convective relations (1 = connected)
	MatrixXi RadiativeConn;			   //!< Matrix for radiative relations (1 = connected)

	VectorXd Tinit;					   //!< Initial temperatures of the matel nodes [K]
	VectorXd T1;					   //!< Temperature of the metal nodes at current time step [K]
	VectorXd T0;					   //!< Temperature of the metal nodes at previous time step [K]
	VectorXd Tb;					   //!< Temperature of the boundaries [K]

	ArrayXd T_oil;					   //!< Oil temperature (0=inlet 2=outlet) [K]
	double m_oil;					   //!< Oil mass flow [kg / s]
	double T_oilH2;					   //!< Oil temperature at H2 position [K]
	double T_oilMech;				   //!< Oil temperature for mechanical losses [K]

	double m_cool;					   //!< Coolant mass flow [kg / s]
	ArrayXd T_cool;					   //!< Coolant temperature (0=inlet, 1=outlet)

	TurboBearings *FMechLosses;		   //!< The mechanical losses object

	VectorXd D_in;					   //!< Inlet ports diameter [m]
	VectorXd D_out;					   //!< Outlet ports diameter [m]

	/*!
	 * \fn	double NusseltFreeConv(double Gr, double Pr);
	 *
	 * \brief	Nusselt number for free convection arround an isotherm cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	Gr	The Grashof number [-].
	 * \param	Pr	The Prandtl number [-].
	 *
	 * \return	The Nusselt number [-].
	 */

	double NusseltFreeConv(double Gr, double Pr);

	/*!
	 * \fn	double NusseltForcConv(double Re, double Pr);
	 *
	 * \brief	Nusselt number for force convection arround an isotherm cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	Re	The Reynolds number [-].
	 * \param	Pr	The Prandtl number [-].
	 *
	 * \return	The Nusselt number [-].
	 */

	double NusseltForcConv(double Re, double Pr);

	/*!
	 * \fn	double ViewFactorDiskDisk(double r1, double r2, double rc, double h);
	 *
	 * \brief	View factor between two rings.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	r1	External radio of the ring 1 [m].
	 * \param	r2	External radio of the ring 2 [m].
	 * \param	rc	Internal radio of both rings [m].
	 * \param	h 	Distance between both rings [m].
	 *
	 * \return	View factor [-].
	 */

	double ViewFactorDiskDisk(double r1, double r2, double rc, double h);

	/*!
	 * \fn	double ViewFactorDiskCylinder(double r1, double r2, double h);
	 *
	 * \brief	View factor between a ring and a cylinder.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	r1	The enternal radius of the ring [m].
	 * \param	r2	The cylinder radius [m].
	 * \param	h 	The height of the cylinder [m].
	 *
	 * \return	View factor [-].
	 */

	double ViewFactorDiskCylinder(double r1, double r2, double h);

	/*!
	 * \fn	double KRadiationExternal(double L, double D, double e, double F, double T, double Tamb);
	 *
	 * \brief	Radiation conductance from a cylinder to the ambient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	L   	Lenght of the cylinder [m].
	 * \param	D   	Diameter of the cylinder [m].
	 * \param	e   	Emissivity [-].
	 * \param	F   	View factor [-].
	 * \param	T   	Surface temperature [K].
	 * \param	Tamb	Ambient temperature [K].
	 *
	 * \return	Radiation conductance [W / K].
	 */

	double KRadiationExternal(double L, double D, double e, double F, double T, double Tamb);

	/*!
	 * \fn	double KRadiationExternalHous(double A, double e, double F, double T, double Tamb);
	 *
	 * \brief	Radiation conductance from turbocharger housing to the ambient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	A   	Area of the surface [m^2].
	 * \param	e   	Emissivity [-].
	 * \param	F   	View factor [-].
	 * \param	T   	Surface temperature [K].
	 * \param	Tamb	Ambient temperature [K].
	 *
	 * \return	Radiation conductance [W / K].
	 */

	double KRadiationExternalHous(double A, double e, double F, double T, double Tamb);

	/*!
	 * \fn	double KRadiationSurfaces(double A1, double A2, double e1, double e2, double F12,
	 * 		double T1, double T2);
	 *
	 * \brief	Calculate radiation conductance between two surfaces.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	A1 	Area of surface 1 [m^2].
	 * \param	A2 	Area of surface 2 [m^2].
	 * \param	e1 	Emissivity surface 1 [-].
	 * \param	e2 	Emissivity surface 2 [-].
	 * \param	F12	View factor from surface 1 to surface 2 [-].
	 * \param	T1 	Temperature surface 1 [K].
	 * \param	T2 	Temperature surface 2 [K].
	 *
	 * \return	Radiation conductance [W / K].
	 */

	double KRadiationSurfaces(double A1, double A2, double e1, double e2, double F12, double T1, double T2);

	/*!
	 * \fn	void PrintTurbineAdiMap(TurboMachine T, TurboMachine C);
	 *
	 * \brief	Print turbine adiabatic map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	The turbine data structure.
	 * \param	C	The compressor data structure.
	 */

	void PrintTurbineAdiMap(TurboMachine T, TurboMachine C);

	/*!
	 * \fn	void PrintCompressorAdiMap(TurboMachine T, TurboMachine C);
	 *
	 * \brief	Print compressor adiabatic map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	The turbine data structure.
	 * \param	C	The compressor data structure.
	 */

	void PrintCompressorAdiMap(TurboMachine T, TurboMachine C);

  public:

	/*!
	 * \fn	TTC_HTM2(TFluid *Oil);
	 *
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param [in,out]	Oil	fluid object.
	 */

	TTC_HTM2(TFluid *Oil);

	/*!
	 * \fn	virtual ~TTC_HTM2();
	 *
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 */

	virtual ~TTC_HTM2();

	/*!
	 * \fn	void Read_HTMXML(xml_node node_ht);
	 *
	 * \brief	Read the input data from a xml file.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	node_ht	xml node that contains the input data of the heat transfer model.
	 */

	void Read_HTMXML(xml_node node_ht);

	/*!
	 * \fn	void SolveNodeTemperatures(double dt, int mode);
	 *
	 * \brief	Solve node temperatures during a time step dt.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt  	The time step [s].
	 * \param	mode	The mode (0 = steady, 1 = transient).
	 */

	void SolveNodeTemperatures(double dt, int mode);

	/*!
	 * \fn	void SolveExplicit(double dt, int mode);
	 *
	 * \brief	Solve explicitly the node temperature at time step dt.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	dt  	Time step [s].
	 * \param	mode	The mode (0 = steady, 1 = transient).
	 */

	void SolveExplicit(double dt, int mode);

	/*!
	 * \fn	void SolveImplicit(double dt, int mode);
	 *
	 * \brief	Solve implicitly the node temperature at time step dt.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	dt  	Time step [s].
	 * \param	mode	The mode (0 = steady, 1 = transient).
	 */

	void SolveImplicit(double dt, int mode);

	/*!
	 * \fn	void setHeatFromGas(double dt, double T, double m);
	 *
	 * \brief	Calculate the heat exchange between the turbine gas and the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt	The time step [s].
	 * \param	T 	The gas temperature [K].
	 * \param	m 	The gas mass flow [kg/s].
	 */

	void setHeatFromGas(double dt, double T, double m);

	/*!
	 * \fn	void setKFromGas(double T, double m);
	 *
	 * \brief	Sets the convective conductance from gas.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Gas temperature [K].
	 * \param	m	Gas mass flow [kg / s].
	 */

	void setKFromGas(double T, double m);

	/*!
	 * \fn	void setHeatFromAir(double dt, double T, double m);
	 *
	 * \brief	Calculate the heat exchange between the compressor air and the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt	The time step [s].
	 * \param	T 	The air temperature [K].
	 * \param	m 	The air mass flow [kg/s].
	 */

	void setHeatFromAir(double dt, double T, double m);

	/*!
	 * \fn	void setKFromAir(double T, double m);
	 *
	 * \brief	Sets the convective conductance from air.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T	Air temperature [K].
	 * \param	m	Air mass flow [kg / s]
	 */

	void setKFromAir(double T, double m);

	/*!
	 * \fn	void setHeatFromOil(double dt, double T, double m, double mech_pow);
	 *
	 * \brief	Calculate the heat exchange between the oil and the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt			The time step [s].
	 * \param	T			The oil inlet temperature [K].
	 * \param	m			The oil mass flow [kg/s].
	 * \param	mech_pow	The mechanical losses in the shaft [W].
	 */

	void setHeatFromOil(double dt, double T, double m, double mech_pow);

	/*!
	 * \fn	void setKFromOil(double T, double m, double mech_pow);
	 *
	 * \brief	Sets convectiv conductance from oil.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T			The oil inlet temperature [K].
	 * \param	m			The oil mass flow [kg/s].
	 * \param	mech_pow	The mechanical losses in the shaft [W].
	 */

	void setKFromOil(double T, double m, double mech_pow);

	/*!
	 * \fn	void setHeatFromCoolant(double dt, double T, double m);
	 *
	 * \brief	Calculate the heat exchange between the coolant and the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt	The time step [s].
	 * \param	T 	The coolant inlet temperature [K].
	 * \param	m 	The coolant inlet mass flow [kg/s].
	 */

	void setHeatFromCoolant(double dt, double T, double m);

	/*!
	 * \fn	void setKFromCoolant(double T, double m);
	 *
	 * \brief	Sets convective conductance from coolant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T 	The coolant inlet temperature [K].
	 * \param	m 	The coolant inlet mass flow [kg/s].
	 */

	void setKFromCoolant(double T, double m);

	/*!
	 * \fn	void setHeatFromAmbient(double dt, double T, double p);
	 *
	 * \brief	Calculate the heat exchange between the ambient and the metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	24/02/2016
	 *
	 * \param	dt	The time step [s].
	 * \param	T 	The ambient temperature [K].
	 * \param	p 	The ambient pressure [Pa].
	 */

	void setHeatFromAmbient(double dt, double T, double p);

	/*!
	 * \fn	void setKFromAmbient(double T, double p);
	 *
	 * \brief	Sets convective conductances from ambient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T 	The ambient temperature [K].
	 * \param	p 	The ambient pressure [Pa].
	 */

	void setKFromAmbient(double T, double p);

	/*!
	 * \fn	void setRadiatioFromAmbient(double dt, double T);
	 *
	 * \brief	Sets radiative heat from ambient to metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	dt	Time step [s].
	 * \param	T 	Ambient temperature [K].
	 */

	void setRadiatioFromAmbient(double dt, double T);

	/*!
	 * \fn	void setKRadiatioFromAmbient(double T);
	 *
	 * \brief	Sets radiative conductances from ambient.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T 	Ambient temperature [K].
	 */

	void setKRadiatioFromAmbient(double T);

	/*!
	 * \fn	void setViewFactorBetweenNodes(double Dt, double Dc, double Dh, double L);
	 *
	 * \brief	Sets view factor between nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	Dt	Turbine external diameter [m].
	 * \param	Dc	Compressor external diameter [m].
	 * \param	Dh	Housing external diameter [m].
	 * \param	L 	Housing external length [m].
	 */

	void setViewFactorBetweenNodes(double Dt, double Dc, double Dh, double L);

	/*!
	 * \fn	void setViewFactorToBoundaries(double Dt, double Dc, double Dh, double L);
	 *
	 * \brief	Sets view factor to boundaries.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	Dt	Turbine external diameter [m].
	 * \param	Dc	Compressor external diameter [m].
	 * \param	Dh	Housing external diameter [m].
	 * \param	L 	Housing external length [m].
	 */

	void setViewFactorToBoundaries(double Dt, double Dc, double Dh, double L);

	/*!
	 * \fn	void setKCondMatrix();
	 *
	 * \brief	Sets conductive conductances matrix between metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 */

	void setKCondMatrix();

	/*!
	 * \fn	void setCapMatrix();
	 *
	 * \brief	Sets capacitances matrix of metal nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 */

	void setCapMatrix();

	/*!
	 * \fn	void setKRadMatrix();
	 *
	 * \brief	Sets radiative conductances matrix between nodes.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 */

	void setKRadMatrix();

	/*!
	 * \brief	Gets the ambient heat flow by convection.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientHeatFlowConv();

	/*!
	 * \brief	Gets the ambient to turbine heat flow by convection.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToTurbineHeatFlowConv();

	/*!
	 * \brief	Gets the ambient to housing heat flow by convection.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToHousingHeatFlowConv();

	/*!
	 * \brief	Gets the ambient to compressor heat flow by convection.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToCompressorHeatFlowConv();

	/*!
	 * \brief	Gets the ambient heat flow by radiation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientHeatFlowRad();

	/*!
	 * \brief	Gets the ambient to turbine heat flow by radiation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToTurbineHeatFlowRad();

	/*!
	 * \brief	Gets the ambient to housing heat flow by radiation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToHousingHeatFlowRad();

	/*!
	 * \brief	Gets the ambient to compressor heat flow by radiation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	29/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientToCompressorHeatFlowRad();

	/*!
	 * \brief	Gets the total ambient heat flow.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double AmbientHeatFlowTotal();

	/*!
	 * \brief	Gets the heat flow from the coolant.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double WaterHeatFlow();

	/*!
	 * \fn	double TurbineHeatFlow();
	 *
	 * \brief	Gets heat flow from the turbine gas.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double TurbineHeatFlow();

	/*!
	 * \fn	double CompressorHeatFlow();
	 *
	 * \brief	Gets heat flow from the compressor air.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double CompressorHeatFlow();

	/*!
	 * \fn	double OilHeatFlow();
	 *
	 * \brief	Gets heat flow from the oil.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \return	Heat flow [W].
	 */

	double OilHeatFlow();

	/*!
	 * \fn	double getAdiabaticEfficiency(string Case, TurboMachine T, TurboMachine C);
	 *
	 * \brief	Gets the adiabatic efficiency.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	Case	The case ("Turbine" or "Compressor").
	 * \param	T   	The turbine data structure.
	 * \param	C   	The compressor data structure.
	 *
	 * \return	The adiabatic efficiency [-].
	 */

	double getAdiabaticEfficiency(string Case, TurboMachine T, TurboMachine C);

	/*!
	 * \fn	void setOilConditions(double m, double T);
	 *
	 * \brief	Sets oil conditions.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	m	Mass flow [kg / s].
	 * \param	T	Temperature [K].
	 */

	void setOilConditions(double m, double T);

	/*!
	 * \fn	void setCoolantConditions(double m, double T);
	 *
	 * \brief	Sets coolant conditions.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	m	Mass flow [kg / s].
	 * \param	T	Temperature [K].
	 */

	void setCoolantConditions(double m, double T);

	/*!
	 * \fn	void TurbochargerData(double Doil, double DWater, double DT,
	 * 		double LT, double DC, double LC, double DH, double LH);
	 *
	 * \brief	Add turbocharger data.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	Doil  	The oil port diameter [m].
	 * \param	DWater	The water port diameter [m].
	 * \param	DT	  	The turbine external diameter [m].
	 * \param	LT	  	The turbine external length [m].
	 * \param	DC	  	The compressor external diameter [m].
	 * \param	LC	  	The compressor external length [m].
	 * \param	DH	  	The housing external diameter [m].
	 * \param	LH	  	The housing external length [m].
	 */

	void TurbochargerData(double Doil, double DWater, double DT, double LT, double DC, double LC,
						  double DH, double LH);

	/*!
	 * \fn	void AsignTCMechLosses(TurboBearings *MechLosses)
	 *
	 * \brief	Asign the mechanical losses object.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param [in,out]	MechLosses	If non-null, the mechanical losses object.
	 */

	void AsignTCMechLosses(TurboBearings *MechLosses) {
		FMechLosses = MechLosses;
	}

	/*!
	 * \fn	void PutEffTurbMax(double EffMax)
	 *
	 * \brief	Puts the turbine maximum efficiency.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	EffMax	The maximum efficiency [-].
	 */

	void PutEffTurbMax(double EffMax) {
		TurbMaxEfficiency = EffMax;
	}

	/*!
	 * \fn	void setDiameters(string boundary, double din, double dout)
	 *
	 * \brief	Sets ports diameters.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	boundary	The boundary label.
	 * \param	din			Inlet diameter [m].
	 * \param	dout		Outlet diameter [m].
	 */

	void setDiameters(string boundary, double din, double dout) {
		D_in(Boundary[boundary]) = din;
		D_out(Boundary[boundary]) = dout;
	}

	/*!
	 * \fn	double getTemperatureNode(string node)
	 *
	 * \brief	Gets the temperature of the node.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	node	The node label.
	 *
	 * \return	The temperature [K].
	 */

	double getTemperatureNode(string node) {
		return T1(MetalNode[node]);
	}

	/*!
	 * \fn	double getTemperatureNode(int index)
	 *
	 * \brief	Gets the temperature of the node.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	index	The node index.
	 *
	 * \return	The temperature [K].
	 */

	double getTemperatureNode(int index) {
		return T1(index);
	}

	/*!
	 * \fn	double getTemperatureOil(int i)
	 *
	 * \brief	Gets the oil temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	i	Index of the step (0 = inlet, 2 = outlet).
	 *
	 * \return	The temperature [K].
	 */

	double getTemperatureOil(int i) {
		return T_oil(i);
	}

	/*!
	 * \fn	double getTemperatureCoolant(int i)
	 *
	 * \brief	Gets the coolant temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	i	Index of the step (0 = inlet, 1 = outlet).
	 *
	 * \return	The temperature [K].
	 */

	double getTemperatureCoolant(int i) {
		return T_cool(i);
	}

	/*!
	 * \fn	double getT2avg(double T2, double m, double heat)
	 *
	 * \brief	Gets compressor temperature average between outlet and diffusor inlet.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param	T2  	Compressor outlet temperature [K].
	 * \param	m   	Air mass flow [kg / s].
	 * \param	heat	Heat power from the air to the metal node [W].
	 *
	 * \return	The temperature [K].
	 */

	double getT2avg(double T2, double m, double heat) {
		return T2 - heat / m / Fluid["AIR"]->FunCp(T2) / 2;
	}

	/*!
	 * \fn	void HeaderInsTemperatures(stringstream & insoutput, int i);
	 *
	 * \brief	Header for output temperatures.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param [in,out]	insoutput	Container for the instantaneous results.
	 * \param	i				 	Turbocharger index.
	 */

	void HeaderInsTemperatures(stringstream & insoutput, int i);

	/*!
	 * \fn	void PrintInsTemperatures(stringstream & insoutput);
	 *
	 * \brief	Output for instantaneous metal node temperatures.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param [in,out]	insoutput	Container for the instantaneous results.
	 */

	void PrintInsTemperatures(stringstream & insoutput);

	/*!
	 * \fn	void PrintInsHeatFlow(stringstream & insoutput);
	 *
	 * \brief	Output for instantaneous heat flows.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param [in,out]	insoutput	Container for the instantaneous results.
	 */

	void PrintInsHeatFlow(stringstream & insoutput);

	/*!
	 * \fn	void HeaderInsHeatFlow(stringstream & insoutput, int num);
	 *
	 * \brief	Header for the heat flows output.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	03/03/2016
	 *
	 * \param [in,out]	insoutput	Container for the instantaneous results.
	 * \param	num				 	Turbocharger index.
	 */

	void HeaderInsHeatFlow(stringstream & insoutput, int num);

	/*!
	 * \fn	void Put_ExternalVelocity(double val)
	 *
	 * \brief	Impose externally the external velocidy.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	val	Air velocity [m / s].
	 */

	void Put_ExternalVelocity(double val) {
		FExtVelocity = val;
	}

	/*!
	 * \fn	void Put_ExternalHeat(double val)
	 *
	 * \brief	Impose externally the multiplier for radiation and external convection heat transfer.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	07/12/2015
	 *
	 * \param	val	The multiplier [-].
	 */

	void Put_ExternalHeat(double val) {
		FFitRadiation = val;
		FFitConvection = val;
	}

	/*!
	 * \brief	Initializes the metal temperatures.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	16/03/2016
	 *
	 * \param	T3	  	Turbine inlet temperature [K].
	 * \param	T2	  	Compressor outlet temperature [K].
	 * \param	Toil  	Oil temperature [K].
	 * \param	Twater	Water temperature [K].
	 */

	void InitializeTemp(double T3, double T2, double Toil, double Twater);

};

#endif /* SOURCE_TURBOCOMPRESSOR_TTCHTM2_H_ */
