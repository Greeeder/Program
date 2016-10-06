/**
* @file TSAEMap.h
* @author Francisco Jose Arnau <farnau@mot.upv.es>
* @date 19 de mar. de 2016
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
* Compressor map in SAE format.
*/

#ifndef TSAEMapH
#define TSAEMapH

#include <vector>
#include <cstdio>

#include "TCompressorMap.h"
#include <ExtrapMapFunctions.h>

// ---------------------------------------------------------------------------

/*!
 * \brief	Compressor maps in SAE format. Includes different methods to extrapolate, interpolate and obtain 
 * 			the adiabatic efficiency.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	19/03/2016
 */

class TSAEMap: public TCompressorMap {
  private:
	int FNumeroCompresor;			   //!< The number of the compressor

	double FTempRef;				   //!< The reference temperature [K]
	double FPresionRef;				   //!< The reference pressure [bar]

	bool FIsAdiabatic;				   //!< true if the map is adiabatic
	double FTempMeasure;			   //!< The compressor inlet temperature used to measure the map [K].

	dMatrix FSpeed;					   //!< The compressor corrected speed for the different points of the map [rpm].
	dVector FSpeedVec;				   //!< The compressor corrected speed for the different iso-speed lines of the map [rpm]
	dVector FSpeedAdim;				   //!< The normalized compressor corrected speed for the different iso-speed lines of the map [-]

//	dVector FSpeedIndex;			   //!< Zero-based index of the speed
	double FSpeedMAX;				   //!< The maximum compressor speed [rpm]

	int FNumLines;					   //!< Number of iso-speed lines

	dMatrix FMass;					   //!< The corrected mass flow for the different points [kg / s]
	dMatrix FMassAdim;				   //!< The normalized corrected mass flow for the different points [-]
	Base_interp* FMassMAX_int;		   //!< The interpolator object for the maximum corrected mass.
	dVector FMassMAX;				   //!< The maximum corrected mass flow for each iso-speed line [kg / s]
	dVector FMassMAXAdim;			   //!< The normalized maximum corrected mass flow for each iso-speed line [-]
	double FMassMAXMAX;				   //!< The maximum corrected mass flow of the map [kg / s]

	dMatrix FPres;					   //!< The compressor ratio for the different points [-]
	dMatrix FPresAdim;				   //!< The normalized compressor ratio for the different points [-]
	Base_interp* FPresMAX_int;		   //!< The interpolator object for the maximum compressor ratio.
	dVector FPresMAX;				   //!< The maximum compressor ratio for each iso-speed line [-]
	dVector FPresMAXAdim;			   //!< The normalized maximum compressor ratio for each iso-speed line [-]
	double FPresMAXMAX;				   //!< The maximum compressor ratio of the map [-]

	dMatrix FEff;					   //!< The efficiency for the different points [-]
	dMatrix FEffAdim;				   //!< The normalized efficiency for the different points [-]
	Base_interp* FEffMAX_int;		   //!< The interpolator object for the maximum efficiency [-]
	dVector FEffMAX;				   //!< The maximum efficiency for each iso-speed line [-]
	dVector FEffMAXAdim;			   //!< The normalized maximum efficiency for each iso-speed line [-]
	double FEffMAXMAX;				   //!< The maximum efficiency of the map [-]

//	dMatrix FCoefCR;				   //!< The coef carriage return
//	dMatrix FCoefEff;				   //!< The coef eff

	std::vector<Base_interp*> FPre_MassCurve;   //!< The interpolator object for the normalized compressor ratio
	std::vector<Base_interp*> FEff_MassCurve;   //!< The interpolator object for the normalized efficiency

	int FCurrentIND;				   //!< The current iso-speed line index
	double FRTCAdim;				   //!< The current normalized corrected speed [-]
	double FCurrentMassMAX;			   //!< The current maximum corrected mass flow [kg / s]
	double FCurrentPresMAX;			   //!< The current maximum compressor ratio [-]
	double FCurrentEffMAX;			   //!< The current maximum efficiency [-]
	double FDeltaLow;				   //!< The normalized difference to the low iso-speed of the map.

	/*!
	 * \brief	Gets the vd eff adim.
	 *
	 * \return	The vd eff adim.
	 */

	//dVector *vdMassAdim;
	//dVector *vdPresAdim;
	//dVector *vdEffAdim;

	iVector FMaxIndex;				   //!< Index of the maximum compressor ratio for each iso-speed line

	dMatrix FPres_ex;				   //!< The extrapolated compressor ratio for the different points [-]
	dMatrix FSpeed_ex;				   //!< The extrapolated compressor corrected speed for the different points of the map [rpm].
	dMatrix FEff_ex;				   //!< The extrapolated efficiency for the different points [-]
	dMatrix FMass_ex;				   //!< The extrapolated corrected mass flow for the different points [kg / s]

	dMatrix FPres_ns;				   //!< The compression ratio for negative slope zone of the map [-]
	dMatrix FSpeed_ns;				   //!< The corrected speed for negative slope zone of the map [rpm]
	dMatrix FEff_ns;				   //!< The efficiency for negative slope zone of the map [-]
	dMatrix FMass_ns;				   //!< The corrected mass flow for negative slope of the map [-]

	dMatrix FPres_LPR;				   //!< The extrapolated compression ratio for low pressure ratio zone of the map [-]
	dMatrix FSpeed_LPR;				   //!< The extrapolated corrected speed for low pressure ratio zone of the map [rpm]
	dMatrix FEff_LPR;				   //!< The extrapolated efficiency for low pressure ratio zone of the map [-]
	dMatrix FMass_LPR;				   //!< The extrapolated corrected mass flow for low pressure ratio zone of the map [kg / s]

	dMatrix FPres_LS;				   //!< The extrapolated compression ratio for low speed zone of the map [-]
	dMatrix FSpeed_LS;				   //!< The extrapolated corrected speed for low speed zone of the map [rpm]
	dMatrix FEff_LS;				   //!< The extrapolated efficiency for low speed zone of the map [-]
	dMatrix FMass_LS;				   //!< The extrapolated corrected mass flow for low speed zone of the map [kg / s]

	dMatrix FPres_HS;				   //!< The extrapolated compression ratio for high speed zone of the map [-]
	dMatrix FSpeed_HS;				   //!< The extrapolated corrected speed for high speed zone of the map [rpm]
	dMatrix FEff_HS;				   //!< The extrapolated efficiency for high speed zone of the map [-]
	dMatrix FMass_HS;				   //!< The extrapolated corrected mass flow for high speed zone of the map [kg / s]

	dVector FSpeedEli;				   //!< Corrected speed array for fitting the elipse (iso-speed lines with more than 3 points)
	dVector FWmax;					   //!< The maximum corrected mass flow for the elipse [kg / s]
	dVector FWzs1;					   //!< The corrected mass flow for maximum compression ratio [kg / s]
	dVector FRCzs1;					   //!< The maximum compression ratio [-]
	dVector FWrate;					   //!< The ration betwen FWzs1 and FWmax [-]

	dVector FEffMax;				   //!< The maximum efficiency of the curves [-]
	dVector FMassEffMax;			   //!< The corrected mass flow for the maximum efficiency [kg / s]

	stCorrelation *FLeufven;		   //!< The leufven correlation
	stCorrelation *FSurge;			   //!< The surge correlation
	stCorrelation *FMaxEff;			   //!< The maximum efficiency correlation

	/*!
	 * \brief	Fit elipse curve for each iso-speed line.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void FitElipse();

	/*!
	 * \brief	Fit leufven model to the compressor map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void FitLeufvenModel();

	/*!
	 * \brief	Extrapolate low pressure ratio points.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void ExtrapolateLPR();

	/*!
	 * \brief	Extrapolate low speed zone of the map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	Dwheel	The compressor wheel diameter [m].
	 */

	void ExtrapolateLSpeed(double Dwheel);

	/*!
	 * \brief	Extrapolate high speed zone of the map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void ExtrapolateHSpeed();

	/*!
	 * \brief	Fit maximum efficiency curve.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void MaximumEfficiencyCurve();

	/*!
	 * \brief	Extrapolate efficiency.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void ExtrapolateEfficiency();

	/*!
	 * \brief	Buid the final extrapolated map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void BuidExtrapolatedMap();

  public:

	/*!
	 * \brief	Constructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	i	The compressor number.
	 */

	TSAEMap(int i);

	/*!
	 * \brief	Destructor.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	~TSAEMap();

	/*!
	 * \brief	Reads SAE compressor map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	fich	Pointer to FILE.
	 */

	void ReadSAECompressorMap(FILE *fich);

	/*!
	 * \brief	Reads SAE compressor map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	node_map	The xml node that includes the map information.
	 */

	void ReadSAECompressorMapXML(xml_node node_map);

	/*!
	 * \brief	Adimensionalize the map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 */

	void AdimensionalizeMap();

	/*!
	 * \brief	Gets the current compression ratio.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	Mass	The corrected mass flow [kg / s].
	 *
	 * \return	The compression ratio [-].
	 */

	double GetCurrentRC(double Mass);

	/*!
	 * \brief	Gets the maximum compression ratio.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The maximum compression ratio [-].
	 */

	double getMaxCompRatio() {
		return FCurrentPresMAX;
	}


	/*!
	 * \brief	Gets the current efficiency.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	Mass	The corrected mass flow [kg / s].
	 *
	 * \return	The efficiency [-].
	 */

	double GetCurrentEff(double Mass);

	/*!
	 * \brief	Interpolate the map for a value of corrected speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	RTC	The corrected speed [rpm].
	 */

	void InterpolateMAP(double RTC);

	/*!
	 * \brief	Reads the map information.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	fich	The pointer to the FILE.
	 */

	void LeeMapa(FILE *fich);

	/*!
	 * \brief	Reads the map information.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	node_compressor	The xml node of the compressor.
	 */

	void LeeMapaXML(xml_node node_compressor);

	/*!
	 * \brief	Gets the compression ratio by the Hermite interpolation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	mass	The corrected mass flow [kg / s].
	 *
	 * \return	Compression ratio [-].
	 */

	double EvaluaRCHermite(double mass);

	/*!
	 * \brief	Gets the efficiency by the spline interpolation.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	mass	The corrected mass flow [kg / s].
	 *
	 * \return	The efficiency [-].
	 */

	double EvaluaRendSplines(double mass);

	/*!
	 * \brief	Gets the reference temperature.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The reference temperature [K].
	 */

	double getTempRef() {
		return FTempRef;
	}


	/*!
	 * \brief	Gets reference pressure.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The reference pressure [Pa].
	 */

	double getPresionRef() {
		return FPresionRef;
	}

	/*!
	 * \brief	Gets the temperature used to measure the map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The temperature [K].
	 */

	double getTempMeasure() {
		return FTempMeasure;
	}

	/*!
	 * \brief	Gets the corrected mass flow for compression ratio equal to 1.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The corrected mass flow [kg / s].
	 */

	double getGastoRelComp1();

	/*!
	 * \brief	Gets the compression ratio for surge.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The compression ratio [-].
	 */

	double getRelCompBombeo();

	/*!
	 * \brief	Interpolate the mapa.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	rtc	The turbocharger speed [rpm].
	 * \param	T10	The compressor inlet temperature [K].
	 */

	void InterpolaMapa(double rtc, double T10);

	/*!
	 * \brief	Gets corrected mass flow for surge.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The corrected mass flow [-].
	 */

	double getGastoBombeo();

	/*!
	 * \brief	Gets the corrected speed.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \return	The corrected speed [rpm].
	 */

	double getRegimenCorregido();

	/*!
	 * \brief	Sets the reference conditions.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	pref	The reference pressure [Pa].
	 * \param	tref	The reference temperature [K].
	 */

	void PutReference(double pref, double tref);

	/*!
	 * \brief	Gets the adiabatic compressor map
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param [in,out]	HTM	The object for the heat transfer in the turbocharger.
	 * \param	TinT	   	The temperature used to measure the turbine [K].
	 */

	void CalculateAdiabaticEfficiency(TTC_HTM2 *HTM, double TinT);

	/*!
	 * \brief	Extrapolate the compressor map.
	 *
	 * \author	F.J. Arnau (farnau@mot.upv.es)
	 * \date	19/03/2016
	 *
	 * \param	Dw	The compressor wheel diameter [m].
	 */

	virtual void ExtrapolateMap(double Dw);

};

#endif
