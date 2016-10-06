/*--------------------------------------------------------------------------------*\
|==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 | License
 |
 |	This file is part of OpenWAM.
 |
 |	OpenWAM is free software: you can redistribute it and/or modify
 |	it under the terms of the GNU General Public License as published by
 |	the Free Software Foundation, either version 3 of the License, or
 |	(at your option) any later version.
 |
 |	OpenWAM is distributed in the hope that it will be useful,
 |	but WITHOUT ANY WARRANTY; without even the implied warranty of
 |	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 |	GNU General Public License for more details.
 |
 |	You should have received a copy of the GNU General Public License
 |	along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 |
 \*--------------------------------------------------------------------------------*/

// ---------------------------------------------------------------------------
#ifndef TOpenWAMH
#define TOpenWAMH

#ifdef __BORLANDC__
#include <vcl.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <regex>
#include <fstream>
#include <sstream>
#include <memory>
#pragma hdrstop
#include "Globales.h"

#include "TTimeControl.h"
#include "pugixml.hpp"

// ENGINE BLOCK AND CYLINDERS
#include "TBloqueMotor.h"
#include "TEngineBlock.h"
#include "TCilindro4T.h"

// COMPRESSOR
#include "TCompresorDep.h"
#include "TCompTubDep.h"
#include "TCompTubos.h"

// EXTERNAL CALCULATIONS
#include "TCalculoExtern.h"
#include "TRemansoMatlab.h"
#include "TCoefDescarga.h"
#include "TControlFuel.h"

// VALVES
#include "TCDExterno.h"
#include "TEstatorTurbina.h"
#include "TRotorTurbina.h"
#include "TWasteGate.h"
#include "TValvulaContr.h"
#include "TDiscoRotativo.h"
#include "TLumbrera.h"
#include "TCDFijo.h"
#include "TValvula4T.h"
#include "TLamina.h"
#include "TMariposa.h"

// PIPES
#include "TTubo.h"
#include "TPipe.hpp"

// PLENUMS
#include "TDepVolVariable.h"
#include "TDepVolCte.h"
#include "TTurbinaSimple.h"
#include "TTurbinaTwin.h"
#include "TVenturi.h"
#include "TUnionDireccional.h"
#include "TBasicPlenum.h"

// BOUNDARY CONDITIONS
#include "TCCDescargaExtremoAbierto.h"
#include "TCCExtremoAnecoico.h"
#include "TCCExtremoCerrado.h"
#include "TCCPulso.h"
#include "TCCCilindro.h"
#include "TCCUnionEntreTubos.h"
#include "TCCPerdidadePresion.h"
#include "TCCDeposito.h"
#include "TCCRamificacion.h"
#include "TCCExtremoInyeccion.h"
#include "TCCEntradaCompresor.h"
#include "TCCUnionEntreDepositos.h"
#include "TCCCompresorVolumetrico.h"
#include "TCCCompresor.h"
#include "TCCPreVble.h"
#include "TCFDConnection.h"
#include "TCCExternalConnection.h"
#include "TCCExternalConnectionVol.h"
#include "AllBoundaryConditions.hpp"

// FINITE-VOLUME TURBINES
#include "TQ2DTurbine.h"

// TURBOCHARGER AXIS
#include "TEjeTurbogrupo.h"

// DIESEL PARTICULATE FILTER
#ifdef ParticulateFilter
#include "TDPF.h"
#include "TCanalDPF.h"
#endif

// CONCENTRIC 1D ELEMENTS
#ifdef ConcentricElement
#include "TConcentricoTubos.h"
#include "TConcentricoDPF.h"
#endif

// CONTROL DEVICES
#include "TSensor.h"
#include "TPIDController.h"
#include "TTable1D.h"
#include "TDecisor.h"
#include "TGain.h"
#include "TExternalController.h"
#include "TControlUnit.h"

#include "TIntegrable.hpp"

// OUTPUT RESULTS
#include "TOutputResults.h"
#define completo 1

/* ! \def gestorcom
 Allow the communication with WAMer
 */

#include <sys/timeb.h>

#ifdef __BORLANDC__
#define gestorcom true
#define graphicalout true
#else
//#define gestorcom 0
//#define graphicalout 0
#endif

#include "Math_wam.h"

#ifdef gestorcom
#include "TCGestorWAM.h"
#endif

#define xml_input true

class TOpenWAM {
  private:

#ifdef gestorcom

	TCGestorWAM *GestorWAM;
#endif

	std::string tzstr;
	struct timeb begining, final, current;

	stRun Run;

	stDatosTGV *DatosTGV;
	std::string fileinput;

	FILE *FileInput;
	// !< Pointers to input and output files.
	FILE *fc; // !< Pointers to input and output files.

	// XML Input
	xml_document FileInputXML;

	//char fileinput[8];

	TBloqueMotor** Engine;
	EngineBlock_ptr NewEngine;
	TCompresor** Compressor;
	TCalculoExtern* EXTERN;
	TEjeTurbogrupo** Axis;

	// ! ARRAY OF TYPES OF VALVES
	TTipoValvula** TypeOfValve;

	// ! POINTERS ARRAY TO VALVES TYPE TURBINE STATOR
	TEstatorTurbina*** StatorTurbine;
	// ! POINTERS ARRAY TO VALVES TYPE TURBINE ROTOR
	TRotorTurbina** RotorTurbine;
	// ! POINTERS ARRAY TO EXTERNAL CONNECTIONS
	TTipoValvula** CCCalcExtern;
	TTipoValvula** BCButerflyValve;

	// ! ARRAY OF PIPES
	TTubo** Pipe;

	std::vector<Pipe_ptr> NewPipes;

	// ! ARRAY OF CONCENTRIC ELEMENTS
#ifdef ConcentricElement
	TConcentrico** Concentric;
#endif

	// ! ARRAY OF DPFs
#ifdef ParticulateFilter
	TDPF** DPF;
#endif

	// ! ARRAYS OF PLENUMS
	TDeposito** Plenum;
	std::vector<TTurbina*> Turbine;
	TVenturi** Venturi;

	std::vector<BasicPlenum_ptr> NewPlenum;

	vector<Q2DTurbine_ptr> Q2DTurbine; ///< Vector of TQ2DTurbine objects.

	// ! ARRAYS OF BOUNDARY CONDITIONS
	TCondicionContorno** BC;
	TCondicionContorno** BCIntakeValve;
	TCondicionContorno** BCExhaustValve;
	TCondicionContorno** BCReedValve;
	TCondicionContorno** BCWasteGate;

	TCCExternalConnection** BCExtConnection;
	std::vector<TCCExternalConnectionVol*> BCExtConnectionVol;

	TCCCompresorVolumetrico** VolumetricCompressor;
	TCCDescargaExtremoAbierto** MatlabDischarge;
	TCCExtremoInyeccion** InjectionEnd;
	TCCPerdidadePresion **PerdidaPresion;

	std::vector<Integrable_ptr> IntegrableObjects; //!< Vector of integrable objects.

	// !OUTPUT OBJECT
	TOutputResults* Output;
	TOutputResults* NewOutput;

	// ! CONTROL PARAMETERS
	bool FirstIteration;
	int JStepMax;
	int JStepMaxDPF;
	int JCurrent;
	int JCurrentDPF;
	double TimeEndStep;
	double DeltaTPlenums;
	bool Independent;
	bool Is_EndStep;
	bool PipeStepMax;
	bool DPFStepMax;
	bool TimeMinPipe;
	bool TimeMinDPF;

	double CrankAngle;
	double AcumulatedTime;
	double Theta;
	double Theta0;

	// ! SPECIES MODEL PARAMETERS

	stEspecies* SpeciesName;
	int SpeciesNumber;

	nmTipoCalculoEspecies SpeciesModel;

	double* AtmosphericComposition;

	nmTipoCombustible FuelType;
	nmCalculoGamma GammaCalculation;

	std::map<string, TSolid*> MaterialsDB;

	TComponentArray_obj WorkingFluid;
	TComponentArray_ptr WorkingFluid_ptr;
	Eigen::ArrayXd Y_Ambient;

	// ! GENERAL PARAMETERS
	nmTipoMotor EngineType;

	nmTipoModelado SimulationType;
	bool ThereIsEGR;
	bool ThereIsFuel;
	int OpenWAMVersion;
	int Steps;
	int Increment;
	float Percentage;
	double ThetaIni;
	double ene;
	double agincr;
	double thmax;
	double grmax;
	double SimulationDuration;
	int CyclesWithoutThemalInertia;
	double AmbientPressure;
	double AmbientTemperature;
	bool ConvergenceFirstTime;

	// ! DOES THE ENGINE BLOCK EXIST?
	bool EngineBlock;
	bool NewEngineBlock;

	// ! NUMBER OF PIPES
	int NumberOfPipes;

	// ! NUMBER OF CONCENTRIC ELEMENTS
	int NumberOfConcentrics;

	// ! NUMBER OF DIESEL PARTICULATE FILTERS
	int NumberOfDPF;

	// ! VALVES PARAMETERS
	int NumberOfValves;
	int NumberOfReedValves;
	int NumberOfWasteGates;
	int NumberOfExternalCalculatedValves;

	// ! CONNECTIONS PARAMETERS
	int NumberOfConnections;
	int NumberOfVolumetricCompressors;
	int NumberOfExhaustValves;
	int NumberOfIntakeValves;
	int NumberOfCompressorsConnections;
	int NumberOfInjectionEnds;
	int NumberOfConectionsBetweenPlenums;
	int NumberOfButerflyValves;

	// ! NUMBER OF PLENUMS
	int NumberOfPlenums;

	// ! NUMBER OF VENTURIS
	int NumberOfVenturis;

	// ! NUMBER OF DIRECTIONAL JUNCIONS
	int NumberOfDirectionalJunctions;

	// ! PARAMETER FOR THE CONTROL UNIT
	int NumberOfSensors;

	TSensor **Sensor;

	int NumberOfControllers;

	TController **Controller;
	std::vector<TExternalController *> ExtCtrl;

	ControlUnit_ptr ControlUnit;

	// ! EXTERNAL CALCULATION PARAMETERS
	bool ThereIsDLL;
	int controlvalv;
	int nematlab;

	// ! TURBINE PARAMETERS
	int NumberOfTurbines;
	int CountVGT;
	int NumberOfFVTurbines; //!< Number of FV turbines.

	// ! NUMBER OF TURBOCHARGER AXIS
	int NumberOfAxis;

	// ! NUMBER OF COMPRESSORS
	int NumberOfCompressors;

	// ! NUMBER OF PRESSURE LOSSES
	int NumTCCPerdidaPresion;

	int fi_num_threads; ///< Available threads for CalculateFlowIndependent.

	/**
	 * @brief Assigns the number of threads for CalculateFlowIndependent.
	 *
	 * As CalculateFlowFlowIndependent can use up to 3 threads, it counts
	 * the number of available processors and sets fi_num_threads to 1, 2
	 * or 3 accordingly.  Also, if OMP_NUM_THREADS is set to 2 or 1, it
	 * observes it.
	 */
	void InitFlowIndependentNumThreads();

	void CleanLabelsX();

	void CleanLabels();

	void ReadGeneralData();

	void ReadGeneralDataXML();

	void ReadWorkingFluidXML();

	void ReadEngine();

	void ReadEngineXML();

	void ReadNewEngineXML();

	void ReadPipes();

	void ReadPipesXML();

	/*!
	 * \brief Reads and creates TPipe objects using XML data.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 */
	void ReadNewPipesXML();

	void ReadDPF();

	void ReadConcentric();

	void ReadValves();

	void ReadValvesXML();

	void ReadPlenums();

	void ReadPlenumsXML();

	void ReadNewPlenumsXML();

	void ReadCompressors();

	void ReadCompressorsXML();

	void ReadConnections();

	void ReadConnectionsXML();

	/*!
	 * \brief Read turbine data from the XML file.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/15
	 */
	void ReadTurbinesXML();

	/*!
	 * \brief Reads and creates TBoundaryCondition objects using XML data.
	 * 
	 * \author L.M. García-Cuevas González <luiga12@mot.upv.es>
	 * \date 2016/04/07
	 */
	void ReadNewConnectionsXML();

	void ReadTurbochargerAxis();

	void ReadTurbochargerAxisXML();

	/*!
	* \brief	Reads and creates Heat Exchanger objects using XML data.
	*
	* \author	Josep Salvador <josalib@mot.upv.es>
	* \date		2016/07/26
	*/
	void ReadHeatExchangerXML();

	void ReadSensors();

	void ReadSensorsXML();

	void ReadControllers();

	void ReadControllersXML();

	void ReadOutput(char* FileName);

	void ReadOutputXML(char* FileName);

	void ReadDataDLL();

	void BuildMaterialsDB(xml_node node_mat, std::map<string, TSolid*> &MaterialsDB);

	void RunningControl();

	void InitializeRunningAngles();

	void AllocateVGTData();

	void CalculateNewHeatPositions();

	void CalculateDistance(int NodoOrigen, int NodoFin, double Longitud, int NumberOfPlenums, int NumberOfPipes,
						   int NumberOfConnections, TTubo **Pipe, TCondicionContorno **BC);

	int SelectPipe(TTubo **Pipe, int NumberOfPipes, int nodo1, int nodo2);

	void MethodStability();

	void SearchMinimumTimeStep();

	void StudyInflowOutflowMass();

	void SearchMinimumTime(int LNumDepInicial, double* LTMinimo, TDeposito **LPlenum);

	void SearchMinimumTimeGroup(double *LTMinimo, int LNumDeposito, TDeposito **LPlenum);

	void FixTimeStep();

	void FixTimeStepExternal(double deltat);

	void RecalculateStability();

	void SolveAdjacentElements(int PipeEnd, double TiempoActual);

	void SolveBranch(int NumDeposito, double TiempoActual);

	void UpdateEngine();

	void SolveRoadLoadModel();

	void RecalculateStabilitySolver();

	void UpdateTurbocharger();

	void comunica_wam_dll();

	void ModificacionControlEjecucion();

	void Actuadores();

  public:

	TOpenWAM();

	~TOpenWAM();

	void ReadInputData(char* FileName);

	void ReadInputDataXML(char* FileName);

	void InitializeParameters();

	void ConnectFlowElements();

	void ConnectControlElements();

	void InitialHeatTransferParameters();

	void DetermineTimeStepIndependent();

	void DetermineTimeStepCommon();

	void DetermineTimeStep(double t);

	void InitializeOutput();

	void CalculateFlowIndependent();

	void CalculateFlowCommon();

	void ManageOutput();

	bool CalculationEnd();

	void Progress();

	void ProgressBegin();

	void ProgressEnd();

	void NewEngineCycle();

	void GeneralOutput();

	bool IsIndependent() {
		return Independent;
	}
	;

	void UpdateExternalBoundary(int i, double U0, double U1, double T0, double T1, double P0, double P1, double t);

	void InitialParameterTCM();

	void UpdateExternalBoundary(int i, double U0, double T0, double P0, double t);

	void UpdateControllers(int i, double val, double t);

	void UpdateHTMData(int TC, int ID, double val);

	void InitiateExternalBoundary(int i, double D0, double D1, double dX);

	void InitiateExternalBoundary(int i, double A);

	/*!
	 * \brief Gets the boundary pressure, temperature and speed.
	 * 
	 * u is negative it the speed goes out of the OpenWAM domain.
	 * If this BC is attached to the left end of a pipe, the speed is negative if it
	 * goes to the left. If this BC is attached to the right end of a pipe, the speed
	 * is negative if it goes to the right.
	 * 
	 * \param p Pressure at the BC. [bar]
	 * \param T Temperature at the BC. [K]
	 * \param u Speed at the BC. [m / s]
	 */
	void LoadNewData(int i, double* p, double* T, double* u);

	bool GetIs_EndStep();

	double Get_Output(int ID);

};
// ---------------------------------------------------------------------------
#endif
