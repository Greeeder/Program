// ---------------------------------------------------------------------------

#ifndef TOutputResultsH
#define TOutputResultsH

#include "Globales.h"
#include <sstream>

#include "TTubo.h"
#include "TBloqueMotor.h"
#include "TDeposito.h"
#include "TEjeTurbogrupo.h"
#include "TCompresor.h"
#include "TTurbina.h"
#include "TTipoValvula.h"
#include "TCCCilindro.h"
#include "TCCDeposito.h"
#include "TCCUnionEntreDepositos.h"
#include "TCCCompresorVolumetrico.h"
#include "TVenturi.h"
#include "TSensor.h"
#include "TController.h"
#include "TCalculoExtern.h"
#include "TWasteGate.h"
#include "TLamina.h"
#include "TDPF.h"
#include "TPipe.hpp"

#include "TGraphicable.h"

// ---------------------------------------------------------------------------

enum nmTypeOfResults {
	nmLastCyle = 0, nmAllCyclesIndependent = 1, nmAllCyclesConcatenated = 2, nmEveryNCycles = 3
};

class TOutputResults {
  private:

	std::string PlotsFileName;							//!< File name to store the plots
	std::string CycleOutFileName;						//!< File name to store the cycle averaged results

	std::vector<Graphicable_ptr> ObjectForPlots;		//!< Graphicable object array for plots
	std::vector<Graphicable_ptr> ObjectForCycleOut;		//!< Graphicable object array for averaged results

	std::vector<std::string> HeaderForPlots;			//!< String array to store the plots header
	std::vector<std::string> HeaderForCycleOut;			//!< String array to store the cycle averaged results

	std::vector<std::vector<float>> Plots;				//!< Matrix data with plots results
	std::vector<std::vector<float>> CycleOutput;		//!< Matrix data with cycle averaged results

	double sampling_interval;							//!< Sampling interval to store the plots [s]
	double next_t;										//!< Next time when the plots will be stored [s]
	double previous_t;									//!< Last time when the plots were stored [s]

	nmTypeOfResults FTypeOfInsResults;

	int FCyclePeriod;

	double FInsPeriod;

	bool FPlotThisCycle;
	bool FMultipleFiles;
	bool FPlotIns;

	bool FFirstTime;

	bool InsHeaderCreated;
	bool WriteInsHeader;

	bool FWriteSpaceTime;

	double FControlAngle0;
	double FControlAngle1;

	vector<TTubo*> AvgPipe;
	vector<TCilindro*> AvgCylinder;
	TBloqueMotor* AvgEngine;
	vector<TDeposito*> AvgPlenum;
	vector<TEjeTurbogrupo*> AvgAxis;
	vector<TCompresor*> AvgCompressor;
	vector<TTurbina*> AvgTurbine;
	vector<TTipoValvula*> AvgValve;
	iVector AvgValveNode;
	vector<TCCCompresorVolumetrico*> AvgRoot;
	vector<TVenturi*> AvgVenturi;
	vector<TCCUnionEntreDepositos*> AvgConnection;
	vector<TSensor*> AvgSensor;
	vector<TController*> AvgController;
	vector<TDPF*> AvgDPF;

	stringstream FAvgOutput;
	fstream FFileAvg;
	char FAvgFilename[300];

	vector<TCilindro*> InsCylinder;
	vector<TDeposito*> InsPlenum;
	vector<TTubo*> InsPipe;
	vector<TVenturi*> InsVenturi;
	vector<TTipoValvula*> InsValve;
	iVector InsValveNode;
	vector<TEjeTurbogrupo*> InsTurbo;
	vector<TCompresor*> InsCompressor;
	vector<TTurbina*> InsTurbine;
	vector<TCCCompresorVolumetrico*> InsRoot;
	vector<TCCUnionEntreDepositos*> InsConnection;
	vector<TWasteGate*> InsWasteGate;
	vector<TLamina*> InsReedValve;
	vector<TSensor*> InsSensor;
	vector<TController*> InsController;
	vector<TDPF*> InsDPF;
	vector<Integrable_ptr> InsIntegrable;

	stringstream FInsOutput;
	fstream FFileIns;
	char FInsFilename[300];

	char FFileCountC[10];
	int FFileCountI;
	int FCharacters;

	vector<TTubo*> STPipe;
	vector<TDeposito*> STPlenum;
	vector<TCilindro*> STCylinder;

	iVector FParameterSpaceTime;

	FILE *FileOutPressure; // !< Pointers to files for space time results.
	FILE *FileOutTemp; // !< Pointers to files for space time results.
	FILE *FileOutVel; // !< Pointers to files for space time results.
	FILE *FileOutFlow; // !< Pointers to files for space time results.
	FILE *FOutYO2; // !< Pointers to files for space time results.
	FILE *FOutYN2; // !< Pointers to files for space time results.
	FILE *FOutYCO2; // !< Pointers to files for space time results.
	FILE *FOutYH2O; // !< Pointers to files for space time results.
	FILE *FOutYCO; // !< Pointers to files for space time results.
	FILE *FOutYNOx; // !< Pointers to files for space time results.
	FILE *FOutYSoot; // !< Pointers to files for space time results.
	FILE *FOutYHC; // !< Pointers to files for space time results.
	FILE *FOutYFuel; // !< Pointers to files for space time results.
	FILE *FOutYFreshAir; // !< Pointers to files for space time results.
	FILE *FOutYBurntGas; // !< Pointers to files for space time results.
	FILE *FOutFlowO2; // !< Pointers to files for space time results.
	FILE *FOutFlowN2; // !< Pointers to files for space time results.
	FILE *FOutFlowCO2; // !< Pointers to files for space time results.
	FILE *FOutFlowH2O; // !< Pointers to files for space time results.
	FILE *FOutFlowCO; // !< Pointers to files for space time results.
	FILE *FOutFlowNOx; // !< Pointers to files for space time results.
	FILE *FOutFlowSoot; // !< Pointers to files for space time results.
	FILE *FOutFlowHC; // !< Pointers to files for space time results.
	FILE *FOutFlowFuel; // !< Pointers to files for space time results.
	FILE *FOutFlowFreshAir; // !< Pointers to files for space time results.
	FILE *FOutFlowBurntGas; // !< Pointers to files for space time results.

	char salpre[300];
	char saltem[300];
	char salvel[300];
	char salair[300];
	char salYO2[300];
	char salYN2[300];
	char salYCO2[300];
	char salYH2O[300];
	char salYCO[300];
	char salYNOx[300];
	char salYSoot[300];
	char salYHC[300];
	char salYCombustible[300];
	char salYAireFresco[300];
	char salYGasQuemado[300];
	char salGastoO2[300];
	char salGastoN2[300];
	char salGastoCO2[300];
	char salGastoH2O[300];
	char salGastoCO[300];
	char salGastoNOx[300];
	char salGastoSoot[300];
	char salGastoHC[300];
	char salGastoCombustible[300];
	char salGastoAireFresco[300];
	char salGastoGasQuemado[300];

	void ConvertCharacter(int confile, char confile1[], int Characters);

  public:

	/*!
	 * \brief Default constructor.
	 *
	 * \author F. J. Arnau <farnau@mot.upv.es>
	 * \date 2016/06/11
	 */
	TOutputResults();

	/*!
	* \brief Constructor with the case name.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param case_name	String with the name of the case
	*/
	TOutputResults(string case_name);

	/*!
	* \brief Default destructor.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*/
	~TOutputResults();

	/*!
	* \brief Append an object which contains results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param obj	Graphicable object
	*/
	void AppendPlotObj(Graphicable_ptr obj) {
		ObjectForPlots.push_back(obj);
	}

	/*!
	* \brief Append an object which contains cycle averaged results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param obj	Graphicable object
	*/
	void AppendCycleOutObj(Graphicable_ptr obj) {
		ObjectForCycleOut.push_back(obj);
	}

	/*!
	* \brief Create the header for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*/
	void CreatePlotsHeader();

	/*!
	* \brief Create the header cycle averaged results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*/
	void CreateCycleOutHeader();

	/*!
	* \brief Write cycle averaged results to a file.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param separator	String used to separate the data columns
	*/
	void CycleOutToFile(std::string separator);

	/*!
	* \brief Integrate the instantaneous results to get the cycle averaged results.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param t		Current time [s]
	* \param ncycle	Number of the current cycle
	*/
	void IntegrateCycleOutput(double t, int ncycle);

	/*!
	* \brief Store the instantaneous results for plots.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param t		Current time [s]
	*/
	void StorePlots(double t);

	/*!
	* \brief Write plots to a file.
	*
	* \author F. J. Arnau <farnau@mot.upv.es>
	* \date 2016/06/11
	*
	* \param separator	String used to separate the data columns
	*/
	void PlotsToFile(std::string separator);

	double GetFControlAngle1() {
		return FControlAngle1;
	}
	;

	/*!
	 * \brief Initialises a vector of TIntegrble objects with instantaneous results.
	 *
	 * \param integrable Vector of integrable objects.
	 */ 
	void InitIntegrableInstantaneousResults(const vector<Integrable_ptr> & integrables);

	void ReadAverageResults(const char* FileWAM, fpos_t& filepos, TTubo** Pipe, bool EngineBlock, TBloqueMotor** Engine,
							TDeposito **Plenum, TEjeTurbogrupo** Axis, TCompresor** Compressor,
							TTurbina** Turbine, TCondicionContorno** BC, TDPF** DPF, TCCCompresorVolumetrico** Root, TVenturi** Venturi,
							TSensor** Sensor, TController** Controller, int TotalCycles, char* ModelName);

	void ReadAverageResultsXML(xml_node node_openwam, TTubo** Pipe, bool EngineBlock, TBloqueMotor** Engine,
							   TDeposito **Plenum, TEjeTurbogrupo** Axis, TCompresor** Compressor, TTurbina** Turbine,
							   TCondicionContorno** BC, TDPF** DPF, TCCCompresorVolumetrico** Root, TVenturi** Venturi, TSensor** Sensor,
							   TController** Controller, int TotalCycles, char* ModelName);

	void HeaderAverageResults(stEspecies *SpeciesName, TCalculoExtern* EXTERN, bool ThereIsDLL);

	void OutputAverageResults(double AcumulatedTime, TCalculoExtern* EXTERN, bool ThereIsDLL);

	void CopyAverageResultsToFile(int mode);

	void CopyInstananeousResultsToFile(int mode);

	void ReadInstantaneousResults(const char* FileWAM, fpos_t &filepos, TBloqueMotor** Engine, TDeposito** Plenum,
								  TTubo** Pipe, TVenturi** Venturi, TCondicionContorno** BC, TDPF** DPF,
								  TEjeTurbogrupo** Turbo, TCompresor** Compressor, TTurbina** Turbine, TCCCompresorVolumetrico** Root,
								  TCondicionContorno** BCWasteGate, int NumberOfWasteGates, TCondicionContorno** BCReedValve,
								  int NumberOfReedValves, TSensor** Sensor, TController** Controller, char* ModelName);

	void ReadInstantaneousResultsXML(xml_node node_openwam, TBloqueMotor** Engine, TDeposito** Plenum, TTubo** Pipe,
									 TVenturi** Venturi, TCondicionContorno** BC, TDPF** DPF, TEjeTurbogrupo** Turbo,
									 TCompresor** Compressor, TTurbina** Turbine, TCCCompresorVolumetrico** Root, TCondicionContorno** BCWasteGate,
									 int NumberOfWasteGates, TCondicionContorno** BCReedValve, int NumberOfReedValves,
									 TSensor** Sensor, TController** Controller, char* ModelName);

	void ReadSettings(xml_node node);

	void ReadSpaceTimeResults(const char* FileWAM, fpos_t &filepos, TTubo** Pipe, TBloqueMotor** Engine,
							  TDeposito **Plenum);

	void ReadSpaceTimeResultsXML(xml_node node_openwam, TTubo** Pipe, TBloqueMotor** Engine, TDeposito **Plenum);

	void DoSpaceTimeFiles(int SpeciesNumber);

	void HeaderSpaceTimeResults(double thmax, double grmax, double agincr, int SpeciesNumber);

	void PrintSpaceTimeResults(bool EngineBlock, double Theta, double SimulationDuration, TBloqueMotor **Engine,
							   int SpeciesNumber);

	void HeaderInstantaneousResults(TCalculoExtern *EXTERN, bool ThereIsDLL, bool EngineBlock, stEspecies *SpeciesName);

	void PlotThisCycle(TBloqueMotor* Engine, int TotalCycles);

	void OutputInstantaneousResults(TCalculoExtern *EXTERN, bool ThereIsDLL, bool EngineBlock, double Theta,
									TBloqueMotor* Engine, double Time);

	void PlotControl(double Theta0, double Theta, double CycleDuration);

	void WriteInstantaneous(bool EngineBlock, double Angle, double AngStep, TBloqueMotor* Engine, int TotalCycles);

	void WriteSpaceTime(bool EngineBlock, TBloqueMotor* Engine, int TotalCycles);

	void PutInsPeriod(double agincr) {
		FInsPeriod = agincr;
	}
	;

};
#endif
