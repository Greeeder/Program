//---------------------------------------------------------------------------

#ifndef TCCExternalConnectionVolH
#define TCCExternalConnectionVolH
//---------------------------------------------------------------------------

#include "TCondicionContorno.h"

class TCCExternalConnectionVol: public TCondicionContorno {
  private:

	double FDExt;
	double FAExt;

	double FDeltaX;

	double FCurrentTime;

	double FA_AExt;
	double FK_CExt;
	double FA_Dep;

	double FP_Boundary;
	double FT_Boundary;
	double FU_Boundary;

	double FP_BoundarySum;
	double FT_BoundarySum;
	double FU_BoundarySum;

	double *FCC; // Caracteristica conocida del tubo.
	double *FCD; // Caracteristica desconocida del tubo.

	int FNodoFin;
	int FIndiceCC;

	double FTime0;
	double FTimeSum;

  protected:

	double FUExt; //!< External speed. [m / s]
	double FTExt; //!< External temperature. [K]
	double FPExt; //!< External pressure. [Pa]

	int FID; //!< External connection ID.

	/*!
	 * \brief Default constructor.
	 */
	TCCExternalConnectionVol();

  public:
	TCCExternalConnectionVol(nmTypeBC TipoCC, int numCC, nmTipoCalculoEspecies SpeciesModel, int numeroespecies,
							 nmCalculoGamma GammaCalculation, bool ThereIsEGR);

	~TCCExternalConnectionVol();

	virtual void UpdateCurrentExternalProperties(double U0, double T0, double P0, double t);

	void AsignGeometricalData(double A);

	void ExternalCharacteristics(double Time);

	void CalculaCondicionContorno(double Time);

	void ReadBoundaryData(const char *FileWAM, fpos_t &filepos, int NumberOfPipes, TTubo **Pipe, int nDPF, TDPF **DPF);

	void ReadBoundaryDataXML(xml_node node_connect, int NumberOfPipes, TTubo **Pipe, int nDPF, TDPF **DPF);

	int GetID() {
		return FID;
	}
	;

	virtual void LoadNewData(double* p, double* T, double* u);

	double get_T() {
		return FTExt;
	}
	;

	double get_p() {
		return FPExt;
	}
	;

	double get_u() {
		return FUExt;
	}
	;

	double get_A() {
		return FAExt;
	}
	;

};

#endif
