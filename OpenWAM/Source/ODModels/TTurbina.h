/* --------------------------------------------------------------------------------*\
|==========================|
 |\\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
 | \\ |  X  | //  W ave     |
 |  \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
 |   \\/   \//    M odel    |
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


 \*-------------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------
#ifndef TTurbinaH
#define TTurbinaH
#include <iostream>
#ifdef __BORLANDC__
#include <vcl.h>
#endif
#include "TDepVolCteBase.h"
#include "TTurbineMap.h"
#include "TController.h"
//#include "TTC_HTM.h"
#include "TTCHTM2.h"
#include "TAcousticTurbine.h"

class TTurbina: public TDepVolCteBase {

  protected:
	// ---------------------------------------------------------------------------
	// VARIABLES PRIVADAS
	// ---------------------------------------------------------------------------

	TCondicionContorno **FCCEntrada;
	TCondicionContorno **FCCSalida;

//	TTC_HTM *FHTM;
	TTC_HTM2 *FHTM;
	bool FThereIsHTM;

	TTurbineMap *FMapa;

	char* FMapCoefFile;

	double FRack;
	double FMaxEfficiency;
	int FNumControlObject;
	TController *FRackController;
	bool FRackIsControlled;

	nmTipoGas FTipoGas;
	int Fncilin;

	// int FNumeroTurbina;
	int FNumeroEntradas;
	int FNumeroTurbinaTGV;
	int FNodoSalida;
	int FTuboSalida;
	int FExtremoSalida;
	int FSentidoSalida;
	int FTipoAjusteRdTurbina;
	bool asgNumeroTurbina;
	nmRendimientoTurbina FCalRendTurbina;
	double FPresionSalida;
	double FVelocidadSalida;
	double FGastoSalida;
	double FTempSalida;
	double FAsonidoSalida;
	double FRhoSalida;

	double FDiametroRodete;
	double FDiametroTurbinaIn;
	double FDiametroRodeteOut;
	double FDiametroTuerca;
	double FCriticalAngle;
	double FBladeHeight;
	double FLTE;
	double FR2geom;
	int FZ0;
	double VGTtoAlpha1;
	double VGTtoAlpha2;
	double FBeta;
	bool FVaneless;

	double FTrabajoIsen;
	double FTrabajoReal;
	double FTrabajoRealPaso;
	double FTrabajoTotal;
	double FTrabajoFluido;
	double FRendInstantaneo;
	double *FRelacionCinAcum;
	double FAjustRendTurb;
	double FDeltaPaso;
	double FTrabajoPaso;
	double FRegimen;
	stResMediosTurbina FResMediosTurbina;
	stResInstantTurbina FResInstantTurbina;
	int *FNodoEntrada;
	int *FTuboEntrada;
	int *FExtremoEntrada;
	int *FSentidoEntrada;
	double *FPresionEntrada;
	double *FVelocidadEntrada;
	double *FGastoEntrada;
	double *FTempEntrada;
	double *FAsonidoEntrada;
	double *FRhoEntrada;
	double *FTemp0Entrada;
	double *FEntalpia0Entrada;
	double *FPresion0Entrada;
	double *FEntalpiaIsenSalida;
	double *FRelacionCinematica;
	double *FRendTurbina;
	double *FRDTAjuste;
	double FRcoptima;
	double FRcmaxima;
	double FRendmaximo;

	double FPotenciaPaso;
	double FPotencia;
	double FTrabajoIsenInstTotal;
	double FRelacionCinGlobalAcum;
	double *FPonderacionRelacionCinematica;

	double *FGastoCorregido;
	double *FRegimenCorregido;
	double *FRelacionExpansion;

	nmTurbineType FTipoTurbina;
	stDatosTGV *FDatosTGV;

	double FTimeTurbina;
	TAcousticTurbine* FAcTurb;
	bool FIsAcoustic;

	// ---------------------------------------------------------------------------
	// FUNCIONES PRIVADAS
	// ---------------------------------------------------------------------------

	void TransformaContorno(double *L, double *B, double *E, double *a, double *v, double *p, int modo, double Gamma);

	// void PutNumeroTurbina(int NumeroTurbina);

	// void PutRegimen(double valor);

	double CpTurbineSimple(double Temperature, double YBurnt);

	double CpTurbineComplete(double YO2, double YCO2, double YH2O, double Temperature);

	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/17
	 */
	TTurbina();

  public:
	// ---------------------------------------------------------------------------
	// VARIABLES PUBLICAS
	// ---------------------------------------------------------------------------

	void PutNumeroTurbina(int NumeroTurbina) {
		if(!asgNumeroTurbina) {
			FNumeroTurbina = NumeroTurbina;
			asgNumeroTurbina = true;
		} else {
			std::cout << "ERROR: Esta turbina ya tiene numero asignada" << std::endl;
			throw Exception("");
		}
	}

	TTurbineMap* getMap() {
		return FMapa;
	}

	int FNumeroTurbina;

	int getNumeroTurbina() {
		return FNumeroTurbina;
	}

	int getNumeroEntradas() {
		return FNumeroEntradas;
	}

	void PutRegimen(double valor) {
		FRegimen = valor;
	}

	double GetMaxEff() {
		return FMaxEfficiency;
	}

	void PutRack(double valor) {
		FRack = valor;
	}

	TCondicionContorno* GetCCEntrada(int i);

	TCondicionContorno* GetCCSalida(int i);

	double getPotenciaPaso() {
		return FPotenciaPaso;
	}

	double GetRelacionCinematica(int i);

	virtual double GetEfficiency() = 0;

	double getDiametroRodete() const {
		return FDiametroRodete;
	}

	double getDiametroRodeteOut() const {
		return FDiametroRodeteOut;
	}

	double getDiametroTuerca() const {
		return FDiametroTuerca;
	}
	;

	double getBladeHeight() const {
		return FBladeHeight;
	}

	double getDiametroTurbinaIn() const {
		return FDiametroTurbinaIn;
	}

	double getCriticalAngle() const {
		return FCriticalAngle;
	}

	double getFr2geom() const {
		return FR2geom;
	}

	int getFz0() const {
		return FZ0;
	}

	double getVgTtoAlpha1() const {
		return VGTtoAlpha1;
	}

	double getVgTtoAlpha2() const {
		return VGTtoAlpha2;
	}

	double getBeta() const {
		return FBeta;
	}

	double getFlte() const {
		return FLTE;
	}

	bool IsVaneless() const {
		return FVaneless;
	}

	double CalculateL2(double Position);

	// ---------------------------------------------------------------------------
	// FUNCIONES PUBLICAS
	// ---------------------------------------------------------------------------

	TTurbina(int i, nmTipoDeposito TipoDeposito, int nentradas, nmTipoCalculoEspecies SpeciesModel, int numeroespecies,
			 nmCalculoGamma GammaCalculation, bool ThereIsEGR);

	virtual ~TTurbina() = 0;

	void ActualizaPropiedades(double TimeCalculo);

	void LeeTurbina(const char *FileWAM, fpos_t &filepos);

	void LeeTurbinaXML(xml_node node_turb);

	void AsignaDatosSalida(int nodsaltur, int tubsaltur, int extremo, int sentido);

	virtual void ReadAverageResultsTurb(const char *FileWAM, fpos_t &filepos) = 0;

	virtual void ReadAverageResultsTurbXML(xml_node node_trb) = 0;

	virtual void CabeceraResultadosMedTurb(stringstream& medoutput) = 0;

	virtual void ImprimeResultadosMedTurb(stringstream& medoutput) = 0;

	virtual void IniciaMedias() = 0;

	virtual void AcumulaMedias(double Tiempo) = 0;

	virtual void LeeResultadosInstantTurb(const char *FileWAM, fpos_t &filepos) = 0;

	virtual void LeeResultadosInstantTurbXML(xml_node node_turb) = 0;

	virtual void CabeceraResultadosInstantTurb(stringstream& insoutput) = 0;

	virtual void ImprimeResultadosInstantTurb(stringstream& insoutput) = 0;

	virtual void ResultadosInstantTurb() = 0;

	virtual void CalculaResultadosMediosTurb() = 0;

	void AsignaDatosEntrada(int nodentur, int tubsaltur, int extremo, int sentido, int n);

	virtual void CalculaCondicionTurbina(double TimeCalculo) = 0;

	virtual void AsignaEntradaSalidaCC() = 0;

	void CalculoPotenciaPaso();

	virtual void ImprimeResultadosMediosPantalla() = 0;

	void AsignaRackController(TController **Controller);

	void AllocateDatosTGV(stDatosTGV *DatosTGV);

//	void AsignTCHTM(TTC_HTM *HTM) {
//		FHTM = HTM;
//		FThereIsHTM = true;
//	}
//	;

	void AsignTCHTM(TTC_HTM2 *HTM) {
		FHTM = HTM;
		FThereIsHTM = true;
	}


	void PreprocessMap(double TinC);

	double TempMeasure() {
		return FMapa->TempMeasure();
	}
	;

	void AsignAcousticElements(TTubo **Pipe);

	TAcousticTurbine* AcousticT() {
		return FAcTurb;
	}

};
// ---------------------------------------------------------------------------
#endif

