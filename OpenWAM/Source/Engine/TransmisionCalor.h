#include <string>
#include <vector>
#include "Constantes.h"
//#include "Comun.h"

class Ensayo;
class Cilindro;

//namespace CalculoCALMEC
//{
#ifndef __TRANSMISIONCALOR_H
#define __TRANSMISIONCALOR_H

struct stInstantaneosNodal
{
	double ang;
	double h_adm;
	double h_esc;
	double h_gas;
	double hi;
	double area_cil;
	double T;
	double Qrad_cul;
	double Qrad_cil;
	double Qrad_pis;
	double Tgas_media_pipaesc;
	double Tgas_salida_pipaesc;
};

class TransmisionCalor
{
public:

	struct stNodos_Cil {
		int ID_NODO;
		std::string NOMBRE;
		int numero;
		int POS_AX = 0;
		int POS_RAD = 0;
		int POS_CIRC = 0;
		double D_MIN = 0;
		double D_MAX = 0;
		double L_SUP = 0;
		double L_INF = 0;
		double ANG_INI = 0;
		double ANG_FIN = 0;
		double dens = 0;
		double k = 0;
		double c = 0;
		int id_motor;

	};

	struct stNodos_Cul {

		int ID_NODO;
		std::string NOMBRE;
		int numero;
		double x = 0;
		double Y = 0;
		double Z = 0;
		double AREA_GAS=0;
		double AREA_ADM=0;
		double AREA_ESC = 0;
		double AREA_COOL = 0;
		double dens = 0;
		double k = 0;
		double c = 0;
		int id_motor;

	};

	struct stNodos_Pis {

		int ID_NODO;
		std::string NOMBRE;
		int numero;
		double X_II = 0;
		double Y_II = 0;
		double X_SI = 0;
		double Y_SI = 0;
		double X_SD = 0;
		double Y_SD = 0;
		double X_ID = 0;
		double Y_ID = 0;
		double ANCHURA = 0;
		double ALTURA = 0;
		double X_CG = 0;
		double Y_CG = 0;
		double A_S = 0;
		double A_I = 0;
		double A_LI = 0;
		double A_LD = 0;
		double Area = 0;
		double dens = 0;
		double k = 0;
		double c = 0;
		int id_motor;

	};

	struct stConductancias_Conduct {

		int ID_NODO1;
		int ID_NODO2;
		double k = 0;

	};

	struct stNodos_Fluido {

		// Esto formara un matriz si queremos
		// los vectores deben ser de la misma longitud
		//
		int ID_NODO;
		int numero;
		std::string nombre;
		std::string  tipo_fluido;
		std::string  tipo_nodo;
		int pos_matriz = 0;

		int id_motor;

	};

	struct stNodos {

		int ID_NODO;
		int id_motor;
		int pos_matriz;
		int numero;
		std::string  nombre;
		std::string  tipo_nodo;

	};

	struct stResultadosNodalT { // temperaturas

		int ID_NODO1;
		int Z;
		double valor = 0;

	};

	struct stResultadosNodalK { // k's

		int ID_NODO1;
		int ID_NODO2;
		int Z;
		double valor = 0;

	};

	struct stResultadosNodalF { // flujos

		int ID_NODO1;
		int ID_NODO2;
		int Z;
		double valor = 0;
	};


	struct stInicializar
	{
		double ConfMotor_BASE_D;
		double ConfMotor_BASE_S;
		double ConfMotor_BASE_LM;
		double ConfMotor_BASE_LB;
		double ConfMotor_BASE_E;
		double ConfMotor_BASE_H_CIL;
		double ConfMotor_BASE_VD;
		double ConfMotor_BASE_WR1A;
		double ConfMotor_BASE_WR1B;
		double ConfMotor_BASE_W2;
		int ConfMotor_BASE_NODOS_AXIALES;
		int ConfMotor_BASE_NODOS_RADIALES;
		int ConfMotor_BASE_NODOS_CIRCUNF;
		double ConfMotor_BASE_CONDUCTIVIDAD;
		double ConfMotor_PISTON_HLI;
		double ConfMotor_PISTON_DB;
		double ConfMotor_PISTON_AP;
		double ConfMotor_PISTON_PIS2OIL;
		double ConfMotor_PISTON_EX_RE_PIS2OIL;
		double ConfMotor_PISTON_CONDUCTIVIDAD;
		double ConfMotor_PISTON_DIAM_INT_GAL;
		double ConfMotor_PISTON_EX_PR_PIS2OIL;
		double ConfMotor_CULATA_DVA;
		double ConfMotor_CULATA_DVE;
		double ConfMotor_CULATA_A_PIPA_ADM;
		double ConfMotor_CULATA_A_PIPA_ESC;
		double ConfMotor_CULATA_CONDUCTIVIDAD;
		double ConfMotor_CULATA_H_CUL;
		double ConfMotor_CULATA_ACUL;
		double ConfMotor_CULATA_CTM;
		double ConfInstrumentacion_NPCL;
		double Constante_TRANS_CALOR_C_ECIL;
		double Constante_TRANS_CALOR_C_ECUL;
		double Constante_TRANS_CALOR_C_DIAM_INT_GAL;
		double Constante_TRANS_CALOR_C_F_VAST_ADM;
		double Constante_TRANS_CALOR_C_F_VAST_ESC;
		double Constante_TRANS_CALOR_C_DINY;
		double Constante_TRANS_CALOR_C_ANG_ASIENTO;
		double Constante_TRANS_CALOR_C_K_CONTACT_VALV;
		double Constante_TRANS_CALOR_C_A_WOSCHNI;
		double Constante_TRANS_CALOR_C_B_WOSCHNI;
		double Constante_TRANS_CALOR_C_C_WOSCHNI;
		double Constante_TRANS_CALOR_C_N;
		double Constante_TRANS_CALOR_C_WH1;
		bool OpcionesGenerales_MODELONODALAVANZADO;
		std::string ConfMotor_MOTOR_TIPO_ENCENDIDO;
		std::string OpcionesGenerales_FICHERO_COND_COND_MNA;//opcional
		std::string  OpcionesGenerales_FICHERO_AREA_CONV_MNA;//opcional
		std::string  OpcionesGenerales_FICHERO_GEOM_CIL_MNA;//opcional
		double VarMed_TAC;
		double VarMed_TRS;
		double VarMed_TA;
		double VarMed_TE;
		double VarMed_N;
		double W2; // Constante_TRANS_CALOR_C_W2 o si hay radiación ConfMotor_BASE_W2_CON_RAD;
		double DatosCalcGeneral_CM;
		std::string  TipoEnsayo; //C = combustión - A = Arrastre

	};

	struct stInicializarAlRCA
	{
		double TPIS;
		double TCUL;
		double TCIL;
		double PCA;
		double VCA;
		double TCA;
		double CW1;
		double CW2;
	};

	struct stInicializarPorAngulo
	{
		double deltat;
		double ang;
		double areacil;
		bool ciclo_cerrado;
		double VCIL;
	};


	std::vector<stInstantaneosNodal> var_inst_nodal;



	void CalcularTemperaturas(short Z, double *TCUL, double *TCULMAT, double *TCULVALV, double *TPIS, double *TCIL, double *TCULMAX, double *TCULMIN, double *TPISMAX, double *TPISMIN, double *TCILMAX, double *TCILMIN, double *TGASM);

	void calcular_QPIPA_ADM(short Z, double *QPIPA_ADM);
	void calcular_QPIPA_ESC(short Z, double *QPIPA_ESC);
	void calcular_QREF(short Z, double *QREF);
	void calcular_QAC(short Z, double *QAC);

	void Inicializar(stInicializar Inicializar);
	void InicializarAlRCA(stInicializarAlRCA InicializarAlRCA);
	void InicializarPorAngulo(stInicializarPorAngulo InicializarPorAngulo);
	//void Calcula_Calor_Woschni(double ang, double p, double T, double areapis, double areacil, double areacul, double TPIS, double TCUL, double TCIL, double CW1,
	//	double CW2, double C_W2, double COMB, double inc, double *QW, double *h, double *cu, double *ccul, double *cpis, double *ccil,
	//	bool ciclo_cerrado, double c1, std::string MOTOR_TIPO_ENCENDIDO, double CM, double BASE_D, double PISTON_DB, double N, double CTM, double TRANS_CALOR_C_A_WOSCHNI,
	//	double TRANS_CALOR_C_B_WOSCHNI, double TRANS_CALOR_C_C_WOSCHNI, double BASE_LM, double BASE_LB, double BASE_E, double BASE_S, double TRANS_CALOR_C_N, double TRANS_CALOR_C_WH1);
	
	void Calcula_Calor_Woschni(double T, double p);
	void Calcular_Temp_Media_Nodos(short Z, std::string nombre_nodo, double *temp_nodo);

	double getQW() const { return _QW; }
	double geth() const { return _h; }
	double getCU() const { return _CU; }
	double getQCUL() const { return _QCUL; }
	double getQPIS() const { return _QPIS; }
	double getQCIL() const { return _QCIL; }
	double getQRAD() const { return _QRAD; }
	void setQRAD(double value) { _QRAD= value; }
	double getQPIPA_ADM() const { return _QPIPA_ADM; }
	double getQPIPA_ESC() const { return _QPIPA_ESC; }

private:
	Ensayo *_ensayo;
	//Comun *_comun;
	Cilindro *_cil;
	stInicializar _Inicializar;
	stInicializarAlRCA _InicializarAlRCA;
	stInicializarPorAngulo _InicializarPorAngulo;
	//Ensayo::stConfiguracionMotor *_ConfMotor;
	//Ensayo::stVariablesMedias *_VarMed;
	short _pos_matriz;

	//resultados
	double _QW;
	double _h; //coeficiente de película
	double _CU; 
	double _QCUL;
	double _QPIS;
	double _QCIL;
	double _QRAD;
	double _QPIPA_ADM;
	double _QPIPA_ESC;

	std::vector<stNodos_Cil> _NODOS_CIL;
	std::vector<stNodos_Pis> _NODOS_PIS;
	std::vector<stNodos_Cul> _NODOS_CUL;
	std::vector<stNodos_Fluido> _NODOS_FLUIDOS;
	std::vector<stConductancias_Conduct> _CONDUCTANCIAS_CONDUCT;
	std::vector<stResultadosNodalK> _RESULTADOS_K;
	std::vector<stResultadosNodalF> _RESULTADOS_F;
	std::vector<stResultadosNodalT> _RESULTADOS_T;


	void crear_nodos_fluidos();
	void crear_nodos_base();
	void crear_nodos_piston();
	void crear_nodos_culata();
	void calcular_conductancias_conductivas();
	double resuelve_suma_geometrica(double suma, short n);
	double conductancia_radial(double k, double l, double r1, double r2, double FI);
	double distancia_nodos_CH(short nodo1, short nodo2);
	void Promedia_h(double *hmean_adm, double *hmean_exh, double *hmean_gas);
	void tgas(double *hmean_gas, double *Tmean_gas, double *TGASM);
	double rho_oil(double temp_oil);
	double mu_oil(double temp_oil);
	double cp_oil(double temp_oil);
	double k_oil(double temp_oil);
	void calculo_vcil(double Lcil, std::vector< std::vector<double> > * area_inst_nodos_cil);
	void calculo_vcil_avanzado(double Lcil, std::vector< std::vector<double> > * area_inst_nodos_cil, std::vector< std::vector<std::string> > * matriz_geom_cil, short naxiales, short nradiales, short ncircunf);
	double ha_mean(std::vector<double>* A);
	double mult_vector(std::vector<double>* m1, std::vector<double>* m2);
	double T_mean(std::vector<double>* A);
	double Qrad_mean(std::vector<double>* A_nodo_inst);
	double hoil_splash(double n, double toil_K);
	void diagonal_suma(std::vector< std::vector<double> > *matriz);
	void invNxN(std::vector< std::vector<double> > *x1, std::vector< std::vector<double> > *matriz_salida);
	void complete_triangular_matrix(std::vector< std::vector<double> >* A, std::vector< std::vector<double> >* matriz_salida);
	void pivot_triangular_matrix(std::vector< std::vector<double> >* A, std::vector<double>* matriz_salida);
	void mult_matriz(std::vector< std::vector<double> >* m1, std::vector<double>* m2, std::vector<double>* matriz_salida);
	void calcula_flujos(std::vector< std::vector<double> >* matriz_conductances, std::vector<double>* temperatura_nodos, std::vector< std::vector<double> >* mat_flux);
	void guardar_resultado_F_K(short num_cil, std::string tipo_calculo, std::vector< std::vector<double> >* matriz);
	void guardar_resultado_T(short num_cil, std::vector<double>* matriz);
	double calcular_Tg_media();
	double calcular_hA_media();
	void GenerarMatrices(std::vector< std::vector<double> >  *k_cil,
		std::vector< std::vector<double> >  *k_pis,
		std::vector< std::vector<double> >  *k_cul,
		std::vector< std::vector<double> >  *k_cil_pis,
		std::vector< std::vector<double> >  *k_cil_cul,
		std::vector< std::vector<double> >  *k_pis_cul,
		std::vector< std::vector<double> >  *k_pis_cil,
		std::vector< std::vector<double> >  *k_cul_cil,
		std::vector< std::vector<double> >  *k_cul_pis,
		short *npis, short *ncul, short *naxiales, short *nradiales, short *ncircunf);
	void GenerarMatricesDeFicheros(std::vector< std::vector<double> >  *k_cil,
		std::vector< std::vector<double> >  *k_pis,
		std::vector< std::vector<double> >  *k_cul,
		std::vector< std::vector<double> >  *k_cil_pis,
		std::vector< std::vector<double> >  *k_cil_cul,
		std::vector< std::vector<double> >  *k_pis_cul,
		std::vector< std::vector<double> >  *k_pis_cil,
		std::vector< std::vector<double> >  *k_cul_cil,
		std::vector< std::vector<double> >  *k_cul_pis,
		std::vector<std::vector<std::string> > *matriz_areas_cil,
		std::vector<std::vector<std::string> > *matriz_areas_pis,
		std::vector<std::vector<std::string> > *matriz_areas_cul,
		std::vector<std::vector<std::string> > *matriz_geom_cil,
		short *npis, short *ncul, short *naxiales, short *nradiales, short *ncircunf);
	void concatena_4_matrices(std::vector< std::vector<double> >* SI, std::vector< std::vector<double> >* SD, std::vector< std::vector<double> >* II, std::vector< std::vector<double> >* Id, std::vector< std::vector<double> >*matriz);
	void concatena_matrices_up_down(std::vector< std::vector<double> >* m1, std::vector< std::vector<double> >* m2, std::vector< std::vector<double> >*matriz);
	void concatena_matrices_left_right(std::vector< std::vector<double> >* m1, std::vector< std::vector<double> >* m2, std::vector< std::vector<double> >*matriz);
	int ObtenerPosicionNodoCilindro(std::string nodo, int naxiales, int nradiales, int ncircunf);
	void calcular_TCUL(short Z, double *TCUL, double *TCULMAX, double *TCULMIN);
	void calcular_TCIL(short Z, double *TCIL, double *TCILMAX, double *TCILMIN);
	void calcular_TCULMAT(short Z, double *TCULMAT);
	void calcular_TCULVALV(short Z, double *TCULVALV);
	void calcular_TPIS(short Z, double *TPIS, double *TPISMAX, double *TPISMIN);
	void ImprimirMatriz(std::vector< std::vector<double> >* m1, std::string nombre);
	void ImprimirVector(std::vector<double> * m1, std::string nombre);

};
#endif
//}