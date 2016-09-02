#include <string>
#include <vector>
#include "Comun.h"

class Ensayo;
class Cilindro;

//namespace CalculoCALMEC
//{
#ifndef __TRANSMISIONCALOR_H
#define __TRANSMISIONCALOR_H



class TransmisionCalor
{
public:

	struct stNodos_Cil {	//!< Liner nodes data (geometry, etc.)
		int ID_NODO;	//!< Node ID (number among all nodes of the model)	//	[CALMEC]
		int numero;	//!< Node number (inside liner)	REVISAR	//	[CALMEC]
		std::string NOMBRE;	//!< Node name: "CIL" + Axial position (int) + Radial position (int) + Circumferential position (int). Example: CIL123	REVISAR
		int POS_AX = 0;	//!< Axial position (integer). From top to bottom.
		int POS_RAD = 0;	//!< Radial position (integer). From inner to outer.
		int POS_CIRC = 0;	//!< Circumferential position (integer). Counter-clockwise, starting from the uncooled node.
		double D_MIN = 0;	//!< Internal diameter (minor) [m]
		double D_MAX = 0;	//!< External diameter (major) [m]
		double L_SUP = 0;	//!< height of the top border of the node, from the liner bottom [m]
		double L_INF = 0;	//!< height of the bottom border of the node, from the liner bottom [m]
		double ANG_INI = 0;	//!< Node start angle [deg]
		double ANG_FIN = 0;	//!< Node end angle [deg]
		double volume = 0;	//!< Node volume [m3]
		double k = 0;	//!< Node conductivity based on its temperature (not used) [W/m/K]	REVISAR
		int id_motor;	//!< Engine ID	REVISAR	//	[CALMEC]
		//(JS) removed unused fields "dens" and "c" (also in stNodos_Pis & stNodos_Cul)
	};

	struct stNodos_Pis {	//!< Piston nodes data (geometry, etc.)

		int ID_NODO;	//!< Node ID (number among all nodes of the model)	//	[CALMEC]
		int numero;	//!< Node number (inside piston)	REVISAR
		std::string NOMBRE;	//!< Node name	REVISAR
		//!< X axis: radial. From outer piston (0) to the central axis (dp/2).
		//!< Y axis: vertical. From bottom (0) to top.
		double X_II = 0;	//!< X coordinate of lower left corner [m]
		double Y_II = 0;	//!< Y coordinate of lower left corner [m]
		double X_SI = 0;	//!< X coordinate of upper left corner [m]
		double Y_SI = 0;	//!< Y coordinate of upper left corner [m]
		double X_SD = 0;	//!< X coordinate of upper right corner [m]
		double Y_SD = 0;	//!< Y coordinate of upper right corner [m]
		double X_ID = 0;	//!< X coordinate of lower right corner [m]
		double Y_ID = 0;	//!< Y coordinate of lower right corner [m]
		double ANCHURA = 0;	//!< Node width [m]
		double ALTURA = 0;	//!< Node height [m]
		double X_CG = 0;	//!< X coordinate of the geometric center of the node [m]
		double Y_CG = 0;	//!< Y coordinate of the geometric center of the node [m]
		double A_S = 0;	//!< Area of the node top face [m2]
		double A_I = 0;	//!< Area of the node bottom face [m2]
		double A_LI = 0;	//!< Area of the left side of the node [m2]
		double A_LD = 0;	//!< Area of the right side of the node [m2]
		double Area = 0;	//!< Node internal area on the XY section [m2]
		double volume = 0;	//!< Node volume [m3]
		double k = 0;	//!< Node conductivity based on its temperature (not used) [W/m/K]	REVISAR
		int id_motor;	//!< Engine ID	REVISAR	//	[CALMEC]

	};

	struct stNodos_Cul {	//!< Cylinder-head nodes data (geometry, etc.)

		int ID_NODO;	//!< Node ID (number among all nodes of the model)	//	[CALMEC]
		int numero;	//!< Node number (inside cylinder-head)	REVISAR
		std::string NOMBRE;	//!< Node name	REVISAR
		double x = 0;	//!< X coordinate of the geometric center of the node. X axis is transverse as intake and exhaust ports. [m]
		double Y = 0;	//!< Y coordinate of the geometric center of the node. Y axis is longitudinal as crankshaft axis. [m]
		double Z = 0;	//!< Z coordinate of the geometric center of the node. Z axis is vertical as cylinder axis. [m]
		double AREA_GAS=0;	//!< Contact area between the node and the combustion chamber gas [m2]
		double AREA_ADM = 0;	//!< Contact area between the node and the intake gas (in ports) [m2]
		double AREA_ESC = 0;	//!< Contact area between the node and the exhaust gas (in ports) [m2]
		double AREA_COOL = 0;	//!< Contact area between the node and the cooling fluid [m2]
		double volume = 0;	//!< Node volume [m3]
		double k = 0;	//!< Node conductivity based on its temperature (not used) [W/m/K]	REVISAR
		int id_motor;	//!< Engine ID	REVISAR	//	[CALMEC]

	};

	struct stConductancias_Conduct {	//!< Conductive conductances among solid nodes

		int ID_NODO1;	//!< Node ID (number among all nodes of the model)
		int ID_NODO2;	//!< Node ID (number among all nodes of the model)
		double k = 0;	//!< Conductance [W/K]

	};

	struct stNodos_Fluido {	//!< Fluid nodes data	//	[CALMEC]

		// Esto formara una matriz si queremos
		// los vectores deben ser de la misma longitud
		
		int ID_NODO;	//!< Node ID (number among all nodes of the model)
		int numero;	//!< Node number (inside fluid nodes)	REVISAR
		std::string nombre;	//!< Node name	REVISAR
		std::string  tipo_fluido;	//!< Chamber gas, intake gas, exhaust gas, water, oil, ...	REVISAR
		std::string  tipo_nodo;	//!< General convective nodes, gas in contact with liner, ...	REVISAR
		int pos_matriz = 0;	//!< Node number among all nodes of the model	REVISAR
		int id_motor;	//!< Engine ID	REVISAR

	};

	struct stNodos {	//!< Nodes data	//	[CALMEC]

		int ID_NODO;	//!< Node ID (number among all nodes of the model)
		int numero;	//!< Node number (inside fluid nodes)	REVISAR
		std::string  nombre;	//!< Node name	REVISAR
		std::string  tipo_nodo;	//!< ?????	REVISAR
		int pos_matriz;	//!< ?????	REVISAR
		int id_motor;	//!< Engine ID	REVISAR

	};

	struct stResultadosNodalT { // temperaturas	//!< Temperature of nodes (result)

		int ID_NODO1;	//!< Node ID (number among all nodes of the model)
		int Z;	//!< Cylinder number	REVISAR
		double valor = 0;	//!< Node temperature [K]

	};

	struct stResultadosNodalK { // k's	//!< Conductances among nodes

		int ID_NODO1;	//!< Node ID (number among all nodes of the model)
		int ID_NODO2;	//!< Node ID (number among all nodes of the model)
		int Z;	//!< Cylinder number	REVISAR
		double valor = 0;	//!< Conductance between both nodes [W/K]

	};

	struct stResultadosNodalF { // flujos	//!< Heat fluxes among nodes (result)

		int ID_NODO1;	//!< Node ID (number among all nodes of the model)
		int ID_NODO2;	//!< Node ID (number among all nodes of the model)
		int Z;	//!< Cylinder number	REVISAR
		double valor = 0;	//!< Heat flux between both nodes [W]
	};


	struct stInicializar	//!< Several variables and constants
	{
		double ConfMotor_BASE_D;	//!< Bore [m] 
		double ConfMotor_BASE_S;	//!< Stroke [m] 
		double ConfMotor_BASE_LM;	//!< Crankshaft radius [m]
		double ConfMotor_BASE_LB;	//!< Connecting rod length [m]
		double ConfMotor_BASE_E;	//!< Piston pin decentering [m]
		double ConfMotor_BASE_H_CIL;	// Heat transfer coefficient at liner wall [W/m2/K]
		double ConfMotor_BASE_VD;	//!< Displaced volume [m3]
		double ConfMotor_BASE_WR1A;	// constantes Woschni??
		double ConfMotor_BASE_WR1B;	// constantes Woschni??
		double ConfMotor_BASE_W2;	// constantes Woschni??
		int ConfMotor_BASE_NODOS_AXIALES;	//!< Number of nodes in axial direction
		int ConfMotor_BASE_NODOS_RADIALES;	//!< Number of nodes in radial direction
		int ConfMotor_BASE_NODOS_CIRCUNF;	//!< Number of nodes in circumferential direction
		double ConfMotor_BASE_CONDUCTIVIDAD;	//!< Conductivity of liner material [W/m/K]
		double ConfMotor_BASE_ANG_NO_REFR_CIL;	//!< Opening angle of the node column in the liner that is not in contact with coolant [deg] //(JS) ADDED!!!!
		double ConfMotor_PISTON_HLI;	//!< Clearance distance between piston top and firedeck at TDC [m]
		double ConfMotor_PISTON_DB;	//!< Piston bowl diameter [m]
		double ConfMotor_PISTON_AP;	//!< Piston area exposed to chamber gas, optional [m2]	//	[CALMEC]
		double ConfMotor_PISTON_PIS2OIL;
		double ConfMotor_PISTON_EX_RE_PIS2OIL;
		double ConfMotor_PISTON_CONDUCTIVIDAD;	//!< Conductivity of piston material [W/m/K] (not used)
		double ConfMotor_PISTON_DIAM_INT_GAL;	//!< Internal diameter of oil gallery in the piston  [m]
		double ConfMotor_PISTON_EX_PR_PIS2OIL;
		double ConfMotor_CULATA_DVA;	//!< Diameter of intake valves [m]
		double ConfMotor_CULATA_DVE;	//!< Diameter of exhaust valves [m]
		double ConfMotor_CULATA_A_PIPA_ADM;	//!< Surface area of one intake port [m2]	//	[CALMEC]
		double ConfMotor_CULATA_A_PIPA_ESC;	//!< Surface area of one exhaust port [m2]	//	[CALMEC]
		double ConfMotor_CULATA_CONDUCTIVIDAD;	//!< Conductivity of cylinder-head material [W/m/K]
		double ConfMotor_CULATA_H_CUL;
		double ConfMotor_CULATA_ACUL;
		double ConfMotor_CULATA_CTM;	//!< Swirl ratio at IVC
		double ConfInstrumentacion_NPCL;
		double Constante_TRANS_CALOR_C_ECIL;	//!< Factor that multiplies bore to obtain liner thickness
		double Constante_TRANS_CALOR_C_ECUL;	//!< Factor that multiplies bore to obtain firedeck thickness
		double Constante_TRANS_CALOR_C_DIAM_INT_GAL;	//!< Factor that multiplies bore to obtain diameter of the oil gallery in the piston (not used)
		double Constante_TRANS_CALOR_C_F_VAST_ADM;	//!< Factor that multiplies intake valve diameter to obtain stem diameter
		double Constante_TRANS_CALOR_C_F_VAST_ESC;	//!< Factor that multiplies exhaust valve diameter to obtain stem diameter
		double Constante_TRANS_CALOR_C_DINY;	//!< Factor that multiplies bore to obtain injector diameter
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
		double VarMed_TAC;	//!< Oil temperature [K]
		double VarMed_TRS;	// Coolant temperature [K]??
		double VarMed_TA;	//!< Intake gas temperature [K]
		double VarMed_TE;	//!< Exhaust gas temperature [K]
		double VarMed_N;	//!< Engine speed [rps]
		double W2; // Constante_TRANS_CALOR_C_W2 o si hay radiación ConfMotor_BASE_W2_CON_RAD;
		double DatosCalcGeneral_CM;
		std::string  TipoEnsayo; //C = combustión - A = Arrastre

	};

	struct stInicializarAlRCA	//!< Conditions at IVC
	{
		double TPIS;	//!< Mean piston temperature [K]
		double TCUL;	//!< Mean cylinder-head temperature [K]
		double TCIL;	//!< Mean liner temperature [K]
		double PCA;	//!< Chamber pressure [??]
		double VCA;	//!< Chamber volume [m3]
		double TCA;	//!< Chamber gas temperature [K]
		double CW1;	//!< Coefficient for mean piston speed in Woschni correlation
		double CW2;	//!< Coefficient for tangential speed in Woschni correlation
	};

	struct stInicializarPorAngulo
	{
		double deltat;	//!< Time step
		double ang;	//!< Crank angle degree
		double areacil;	//!< Instantaneous exposed liner area
		bool ciclo_cerrado;	//!< Flag for closed cycle
		double VCIL;	//!< Instantaneous cylinder volume
	};


	std::vector<Comun::stInstantaneosNodal> var_inst_nodal;



	void CalcularTemperaturas(short Z, double *TCUL, double *TCULMAT, double *TCULVALV, double *TPIS, double *TCIL, double *TCULMAX, double *TCULMIN, double *TPISMAX, double *TPISMIN, double *TCILMAX, double *TCILMIN, double *TGASM);

	void calcular_QPIPA_ADM(short Z, double *QPIPA_ADM);
	void calcular_QPIPA_ESC(short Z, double *QPIPA_ESC);
	void calcular_QREF(short Z, double *QREF);
	void calcular_QAC(short Z, double *QAC);

	void Inicializar(Comun *comun, stInicializar Inicializar);
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
	Comun *_comun;
	Cilindro *_cil;
	stInicializar _Inicializar;
	stInicializarAlRCA _InicializarAlRCA;
	stInicializarPorAngulo _InicializarPorAngulo;
	//Ensayo::stConfiguracionMotor *_ConfMotor;
	//Ensayo::stVariablesMedias *_VarMed;
	short _pos_matriz;

	//resultados
	double _QW;
	double _h; //coeficiente de película	//!< Heat transfer coefficient [W/m2/k]
	double _CU; 
	double _QCUL;
	double _QPIS;
	double _QCIL;
	double _QRAD;
	double _QPIPA_ADM;
	double _QPIPA_ESC;


	std::vector<stNodos_Cil> _NODOS_CIL;	//!< Liner nodes data (geometry, etc.)
	std::vector<stNodos_Pis> _NODOS_PIS;	//!< Piston nodes data (geometry, etc.)
	std::vector<stNodos_Cul> _NODOS_CUL;	//!< Cylinder-head nodes data (geometry, etc.)
	std::vector<stNodos_Fluido> _NODOS_FLUIDOS;	//!< Fluid nodes data
	std::vector<stConductancias_Conduct> _CONDUCTANCIAS_CONDUCT;	//!< Conductive conductances among solid nodes
	std::vector<stResultadosNodalT> _RESULTADOS_T;	//!< Temperature of nodes (result)
	std::vector<stResultadosNodalK> _RESULTADOS_K;	//!< Conductances among nodes
	std::vector<stResultadosNodalF> _RESULTADOS_F;	//!< Heat fluxes among nodes (result)


	/*!
	* \brief Creates cylinder-head nodes and stores identification data
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 03/06/2016
	*/
	void crear_nodos_fluidos();

	/*!
	* \brief Creates liner nodes and stores geometry and identification data
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 01/06/2016
	*/
	void crear_nodos_base();

	/*!
	* \brief Looks for the factor of geometric growth used to calculate the heights of the liner nodes
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 02/06/2016
	*
	* \param suma	Liner total lenght divided by the height of the first node
	* \param n   	Number of axial nodes
	*
	* \return Factor of geometric growth
	*/
	double resuelve_suma_geometrica(double suma, short n);

	/*!
	* \brief Creates piston nodes and stores geometry and identification data
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 02/06/2016
	*/
	void crear_nodos_piston();

	/*!
	* \brief Creates cylinder-head nodes and stores geometry and identification data
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 13/06/2016
	*/
	void crear_nodos_culata();

	void calcular_conductancias_conductivas();

	/*!
	* \brief Calculates conductance between two nodes that are in contact in radial direction
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 14/06/2016
	*
	* \param k	Thermal conductivity [W/m/K]
	* \param l  Height of nodes [m]
	* \param r1	Radius of the geometric center of the first node [m]
	* \param r2	Radius of the geometric center of the second node [m]
	* \param FI	Aperture angle of the nodes [rad]
	*
	* \return Radial conductance
	*/
	double conductancia_radial(double k, double l, double r1, double r2, double FI);

	/*!
	* \brief Calculates distance between geometric centers of two nodes of the cylinder-head
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 14/06/2016
	*
	* \param nodo1	Number of first node among the nodes of the cylinder-head
	* \param nodo2  Number of second node among the nodes of the cylinder-head
	*
	* \return Distance between geometric centers of the two nodes
	*/
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


	
	/*!
	* \brief Extracts conductive conductances from the common list and stores them in 2-dim conductance vectors of each element
	*
	* \author P. Olmeda, V. Méndez, N. Molina, J. Salvador
	* \date 27/06/2016
	*
	* \param *k_cil		By reference: matrix of conductances among the liner nodes
	* \param *k_pis		By reference: matrix of conductances among the piston nodes
	* \param *k_cul		By reference: matrix of conductances among the cyl.-head nodes
	* \param *k_cil_pis	By reference: matrix of conductances between liner and piston nodes
	* \param *k_cil_cul	By reference: matrix of conductances between liner and cyl.-head nodes
	* \param *k_pis_cul	By reference: matrix of conductances between piston and cyl.-head nodes
	* \param *k_pis_cil	By reference: matrix of conductances between piston and liner nodes
	* \param *k_cul_cil	By reference: matrix of conductances between cyl.-head and liner nodes
	* \param *k_cul_pis	By reference: matrix of conductances between cyl.-head and piston nodes
	* \param *npis		By reference: number of piston nodes
	* \param *ncul		By reference: number of cyl.-head nodes
	* \param *naxiales	By reference: number of liner nodes in axial direction
	* \param *nradiales	By reference: number of liner nodes in radial direction
	* \param *ncircunf	By reference: number of liner nodes in circumferential direction
	*
	*/
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