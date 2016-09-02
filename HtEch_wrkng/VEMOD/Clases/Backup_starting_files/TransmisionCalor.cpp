#include <iostream>
#include <fstream>
#include "TransmisionCalor.h"
#include <fstream>
#include <sstream>
#include "Constantes.h"
#include "Math.h"

using namespace std;
//namespace CalculoCALMEC
//{


void TransmisionCalor::Inicializar(Comun *comun, stInicializar Inicializar)
{

	_comun = comun;

	_Inicializar = Inicializar;

	_QRAD = 0.;

	var_inst_nodal.resize(_Inicializar.ConfInstrumentacion_NPCL + 1);

	if (!_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		crear_nodos_fluidos();
		crear_nodos_base();
		crear_nodos_piston();
		crear_nodos_culata();
		calcular_conductancias_conductivas();
	}
}

void TransmisionCalor::InicializarAlRCA(stInicializarAlRCA InicializarAlRCA)
{
	_InicializarAlRCA = InicializarAlRCA;
}

void TransmisionCalor::InicializarPorAngulo(stInicializarPorAngulo InicializarPorAngulo)
{
	_InicializarPorAngulo = InicializarPorAngulo;
}


void TransmisionCalor::crear_nodos_base()
{
	//h0 altura cámara combustión
	//kcil conductividad del material en este caso el cilindro
	//S carrera del cilindro
	//naxiales número de nodos axiales
	//nradiales número de nodos radiales
	//ncircunf número de nodos circunferenciales
	//espesor espesor del cilindro
	//Dp diámetro pistón

	short cont1;
	short cont2;
	short cont3;
	short cont4;
	short fila;
	double A;
	double delta_fi;
	double ymedio;
	double y0;
	double FI;
	double k_cil;
	double Lapoyo;
	short ncilin;
	double Lcil;
	double ngases;
	double esp_agua;
	double espesor;
	double C_ECIL;
	short naxiales;
	short nradiales;
	short ncircunf;
	double dp;
	short numero;
	double s;
	double h0;

	std::vector<double> L_nodos;

	numero = 1;
	dp = _Inicializar.ConfMotor_BASE_D;
	s = _Inicializar.ConfMotor_BASE_S;
	h0 = _Inicializar.ConfMotor_PISTON_HLI;


	naxiales = _Inicializar.ConfMotor_BASE_NODOS_AXIALES;
	nradiales = _Inicializar.ConfMotor_BASE_NODOS_RADIALES;
	ncircunf = _Inicializar.ConfMotor_BASE_NODOS_CIRCUNF;


	espesor = _Inicializar.Constante_TRANS_CALOR_C_ECIL * dp;

	esp_agua = espesor / 2;

	ncilin = naxiales * nradiales * ncircunf + ncircunf;


	L_nodos.resize(naxiales + 1); //Aqui almaceno las longitudes de cada uno de los nodos axiales
	_NODOS_CIL.resize(ncilin + 1);

	k_cil = _Inicializar.ConfMotor_BASE_CONDUCTIVIDAD;

	//la longitud total es: la carrera + h cámara combustión + falda (~0.10 carrera)
	Lcil = 1.1 * s + h0;

	// la longitud axial de los nodos es variable: más pequeños en la parte superior, mayores en la inferior.
	// He supuesto que la longitud del nodo inmediatamente inferior es igual al superior por un factor de crecimiento (A)
	// Calculo el factor de crecimiento
	A = resuelve_suma_geometrica(Lcil / h0, naxiales);

	// Los nodos circunferenciales son todos del mismo tamaño:
	delta_fi = 360. / ncircunf;

	//Escribo la posición de los nodos (axial radial circunferencial) así como sus posiciones
	L_nodos.at(1) = h0;
	fila = 2;

	for (cont1 = 2; cont1 <= naxiales; cont1++) {
		L_nodos.at(cont1) = L_nodos.at(cont1 - 1) * A;
	}

	for (cont1 = 1; cont1 <= naxiales; cont1++) {
		ymedio = L_nodos.at(cont1) / 2;
		if (cont1 == 1) {
			y0 = ymedio;
		}
		else {
			y0 = ymedio + L_nodos.at(cont1 - 1);
		}

		for (cont2 = 1; cont2 <= nradiales; cont2++) {
			for (cont3 = 1; cont3 <= ncircunf; cont3++) {
				FI = (cont3 - 1) * delta_fi;
				_pos_matriz++;
				_NODOS_CIL.at(fila - 1).ID_NODO = _pos_matriz;
				_NODOS_CIL.at(fila - 1).NOMBRE = "CIL" + std::to_string(cont1) + std::to_string(cont2) + std::to_string(cont3);
				_NODOS_CIL.at(fila - 1).numero = numero;
				_NODOS_CIL.at(fila - 1).POS_AX = cont1;
				_NODOS_CIL.at(fila - 1).POS_RAD = cont2;
				_NODOS_CIL.at(fila - 1).POS_CIRC = cont3;
				_NODOS_CIL.at(fila - 1).D_MIN = dp + 2 * espesor / nradiales * (cont2 - 1);
				_NODOS_CIL.at(fila - 1).D_MAX = dp + 2 * espesor / nradiales * (cont2);
				if (cont1 == 1) {
					_NODOS_CIL.at(fila - 1).L_SUP = Lcil;
				}
				else {
					Lapoyo = 0;
					for (cont4 = 1; cont4 <= cont1 - 1; cont4++) {
						Lapoyo = Lapoyo + L_nodos.at(cont4);
					}
					_NODOS_CIL.at(fila - 1).L_SUP = Lcil - Lapoyo;
				}
				_NODOS_CIL.at(fila - 1).L_INF = _NODOS_CIL.at(fila - 1).L_SUP - L_nodos.at(cont1);
				_NODOS_CIL.at(fila - 1).ANG_INI = delta_fi * (cont3 - 1);
				_NODOS_CIL.at(fila - 1).ANG_FIN = delta_fi * (cont3);
				_NODOS_CIL.at(fila - 1).k = k_cil;
				numero = numero + 1;
				fila = fila + 1;
			}
		}
	}

	//Nodo del cilindro  justo encima del refrigerante y de longitud h0
	for (cont1 = (1 + naxiales * nradiales * ncircunf); cont1 <= ncilin; cont1++) {
		//cont2 = cont1 - naxiales * nradiales * ncircunf + ncircunf 'identifica al nodo del cilindro de mayor radio y posicion superior
		cont2 = cont1 - naxiales * nradiales * ncircunf + (nradiales - 1) * (ncircunf);
		_pos_matriz++;
		_NODOS_CIL.at(cont1).ID_NODO = _pos_matriz;
		_NODOS_CIL.at(cont1).NOMBRE = "CIL" + std::to_string(cont1) + std::to_string(cont2) + std::to_string(cont3);
		_NODOS_CIL.at(cont1).id_motor = _NODOS_CIL.at(cont2).id_motor;
		_NODOS_CIL.at(cont1).numero = numero;
		_NODOS_CIL.at(cont1).POS_AX = _NODOS_CIL.at(cont2).POS_AX;
		_NODOS_CIL.at(cont1).POS_RAD = _NODOS_CIL.at(cont2).POS_RAD + 1;
		_NODOS_CIL.at(cont1).POS_CIRC = _NODOS_CIL.at(cont2).POS_CIRC;
		_NODOS_CIL.at(cont1).D_MIN = _NODOS_CIL.at(cont2).D_MAX;
		_NODOS_CIL.at(cont1).D_MAX = _NODOS_CIL.at(cont2).D_MAX + 2 * esp_agua;
		_NODOS_CIL.at(cont1).L_SUP = _NODOS_CIL.at(cont2).L_SUP;
		_NODOS_CIL.at(cont1).L_INF = _NODOS_CIL.at(cont2).L_INF;
		_NODOS_CIL.at(cont1).ANG_INI = _NODOS_CIL.at(cont2).ANG_INI;
		_NODOS_CIL.at(cont1).ANG_FIN = _NODOS_CIL.at(cont2).ANG_FIN;
		_NODOS_CIL.at(cont1).k = _NODOS_CIL.at(cont2).k;
		numero = numero + 1;
	}



}


void TransmisionCalor::crear_nodos_piston()
{
	//d -> diámetro del piston y del cilindro
	short cont1;
	short cont2;
	short f;
	short c;
	short n1;
	short n2;
	double Area_cont;
	double k_pis;
	double dp;

	const double f_Lp = 1.1;
	const double f_L1 = 0.6 * f_Lp;
	const double f_L2 = 0.65 * f_Lp;
	const double f_L3 = 0.8 * f_Lp;
	const double f_L4 = 0.85 * f_Lp;
	const double f_L5 = 1 * f_Lp;
	const double f_D1 = 0.15;
	//Const f_Dgal As Double = 0.05 'de Dp

	double Lp;
	double L1p;
	double L2p;
	double L3p;
	double L4p;
	double L5p;
	double L6p;

	double D1p;
	double D2p;
	double D3p;
	double D4p;
	double Dgal;

	short cont0;
	short cont;
	double Aux;
	double l;
	double fi01;
	double xg1;
	double xg2;
	double xg3;
	double Sup1;
	double Sup2;
	double Sup3;
	double r01;
	double r02;
	short npis;
	double C_DIAM_INT_GAL;
	double DB;

	_NODOS_PIS.resize(11);

	k_pis = _Inicializar.ConfMotor_PISTON_CONDUCTIVIDAD;
	dp = _Inicializar.ConfMotor_BASE_D;

	C_DIAM_INT_GAL = _Inicializar.Constante_TRANS_CALOR_C_DIAM_INT_GAL;
	DB = _Inicializar.ConfMotor_PISTON_DB;
	//Define geometria de nodos del piston (L exialmente y D radialmente)

	if (DB <= 0)
		DB = 0.9 * dp;

	Lp = f_Lp * dp;
	L1p = f_L1 * dp;
	L2p = f_L2 * dp;
	L3p = f_L3 * dp;
	L4p = f_L4 * dp;
	L5p = f_L5 * dp;
	L6p = 0.5 * (L5p + L4p);
	D2p = (dp - DB) / 2;
	D4p = dp / 2;
	D3p = (D2p + D4p) / 2;
	D1p = f_D1 * D3p;
	Dgal = C_DIAM_INT_GAL * dp;

	npis = 10;

	//UPGRADE_WARNING: El límite inferior de la matriz _NODOS_PIS ha cambiado de 1 a 0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	// ERROR: Not supported in C#: ReDimStatement


	cont = 1;

	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = 0;
	_NODOS_PIS.at(cont).Y_II = L4p;
	_NODOS_PIS.at(cont).X_SI = 0;
	_NODOS_PIS.at(cont).Y_SI = L5p;
	_NODOS_PIS.at(cont).X_SD = D2p;
	_NODOS_PIS.at(cont).Y_SD = L5p;
	_NODOS_PIS.at(cont).X_ID = D2p;
	_NODOS_PIS.at(cont).Y_ID = L4p;


	cont = 2;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_CIL.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = 0;
	_NODOS_PIS.at(cont).Y_II = L3p;
	_NODOS_PIS.at(cont).X_SI = 0;
	_NODOS_PIS.at(cont).Y_SI = L4p;
	_NODOS_PIS.at(cont).X_SD = D2p;
	_NODOS_PIS.at(cont).Y_SD = L4p;
	_NODOS_PIS.at(cont).X_ID = D2p;
	_NODOS_PIS.at(cont).Y_ID = L3p;


	cont = 3;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = 0;
	_NODOS_PIS.at(cont).Y_II = L2p;
	_NODOS_PIS.at(cont).X_SI = 0;
	_NODOS_PIS.at(cont).Y_SI = L3p;
	_NODOS_PIS.at(cont).X_SD = D2p;
	_NODOS_PIS.at(cont).Y_SD = L3p;
	_NODOS_PIS.at(cont).X_ID = D2p;
	_NODOS_PIS.at(cont).Y_ID = L2p;


	cont = 4;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = 0;
	_NODOS_PIS.at(cont).Y_II = L1p;
	_NODOS_PIS.at(cont).X_SI = 0;
	_NODOS_PIS.at(cont).Y_SI = L2p;
	_NODOS_PIS.at(cont).X_SD = D1p;
	_NODOS_PIS.at(cont).Y_SD = L2p;
	_NODOS_PIS.at(cont).X_ID = D1p;
	_NODOS_PIS.at(cont).Y_ID = L1p;



	cont = 5;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = 0;
	_NODOS_PIS.at(cont).Y_II = 0;
	_NODOS_PIS.at(cont).X_SI = 0;
	_NODOS_PIS.at(cont).Y_SI = L1p;
	_NODOS_PIS.at(cont).X_SD = D1p;
	_NODOS_PIS.at(cont).Y_SD = L1p;
	_NODOS_PIS.at(cont).X_ID = D1p;
	_NODOS_PIS.at(cont).Y_ID = 0;



	cont = 6;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = D1p;
	_NODOS_PIS.at(cont).Y_II = L1p;
	_NODOS_PIS.at(cont).X_SI = D1p;
	_NODOS_PIS.at(cont).Y_SI = L2p;
	_NODOS_PIS.at(cont).X_SD = D2p;
	_NODOS_PIS.at(cont).Y_SD = L2p;
	_NODOS_PIS.at(cont).X_ID = D1p;
	_NODOS_PIS.at(cont).Y_ID = L1p;


	cont = 7;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = D2p;
	_NODOS_PIS.at(cont).Y_II = L3p;
	_NODOS_PIS.at(cont).X_SI = D2p;
	_NODOS_PIS.at(cont).Y_SI = L4p;
	_NODOS_PIS.at(cont).X_SD = D3p;
	_NODOS_PIS.at(cont).Y_SD = L4p;
	_NODOS_PIS.at(cont).X_ID = D3p;
	_NODOS_PIS.at(cont).Y_ID = L3p;

	cont = 8;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = D2p;
	_NODOS_PIS.at(cont).Y_II = L2p;
	_NODOS_PIS.at(cont).X_SI = D2p;
	_NODOS_PIS.at(cont).Y_SI = L3p;
	_NODOS_PIS.at(cont).X_SD = D3p;
	_NODOS_PIS.at(cont).Y_SD = L3p;
	_NODOS_PIS.at(cont).X_ID = D2p;
	_NODOS_PIS.at(cont).Y_ID = L2p;


	cont = 9;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = D3p;
	_NODOS_PIS.at(cont).Y_II = L4p;
	_NODOS_PIS.at(cont).X_SI = D3p;
	_NODOS_PIS.at(cont).Y_SI = L6p;
	_NODOS_PIS.at(cont).X_SD = D4p;
	_NODOS_PIS.at(cont).Y_SD = L6p;
	_NODOS_PIS.at(cont).X_ID = D4p;
	_NODOS_PIS.at(cont).Y_ID = L4p;


	cont = 10;
	_pos_matriz++;
	_NODOS_PIS.at(cont).ID_NODO = _pos_matriz;
	_NODOS_PIS.at(cont).NOMBRE = "PIS" + std::to_string(cont);
	_NODOS_PIS.at(cont).numero = cont;
	_NODOS_PIS.at(cont).X_II = D3p;
	_NODOS_PIS.at(cont).Y_II = L3p;
	_NODOS_PIS.at(cont).X_SI = D3p;
	_NODOS_PIS.at(cont).Y_SI = L4p;
	_NODOS_PIS.at(cont).X_SD = D4p;
	_NODOS_PIS.at(cont).Y_SD = L4p;
	_NODOS_PIS.at(cont).X_ID = D4p;
	_NODOS_PIS.at(cont).Y_ID = L3p;


	for (cont0 = 1; cont0 <= 10; cont0++) {

		_NODOS_PIS.at(cont0).ANCHURA = _NODOS_PIS.at(cont0).X_SD - _NODOS_PIS.at(cont0).X_II;
		_NODOS_PIS.at(cont0).ALTURA = _NODOS_PIS.at(cont0).Y_SI - _NODOS_PIS.at(cont0).Y_II;

		if ((cont0 == 6) | (cont0 == 8)) {
			_NODOS_PIS.at(cont0).X_CG = _NODOS_PIS.at(cont0).X_II + _NODOS_PIS.at(cont0).ANCHURA / 3;
			_NODOS_PIS.at(cont0).Y_CG = _NODOS_PIS.at(cont0).Y_SI - _NODOS_PIS.at(cont0).ALTURA / 3;
		}
		else if ((cont0 == 9) | (cont0 == 10)) {
			_NODOS_PIS.at(cont0).X_CG = _NODOS_PIS.at(cont0).X_II + _NODOS_PIS.at(cont0).ANCHURA;
			_NODOS_PIS.at(cont0).Y_CG = _NODOS_PIS.at(cont0).Y_II + _NODOS_PIS.at(cont0).ALTURA / 2;
		}
		else {
			_NODOS_PIS.at(cont0).X_CG = _NODOS_PIS.at(cont0).X_II + _NODOS_PIS.at(cont0).ANCHURA / 2;
			_NODOS_PIS.at(cont0).Y_CG = _NODOS_PIS.at(cont0).Y_II + _NODOS_PIS.at(cont0).ALTURA / 2;
		}

		Aux = 0;

		if (((cont0 == 6) | (cont0 == 8))) {
			double pot1 = pow(_NODOS_PIS.at(cont0).ANCHURA, 2);
			double pot2 = pow(_NODOS_PIS.at(cont0).ALTURA, 2);
			double pot3 = pow((pot1 + pot2), 0.5);
			Aux = pot3 - _NODOS_PIS.at(cont0).ALTURA;
		}

		_NODOS_PIS.at(cont0).A_S = Constantes::pi * (pow((dp / 2 - _NODOS_PIS.at(cont0).X_II), 2) - pow((dp / 2 - _NODOS_PIS.at(cont0).X_SD), 2));
		_NODOS_PIS.at(cont0).A_I = Constantes::pi * (pow((dp / 2 - _NODOS_PIS.at(cont0).X_II + Aux), 2) - pow((dp / 2 - _NODOS_PIS.at(cont0).X_SD), 2));
		_NODOS_PIS.at(cont0).A_LI = 2 * Constantes::pi * _NODOS_PIS.at(cont0).ALTURA * (dp / 2 - _NODOS_PIS.at(cont0).X_II);
		if (((cont0 == 6) | (cont0 == 8))) {
			_NODOS_PIS.at(cont0).A_LD = 0;
			_NODOS_PIS.at(cont0).Area = _NODOS_PIS.at(cont0).ANCHURA * _NODOS_PIS.at(cont0).ALTURA / 2;
		}
		else {
			_NODOS_PIS.at(cont0).A_LD = 2 * Constantes::pi * _NODOS_PIS.at(cont0).ALTURA * (dp / 2 - _NODOS_PIS.at(cont0).X_SD);
			_NODOS_PIS.at(cont0).Area = _NODOS_PIS.at(cont0).ANCHURA * _NODOS_PIS.at(cont0).ALTURA;
		}
		_NODOS_PIS.at(cont0).dens = 0;
		_NODOS_PIS.at(cont0).k = k_pis;
		_NODOS_PIS.at(cont0).c = 0;


	}

	//Corregimos el área de los nodos del piston expuestos al gas teniendo en cuenta el área del pistón introducida por el usuario
	double area_nodal;
	double factor_area;
	area_nodal = _NODOS_PIS.at(1).A_S + _NODOS_PIS.at(1).A_LD + _NODOS_PIS.at(7).A_S + _NODOS_PIS.at(9).A_S + _NODOS_PIS.at(9).A_LI;
	factor_area = _Inicializar.ConfMotor_PISTON_AP / area_nodal;
	_NODOS_PIS.at(1).A_S = _NODOS_PIS.at(1).A_S * factor_area;
	_NODOS_PIS.at(1).A_LD = _NODOS_PIS.at(1).A_LD * factor_area;
	_NODOS_PIS.at(7).A_S = _NODOS_PIS.at(7).A_S * factor_area;
	_NODOS_PIS.at(9).A_S = _NODOS_PIS.at(9).A_S * factor_area;
	_NODOS_PIS.at(9).A_LI = _NODOS_PIS.at(9).A_LI * factor_area;



}


void TransmisionCalor::crear_nodos_culata()
{

	double Lcul;
	double L_per;
	double centro_v;
	double c_val;

	//Proporcion entre el diámetro de la válvula y el diamétro del vástago
	double F_VAST_ADM;
	double F_VAST_ESC;

	double Espesor_zona1_2;

	double STotal;
	double xtotal;
	double Ytotal;
	double Siny;
	double xiny;
	double yiny;
	double SVE;
	double xVE;
	double yVE;
	double SVA;
	double xVA;
	double yVA;
	double SCE;
	double xCE;
	double yCE;
	double SCA;
	double xCA;
	double yCA;
	double SCentral1;
	double xCentral1;
	double yCentral1;

	double SCentral;
	double xCentral;
	double yCentral;

	double ST1;
	double ST2;
	double SC3E;
	double SC3A;
	double SEE;
	double SAA;
	//Dim SEA As Double
	//Dim SAE As Double
	double xT1;
	double xT2;
	double xC3E;
	double xC3A;
	double xEE;
	double xAA;
	double xEA;
	double xAE;

	double yT1;
	double yT2;
	double yC3E;
	double yC3A;
	double yEE;
	double yAA;
	double yEA;
	double yAE;

	double SP1_A;
	double SP2_A;
	double SP_A;
	double SP1_E;
	double SP2_E;
	double SP_E;

	double xP1_A;
	double xP2_A;
	double xP_A;
	double xP1_E;
	double xP2_E;
	double xP_E;

	double yP1_A;
	double yP2_A;
	double yP_A;
	double yP1_E;
	double yP2_E;
	double yP_E;

	double zP1_A;
	double zP2_A;
	double zP_A;
	double zP1_E;
	double zP2_E;
	double zP_E;

	double dp;
	double espesor;
	double DVA;
	double DVE;
	double DINY;
	double espesor_FD;
	double espesor_pipas;
	double C_ECIL;
	double C_ECUL;
	double C_DINY;

	short cont;

	_NODOS_CUL.resize(36);

	C_ECIL = _Inicializar.Constante_TRANS_CALOR_C_ECIL;
	C_ECUL = _Inicializar.Constante_TRANS_CALOR_C_ECUL;
	F_VAST_ADM = _Inicializar.Constante_TRANS_CALOR_C_F_VAST_ADM;
	F_VAST_ESC = _Inicializar.Constante_TRANS_CALOR_C_F_VAST_ESC;
	C_DINY = _Inicializar.Constante_TRANS_CALOR_C_DINY;


	//FACTORES
	//F_VAST_ADM = 0.2
	//F_VAST_ESC = 0.25

	//VALORES NECESARIOS
	dp = _Inicializar.ConfMotor_BASE_D;
	espesor = C_ECIL * dp;
	DVA = _Inicializar.ConfMotor_CULATA_DVA;
	DVE = _Inicializar.ConfMotor_CULATA_DVE;

	DINY = dp * C_DINY;

	espesor_FD = C_ECUL * dp;

	espesor_pipas = C_ECIL * dp;


	//'CALCULOS
	Lcul = dp + 2 * espesor;
	L_per = Lcul / 2;
	centro_v = L_per / 2;

	//Posicion del centro de la válvulas tanto en x como en y
	//c_val = (dp / 2 - DINY / 2) / 2 ^ (1 / 2) 'El cateto sabiendo la hipotenusa
	double Raiz2 = pow(2, 0.5);
	c_val = (dp / 4 + DINY / 4) / (Raiz2);

	//Alto de la culata
	if ((dp - c_val + DVA / 2 + 2 * espesor_pipas / 2) > dp) {
		Espesor_zona1_2 = dp - c_val + DVA / 2 + 2 * espesor_pipas / 2;
	}
	else {
		Espesor_zona1_2 = dp;
	}


	///'Posiciones cdg de valvulas, inyector y zona central
	//sector de culata cuadrado
	STotal = c_val * c_val;
	xtotal = c_val / 2;
	Ytotal = c_val / 2;
	//inyector
	Siny = Constantes::pi / 4 * pow((DINY / 2), 2);
	xiny = 4. / 3. * (Raiz2) / Constantes::pi  * (DINY / 2) * cos(Constantes::pi / 4);
	yiny = 4. / 3. * (Raiz2) / Constantes::pi  * (DINY / 2) * sin(Constantes::pi / 4);
	//valvula de escape
	SVE = Constantes::pi / 4 * pow((DVE / 2), 2);
	xVE = (Raiz2 * c_val - 4. / 3. * (Raiz2) / Constantes::pi  * (DVE / 2)) * cos(Constantes::pi / 4);
	yVE = (Raiz2 * c_val - 4. / 3. * (Raiz2) / Constantes::pi * (DVE / 2)) * sin(Constantes::pi / 4);
	//valvula de admision
	SVA = Constantes::pi / 4 * pow((DVA / 2), 2);
	xVA = (Raiz2 * c_val - 4. / 3. * (Raiz2) / Constantes::pi  * (DVA / 2)) * cos(Constantes::pi / 4);
	yVA = (Raiz2 * c_val - 4. / 3. * (Raiz2) / Constantes::pi  * (DVA / 2)) * sin(Constantes::pi / 4);
	//resto del sector cuadrado central (escape)
	SCE = STotal - Siny - SVE;
	xCE = (xtotal * STotal - xiny * Siny - xVE * SVE) / SCE;
	yCE = (Ytotal * STotal - yiny * Siny - yVE * SVE) / SCE;
	//resto del sector cuadrado central (admision)
	SCA = STotal - Siny - SVA;
	xCA = (xtotal * STotal - xiny * Siny - xVA * SVA) / SCA;
	yCA = (Ytotal * STotal - yiny * Siny - yVA * SVA) / SCA;
	//sector rectangular central entre valvulas (2 valvulas - adm+esc -)
	SCentral1 = SCE + SCA;
	//sector rectangular central entre valvulas  (4 valvulas - 2adm+2esc -)
	xCentral = (SCE * yCE - SCA * yCA) / SCentral1;
	SCentral = 2 * SCentral1;
	yCentral = 0;


	//Centros de gravedad material culata (exterior) contacto gases EAS,EES,EAI,EEI;AES,AAS,AAI
	//La primera letra: Escape-Admisión
	//La segunda: con quien contacta Escape-Admisión
	//La tercera: Superior e Inferior

	ST1 = L_per * L_per / 2;
	ST2 = c_val * c_val / 2;
	SC3E = (Constantes::pi / 2 + Constantes::pi / 4) * 1. / 2. * pow((DVE / 2), 2);
	SC3A = (Constantes::pi / 2 + Constantes::pi / 4) * 1. / 2. * pow((DVA / 2), 2);
	SEE = ST1 - ST2 - SC3E;
	SAA = ST1 - ST2 - SC3A;

	xT1 = 2. / 3. * L_per;
	xT2 = 2. / 3. * c_val;
	xC3E = c_val + 4. / 3. * sin(3 * Constantes::pi / 4) / (3 * Constantes::pi / 4) * (DVE / 2) * sin((Constantes::pi / 2 + Constantes::pi / 4) / 2);
	xC3A = c_val + 4. / 3. * sin(3 * Constantes::pi / 4) / (3 * Constantes::pi / 4) * (DVA / 2) * sin((Constantes::pi / 2 + Constantes::pi / 4) / 2);
	xEE = (xT1 * ST1 - xT2 * ST2 - xC3E * SC3E) / SEE;
	xAA = (xT1 * ST1 - xT2 * ST2 - xC3A * SC3A) / SAA;

	yT1 = 1. / 3. * L_per;
	yT2 = 1. / 3. * c_val;

	yC3E = c_val - 4. / 3. * sin(3 * Constantes::pi / 4) / (3 * Constantes::pi / 4) * (DVE / 2) * cos((Constantes::pi / 2 + Constantes::pi / 4) / 2);
	yC3A = c_val - 4. / 3. * sin(3 * Constantes::pi / 4) / (3 * Constantes::pi / 4) * (DVA / 2) * cos((Constantes::pi / 2 + Constantes::pi / 4) / 2);
	yEE = (yT1 * ST1 - yT2 * ST2 - yC3E * SC3E) / SEE;
	yAA = (yT1 * ST1 - yT2 * ST2 - yC3A * SC3A) / SAA;

	xEA = yEE;
	xAE = yAA;

	yEA = xEE;
	yAE = xAA;


	//cambios 20/04/2012
	if (_Inicializar.ConfMotor_CULATA_A_PIPA_ADM == 0) {
		SP1_A = Constantes::pi * DVA * (dp / 2 - c_val);
		SP2_A = Constantes::pi * DVA * (dp / 2 - c_val);
		SP_A = SP1_A + SP2_A;
	}
	else {
		SP_A = _Inicializar.ConfMotor_CULATA_A_PIPA_ADM;
	}
	_Inicializar.ConfMotor_CULATA_A_PIPA_ADM = SP_A;

	if (_Inicializar.ConfMotor_CULATA_A_PIPA_ESC == 0) {
		SP1_E = Constantes::pi * DVE * (dp / 2 - c_val);
		SP2_E = Constantes::pi * DVE * (dp / 2 - c_val);
		SP_E = SP1_E + SP2_E;
	}
	else {
		SP_E = _Inicializar.ConfMotor_CULATA_A_PIPA_ESC;
	}
	_Inicializar.ConfMotor_CULATA_A_PIPA_ESC = SP_E;


	xP1_A = -c_val;
	xP2_A = -(L_per - (dp / 2 - c_val) / 2);
	xP_A = (xP1_A * SP1_A + xP2_A * SP2_A) / SP_A;
	xP1_E = c_val;
	xP2_E = L_per - (dp / 2 - c_val) / 2;
	xP_E = (xP1_E * SP1_E + xP2_E * SP2_E) / SP_E;

	yP1_A = c_val;
	yP2_A = c_val;
	yP_A = (yP1_A * SP1_A + yP2_A * SP2_A) / SP_A;
	yP1_E = c_val;
	yP2_E = c_val;
	yP_E = (yP1_E * SP1_E + yP2_E * SP2_E) / SP_E;

	zP1_A = (dp / 2 - c_val) / 2 + espesor_FD;
	zP2_A = dp / 2 - c_val + espesor_FD;
	zP_A = (zP1_A * SP1_A + zP2_A * SP2_A) / SP_A;

	zP1_E = (dp / 2 - c_val) / 2 + espesor_FD;
	zP2_E = dp / 2 - c_val + espesor_FD;
	zP_E = (zP1_E * SP1_E + zP2_E * SP2_E) / SP_E;

	//UPGRADE_WARNING: El límite inferior de la matriz _NODOS_CUL ha cambiado de 1 a 0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	// ERROR: Not supported in C#: ReDimStatement


	for (cont = 1; cont <= 35; cont++) {
		_pos_matriz++;
		_NODOS_CUL.at(cont).ID_NODO = _pos_matriz;
		_NODOS_CUL.at(cont).numero = cont;

		switch (cont)
		{
		case 1:
			_NODOS_CUL.at(cont).NOMBRE = "CUL1_EAS_A";
			break;
		case 2:
			_NODOS_CUL.at(cont).NOMBRE = "CUL2_EES_A";
			break;
		case 3:
			_NODOS_CUL.at(cont).NOMBRE = "CUL3_EEI_A";
			break;
		case 4:
			_NODOS_CUL.at(cont).NOMBRE = "CUL4_EAI_A";
			break;
		case 5:
			_NODOS_CUL.at(cont).NOMBRE = "CUL5_AEI_A";
			break;
		case 6:
			_NODOS_CUL.at(cont).NOMBRE = "CUL6_AAI_A";
			break;
		case 7:
			_NODOS_CUL.at(cont).NOMBRE = "CUL7_AAS_A";
			break;
		case 8:
			_NODOS_CUL.at(cont).NOMBRE = "CUL8_AES_A";
			break;
		case 9:
			_NODOS_CUL.at(cont).NOMBRE = "CUL9_CENTRAL_A";
			break;
		case 10:
			_NODOS_CUL.at(cont).NOMBRE = "CUL10_INY_A";
			break;
		case 11:
			_NODOS_CUL.at(cont).NOMBRE = "CUL11_EAS_B";
			break;
		case 12:
			_NODOS_CUL.at(cont).NOMBRE = "CUL12_EES_B";
			break;
		case 13:
			_NODOS_CUL.at(cont).NOMBRE = "CUL13_EEI_B";
			break;
		case 14:
			_NODOS_CUL.at(cont).NOMBRE = "CUL14_EAI_B";
			break;
		case 15:
			_NODOS_CUL.at(cont).NOMBRE = "CUL15_AEI_B";
			break;
		case 16:
			_NODOS_CUL.at(cont).NOMBRE = "CUL16_AAI_B";
			break;
		case 17:
			_NODOS_CUL.at(cont).NOMBRE = "CUL17_AAS_B";
			break;
		case 18:
			_NODOS_CUL.at(cont).NOMBRE = "CUL18_AES_B";
			break;
		case 19:
			_NODOS_CUL.at(cont).NOMBRE = "CUL19_CENTRAL_B";
			break;
		case 20:
			_NODOS_CUL.at(cont).NOMBRE = "CUL20_INY_B";
			break;
		case 21:
			_NODOS_CUL.at(cont).NOMBRE = "CUL21_VLV_AS";
			break;
		case 22:
			_NODOS_CUL.at(cont).NOMBRE = "CUL22_VLV_AI";
			break;
		case 23:
			_NODOS_CUL.at(cont).NOMBRE = "CUL23_VLV_ES";
			break;
		case 24:
			_NODOS_CUL.at(cont).NOMBRE = "CUL24_VLV_EI";
			break;
		case 25:
			_NODOS_CUL.at(cont).NOMBRE = "CUL25_VAST_AS";
			break;
		case 26:
			_NODOS_CUL.at(cont).NOMBRE = "CUL26_VAST_AI";
			break;
		case 27:
			_NODOS_CUL.at(cont).NOMBRE = "CUL27_VAST_ES";
			break;
		case 28:
			_NODOS_CUL.at(cont).NOMBRE = "CUL28_VAST_EI";
			break;
		case 29:
			_NODOS_CUL.at(cont).NOMBRE = "CUL29_INY";
			break;
		case 30:
			_NODOS_CUL.at(cont).NOMBRE = "CUL30_PIP_AS";
			break;
		case 31:
			_NODOS_CUL.at(cont).NOMBRE = "CUL31_PIP_AI";
			break;
		case 32:
			_NODOS_CUL.at(cont).NOMBRE = "CUL32_PIP_ES";
			break;
		case 33:
			_NODOS_CUL.at(cont).NOMBRE = "CUL33_PIP_EI";
			break;
		case 34:
			_NODOS_CUL.at(cont).NOMBRE = "CUL34_MA";
			break;
		case 35:
			_NODOS_CUL.at(cont).NOMBRE = "CUL35_ME";
			break;

		}
	}

	for (cont = 1; cont <= 2; cont++) {
		_NODOS_CUL.at(10 * (cont - 1) + 1).x = xEA;
		_NODOS_CUL.at(10 * (cont - 1) + 2).x = xEE;
		_NODOS_CUL.at(10 * (cont - 1) + 3).x = xEE;
		_NODOS_CUL.at(10 * (cont - 1) + 4).x = xEA;
		_NODOS_CUL.at(10 * (cont - 1) + 5).x = -xAE;
		_NODOS_CUL.at(10 * (cont - 1) + 6).x = -xAA;
		_NODOS_CUL.at(10 * (cont - 1) + 7).x = -xAA;
		_NODOS_CUL.at(10 * (cont - 1) + 8).x = -xAE;
		_NODOS_CUL.at(10 * (cont - 1) + 9).x = xCentral;
		_NODOS_CUL.at(10 * (cont - 1) + 10).x = 0;


		_NODOS_CUL.at(10 * (cont - 1) + 1).Y = yEA;
		_NODOS_CUL.at(10 * (cont - 1) + 2).Y = yEE;
		_NODOS_CUL.at(10 * (cont - 1) + 3).Y = -yEE;
		_NODOS_CUL.at(10 * (cont - 1) + 4).Y = -yEA;
		_NODOS_CUL.at(10 * (cont - 1) + 5).Y = -yAE;
		_NODOS_CUL.at(10 * (cont - 1) + 6).Y = -yAA;
		_NODOS_CUL.at(10 * (cont - 1) + 7).Y = yAA;
		_NODOS_CUL.at(10 * (cont - 1) + 8).Y = yAE;
		_NODOS_CUL.at(10 * (cont - 1) + 9).Y = yCentral;
		_NODOS_CUL.at(10 * (cont - 1) + 10).Y = 0;
	}


	_NODOS_CUL.at(21).x = -c_val;
	_NODOS_CUL.at(22).x = -c_val;
	_NODOS_CUL.at(23).x = c_val;
	_NODOS_CUL.at(24).x = c_val;
	_NODOS_CUL.at(25).x = -c_val;
	_NODOS_CUL.at(26).x = -c_val;
	_NODOS_CUL.at(27).x = c_val;
	_NODOS_CUL.at(28).x = c_val;
	_NODOS_CUL.at(29).x = 0;

	_NODOS_CUL.at(21).Y = c_val;
	_NODOS_CUL.at(22).Y = -c_val;
	_NODOS_CUL.at(23).Y = c_val;
	_NODOS_CUL.at(24).Y = -c_val;
	_NODOS_CUL.at(25).Y = c_val;
	_NODOS_CUL.at(26).Y = -c_val;
	_NODOS_CUL.at(27).Y = c_val;
	_NODOS_CUL.at(28).Y = -c_val;
	_NODOS_CUL.at(29).Y = 0;



	_NODOS_CUL.at(30).x = xP_A;
	_NODOS_CUL.at(31).x = xP_A;
	_NODOS_CUL.at(32).x = xP_E;
	_NODOS_CUL.at(33).x = xP_E;

	_NODOS_CUL.at(30).Y = yP_A;
	_NODOS_CUL.at(31).Y = -yP_A;
	_NODOS_CUL.at(32).Y = yP_E;
	_NODOS_CUL.at(33).Y = -yP_E;

	_NODOS_CUL.at(30).Z = zP_A;
	_NODOS_CUL.at(31).Z = zP_A;
	_NODOS_CUL.at(32).Z = zP_E;
	_NODOS_CUL.at(33).Z = zP_E;


	_NODOS_CUL.at(34).x = -(L_per / 2 - DINY / 2);
	_NODOS_CUL.at(35).x = (L_per / 2 - DINY / 2);

	_NODOS_CUL.at(34).Y = 0;
	_NODOS_CUL.at(35).Y = 0;


	_NODOS_CUL.at(34).Z = Espesor_zona1_2 / 2 + espesor_FD;
	_NODOS_CUL.at(35).Z = Espesor_zona1_2 / 2 + espesor_FD;


	///'Zs
	for (cont = 1; cont <= 10; cont++) {
		_NODOS_CUL.at(cont).Z = espesor_FD / 2 / 2;
	}

	for (cont = 11; cont <= 20; cont++) {
		_NODOS_CUL.at(cont).Z = espesor_FD / 2 + espesor_FD / 2 / 2;
	}

	for (cont = 21; cont <= 24; cont++) {
		_NODOS_CUL.at(cont).Z = espesor_FD / 2 / 2;
	}



	//'Febrero 2012: dia 17: signos

	for (cont = 25; cont <= 26; cont++) {
		_NODOS_CUL.at(cont).Z = (dp / 2 - c_val - DVA / 2) / 2 + espesor_FD;
	}

	for (cont = 27; cont <= 28; cont++) {
		_NODOS_CUL.at(cont).Z = (dp / 2 - c_val - DVE / 2) / 2 + espesor_FD;
	}

	_NODOS_CUL.at(29).Z = (dp / 2 - c_val - DVA / 2) / 2 + espesor_FD;

	//For cont = 25 To 28
	//    _NODOS_CUL.at(cont).Z = espesor_FD / 2 + Espesor_zona1_2 / 2
	//Next cont
	//_NODOS_CUL.at(29).Z = espesor_FD + Espesor_zona1_2 / 2



	//Superficies contacto GAS cilindro
	double pot1 = pow((dp / 2), 2);
	for (cont = 1; cont <= 4; cont++) {
		_NODOS_CUL.at(cont).AREA_GAS = SEE - 1. / 2. * (L_per * L_per - Constantes::pi / 4 * pot1);
	}

	for (cont = 5; cont <= 8; cont++) {
		_NODOS_CUL.at(cont).AREA_GAS = SAA - 1. / 2. * (L_per * L_per - Constantes::pi / 4 * pot1);
	}

	_NODOS_CUL.at(9).AREA_GAS = SCentral;
	double pot2 = pow((DINY / 2), 2);
	_NODOS_CUL.at(10).AREA_GAS = Constantes::pi * pot2;

	//'valvulas ADMISION
	double pot3 = pow((DVA / 2), 2);
	for (cont = 21; cont <= 22; cont++) {
		_NODOS_CUL.at(cont).AREA_GAS = Constantes::pi * pot3;
	}

	//'valvulas Escape
	double pot4 = pow((DVE / 2), 2);
	for (cont = 23; cont <= 24; cont++) {
		_NODOS_CUL.at(cont).AREA_GAS = Constantes::pi * pot4;
	}

	//Superficies contacto GAS admisión

	//Valvula
	double pot5 = pow((F_VAST_ADM * DVA), 2);
	for (cont = 21; cont <= 22; cont++) {
		_NODOS_CUL.at(cont).AREA_ADM = Constantes::pi / 4 * (pot3 - pot5);
	}


	//Vastago
	//For cont = 25 To 26
	//    _NODOS_CUL.at(cont).AREA_ADM = Constantes::pi * (F_VAST_ADM * DVA) * (dp - c_val + espesor_FD / 2 + DVA / 2)
	//Next cont
	for (cont = 25; cont <= 26; cont++) {
		_NODOS_CUL.at(cont).AREA_ADM = Constantes::pi * (F_VAST_ADM * DVA) * (dp / 2 - c_val + espesor_FD / 2 - DVA / 2);
	}
	//Pipas
	for (cont = 30; cont <= 31; cont++) {
		_NODOS_CUL.at(cont).AREA_ADM = SP_A;
	}


	//Superficies contacto GAS escape

	//Valvula
	double pot6 = pow(DVE, 2);
	double pot7 = pow((F_VAST_ESC * DVE), 2);
	for (cont = 23; cont <= 24; cont++) {
		_NODOS_CUL.at(cont).AREA_ESC = Constantes::pi / 4 * (pot6 - pot7);
	}

	//Vastago
	//For cont = 27 To 28
	//    _NODOS_CUL.at(cont).AREA_ESC = Constantes::pi * (F_VAST_ESC * DVE) * (dp - c_val + espesor_FD / 2 + DVE / 2)
	//Next cont

	//Febrero 2012
	for (cont = 27; cont <= 28; cont++) {
		_NODOS_CUL.at(cont).AREA_ESC = Constantes::pi * (F_VAST_ESC * DVE) * (dp / 2 - c_val + espesor_FD / 2 - DVE / 2);
	}



	//Pipas
	for (cont = 32; cont <= 33; cont++) {
		_NODOS_CUL.at(cont).AREA_ESC = SP_E;
	}


	//'Superficies contacto COOLANT


	//Admision
	//
	//_NODOS_CUL.at(34).AREA_COOL = L_per * L_per - 1. / 2. * (Constantes::pi * (DINY / 2) ^ 2) - 2 * (Constantes::pi * (DVA / 2) ^ 2)
	//'Escape
	//_NODOS_CUL.at(35).AREA_COOL = L_per * L_per - 1. / 2. * (Constantes::pi * (DINY / 2) ^ 2) - 2 * (Constantes::pi * (DVE / 2) ^ 2)

	//Febrero 2012
	//_NODOS_CUL.at(34).AREA_COOL = Lcul * Lcul + 2 * L_per * Espesor_zona1_2 + Lcul * Espesor_zona1_2 - 2 * Constantes::pi / 4 * DVA ^ 2 'Febrero 2012
	//Escape
	//_NODOS_CUL.at(35).AREA_COOL = Lcul * Lcul + 2 * L_per * Espesor_zona1_2 + Lcul * Espesor_zona1_2 - 2 * Constantes::pi / 4 * DVE ^ 2 'Febrero 2012


	//cambios 27/01/2014

	double pot8 = pow((L_per - c_val), 2);
	double pot9 = pow((DVA / 2 + espesor_pipas), 2);
	double pot10 = pow((Lcul / 2), 2);
	double pot11 = pow(c_val, 2);
	double pot12 = pow((DVE / 2 + espesor_pipas), 2);
	double pot13 = pow((DVA / 2), 2);
	double pot14 = pow((DVE / 2), 2);
	//Admision 'Modificadas para nuevo conducto en la culata:
	_NODOS_CUL.at(15).AREA_COOL = 0.5 * pot8 - Constantes::pi * pot9 / 8;
	_NODOS_CUL.at(16).AREA_COOL = (pot10 - pot11 - 3. / 4. * Constantes::pi * pot9) / 2;
	_NODOS_CUL.at(17).AREA_COOL = (pot10 - pot11 - 3. / 4. * Constantes::pi * pot9) / 2;
	_NODOS_CUL.at(18).AREA_COOL = 0.5 * pot8 - Constantes::pi * pot9 / 8;
	_NODOS_CUL.at(30).AREA_COOL = Constantes::pi * (DVA / 2 + espesor_pipas) * 2. / 3. * (dp / 2 - c_val - DVA / 2 - espesor_pipas);
	_NODOS_CUL.at(31).AREA_COOL = Constantes::pi * (DVA / 2 + espesor_pipas) * 2. / 3. * (dp / 2 - c_val - DVA / 2 - espesor_pipas);
	_NODOS_CUL.at(34).AREA_COOL = Lcul * L_per + 2 * (L_per * Espesor_zona1_2 - (L_per - c_val) * 2. / 3. * (dp / 2 - c_val - DVA / 2 - espesor_pipas)) + Lcul * Espesor_zona1_2 - 2. / 3. * (dp / 2 - c_val - DVA / 2 - espesor_pipas) * Lcul - 2 * Constantes::pi * pot13 + (L_per - c_val) * Lcul - Constantes::pi * pot9 + (Lcul - 2 * DVA) * 2. / 3. * (dp / 2 - c_val - DVA / 2 - espesor_pipas);

	//Escape 'Modificadas para nuevo conducto en la culata:
	_NODOS_CUL.at(11).AREA_COOL = 0.5 * pot8 - Constantes::pi * pot12 / 8;
	_NODOS_CUL.at(12).AREA_COOL = (pot10 - pot11 - 3. / 4. * Constantes::pi * pot12) / 2;
	_NODOS_CUL.at(13).AREA_COOL = (pot10 - pot11 - 3. / 4. * Constantes::pi * pot12) / 2;
	_NODOS_CUL.at(14).AREA_COOL = 0.5 * pot8 - Constantes::pi * pot12 / 8;
	_NODOS_CUL.at(32).AREA_COOL = Constantes::pi * (DVE / 2 + espesor_pipas) * 2. / 3. * (dp / 2 - c_val - DVE / 2 - espesor_pipas);
	_NODOS_CUL.at(33).AREA_COOL = Constantes::pi * (DVE / 2 + espesor_pipas) * 2. / 3. * (dp / 2 - c_val - DVE / 2 - espesor_pipas);
	_NODOS_CUL.at(35).AREA_COOL = Lcul * L_per + 2 * (L_per * Espesor_zona1_2 - (L_per - c_val) * 2. / 3. * (dp / 2 - c_val - DVE / 2 - espesor_pipas)) + Lcul * Espesor_zona1_2 - 2. / 3. * (dp / 2 - c_val - DVE / 2 - espesor_pipas) * Lcul - 2 * Constantes::pi * pot14 + (L_per - c_val) * Lcul - Constantes::pi * pot12 + (Lcul - 2 * DVE) * 2. / 3. * (dp / 2 - c_val - DVE / 2 - espesor_pipas);



	for (cont = 1; cont <= 35; cont++) {
		_NODOS_CUL.at(cont).k = _Inicializar.ConfMotor_CULATA_CONDUCTIVIDAD;
		_NODOS_CUL.at(cont).c = 0;
		_NODOS_CUL.at(cont).dens = 0;
		if (_NODOS_CUL.at(cont).AREA_COOL < 0)
		{
			_NODOS_CUL.at(cont).AREA_COOL = 0;
		}
	}

}

double TransmisionCalor::resuelve_suma_geometrica(double suma, short n)
{
	///''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	//Se utiliza para obtener el factor de crecimiento geométrico
	//suma=[A^n - 1]/[A-1]
	// suma es la longitud total y n es el núemro de datos
	///''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	double A0;
	double a1;
	double A2;
	double suma_inf;
	double suma_sup;
	double suma_mean;

	//Valores iniciales
	A0 = 1.000000001;
	a1 = 100;
	suma_inf = (pow(A0, n) - 1) / (A0 - 1);
	suma_sup = (pow(a1, n) - 1) / (a1 - 1);

	//'Procedimiento -> busco el valor a la mitad de los dos anteriores que dejen en medio el buscado

	while (((suma - suma_inf) / suma) > 1E-07) {
		A2 = (A0 + a1) / 2;
		suma_mean = (pow(A2, n) - 1) / (A2 - 1);
		if (suma > suma_mean) {
			suma_inf = suma_mean;
			A0 = A2;
		}
		else if (suma < suma_mean) {
			suma_sup = suma_mean;
			a1 = A2;
		}
		else if (suma == suma_mean) {
			suma_inf = suma_mean;
			suma_sup = suma_mean;
			A0 = A2;
			a1 = A2;
		}
	}

	return A0;

}

void TransmisionCalor::calcular_conductancias_conductivas()
{
	double c_val;
	double DINY;
	double DVE;
	double DVA;
	double C_ANG_ASIENTO;
	double C_ECUL;
	double C_ECIL;
	int cont3;
	double filas_total;
	filas_total = 3600;

	short cont1;
	short cont2;
	short n1;
	short n2;
	double Area_cont;
	double K_piston;
	short ncilin;
	short npis;
	double d;

	short cont0;
	double l;
	double fi01;
	double xg1;
	double xg2;
	double xg3;
	double Sup1;
	double Sup2;
	double Sup3;
	double r01;
	double r02;

	stConductancias_Conduct reg_cond_cond;

	std::vector<double> spiston(filas_total + 1);

	std::vector<double> vec_cp_agua(7);
	std::vector<double> vec_rho_agua(7);
	std::vector<double> vec_mu_agua(7);
	std::vector<double> vec_k_agua(7);
	std::vector<double> vec_cp_oil(7);
	std::vector<double> vec_rho_oil(7);
	std::vector<double> vec_mu_oil(7);
	std::vector<double> vec_k_oil(7);
	double k;
	double espesor;

	double dp;
	int i;

	npis = (_NODOS_PIS).size() - 1;


	ncilin = (_NODOS_CIL).size() - 1;



	K_piston = _Inicializar.ConfMotor_PISTON_CONDUCTIVIDAD;
	dp = _Inicializar.ConfMotor_BASE_D;

	//UPGRADE_WARNING: El límite inferior de la matriz k_pis ha cambiado de 1,1 a 0,0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'



	//**********************PISTON - PISTON***************************
	n1 = 1;
	n2 = 2;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_I;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	n1 = 2;
	n2 = 3;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_I;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	i = 4;
	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	n1 = 2;
	n2 = 7;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_LD;

	l = _NODOS_PIS.at(n1).ALTURA;

	fi01 = 2 * Constantes::pi;

	xg1 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n2).X_SD);
	xg2 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_SD);
	xg3 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_II);

	Sup1 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n2).X_SD), 2);
	Sup2 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_SD), 2);
	Sup3 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_II), 2);
	r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1);
	r02 = (xg3 * Sup3 - xg2 * Sup2) / (Sup3 - Sup2);

	k = -conductancia_radial(K_piston, l, r01, r02, fi01);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	n1 = 3;
	n2 = 4;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n2).A_S;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	n1 = 3;
	n2 = 6;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n2).A_S;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	n1 = 3;
	n2 = 8;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_LD;

	l = _NODOS_PIS.at(n1).ALTURA;
	fi01 = 2 * Constantes::pi;

	xg1 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n2).X_SD);
	xg2 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_SD);
	xg3 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_II);

	Sup1 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n2).X_SD), 2);
	Sup2 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_SD), 2);
	Sup3 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_II), 2);
	r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1);
	r02 = (xg3 * Sup3 - xg2 * Sup2) / (Sup3 - Sup2);

	k = -conductancia_radial(K_piston, l, r01, r02, fi01);


	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	n1 = 4;
	n2 = 5;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_I;
	k = -K_piston * Area_cont / d;


	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	n1 = 4;
	n2 = 6;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_LD;

	l = _NODOS_PIS.at(n1).ALTURA;
	fi01 = 2 * Constantes::pi;

	xg1 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n2).X_SD);
	xg2 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_SD);
	xg3 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_II);

	Sup1 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n2).X_SD), 2);
	Sup2 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_SD), 2);
	Sup3 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_II), 2);
	r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1);
	r02 = (xg3 * Sup3 - xg2 * Sup2) / (Sup3 - Sup2);

	k = -conductancia_radial(K_piston, l, r01, r02, fi01);


	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	n1 = 7;
	n2 = 8;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_I;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	n1 = 7;
	n2 = 10;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_LD;

	l = _NODOS_PIS.at(n1).ALTURA;
	fi01 = 2 * Constantes::pi;

	xg1 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n2).X_SD);
	xg2 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_SD);
	xg3 = 4. / 3. * (sin(fi01 / 2) / fi01) * (dp / 2 - _NODOS_PIS.at(n1).X_II);

	Sup1 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n2).X_SD), 2);
	Sup2 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_SD), 2);
	Sup3 = fi01 / 2 * pow((dp / 2 - _NODOS_PIS.at(n1).X_II), 2);
	r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1);
	r02 = (xg3 * Sup3 - xg2 * Sup2) / (Sup3 - Sup2);

	k = -conductancia_radial(K_piston, l, r01, r02, fi01);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	n1 = 9;
	n2 = 10;

	d = sqrt(pow((_NODOS_PIS.at(n1).X_CG - _NODOS_PIS.at(n2).X_CG), 2) + pow((_NODOS_PIS.at(n1).Y_CG - _NODOS_PIS.at(n2).Y_CG), 2));
	Area_cont = _NODOS_PIS.at(n1).A_I;
	k = -K_piston * Area_cont / d;

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(n2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(n1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	//**********************PISTON - CILINDRO***************************

	double Lmaxim;
	Lmaxim = _NODOS_CIL.at(1).L_SUP;

	short nodo_cil_seg;
	std::vector<double> pos_seg(filas_total + 1);
	short cont4;
	double acum;
	double Kseg;
	short naxiales;
	short nradiales;
	short ncircunf;

	naxiales = _Inicializar.ConfMotor_BASE_NODOS_AXIALES;
	nradiales = _Inicializar.ConfMotor_BASE_NODOS_RADIALES;
	ncircunf = _Inicializar.ConfMotor_BASE_NODOS_CIRCUNF;

	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto filas_total. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'


	// ERROR: Not supported in C#: ReDimStatement

	// ERROR: Not supported in C#: ReDimStatement


	double ang;
	double dist1;
	ang = 0;
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto filas_total. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	for (cont1 = 1; cont1 <= filas_total; cont1++) {
		dist1 = _comun->Dist_Despl_inst(ang, (_Inicializar.ConfMotor_BASE_LM), (_Inicializar.ConfMotor_BASE_LB), (_Inicializar.ConfMotor_BASE_E));
		spiston.at(cont1) = _Inicializar.ConfMotor_PISTON_HLI + dist1;
		pos_seg.at(cont1) = spiston.at(cont1) + _NODOS_PIS.at(1).ALTURA + _NODOS_PIS.at(2).ALTURA + _NODOS_PIS.at(3).ALTURA / 2;
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto filas_total. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		ang = ang + 360 / filas_total;
	}

	acum = 0;

	for (cont1 = 1; cont1 <= naxiales; cont1++) {
		for (cont2 = 1; cont2 <= nradiales; cont2++) {
			for (cont3 = 1; cont3 <= ncircunf; cont3++) {
				if (cont2 == 1) {
					//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto cont3. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
					nodo_cil_seg = (cont1 - 1) * nradiales * ncircunf + (cont2 - 1) * ncircunf + cont3;
					//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto filas_total. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
					for (cont4 = 1; cont4 <= filas_total; cont4++) {
						if (pos_seg.at(cont4) > (Lmaxim - _NODOS_CIL.at(nodo_cil_seg).L_INF)) {
							acum = acum + 0;
						}
						else if (pos_seg.at(cont4) < (Lmaxim - _NODOS_CIL.at(nodo_cil_seg).L_SUP)) {
							acum = acum + 0;
						}
						else {
							//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto filas_total. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
							acum = acum + (_NODOS_PIS.at(3).ALTURA) / filas_total;
						}
					}

					Kseg = (_NODOS_CIL.at(nodo_cil_seg).ANG_FIN - _NODOS_CIL.at(nodo_cil_seg).ANG_INI) * Constantes::pi / 180 * 4 * acum / log((1 + _NODOS_CIL.at(nodo_cil_seg).D_MAX / _NODOS_CIL.at(nodo_cil_seg).D_MIN) / 2);


					if (Kseg != 0) {
						k = -Kseg;

						reg_cond_cond.ID_NODO1 = _NODOS_PIS.at(3).ID_NODO;
						reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(nodo_cil_seg).ID_NODO;
						reg_cond_cond.k = k;
						_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

						reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(nodo_cil_seg).ID_NODO;
						reg_cond_cond.ID_NODO2 = _NODOS_PIS.at(3).ID_NODO;
						reg_cond_cond.k = k;
						_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

					}
					acum = 0;
				}
			}
		}
	}



	//**********************CILINDRO - CILINDRO***************************


	short c1;
	short r1;
	short ax1;
	short ax2;
	short r2;
	short c2;
	double k_cil;
	double L1;
	double fi02;
	double L2;
	double resta;
	double re;
	double e_n1_n2;
	double rr;
	double Area;


	//Busco conexiones entre nodos

	for (cont1 = 1; cont1 <= ncilin; cont1++) {
		ax1 = _NODOS_CIL.at(cont1).POS_AX;
		r1 = _NODOS_CIL.at(cont1).POS_RAD;
		c1 = _NODOS_CIL.at(cont1).POS_CIRC;

		for (cont2 = cont1; cont2 <= ncilin; cont2++) {
			ax2 = _NODOS_CIL.at(cont2).POS_AX;
			r2 = _NODOS_CIL.at(cont2).POS_RAD;
			c2 = _NODOS_CIL.at(cont2).POS_CIRC;

			k_cil = (_NODOS_CIL.at(cont1).k + _NODOS_CIL.at(cont2).k) / 2;

			if ((ax1 == ax2)) {
				if (((c1 == c2) & (abs(r2 - r1) == 1)))
				{
					l = _NODOS_CIL.at(cont1).L_SUP - _NODOS_CIL.at(cont1).L_INF;
					fi01 = (_NODOS_CIL.at(cont1).ANG_INI - _NODOS_CIL.at(cont1).ANG_FIN) * Constantes::pi / 180;
					//cdg de los nodos
					xg1 = 4. / 3. * (sin(fi01 / 2) / fi01) * _NODOS_CIL.at(cont1).D_MIN / 2;
					xg2 = 4. / 3. * (sin(fi01 / 2) / fi01) * _NODOS_CIL.at(cont1).D_MAX / 2;
					xg3 = 4. / 3. * (sin(fi01 / 2) / fi01) * _NODOS_CIL.at(cont2).D_MAX / 2;
					Sup1 = fi01 / 2 * pow((_NODOS_CIL.at(cont1).D_MIN / 2), 2);
					Sup2 = fi01 / 2 * pow((_NODOS_CIL.at(cont1).D_MAX / 2), 2);
					Sup3 = fi01 / 2 * pow((_NODOS_CIL.at(cont2).D_MAX / 2), 2);
					r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1);
					r02 = (xg3 * Sup3 - xg2 * Sup2) / (Sup3 - Sup2);

					k = -conductancia_radial(k_cil, l, r01, r02, fi01);

					reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont1).ID_NODO;
					reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont2).ID_NODO;
					reg_cond_cond.k = k;
					_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


					reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont2).ID_NODO;
					reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont1).ID_NODO;
					reg_cond_cond.k = k;
					_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


				}

				if ((r1 == r2)) {
					if (c2 >= c1) {
						resta = c2 - c1;
					}
					else {
						resta = c1 - c2;
					}

					if (((resta == 1) | ((ncircunf - resta) == 1))) {
						l = _NODOS_CIL.at(cont1).L_SUP - _NODOS_CIL.at(cont1).L_INF;
						r01 = _NODOS_CIL.at(cont1).D_MIN;
						r02 = _NODOS_CIL.at(cont1).D_MAX;


						fi01 = abs((_NODOS_CIL.at(cont1).ANG_INI - _NODOS_CIL.at(cont1).ANG_FIN) * Constantes::pi / 180);

						//27/01/2014
						//xg1 = 4. / 3. * (Sin(fi01 / 2) / fi01) * _NODOS_CIL.at(cont1).D_MIN / 2
						//xg2 = 4. / 3. * (Sin(fi01 / 2) / fi01) * _NODOS_CIL.at(cont1).D_MAX / 2
						//Sup1 = fi01 / 2 * (_NODOS_CIL.at(cont1).D_MIN / 2) ^ 2
						//Sup2 = fi01 / 2 * (_NODOS_CIL.at(cont1).D_MAX / 2) ^ 2
						//r01 = (xg2 * Sup2 - xg1 * Sup1) / (Sup2 - Sup1) '27/01/2014
						//e_n1_n2 = r01 * Abs((2 * (1 - Cos(fi01))) ^ (1 / 2))

						rr = abs(r02 - r01) / 2;

						//27/01/2014
						//fi02 = ((_NODOS_CIL.at(cont1).ANG_INI + _NODOS_CIL.at(cont1).ANG_FIN) / 2) - ((_NODOS_CIL.at(cont2).ANG_INI + _NODOS_CIL.at(cont2).ANG_FIN) / 2) * Constantes::pi / 180

						re = (r01 / 2 + r02 / 2) / 2;
						e_n1_n2 = fi01 * re;

						//27/01/2014
						//e_n1_n2 = r01 * Abs((2 * (1 - Cos(fi01))) ^ (1 / 2))
						k = -k_cil * (l * rr) / e_n1_n2;


						reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont1).ID_NODO;
						reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont2).ID_NODO;
						reg_cond_cond.k = k;
						_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


						reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont2).ID_NODO;
						reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont1).ID_NODO;
						reg_cond_cond.k = k;
						_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


					}
				}
			}
			else {
				if (((c1 == c2) & (r1 == r2) & (abs(ax1 - ax2) == 1))) {
					L1 = (_NODOS_CIL.at(cont1).L_SUP + _NODOS_CIL.at(cont1).L_INF) / 2;
					L2 = (_NODOS_CIL.at(cont2).L_SUP + _NODOS_CIL.at(cont2).L_INF) / 2;
					l = abs(L1 - L2);

					r01 = _NODOS_CIL.at(cont1).D_MIN / 2;
					r02 = _NODOS_CIL.at(cont1).D_MAX / 2;
					fi01 = _NODOS_CIL.at(cont1).ANG_INI * Constantes::pi / 180;
					fi02 = _NODOS_CIL.at(cont1).ANG_FIN * Constantes::pi / 180;

					Area = 1. / 2. * abs((pow(r01, 2) - pow(r02, 2)) * (fi02 - fi01));
					k = -k_cil * Area / l;

					reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont1).ID_NODO;
					reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont2).ID_NODO;
					reg_cond_cond.k = k;
					_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


					reg_cond_cond.ID_NODO1 = _NODOS_CIL.at(cont2).ID_NODO;
					reg_cond_cond.ID_NODO2 = _NODOS_CIL.at(cont1).ID_NODO;
					reg_cond_cond.k = k;
					_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


				}
			}

		}
	}



	//**********************CULATA - CULATA***************************


	std::vector<double> K_CH_mat;
	double k_cul2;
	double kval;
	double dist;

	short nodo1;
	short nodo2;
	double L_per;

	double raiz2;
	double Kseat;
	//Proporcion entre el diámetor de la válvula y el diamétro del vástago
	double F_VAST_ADM;
	double F_VAST_ESC;
	double C_DINY;
	double C_K_CONTACT_VALV;

	//Proporción entre el coeficiente de película en la pipa y el del vástago
	double f_h_vast_adm;
	double f_h_vast_esc;

	double espesor_pipas;

	double espesor_FD;

	double Espesor_zona1_2;
	double L_cul;
	double A_seat_adm;
	double A_seat_esc;

	raiz2 = pow(2, 0.5);

	k_cul2 = _Inicializar.ConfMotor_CULATA_CONDUCTIVIDAD;

	C_ECIL = _Inicializar.Constante_TRANS_CALOR_C_ECIL;
	C_ECUL = _Inicializar.Constante_TRANS_CALOR_C_ECUL;
	F_VAST_ADM = _Inicializar.Constante_TRANS_CALOR_C_F_VAST_ADM;
	F_VAST_ESC = _Inicializar.Constante_TRANS_CALOR_C_F_VAST_ESC;
	C_DINY = _Inicializar.Constante_TRANS_CALOR_C_DINY;
	C_ANG_ASIENTO = _Inicializar.Constante_TRANS_CALOR_C_ANG_ASIENTO;
	C_K_CONTACT_VALV = _Inicializar.Constante_TRANS_CALOR_C_K_CONTACT_VALV;

	espesor = C_ECIL * dp;

	espesor_FD = C_ECUL * dp;

	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto C_ECIL. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	espesor_pipas = C_ECIL * dp;

	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	DVA = _Inicializar.ConfMotor_CULATA_DVA;
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	DVE = _Inicializar.ConfMotor_CULATA_DVE;
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	DINY = dp * C_DINY;

	//FACTORES
	//F_VAST_ADM = 0.2
	//F_VAST_ESC = 0.25
	//C_H_ASIENTO
	//
	//A_seat_adm = Constantes::pi * DVA * (espesor_FD / 2) / Sin(C_ANG_ASIENTO * Constantes::pi / 180) 'se asume altura de asiento de la válvula como la mitad del firedeck
	//A_seat_esc = Constantes::pi * DVE * (espesor_FD / 2) / Sin(C_ANG_ASIENTO * Constantes::pi / 180)
	//'17 de febrero 2012: es un tronco de cono y no un cono entero.
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto C_ANG_ASIENTO. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	A_seat_adm = (Constantes::pi * espesor_FD / (2 * sin(2 * C_ANG_ASIENTO * Constantes::pi / 180))) * (2 * DVA * cos(C_ANG_ASIENTO * Constantes::pi / 180) - espesor_FD * sin(C_ANG_ASIENTO * Constantes::pi / 180));
	//se asume altura de asiento de la válvula como la mitad del firedeck
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto C_ANG_ASIENTO. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	A_seat_esc = (Constantes::pi * espesor_FD / (2 * sin(2 * C_ANG_ASIENTO * Constantes::pi / 180))) * (2 * DVE * cos(C_ANG_ASIENTO * Constantes::pi / 180) - espesor_FD * sin(C_ANG_ASIENTO * Constantes::pi / 180));


	//Kseat = 0.7 * 3000 ''Conductancia de contacto válvulas culata, como en CALMEC
	Kseat = 0.7 * (0.002 / 0.0025) * C_K_CONTACT_VALV;

	//UPGRADE_WARNING: El límite inferior de la matriz K_CH_mat ha cambiado de 1,1 a 0,0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	// ERROR: Not supported in C#: ReDimStatement


	//Ya estaba definida en nodos culata
	L_cul = dp + 2 * espesor;
	L_per = (dp + 2 * espesor) / 2;
	//c_val = (dp / 2 - DINY / 2) / raiz2
	//Febrero 2012
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	c_val = (dp / 4 + DINY / 4) / (raiz2);



	//Alto de la culata
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	if ((dp - c_val + DVA / 2 + 2 * espesor_pipas / 2) > dp) {
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		Espesor_zona1_2 = dp - c_val + DVA / 2 + 2 * espesor_pipas / 2;
	}
	else {
		Espesor_zona1_2 = dp;
	}

	nodo1 = 1;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 2;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVE) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 8;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 23;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVE / 2) * espesor_FD / 2 'de pablo
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 11;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 2;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 3;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 23;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVE / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 12;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 3;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 4;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVE) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 24;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVE / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 13;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 4;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 5;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 24;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVE / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 14;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 5;


	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 6;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVA) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 22;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 15;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 6;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 7;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 22;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 16;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo1 = 7;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 8;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVA) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 21;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 17;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 8;

	nodo2 = 9;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 21;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (Constantes::pi / 2 + Constantes::pi / 4) * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2 + Constantes::pi / 4) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 18;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow(L_per, 2)) / 2 - (pow(c_val, 2)) / 2 - 3 * Constantes::pi / 8 * pow((DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 9;

	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = pow((2 * c_val), 2) - 1. / 2. * Constantes::pi * (pow((DVA / 2), 2) + pow((DVE / 2), 2)) - Constantes::pi * pow((DINY / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 10;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DINY / 2) * espesor_FD / 2;
	//=2*Constantes::pi()*($B$3/2)*$B$4/2
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 22;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = Constantes::pi / 2 * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo2 = 21;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = Constantes::pi / 2 * (DVA / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2) / (2 * Constantes::pi) * A_seat_adm;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 24;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = Constantes::pi / 2 * (DVE / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 23;
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = Constantes::pi / 2 * (DVE / 2) * espesor_FD / 2
	Area = (Constantes::pi / 2) / (2 * Constantes::pi) * A_seat_esc;
	k = -Kseat * Area;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 10;

	nodo2 = 20;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((DINY / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 11;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 12;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVE) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 18;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 32;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 12;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 13;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 32;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 13;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 14;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVE) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 33;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);




	nodo1 = 14;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVE / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 15;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 33;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 15;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 16;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVA) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 31;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 16;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 17;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_per - c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 31;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo1 = 17;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 18;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (raiz2 * L_per - (DVA) / 2 - raiz2 * c_val) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 30;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 18;


	nodo2 = 19;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (c_val - DVA / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);




	nodo2 = 30;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 3 / 2 * Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 19;

	nodo2 = 20;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DINY / 2) * espesor_FD / 2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 30;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 31;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi / 4 * (pow((DVA / 2 + espesor_pipas), 2) - pow((DVA / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 32;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 33;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi / 4 * (pow((DVE / 2 + espesor_pipas), 2) - pow((DVE / 2), 2));
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	///febrero 2012'''''

	//nodo2 = 34
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (2 * c_val) ^ 2 - 2 * 1 / 4 * Constantes::pi * ((DVA / 2) ^ 2 + (DVE / 2) ^ 2) - Constantes::pi * (DINY / 2) ^ 2
	//k = -k_cul2 * Area / dist
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar


	nodo2 = 34;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (2 * c_val * c_val) - 1. / 2. * Constantes::pi * (pow((DVA / 2 + espesor_pipas), 2)) - 1. / 2. * Constantes::pi * pow((DINY / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo2 = 35;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (2 * c_val * c_val) - 1. / 2. * Constantes::pi * (pow((DVE / 2 + espesor_pipas), 2)) - 1. / 2. * Constantes::pi * pow((DINY / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);






	nodo1 = 20;

	nodo2 = 29;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((DINY / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 21;

	nodo2 = 25;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((F_VAST_ADM * DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 22;

	nodo2 = 26;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((F_VAST_ADM * DVA / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 23;

	nodo2 = 27;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((F_VAST_ESC * DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 24;

	nodo2 = 28;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * pow((F_VAST_ESC * DVE / 2), 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 25;

	nodo2 = 30;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * F_VAST_ADM * DVA * espesor_pipas;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	//febrero 2012:

	//nodo2 = 34
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (dp - (dp - c_val) - (DVA / 2) - espesor_pipas) * (2 * Constantes::pi * (F_VAST_ADM * DVA / 2))
	//k = -k_cul2 * Area / dist
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar


	nodo1 = 26;

	nodo2 = 31;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * F_VAST_ADM * DVA * espesor_pipas;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	//febrero 2012

	//nodo2 = 34
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (dp - (dp - c_val) - (DVA / 2) - espesor_pipas) * (2 * Constantes::pi * (F_VAST_ADM * DVA / 2))
	//k = -k_cul2 * Area / dist
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar


	nodo1 = 27;

	nodo2 = 32;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * F_VAST_ESC * DVE * espesor_pipas;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	//febrero 2012
	//nodo2 = 35
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (dp - (dp - c_val) - (DVE / 2) - espesor_pipas) * (2 * Constantes::pi * (F_VAST_ESC * DVE / 2))
	//k = -k_cul2 * Area / dist
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar



	nodo1 = 28;

	nodo2 = 33;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = Constantes::pi * F_VAST_ESC * DVE * espesor_pipas;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	//febrero 2012
	//nodo2 = 35
	//dist = distancia_nodos_CH(nodo1, nodo2)
	//Area = (dp - (dp - c_val) - (DVE / 2) - espesor_pipas) * (2 * Constantes::pi * (F_VAST_ESC * DVE / 2))
	//k = -k_cul2 * Area / dist
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar
	//Set reg_cond_cond = New CConductanciaConductiva
	//reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//reg_cond_cond.k = k
	//reg_cond_cond.salvar



	nodo1 = 29;

	nodo2 = 34;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = 1. / 2. * Constantes::pi * DINY * Espesor_zona1_2
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 1. / 2. * Constantes::pi * DINY * (dp / 2 - c_val - DVA / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo2 = 35;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = 1. / 2. * Constantes::pi * DINY * Espesor_zona1_2
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 1. / 2. * Constantes::pi * DINY * (dp / 2 - c_val - DVA / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 30;

	nodo2 = 34;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = Constantes::pi * (DVA + 2 * espesor_pipas) * (2 * (dp - c_val) + espesor_pipas)
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DVA + 2 * espesor_pipas) * (dp / 2 - c_val - DVA / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 31;


	nodo2 = 34;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = Constantes::pi * (DVA + 2 * espesor_pipas) * (2 * (dp - c_val) + espesor_pipas)
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DVA + 2 * espesor_pipas) * (dp / 2 - c_val - DVA / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 32;

	nodo2 = 35;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = Constantes::pi * (DVE + 2 * espesor_pipas) * (2 * (dp - c_val) + espesor_pipas)
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DVE + 2 * espesor_pipas) * (dp / 2 - c_val - DVE / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo1 = 33;

	nodo2 = 35;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//Area = Constantes::pi * (DVE + 2 * espesor_pipas) * (2 * (dp - c_val) + espesor_pipas)
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = 2 * Constantes::pi * (DVE + 2 * espesor_pipas) * (dp / 2 - c_val - DVE / 2);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	nodo1 = 34;

	nodo2 = 35;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DINY. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (L_cul - DINY) * Espesor_zona1_2;
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);




	///''''''''''''''''''''Febrero 2012

	nodo1 = 34;

	//cambios del  27/01/2014

	nodo2 = 15;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow((L_cul / 2), 2) - pow(c_val, 2) - 3. / 4. * Constantes::pi * pow((DVA / 2 + espesor_pipas), 2)) / 2 - (0.5 * pow((L_per - c_val), 2) - Constantes::pi * pow((DVA / 2 + espesor_pipas), 2) / 8);

	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);



	nodo2 = 18;
	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVA. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow((L_cul / 2), 2) - pow(c_val, 2) - 3. / 4. * Constantes::pi * pow((DVA / 2 + espesor_pipas), 2)) / 2 - (0.5 * pow((L_per - c_val), 2) - Constantes::pi * pow((DVA / 2 + espesor_pipas), 2) / 8);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	//For nodo2 = 15 To 18 'AEI_B 'AAI_B 'AAS_B 'AES_B
	//    dist = distancia_nodos_CH(nodo1, nodo2)
	//    Area = ((L_cul / 2) ^ 2) - (c_val ^ 2) - (3. / 4. * Constantes::pi / 4 * (DVA ^ 2)) / 2 'febrero 2012'
	//    k = -k_cul2 * Area / dist
	//    Set reg_cond_cond = New CConductanciaConductiva
	//    reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//    reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//    reg_cond_cond.k = k
	//    reg_cond_cond.salvar
	//    Set reg_cond_cond = New CConductanciaConductiva
	//    reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//    reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//    reg_cond_cond.k = k
	//    reg_cond_cond.salvar
	//Next nodo2



	nodo1 = 35;

	nodo2 = 11;

	//cambios del  27/01/2014

	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow((L_cul / 2), 2) - pow(c_val, 2) - 3. / 4. * Constantes::pi * pow((DVE / 2 + espesor_pipas), 2)) / 2 - (0.5 * pow((L_per - c_val), 2) - Constantes::pi * pow((DVE / 2 + espesor_pipas), 2) / 8);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	nodo2 = 14;

	dist = distancia_nodos_CH(nodo1, nodo2);
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto DVE. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto c_val. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	Area = (pow((L_cul / 2), 2) - pow(c_val, 2) - 3. / 4. * Constantes::pi * pow((DVE / 2 + espesor_pipas), 2)) / 2 - (0.5 * pow((L_per - c_val), 2) - Constantes::pi * pow((DVE / 2 + espesor_pipas), 2) / 8);
	k = -k_cul2 * Area / dist;

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);

	reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO;
	reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO;
	reg_cond_cond.k = k;
	_CONDUCTANCIAS_CONDUCT.push_back(reg_cond_cond);


	//For nodo2 = 11 To 14 'EAS_B 'EES_B 'EEI_B 'EAI_B
	//    dist = distancia_nodos_CH(nodo1, nodo2)
	//    Area = ((L_cul / 2) ^ 2) - (c_val ^ 2) - (3. / 4. * Constantes::pi / 4 * (DVE ^ 2)) / 2 'febrero 2012'
	//    k = -k_cul2 * Area / dist
	//    Set reg_cond_cond = New CConductanciaConductiva
	//    reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo1).ID_NODO
	//    reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo2).ID_NODO
	//    reg_cond_cond.k = k
	//    reg_cond_cond.salvar
	//    Set reg_cond_cond = New CConductanciaConductiva
	//    reg_cond_cond.ID_NODO1 = _NODOS_CUL.at(nodo2).ID_NODO
	//    reg_cond_cond.ID_NODO2 = _NODOS_CUL.at(nodo1).ID_NODO
	//    reg_cond_cond.k = k
	//    reg_cond_cond.salvar
	//Next nodo2




}

double TransmisionCalor::conductancia_radial(double k, double l, double r1, double r2, double FI)
{
	///'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	//Calcula la conductancia radial entre dos nodos a partir de:
	//la conductividad (k)
	//la longitud del cilindro (L)
	//los radios interiores y exteriores (r1 y r2)
	//el ángulo considerado (fi)
	///'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

	return abs(FI * k * l / log(r2 / r1));

}

double TransmisionCalor::distancia_nodos_CH(short nodo1, short nodo2)
{

	return pow(((pow((_NODOS_CUL.at(nodo1).x - _NODOS_CUL.at(nodo2).x), 2)) + (pow((_NODOS_CUL.at(nodo1).Y - _NODOS_CUL.at(nodo2).Y), 2)) + (pow((_NODOS_CUL.at(nodo1).Z - _NODOS_CUL.at(nodo2).Z), 2))), 0.5);
}


void TransmisionCalor::crear_nodos_fluidos()
{
	short i;
	short numero;

	_NODOS_FLUIDOS.resize(7 + _Inicializar.ConfMotor_BASE_NODOS_AXIALES - 1);

	_pos_matriz = 1;
	numero = 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "GAS" + std::to_string(numero);
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "GAS";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	_pos_matriz = _pos_matriz + 1;
	numero = numero + 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "ADM1";
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "GAS_ADM";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	_pos_matriz = _pos_matriz + 1;
	numero = numero + 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "ESC1";
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "GAS_ESC";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	_pos_matriz = _pos_matriz + 1;
	numero = numero + 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "OIL1";
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "ACEITE";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	_pos_matriz = _pos_matriz + 1;
	numero = numero + 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "COOL_CIL1";
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "AGUA";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	_pos_matriz = _pos_matriz + 1;
	numero = numero + 1;

	_NODOS_FLUIDOS.at(numero).ID_NODO = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).numero = numero;
	_NODOS_FLUIDOS.at(numero).pos_matriz = _pos_matriz;
	_NODOS_FLUIDOS.at(numero).nombre = "COOL_CUL1";
	_NODOS_FLUIDOS.at(numero).tipo_fluido = "AGUA";
	_NODOS_FLUIDOS.at(numero).tipo_nodo = "FLUIDO1";

	//****************************NODOS GAS**********************************
	numero = 1;

	for (i = 1; i <= _Inicializar.ConfMotor_BASE_NODOS_AXIALES - 1; i++) {
		_pos_matriz = _pos_matriz + 1;
		_NODOS_FLUIDOS.at(_pos_matriz).ID_NODO = _pos_matriz;
		_NODOS_FLUIDOS.at(_pos_matriz).numero = numero;
		_NODOS_FLUIDOS.at(_pos_matriz).pos_matriz = _pos_matriz;
		_NODOS_FLUIDOS.at(_pos_matriz).nombre = "GAS" + std::to_string(numero + 1);
		_NODOS_FLUIDOS.at(_pos_matriz).tipo_fluido = "GAS";
		_NODOS_FLUIDOS.at(_pos_matriz).tipo_nodo = "FLUIDO2";
		numero = numero + 1;
	}

}



void TransmisionCalor::CalcularTemperaturas(short Z, double *TCUL, double *TCULMAT, double *TCULVALV, double *TPIS, double *TCIL, double *TCULMAX, double *TCULMIN, double *TPISMAX, double *TPISMIN, double *TCILMAX, double *TCILMIN, double *TGASM)
{
	int n;
	double i;
	double cont;
	double hmean_oil;
	double T_Exhaust;
	double T_Intake;
	double twater_cul;
	double twater_cil;
	double toil;
	double f_h_vast_adm;
	double f_h_vast_esc;
	double f_h_vlv_esc;
	double f_h_vlv_adm;
	short filas_total;
	short npis;
	short ncul;
	short nodos_conv;
	short ngases;
	short ncilin;
	short naxiales;
	short nradiales;
	short ncircunf;
	double hmean_adm;
	double hmean_exh;
	double hmean_gas;
	double Tmean_gas;
	double CM;
	double vcc;
	double h0;
	double Lcil;
	double Dgal;
	short NPC;
	int num_nodos;
	double Re_oil;
	double Pr_oil;
	double Dh_gal_oil;

	_RESULTADOS_F.clear();
	_RESULTADOS_T.clear();
	_RESULTADOS_K.clear();

	//FILE *fp;
	//std::string fichero = "C:\\temp\\Var_inst_nodal_cil" + std::to_string(Z) + ".dat";


	//fp = fopen(fichero.c_str(), "w");
	//fprintf(fp, "posicion_ciclo, ang, area_cil, h_adm, h_esc, h_gas, hi, T, Tgas_media_pipaesc, Tgas_salida_pipaesc\n");

	//for (int posicion_ciclo = 1; posicion_ciclo < var_inst_nodal.size(); posicion_ciclo++)
	//{
	//	double ang = var_inst_nodal.at(posicion_ciclo).ang;
	//	double area_cil = var_inst_nodal.at(posicion_ciclo).area_cil;
	//	double h_adm = var_inst_nodal.at(posicion_ciclo).h_adm;
	//	double h_esc = var_inst_nodal.at(posicion_ciclo).h_esc;
	//	double h_gas = var_inst_nodal.at(posicion_ciclo).h_gas;
	//	double hi = var_inst_nodal.at(posicion_ciclo).hi;
	//	double T = var_inst_nodal.at(posicion_ciclo).T;
	//	double Tgas_media_pipaesc = var_inst_nodal.at(posicion_ciclo).Tgas_media_pipaesc;
	//	double Tgas_salida_pipaesc = var_inst_nodal.at(posicion_ciclo).Tgas_salida_pipaesc;
	//	//fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", posicion_ciclo, var_inst_nodal.at(posicion_ciclo).ang, var_inst_nodal.at(posicion_ciclo).area_cil, var_inst_nodal.at(posicion_ciclo).h_adm, var_inst_nodal.at(posicion_ciclo).h_esc, var_inst_nodal.at(posicion_ciclo).h_gas, var_inst_nodal.at(posicion_ciclo).hi, var_inst_nodal.at(posicion_ciclo).Qrad_cil, var_inst_nodal.at(posicion_ciclo).Qrad_cul, var_inst_nodal.at(posicion_ciclo).Qrad_cil, var_inst_nodal.at(posicion_ciclo).Qrad_pis, var_inst_nodal.at(posicion_ciclo).T, var_inst_nodal.at(posicion_ciclo).Tgas_media_pipaesc, var_inst_nodal.at(posicion_ciclo).Tgas_salida_pipaesc);
	//	fprintf(fp, "%d %f %f %f %f %f %f %f %f %f\n", posicion_ciclo, ang, area_cil, h_adm, h_esc, h_gas, hi, T, Tgas_media_pipaesc, Tgas_salida_pipaesc);

	//}

	//fclose(fp);



	std::vector<double> L_nodos;

	std::vector< std::vector<double> >  k_cil;
	std::vector< std::vector<double> >  k_pis;
	std::vector< std::vector<double> >  k_cul;
	std::vector< std::vector<double> >  k_cil_pis;
	std::vector< std::vector<double> >  k_cil_cul;
	std::vector< std::vector<double> >  k_pis_cul;
	std::vector< std::vector<double> >  k_pis_cil;
	std::vector< std::vector<double> >  k_cul_cil;
	std::vector< std::vector<double> >  k_cul_pis;

	std::vector<std::vector<std::string> > matriz_areas_cil;
	std::vector<std::vector<std::string> > matriz_areas_pis;
	std::vector<std::vector<std::string> > matriz_areas_cul;
	std::vector<std::vector<std::string> > matriz_geom_cil;

	std::vector< std::vector<double> > matriz_conductances;


	std::vector<double> matriz_temp;

	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		GenerarMatricesDeFicheros(&k_cil, &k_pis, &k_cul, &k_cil_pis, &k_cil_cul, &k_pis_cul, &k_pis_cil, &k_cul_cil, &k_cul_pis, &matriz_areas_cil, &matriz_areas_pis, &matriz_areas_cul, &matriz_geom_cil, &npis, &ncul, &naxiales, &nradiales, &ncircunf);
	}
	else
	{
		int num_nodos = _NODOS_CIL.size() - 1 + _NODOS_CUL.size() - 1 + _NODOS_PIS.size() - 1 + _NODOS_FLUIDOS.size() - 1;

		//matriz_conductances(num_nodos, std::vector<double>(num_nodos));

		GenerarMatrices(&k_cil, &k_pis, &k_cul, &k_cil_pis, &k_cil_cul, &k_pis_cul, &k_pis_cil, &k_cul_cil, &k_cul_pis, &npis, &ncul, &naxiales, &nradiales, &ncircunf);
	}

	ncilin = naxiales * nradiales * ncircunf + ncircunf;
	ngases = naxiales - 1;
	nodos_conv = 6; //(GAS, ESCAPE, ADMIS, OIL, AGUA_CIL, AGUA_CUL) Estos son Fijos;
	num_nodos = ncilin + npis + ncul + ngases;


	matriz_temp.resize(nodos_conv + ngases + ncilin + npis + ncul + 1);
	std::vector< std::vector<double> > K_fluidos(nodos_conv + 1, std::vector<double>(nodos_conv + 1));

	std::vector< std::vector<double> > K_gas(ngases + 1, std::vector<double>(ngases + 1));
	std::vector< std::vector<double> > K_fluidos_cil(ncilin + 1, std::vector<double>(nodos_conv + 1));
	std::vector< std::vector<double> > K_gas_cil(ncilin + 1, std::vector<double>(ngases + 1));
	std::vector< std::vector<double> > K_fluidos_pis(npis + 1, std::vector<double>(nodos_conv + 1));
	std::vector< std::vector<double> > K_fluidos_cul(ncul + 1, std::vector<double>(nodos_conv + 1));

	std::vector< std::vector<double> > m_up1(nodos_conv + 1, std::vector<double>(nodos_conv + 1));
	std::vector< std::vector<double> > m_up2(nodos_conv + 1, std::vector<double>(ngases + 1));

	std::vector< std::vector<double> > m_inf1(ngases + 1, std::vector<double>(nodos_conv + 1));
	std::vector< std::vector<double> > m_inf2(ngases + 1, std::vector<double>(ngases + 1));


	std::vector< std::vector<double> > area_inst_nodos_cil;

	std::vector< std::vector<double> > matriz_conductances2;

	std::vector< std::vector<double> > m_up3;
	std::vector< std::vector<double> > m_up4(nodos_conv + ngases + 1, std::vector<double>(ncilin + 1));
	std::vector< std::vector<double> > m_inf3;
	std::vector< std::vector<double> > m_inf4;

	std::vector< std::vector<double> > m_up5;
	std::vector< std::vector<double> > m_up6;
	std::vector< std::vector<double> > m_inf5;
	std::vector< std::vector<double> > m_inf6;

	std::vector< std::vector<double> > m_aux1(nodos_conv + ngases + 1, std::vector<double>(npis + 1));
	std::vector< std::vector<double> > m_aux2(npis + 1, std::vector<double>(ngases + 1));
	std::vector< std::vector<double> > m_aux3;

	std::vector< std::vector<double> >  m_up7;
	std::vector< std::vector<double> >  m_up8;
	std::vector< std::vector<double> >  m_inf7;
	std::vector< std::vector<double> >  m_inf8;

	std::vector< std::vector<double> >  m_aux4;
	std::vector< std::vector<double> >  m_aux5(nodos_conv + ngases + 1, std::vector<double>(ncul + 1));
	std::vector< std::vector<double> >  m_aux6;
	std::vector< std::vector<double> >  m_aux7(ncul + 1, std::vector<double>(ngases + 1));

	std::vector< std::vector<double> >  m_aux8;
	std::vector< std::vector<double> >  m_aux9;
	std::vector< std::vector<double> >  m_aux10;




	h0 = _Inicializar.ConfMotor_PISTON_HLI;
	Lcil = 1.1 * _Inicializar.ConfMotor_BASE_S + h0;
	Dgal = _Inicializar.ConfMotor_PISTON_DIAM_INT_GAL * _Inicializar.ConfMotor_BASE_D;


	//factores para calcular los coeficientes de películas en el vástago y en la parte interior de la cabeza de las válvulas
	f_h_vlv_esc = 0.5;
	f_h_vast_esc = 0.09;
	f_h_vlv_adm = 0.5;
	f_h_vast_adm = 0.086;
	toil = _Inicializar.VarMed_TAC;
	twater_cil = _Inicializar.VarMed_TRS;
	twater_cul = _Inicializar.VarMed_TRS;
	T_Intake = _Inicializar.VarMed_TA;
	T_Exhaust = _Inicializar.VarMed_TE;
	NPC = _Inicializar.ConfInstrumentacion_NPCL;

	//Calculo hs y Ts promedios, las fáciles
	Promedia_h(&hmean_adm, &hmean_exh, &hmean_gas);

	//hmean_adm = 1;
	//hmean_exh = 1;
	//hmean_gas = 1;

	tgas(&hmean_gas, &Tmean_gas, TGASM);

	//Tmean_gas = 300;

	CM = 2 * _Inicializar.ConfMotor_BASE_S * _Inicializar.VarMed_N;

	Re_oil = rho_oil(toil) * (2 * _Inicializar.ConfMotor_BASE_S * _Inicializar.VarMed_N) * (4 * _NODOS_PIS.at(3).ANCHURA) / mu_oil(toil);
	Pr_oil = mu_oil(toil) * cp_oil(toil) / k_oil(toil);

	hmean_oil = (k_oil(toil) / Dgal) * _Inicializar.ConfMotor_PISTON_PIS2OIL * (pow(Re_oil, _Inicializar.ConfMotor_PISTON_EX_RE_PIS2OIL)) * (pow(Pr_oil, (0.3 / 0.4 * _Inicializar.ConfMotor_PISTON_EX_PR_PIS2OIL)));

	///' CONVECCION y Relaciones conductivas pis-cil, cil-cul y pis-cul (no hay)
	///' Voy a crear submatrices conveccion:
	///' 1: Identidad (nodos conv) ----> K_fluidos
	///' 2: Gases cil (identidad naxiales*naxiales)   --->  K_gas, Necesaria porque los nodos cilindro ven T y h diferente,
	///' 3: K_fluidos_cil--> Conveccion con el cilindro (agua, ...) y h fija (nodos superiores, siempre en contacto gas)
	///' 4: Piston --> K_fluidos_pis
	///' 5: culata --> K_fluidos_cul
	///' 6: K_gas_cil: Resto nodos cilindro en contacto, en algun momento no toca al gas


	///' 1: Identidad (nodos conv) ----> K_fluidos
	for (cont = 1; cont <= nodos_conv; cont++) {
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto cont. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		K_fluidos.at(cont).at(cont) = 1;
	}


	///' 2: Gases cil (identidad (naxiales-1)*(naxiales-1))   --->  K_gas, Necesaria porque los nodos cilindro ven T y h diferente,
	for (cont = 1; cont <= ngases; cont++) {
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto cont. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		K_gas.at(cont).at(cont) = 1;
	}


	///' 3: K_fluidos_cil--> Conveccion con el cilindro (agua, ...) y h fija (nodos superiores, siempre en contacto con gas)
	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		calculo_vcil_avanzado(Lcil, &area_inst_nodos_cil, &matriz_geom_cil, naxiales, nradiales, ncircunf);
	}
	else
	{
		calculo_vcil(Lcil, &area_inst_nodos_cil);
	}


	short cont2;
	double areac;
	//Agua

	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		for (int i = 0; i < matriz_areas_cil.size(); i++)
		{
			cont2 = ObtenerPosicionNodoCilindro(matriz_areas_cil.at(0).at(i).c_str(), naxiales, nradiales, ncircunf);
			areac = atof(matriz_areas_cil.at(2).at(i).c_str());
			K_fluidos_cil.at(cont2).at(5) = -_Inicializar.ConfMotor_BASE_H_CIL * areac;
		}
	}
	else
	{

		double espesor = _Inicializar.Constante_TRANS_CALOR_C_ECIL * _Inicializar.ConfMotor_BASE_D;

		double esp_agua = espesor / 2;

		double pot1 = pow((_Inicializar.ConfMotor_BASE_D / 2 + espesor + esp_agua), 2);
		double pot2 = pow((_Inicializar.ConfMotor_BASE_D / 2 + espesor), 2);
		double A_h = (pot1 - pot2) * Constantes::pi / ncircunf;


		for (cont2 = 1; cont2 <= ncilin; cont2++)
		{
			//Veo si es radial final
			if (((_NODOS_CIL.at(cont2).POS_RAD == nradiales) & (_NODOS_CIL.at(cont2).POS_AX != 1) & (_NODOS_CIL.at(cont2).POS_CIRC != 1)))
			{
				areac = abs(((_NODOS_CIL.at(cont2).ANG_FIN - _NODOS_CIL.at(cont2).ANG_INI) / 360 * 2 * Constantes::pi) * (_NODOS_CIL.at(cont2).D_MAX / 2) * (_NODOS_CIL.at(cont2).L_INF - _NODOS_CIL.at(cont2).L_SUP));
				K_fluidos_cil.at(cont2).at(5) = -_Inicializar.ConfMotor_BASE_H_CIL * areac;
			}
			else if (((_NODOS_CIL.at(cont2).POS_RAD == nradiales + 1) & (_NODOS_CIL.at(cont2).POS_CIRC != 1))) // 'Josep 10/04/2015 voladizo en axial = 1, radial = nradiales + 1
			{
				K_fluidos_cil.at(cont2).at(5) = -_Inicializar.ConfMotor_BASE_H_CIL * A_h;
			}
		}
	}

	///' 6: K_gas_cil: Resto nodos cilindro en contacto con gas (no siempre)
	short n1;
	short cont1;
	short cont3;
	std::vector<double> Area2(NPC + 1);
	std::vector< std::vector<double> > hAg(ncilin + 1, std::vector<double>(ngases + 2));
	std::vector<double> matriz_temp_cil(ngases + 1);

	double Qrad_cil_med;
	double Qrad_cil_med_nodo;

	Qrad_cil_med = 0;
	for (i = 1; i <= NPC; i++) {
		Qrad_cil_med = Qrad_cil_med + var_inst_nodal.at(i).Qrad_cil;
	}
	Qrad_cil_med = Qrad_cil_med / NPC;


	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		int  POS_CIRC = 0;

		for (i = 0; i < matriz_geom_cil.size(); i++)
		{
			cont1 = atoi(matriz_geom_cil.at(1).at(i).substr(3, 1).c_str()); //número del gas
			cont2 = ObtenerPosicionNodoCilindro(matriz_geom_cil.at(0).at(i).c_str(), naxiales, nradiales, ncircunf);
			POS_CIRC = atoi(matriz_geom_cil.at(0).at(i).substr(5, 1).c_str());

			for (cont3 = 0; cont3 < NPC; cont3++)
			{
				Area2.at(cont3) = area_inst_nodos_cil.at(cont3).at(cont2);
			}
			if (cont1 == 1) // el gas1 es fluido
			{
				K_fluidos_cil.at(cont2).at(cont1) = -ha_mean(&Area2);
			}
			else
			{
				K_gas_cil.at(cont2).at(cont1 - 1) = -ha_mean(&Area2);
				if (POS_CIRC == 1)
				{
					matriz_temp_cil.at(cont1 - 1) = T_mean(&Area2);
				}
			}
			Qrad_cil_med_nodo = Qrad_mean(&Area2);
			matriz_temp.at(nodos_conv + ngases + cont2) = Qrad_cil_med_nodo;
		}
	}
	else
	{
		for (cont1 = 1; cont1 <= ngases + 1; cont1++)
		{
			for (cont2 = 1; cont2 <= ncilin; cont2++)
			{
				if ((_NODOS_CIL.at(cont2).POS_RAD == 1))
				{
					if (cont1 == _NODOS_CIL.at(cont2).POS_AX)
					{
						for (cont3 = 1; cont3 <= NPC; cont3++)
						{
							Area2.at(cont3) = area_inst_nodos_cil.at(cont3).at(cont2);
						}
						if (cont1 == 1) // el gas1 es fluido
						{
							K_fluidos_cil.at(cont2).at(cont1) = -ha_mean(&Area2);
						}
						else
						{
							K_gas_cil.at(cont2).at(cont1 - 1) = -ha_mean(&Area2);
							if (_NODOS_CIL.at(cont2).POS_CIRC == 1)
							{
								matriz_temp_cil.at(cont1 - 1) = T_mean(&Area2);
							}
						}
						Qrad_cil_med_nodo = Qrad_mean(&Area2);
						matriz_temp.at(nodos_conv + ngases + cont2) = Qrad_cil_med_nodo;
					}
				}
			}
		}
	}

	///' 4: Piston --> K_fluidos_pis
	//calcula Qrad_pis medio
	double Qrad_pis_med;
	double Atotal_pis_rad;

	Qrad_pis_med = 0;
	for (i = 1; i <= NPC; i++) {
		Qrad_pis_med = Qrad_pis_med + var_inst_nodal.at(i).Qrad_pis;
	}
	Qrad_pis_med = Qrad_pis_med / NPC;




	double Qrad_pis_med_nodo;

	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		Atotal_pis_rad = 0;

		for (i = 0; i < matriz_areas_pis.size(); i++)
		{
			n1 = atoi(matriz_areas_pis.at(1).at(i).substr(3, 1).c_str());
			if (n1 == 1)
			{
				Atotal_pis_rad = Atotal_pis_rad + atof(matriz_areas_pis.at(2).at(i).c_str());
			}
		}


		for (i = 0; i < matriz_areas_pis.size(); i++)
		{
			n1 = atoi(matriz_areas_pis.at(0).at(i).substr(3, 1).c_str());
			areac = atof(matriz_areas_pis.at(2).at(i).c_str());

			if (((n1 == 1) | (n1 == 7) | (n1 == 9)))
			{
				K_fluidos_pis.at(n1).at(1) = -areac * hmean_gas;
				///RADIACION NODOS CULATA
				Qrad_pis_med_nodo = Qrad_pis_med * areac / Atotal_pis_rad;
				matriz_temp.at(nodos_conv + ngases + ncilin + n1) = Qrad_pis_med_nodo;
			}

			if (((n1 == 1) | (n1 == 9)))
			{
				K_fluidos_pis.at(n1).at(1) = K_fluidos_pis.at(n1).at(1) - areac * hmean_gas;
				Qrad_pis_med_nodo = Qrad_pis_med * areac / Atotal_pis_rad;
				matriz_temp.at(nodos_conv + ngases + ncilin + n1) = matriz_temp.at(nodos_conv + ngases + ncilin + n1) + Qrad_pis_med_nodo;
			}

			if (n1 == 3)
				K_fluidos_pis.at(3).at(4) = -hmean_oil *areac;

			if (((n1 == 5) | (n1 == 6) | (n1 == 8) | (n1 == 10)))
			{
				K_fluidos_pis.at(n1).at(4) = -hoil_splash(_Inicializar.VarMed_N * 60, toil) * areac;
			}
		}
	}
	else
	{
		Atotal_pis_rad = _NODOS_PIS.at(1).A_S + _NODOS_PIS.at(7).A_S + _NODOS_PIS.at(9).A_S + _NODOS_PIS.at(1).A_LD + _NODOS_PIS.at(9).A_LI;

		for (n1 = 1; n1 <= npis; n1++) {
			if (((n1 == 1) | (n1 == 7) | (n1 == 9))) {
				areac = _NODOS_PIS.at(n1).A_S;
				K_fluidos_pis.at(n1).at(1) = -areac * hmean_gas;
				///RADIACION NODOS CULATA
				Qrad_pis_med_nodo = Qrad_pis_med * areac / Atotal_pis_rad;
				matriz_temp.at(nodos_conv + ngases + ncilin + n1) = Qrad_pis_med_nodo;
			}

			if (((n1 == 1) | (n1 == 9))) {
				if (n1 == 1) {
					areac = _NODOS_PIS.at(n1).A_LD;
				}
				else if (n1 == 9) {
					areac = _NODOS_PIS.at(n1).A_LI;
				}
				K_fluidos_pis.at(n1).at(1) = K_fluidos_pis.at(n1).at(1) - areac * hmean_gas;
				Qrad_pis_med_nodo = Qrad_pis_med * areac / Atotal_pis_rad;
				matriz_temp.at(nodos_conv + ngases + ncilin + n1) = matriz_temp.at(nodos_conv + ngases + ncilin + n1) + Qrad_pis_med_nodo;
			}
		}

		K_fluidos_pis.at(3).at(4) = -hmean_oil * (pow(Constantes::pi, 2)) * 2 * Dgal * _NODOS_PIS.at(3).ANCHURA;
		K_fluidos_pis.at(5).at(4) = -hoil_splash(_Inicializar.VarMed_N * 60, toil) * _NODOS_PIS.at(5).A_LD;
		K_fluidos_pis.at(6).at(4) = -hoil_splash(_Inicializar.VarMed_N * 60, toil) * _NODOS_PIS.at(6).A_I;
		K_fluidos_pis.at(8).at(4) = -hoil_splash(_Inicializar.VarMed_N * 60, toil) * _NODOS_PIS.at(8).A_I;
		K_fluidos_pis.at(10).at(4) = -hoil_splash(_Inicializar.VarMed_N * 60, toil) * _NODOS_PIS.at(10).A_I;
	}
	//

	//for (n1 = 1; n1 <= npis; n1++) {
	//	if (((n1 == 1) | (n1 == 9))) {
	//		if (n1 == 1) {
	//			areac = _NODOS_PIS.at(n1).A_LD;
	//		}
	//		else if (n1 == 9) {
	//			areac = _NODOS_PIS.at(n1).A_LI;
	//		}
	//		K_fluidos_pis.at(n1).at(1) = K_fluidos_pis.at(n1).at(1) - areac * hmean_gas;
	//		Qrad_pis_med_nodo = Qrad_pis_med * areac / Atotal_pis_rad;
	//		matriz_temp.at(nodos_conv + ngases + ncilin + n1) = matriz_temp.at(nodos_conv + ngases + ncilin + n1) + Qrad_pis_med_nodo;
	//	}
	//}

	//K_fluidos_pis.at(3).at(4) = -hmean_oil * (pow(Constantes::pi, 2)) * 2 * Dgal * _NODOS_PIS.at(3).ANCHURA;
	//K_fluidos_pis.at(5).at(4) = -hoil_splash(_Inicializar.VarMed.N * 60, toil) * _NODOS_PIS.at(5).A_LD;
	//K_fluidos_pis.at(6).at(4) = -hoil_splash(_Inicializar.VarMed.N * 60, toil) * _NODOS_PIS.at(6).A_I;
	//K_fluidos_pis.at(8).at(4) = -hoil_splash(_Inicializar.VarMed.N * 60, toil) * _NODOS_PIS.at(8).A_I;
	//K_fluidos_pis.at(10).at(4) = -hoil_splash(_Inicializar.VarMed.N * 60, toil) * _NODOS_PIS.at(10).A_I;


	///' 5: culata
	//calcula Qrad_cul medio
	double Qrad_cul_med;
	double Atotal_cul_rad;

	Atotal_cul_rad = 0;

	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		for (i = 0; i < matriz_areas_pis.size(); i++)
		{
			n1 = atoi(matriz_areas_pis.at(0).at(i).substr(3, 1).c_str());
			Atotal_cul_rad = Atotal_cul_rad + atof(matriz_areas_cul.at(2).at(i).c_str());
		}
	}
	else
	{
		for (n1 = 1; n1 <= ncul; n1++) {
			Atotal_cul_rad = Atotal_cul_rad + _NODOS_CUL.at(n1).AREA_GAS;
		}
	}

	Qrad_cul_med = 0;
	for (i = 1; i <= NPC; i++) {
		Qrad_cul_med = Qrad_cul_med + var_inst_nodal.at(i).Qrad_cul;
	}
	Qrad_cul_med = Qrad_cul_med / NPC;


	double Qrad_cul_med_nodo;


	if (_Inicializar.OpcionesGenerales_MODELONODALAVANZADO)
	{
		std::string nombre;
		for (i = 0; i < matriz_areas_cul.size(); i++)
		{
			n1 = atoi(matriz_areas_cul.at(0).at(i).substr(3, matriz_areas_cul.at(0).at(i).find_first_not_of("_")).c_str());
			nombre = matriz_areas_cul.at(0).at(i).c_str();
			areac = atof(matriz_areas_cul.at(2).at(i).c_str());

			K_fluidos_cul.at(n1).at(1) = -hmean_gas * areac;
			///RADIACION NODOS CULATA
			Qrad_cul_med_nodo = Qrad_cul_med * areac / Atotal_cul_rad;
			matriz_temp.at(nodos_conv + ngases + ncilin + npis + n1) = Qrad_cul_med_nodo;

			if (nombre.find("VLV") >= 0)
			{
				K_fluidos_cul.at(n1).at(2) = -f_h_vlv_adm * hmean_adm * areac;
				K_fluidos_cul.at(n1).at(3) = -f_h_vlv_esc * hmean_exh * areac;
			}
			else if (nombre.find("VAST") >= 0)
			{
				K_fluidos_cul.at(n1).at(2) = -f_h_vast_adm * hmean_adm * areac;
				K_fluidos_cul.at(n1).at(3) = -f_h_vast_esc * hmean_exh * areac;
			}
			else
			{
				K_fluidos_cul.at(n1).at(2) = -hmean_adm * areac;
				K_fluidos_cul.at(n1).at(3) = -hmean_exh * areac;
			}

			K_fluidos_cul.at(n1).at(6) = -_Inicializar.ConfMotor_CULATA_H_CUL * areac;
		}
	}
	else
	{

		for (n1 = 1; n1 <= ncul; n1++) {

			K_fluidos_cul.at(n1).at(1) = -hmean_gas * _NODOS_CUL.at(n1).AREA_GAS;
			///RADIACION NODOS CULATA
			Qrad_cul_med_nodo = Qrad_cul_med * _NODOS_CUL.at(n1).AREA_GAS / Atotal_cul_rad;
			matriz_temp.at(nodos_conv + ngases + ncilin + npis + n1) = Qrad_cul_med_nodo;
			int es_vlv = _NODOS_CUL.at(n1).NOMBRE.find("VLV");
			int es_vast = _NODOS_CUL.at(n1).NOMBRE.find("VAST");

			if ((es_vlv >= 0) && ((_NODOS_CUL.at(n1).AREA_ADM != 0) || (_NODOS_CUL.at(n1).AREA_ESC != 0)))
			{
				K_fluidos_cul.at(n1).at(2) = -f_h_vlv_adm * hmean_adm * _NODOS_CUL.at(n1).AREA_ADM;
				K_fluidos_cul.at(n1).at(3) = -f_h_vlv_esc * hmean_exh * _NODOS_CUL.at(n1).AREA_ESC;
			}
			else if (es_vast >= 0)
			{
				K_fluidos_cul.at(n1).at(2) = -f_h_vast_adm * hmean_adm * _NODOS_CUL.at(n1).AREA_ADM;
				K_fluidos_cul.at(n1).at(3) = -f_h_vast_esc * hmean_exh * _NODOS_CUL.at(n1).AREA_ESC;
			}
			else
			{
				K_fluidos_cul.at(n1).at(2) = -hmean_adm * _NODOS_CUL.at(n1).AREA_ADM;
				K_fluidos_cul.at(n1).at(3) = -hmean_exh * _NODOS_CUL.at(n1).AREA_ESC;
			}

			K_fluidos_cul.at(n1).at(6) = -_Inicializar.ConfMotor_CULATA_H_CUL * _NODOS_CUL.at(n1).AREA_COOL;
		}

	}

	//CONCATENA MATRICES
	//CILINDRO
	m_up1 = K_fluidos;
	m_inf2 = K_gas;

	concatena_4_matrices(&m_up1, &m_up2, &m_inf1, &m_inf2, &matriz_conductances2);

	//ImprimirMatriz(&matriz_conductances2, "matriz_conductances1");

	m_up3 = matriz_conductances2;
	matriz_conductances2.clear();

	concatena_matrices_left_right(&K_fluidos_cil, &K_gas_cil, &m_inf3);

	m_inf4 = k_cil;

	concatena_4_matrices(&m_up3, &m_up4, &m_inf3, &m_inf4, &matriz_conductances2);
	//ImprimirMatriz(&matriz_conductances2, "matriz_conductances2");

	//PISTON
	m_up5 = matriz_conductances2;
	matriz_conductances2.clear();

	concatena_matrices_up_down(&m_aux1, &k_cil_pis, &m_up6);

	concatena_matrices_left_right(&K_fluidos_pis, &m_aux2, &m_aux3);

	concatena_matrices_left_right(&m_aux3, &k_pis_cil, &m_inf5);

	m_inf6 = k_pis;

	concatena_4_matrices(&m_up5, &m_up6, &m_inf5, &m_inf6, &matriz_conductances2);
	//ImprimirMatriz(&matriz_conductances2, "matriz_conductances3");

	//CULATA

	m_up7 = matriz_conductances2;

	concatena_matrices_up_down(&m_aux5, &k_cil_cul, &m_aux6);
	concatena_matrices_up_down(&m_aux6, &k_pis_cul, &m_up8);


	concatena_matrices_left_right(&K_fluidos_cul, &m_aux7, &m_aux8);
	concatena_matrices_left_right(&m_aux8, &k_cul_cil, &m_aux9);
	concatena_matrices_left_right(&m_aux9, &k_cul_pis, &m_inf7);

	m_inf8 = k_cul;

	concatena_4_matrices(&m_up7, &m_up8, &m_inf7, &m_inf8, &matriz_conductances);
	//ImprimirMatriz(&m_up7, "m_up7");
	//ImprimirMatriz(&m_up8, "m_up8");
	//ImprimirMatriz(&m_inf7, "m_inf7");
	//ImprimirMatriz(&m_inf8, "m_inf8");

	//ImprimirMatriz(&matriz_conductances, "matriz_conductances4");

	//std::ofstream output("c:\\temp\\Matriz.txt");
	//for (int k = 1; k < matriz_conductances->size(); k++)
	//{
	//	for (int l = 1; l < matriz_conductances->at(0).size(); l++)
	//	{
	//		output << &matriz_conductances[k][l] << ";"; // behaves like cout - cout is also a stream
	//	}
	//	output << ";" << "\n";
	//}


	//*****************VECTOR TEMPERATURAS******************************

	matriz_temp.at(1) = Tmean_gas;
	matriz_temp.at(2) = T_Intake;
	matriz_temp.at(3) = T_Exhaust;
	matriz_temp.at(4) = toil;
	matriz_temp.at(5) = twater_cil;
	matriz_temp.at(6) = twater_cul;

	for (cont = 1; cont <= ngases; cont++) {
		matriz_temp.at(nodos_conv + cont) = matriz_temp_cil.at(cont);
	}

	//Faltaba calcular la suma de todas conductancias para la diagonal
	diagonal_suma(&matriz_conductances);
	ImprimirMatriz(&matriz_conductances, "matriz_conductances5");
	//If Not EscribirMatriz(matriz_conductances, e.TIPEN & "_matriz_conductances.txt") Then Beep

	n = matriz_conductances.size();

	std::vector<double> temperatura_nodos;

	//Resulve el sistema de ecuaciones
	std::vector< std::vector<double> > inversa(n, std::vector<double>(n));


	invNxN(&matriz_conductances, &inversa);
	mult_matriz(&inversa, &matriz_temp, &temperatura_nodos);

	std::vector< std::vector<double> > matriz_flujos(n, std::vector<double>(n));

	calcula_flujos(&matriz_conductances, &temperatura_nodos, &matriz_flujos);


	guardar_resultado_F_K(Z, "K", &matriz_conductances);
	guardar_resultado_F_K(Z, "F", &matriz_flujos);
	guardar_resultado_T(Z, &temperatura_nodos);

	ImprimirMatriz(&matriz_conductances, "matriz_conductances");
	ImprimirMatriz(&matriz_flujos, "matriz_flujos");

	ImprimirVector(&temperatura_nodos, "temperatura_nodos");


	calcular_TCUL(Z, TCUL, TCULMAX, TCULMIN);
	calcular_TCULMAT(Z, TCULMAT);
	calcular_TCULVALV(Z, TCULVALV);
	calcular_TPIS(Z, TPIS, TPISMAX, TPISMIN);

	calcular_TCIL(Z, TCIL, TCILMAX, TCILMIN);





}


void TransmisionCalor::calcular_TPIS(short Z, double *TPIS, double *TPISMAX, double *TPISMIN)
{
	int id_gas1;
	double K;
	double T;
	double suma_K;

	*TPISMAX = 0;
	*TPISMIN = 999999.999;
	suma_K = 0;

	*TPIS = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "GAS1")
		{
			id_gas1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}

	for (int i = 0; i < _NODOS_PIS.size(); i++)
	{
		K = 0;
		T = 0;

		for (int j = 0; j < _RESULTADOS_K.size(); j++)
		{
			if ((_NODOS_PIS[i].ID_NODO == _RESULTADOS_K[j].ID_NODO1) && (_RESULTADOS_K[j].ID_NODO2 == id_gas1) && (_RESULTADOS_K[j].Z == Z))
			{
				K = _RESULTADOS_K[j].valor;
				suma_K = suma_K + K;
				break;
			}
		}

		for (int j = 0; j < _RESULTADOS_T.size(); j++)
		{
			if ((_NODOS_PIS[i].ID_NODO == _RESULTADOS_T[j].ID_NODO1) && (_RESULTADOS_T[j].Z == Z))
			{
				T = _RESULTADOS_T[j].valor;
				if (T> *TPISMAX)
				{
					*TPISMAX = T;
				}

				if (T < *TPISMIN)
				{
					*TPISMIN = T;
				}
				break;
			}
		}



		*TPIS = *TPIS + (K * T);
	}

	*TPIS = *TPIS / suma_K;
}


void TransmisionCalor::calcular_QPIPA_ADM(short Z, double *QPIPA_ADM)
{
	int id_adm1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "ADM1")
		{
			id_adm1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}


	for (int j = 0; j < _RESULTADOS_F.size(); j++)
	{
		if ((_RESULTADOS_F[j].ID_NODO2 == id_adm1) && (_RESULTADOS_F[j].Z == Z))
		{
			K = _RESULTADOS_F[j].valor;
			suma_K = suma_K + K;
			break;
		}
	}


	*QPIPA_ADM = suma_K;
}

void TransmisionCalor::calcular_QAC(short Z, double *QAC)
{
	int id_oil1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "OIL1")
		{
			id_oil1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}


	for (int j = 0; j < _RESULTADOS_F.size(); j++)
	{
		if ((_RESULTADOS_F[j].ID_NODO2 == id_oil1) && (_RESULTADOS_F[j].Z == Z))
		{
			K = _RESULTADOS_F[j].valor;
			suma_K = suma_K + K;
			break;
		}
	}


	*QAC = suma_K;
}


void TransmisionCalor::calcular_QPIPA_ESC(short Z, double *QPIPA_ESC)
{
	int id_esc1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "ESC1")
		{
			id_esc1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}


	for (int j = 0; j < _RESULTADOS_F.size(); j++)
	{
		if ((_RESULTADOS_F[j].ID_NODO2 == id_esc1) && (_RESULTADOS_F[j].Z == Z))
		{
			K = _RESULTADOS_F[j].valor;
			suma_K = suma_K + K;
			break;
		}
	}

	*QPIPA_ESC = suma_K;
}

void TransmisionCalor::calcular_QREF(short Z, double *QREF)
{
	int id_cool_cil1;
	int id_cool_cul1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "COOL_CIL1")
		{
			id_cool_cil1 = _NODOS_FLUIDOS[i].ID_NODO;
		}
		if (_NODOS_FLUIDOS[i].nombre == "COOL_CUL1")
		{
			id_cool_cul1 = _NODOS_FLUIDOS[i].ID_NODO;
		}

	}


	for (int j = 0; j < _RESULTADOS_F.size(); j++)
	{
		if (((_RESULTADOS_F[j].ID_NODO2 == id_cool_cil1) || (_RESULTADOS_F[j].ID_NODO2 == id_cool_cul1)) && (_RESULTADOS_F[j].Z == Z))
		{
			K = _RESULTADOS_F[j].valor;
			suma_K = suma_K + K;
			break;
		}
	}

	*QREF = suma_K;
}



void TransmisionCalor::calcular_TCIL(short Z, double *TCIL, double *TCILMAX, double *TCILMIN)
{
	int id_gas1;
	double K;
	double T;
	double suma_K;
	double Tg;
	double hA;

	*TCILMAX = 0;
	*TCILMIN = 99999.99;
	suma_K = 0;

	*TCIL = 0;

	double GAS_CIL;
	GAS_CIL = 0;
	for (int i = 1; i < _NODOS_CIL.size(); i++)
	{
		for (int j = 1; j < _NODOS_FLUIDOS.size(); j++)
		{
			if (_NODOS_FLUIDOS.at(j).tipo_fluido == "GAS")
			{
				for (int k = 1; k < _RESULTADOS_F.size(); k++)
				{
					if ((_RESULTADOS_F[k].ID_NODO1 == _NODOS_CIL[i].ID_NODO) && (_RESULTADOS_F[k].ID_NODO2 == _NODOS_FLUIDOS[j].ID_NODO))
					{
						GAS_CIL = GAS_CIL + _RESULTADOS_F[k].valor;
					}
				}

			}
		}
	}

	for (int i = 0; i < _NODOS_CIL.size(); i++)
	{
		T = 0;

		for (int j = 0; j < _RESULTADOS_T.size(); j++)
		{
			if ((_NODOS_CIL[i].ID_NODO == _RESULTADOS_T[j].ID_NODO1) && (_RESULTADOS_T[j].Z == Z))
			{
				T = _RESULTADOS_T[j].valor;
				if (T> *TCILMAX)
				{
					*TCILMAX = T;
				}

				if (T < *TCILMIN)
				{
					*TCILMIN = T;
				}
				break;
			}
		}


	}

	Tg = calcular_Tg_media();
	hA = calcular_hA_media();

	*TCIL = Tg - (GAS_CIL / hA);
}

void TransmisionCalor::calcular_TCUL(short Z, double *TCUL, double *TCULMAX, double *TCULMIN)
{
	int id_gas1;
	double K;
	double suma_K;
	double T;


	*TCUL = 0;
	*TCULMAX = 0;
	*TCULMIN = 99999.99;
	suma_K = 0;

	for (int i = 1; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "GAS1")
		{
			id_gas1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}

	for (int i = 1; i < _NODOS_CUL.size(); i++)
	{
		K = 0;
		T = 0;

		for (int j = 1; j < _RESULTADOS_K.size(); j++)
		{
			if ((_NODOS_CUL[i].ID_NODO == _RESULTADOS_K[j].ID_NODO1) && (_RESULTADOS_K[j].ID_NODO2 == id_gas1) && (_RESULTADOS_K[j].Z == Z))
			{
				K = _RESULTADOS_K[j].valor;
				suma_K = suma_K + K;
				break;
			}
		}

		for (int j = 1; j < _RESULTADOS_T.size(); j++)
		{
			if ((_NODOS_CUL[i].ID_NODO == _RESULTADOS_T[j].ID_NODO1) && (_RESULTADOS_T[j].Z == Z))
			{
				T = _RESULTADOS_T[j].valor;
				if (T> *TCULMAX)
				{
					*TCULMAX = T;
				}

				if (T < *TCULMIN)
				{
					*TCULMIN = T;
				}
				break;
			}
		}

		*TCUL = *TCUL + (K * T);

	}
	*TCUL = *TCUL / suma_K;
}


void TransmisionCalor::calcular_TCULMAT(short Z, double *TCULMAT)
{
	int id_gas1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;

	*TCULMAT = 0;

	for (int i = 1; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "GAS1")
		{
			id_gas1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}

	for (int i = 1; i < _NODOS_CUL.size(); i++)
	{
		K = 0;
		T = 0;
		int esta = _NODOS_CUL[i].NOMBRE.find("VLV");

		for (int j = 1; j < _RESULTADOS_K.size(); j++)
		{
			if ((esta == -1) && (_NODOS_CUL[i].ID_NODO == _RESULTADOS_K[j].ID_NODO1) && (_RESULTADOS_K[j].ID_NODO2 == id_gas1) && (_RESULTADOS_K[j].Z == Z))
			{
				K = _RESULTADOS_K[j].valor;
				suma_K = suma_K + K;
				break;
			}
		}

		for (int j = 1; j < _RESULTADOS_T.size(); j++)
		{
			if ((esta == -1) && (_NODOS_CUL[i].ID_NODO == _RESULTADOS_T[j].ID_NODO1) && (_RESULTADOS_T[j].Z == Z))
			{
				T = _RESULTADOS_T[j].valor;
				*TCULMAT = *TCULMAT + (K * T);
				break;
			}
		}

		//*TCULMAT = *TCULMAT + (K * T);
	}

	*TCULMAT = *TCULMAT / suma_K;
}


void TransmisionCalor::calcular_TCULVALV(short Z, double *TCULVALV)
{
	int id_gas1;
	double K;
	double T;
	double suma_K;

	suma_K = 0;
	*TCULVALV = 0;

	for (int i = 0; i < _NODOS_FLUIDOS.size(); i++)
	{
		if (_NODOS_FLUIDOS[i].nombre == "GAS1")
		{
			id_gas1 = _NODOS_FLUIDOS[i].ID_NODO;
			break;
		}
	}

	for (int i = 0; i < _NODOS_CUL.size(); i++)
	{
		K = 0;
		T = 0;
		int esta = _NODOS_CUL[i].NOMBRE.find("VLV");

		for (int j = 0; j < _RESULTADOS_K.size(); j++)
		{

			if ((esta >= 0) && (_NODOS_CUL[i].ID_NODO == _RESULTADOS_K[j].ID_NODO1) && (_RESULTADOS_K[j].ID_NODO2 == id_gas1) && (_RESULTADOS_K[j].Z == Z))
			{
				K = _RESULTADOS_K[j].valor;
				suma_K = suma_K + K;
				break;
			}
		}

		for (int j = 0; j < _RESULTADOS_T.size(); j++)
		{
			if ((esta >= 0) && (_NODOS_CUL[i].ID_NODO == _RESULTADOS_T[j].ID_NODO1) && (_RESULTADOS_T[j].Z == Z))
			{
				T = _RESULTADOS_T[j].valor;
				*TCULVALV = *TCULVALV + (K * T);
				break;
			}
		}
		//*TCULVALV = *TCULVALV + (K * T);

	}

	*TCULVALV = *TCULVALV / suma_K;

}


void TransmisionCalor::GenerarMatricesDeFicheros(std::vector< std::vector<double> >  *k_cil,
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
	short *npis, short *ncul, short *naxiales, short *nradiales, short *ncircunf)
{
	int max_nodo_cil_circ = 0;
	int max_nodo_cil_rad = 0;
	int max_nodo_cil_ax = 0;
	int max_nodo_pis = 0;
	int max_nodo_cul = 0;

	int nodo_cil_circ = 0;
	int nodo_cil_rad = 0;
	int nodo_cil_ax = 0;
	int nodo_pis = 0;
	int nodo_cul = 0;

	int num_nodos = 0;
	int ncilin = 0;
	int ngases = 0;
	int nodos_conv = 0;
	bool primera_linea = true;
	int col = 0;

	std::string tipo_nodo = "";

	std::string nodo = "";
	std::vector<std::vector<double> > values;
	std::vector<std::string> lineatexto;
	std::vector<double> lineavalores;

	std::ifstream fichero1(_Inicializar.OpcionesGenerales_FICHERO_COND_COND_MNA);
	std::string item;

	if (_Inicializar.OpcionesGenerales_FICHERO_COND_COND_MNA != "")
	{

		for (std::string line; getline(fichero1, line);)
		{
			col = 0;

			std::istringstream in(line);

			while (getline(in, item, ';'))
			{
				if (col != 0)
				{
					if (primera_linea)
					{
						lineatexto.push_back(item.c_str());
					}
					else
						lineavalores.push_back(atof(item.c_str()));
				}
				col++;
			}

			if (primera_linea)
			{
				for (int i = 1; i < lineatexto.size(); i++)
				{
					nodo = lineatexto.at(i);
					tipo_nodo = nodo.substr(0, 3);

					if (tipo_nodo == "CIL") // es un nodo de cilindro
					{
						nodo_cil_ax = atoi(nodo.substr(3, 1).c_str()); //convierte a entero
						if (nodo_cil_ax > max_nodo_cil_ax)
							max_nodo_cil_ax = nodo_cil_ax;

						nodo_cil_rad = atoi(nodo.substr(4, 1).c_str()); //convierte a entero
						if (nodo_cil_rad > max_nodo_cil_rad)
							max_nodo_cil_rad = nodo_cil_rad;

						nodo_cil_circ = atoi(nodo.substr(4, 1).c_str()); //convierte a entero
						if (nodo_cil_circ > max_nodo_cil_circ)
							max_nodo_cil_circ = nodo_cil_circ;
					}
					else if (tipo_nodo == "PIS") // es un nodo de piston
					{
						nodo_pis = atoi(nodo.substr(3).c_str()); //convierte a entero
						if (nodo_pis > max_nodo_pis)
							max_nodo_pis = nodo_pis;
					}
					else if (tipo_nodo == "CUL") // es un nodo de culata
					{
						std::size_t pos = nodo.find("_");      // position of "live" in str

						nodo_cul = atoi(nodo.substr(3, pos - 3).c_str()); //convierte a entero
						if (nodo_cul > max_nodo_cul)
							max_nodo_cul = nodo_cul;
					}
				}
				primera_linea = false;
				*npis = max_nodo_pis;
				*ncul = max_nodo_cul;
				ncilin = max_nodo_cil_ax * max_nodo_cil_rad * max_nodo_cil_circ + max_nodo_cil_circ;
				ngases = max_nodo_cil_ax - 1;
				nodos_conv = 6; //(GAS, ESCAPE, ADMIS, OIL, AGUA_CIL, AGUA_CUL) Estos son Fijos;
				num_nodos = ncilin + max_nodo_pis + max_nodo_cul + ngases;
			}
			else
			{
				values.push_back(lineavalores);
				lineavalores.clear();
			}

		}

		k_cil->resize(ncilin + 1, std::vector<double>(ncilin + 1));
		k_pis->resize(*npis + 1, std::vector<double>(*npis + 1));
		k_cul->resize(*ncul + 1, std::vector<double>(*ncul + 1));
		k_cil_pis->resize(ncilin + 1, std::vector<double>(*npis + 1));
		k_cil_cul->resize(ncilin + 1, std::vector<double>(*ncul + 1));
		k_pis_cul->resize(*npis + 1, std::vector<double>(*ncul + 1));
		k_pis_cil->resize(*npis + 1, std::vector<double>(ncilin + 1));
		k_cul_cil->resize(*ncul + 1, std::vector<double>(ncilin + 1));
		k_cul_pis->resize(*ncul + 1, std::vector<double>(*npis + 1));

		for (int i = 1; i <= ncilin; i++)
		{
			for (int j = 1; j <= ncilin; j++)
			{
				k_cil->at(i).at(j) = values.at(i).at(j);
			}
		}

		for (int i = 1; i <= ncilin; i++)
		{
			for (int j = ncilin + 1; j <= ncilin + *npis; j++)
			{
				k_pis_cil->at(i).at(j) = values.at(i).at(j);
			}
		}
		for (int i = ncilin + 1; i <= ncilin + *npis; i++)
		{
			for (int j = 1; j <= ncilin; j++)
			{
				k_cil_pis->at(i).at(j) = values.at(i).at(j);
			}
		}

		for (int i = ncilin + 1; i <= ncilin + *npis; i++)
		{
			for (int j = ncilin + 1; j <= ncilin + *npis; j++)
			{
				k_pis->at(i).at(j) = values.at(i).at(j);
			}
		}
		for (int i = ncilin + *npis + 1; i <= ncilin + *npis + *ncul; i++)
		{
			for (int j = ncilin + *npis + 1; j <= ncilin + *npis + *ncul; j++)
			{
				k_cul->at(i).at(j) = values.at(i).at(j);
			}
		}
	}


	//GENERAR MATRIZ AREAS
	std::ifstream fichero2(_Inicializar.OpcionesGenerales_FICHERO_AREA_CONV_MNA);

	for (std::string line; getline(fichero2, line);)
	{
		col = 0;

		std::istringstream in(line);

		while (getline(in, item, ';'))
		{
			lineatexto.push_back(item.c_str());

			if (col == 0)
			{
				if (primera_linea)
				{

				}
				else
				{
					nodo = lineatexto.at(0);
					tipo_nodo = nodo.substr(0, 3);
				}
			}
			col++;
		}

		if (!primera_linea)
		{

			if (tipo_nodo == "CIL") // es un nodo de cilindro
			{
				matriz_areas_cil->push_back(lineatexto);
				lineatexto.clear();
			}
			else if (tipo_nodo == "PIS") // es un nodo de piston
			{
				matriz_areas_pis->push_back(lineatexto);
				lineatexto.clear();
			}
			else if (tipo_nodo == "CUL") // es un nodo de culata
			{
				matriz_areas_cul->push_back(lineatexto);
				lineatexto.clear();
			}

		}
		else
		{
			primera_linea = false;
		}

	}


	//GENERAR MATRIZ GEOMETRIA CILINDRO
	std::ifstream fichero3(_Inicializar.OpcionesGenerales_FICHERO_GEOM_CIL_MNA);

	primera_linea = true;

	for (std::string line; getline(fichero3, line);)
	{
		col = 0;

		std::istringstream in(line);

		while (getline(in, item, ';'))
		{
			lineatexto.push_back(item.c_str());
		}

		if (!primera_linea)
		{
			matriz_geom_cil->push_back(lineatexto);
			lineatexto.clear();
		}
		else
		{
			primera_linea = false;
		}
	}

}


void TransmisionCalor::GenerarMatrices(std::vector< std::vector<double> >  *k_cil,
	std::vector< std::vector<double> >  *k_pis,
	std::vector< std::vector<double> >  *k_cul,
	std::vector< std::vector<double> >  *k_cil_pis,
	std::vector< std::vector<double> >  *k_cil_cul,
	std::vector< std::vector<double> >  *k_pis_cul,
	std::vector< std::vector<double> >  *k_pis_cil,
	std::vector< std::vector<double> >  *k_cul_cil,
	std::vector< std::vector<double> >  *k_cul_pis, short *npis, short *ncul, short *naxiales, short *nradiales, short *ncircunf)

{
	short nodos_conv;
	short ngases;
	short ncilin;
	short num_nodos;
	std::vector<double> L_nodos;


	*naxiales = _Inicializar.ConfMotor_BASE_NODOS_AXIALES;
	*nradiales = _Inicializar.ConfMotor_BASE_NODOS_RADIALES;
	*ncircunf = _Inicializar.ConfMotor_BASE_NODOS_CIRCUNF;

	*npis = 10;
	*ncul = _NODOS_CUL.size() - 1;
	nodos_conv = 6;
	ngases = *naxiales - 1;
	ncilin = *naxiales * *nradiales * *ncircunf + *ncircunf;

	num_nodos = _NODOS_CIL.size() - 1 + _NODOS_CUL.size() - 1 + _NODOS_PIS.size() - 1 + _NODOS_FLUIDOS.size() - 1;

	L_nodos.resize(*naxiales + 1); //Aqui almaceno las longitudes de cada uno de los nodos axiales

	k_cil->resize(ncilin + 1, std::vector<double>(ncilin + 1));
	k_pis->resize(*npis + 1, std::vector<double>(*npis + 1));
	k_cul->resize(*ncul + 1, std::vector<double>(*ncul + 1));
	k_cil_pis->resize(ncilin + 1, std::vector<double>(*npis + 1));
	k_cil_cul->resize(ncilin + 1, std::vector<double>(*ncul + 1));
	k_pis_cul->resize(*npis + 1, std::vector<double>(*ncul + 1));
	k_pis_cil->resize(*npis + 1, std::vector<double>(ncilin + 1));
	k_cul_cil->resize(*ncul + 1, std::vector<double>(ncilin + 1));
	k_cul_pis->resize(*ncul + 1, std::vector<double>(*npis + 1));


	int num_nodo1;
	int num_nodo2;
	for (int i = 0; i < _CONDUCTANCIAS_CONDUCT.size(); i++)
	{
		num_nodo1 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO1 - (nodos_conv + ngases);
		num_nodo2 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO2 - (nodos_conv + ngases);

		if ((num_nodo1 > 0) && (num_nodo1 <= ncilin)) //nodo del cilindro
		{
			if ((num_nodo2 > 0) && (num_nodo2 <= ncilin)) //relacionado con nodo del cilindro
			{
				k_cil->at(num_nodo1).at(num_nodo2) = _CONDUCTANCIAS_CONDUCT.at(i).k;
			}
			if (((num_nodo2 - ncilin) > 0) && ((num_nodo2 - ncilin) <= *npis)) //relacionado con nodo el piston
			{

				k_cil_pis->at(num_nodo1).at(num_nodo2 - ncilin) = _CONDUCTANCIAS_CONDUCT.at(i).k;
				k_pis_cil->at(num_nodo2 - ncilin).at(num_nodo1) = _CONDUCTANCIAS_CONDUCT.at(i).k;
			}
			//if (((num_nodo2 - ncilin - *npis) > 0) && ((num_nodo2 - ncilin - *npis) <= *ncul)) //relacionado con nodos culata
			//{

			//	k_cil_cul->at(num_nodo1).at(num_nodo2 - ncilin - *npis) = _CONDUCTANCIAS_CONDUCT.at(i).k;
			//	k_cul_cil->at(num_nodo2 - ncilin - *npis).at(num_nodo1) = _CONDUCTANCIAS_CONDUCT.at(i).k;
			//}

		}
		else
		{
			num_nodo1 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO1 - (nodos_conv + ngases + ncilin);
			num_nodo2 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO2 - (nodos_conv + ngases + ncilin);

			if ((num_nodo1 > 0) && (num_nodo1 <= *npis)) //nodo del piston
			{
				if ((num_nodo2 > 0) && (num_nodo2 <= *npis)) //relacionado con nodo del piston
				{
					k_pis->at(num_nodo1).at(num_nodo2) = _CONDUCTANCIAS_CONDUCT.at(i).k;
				}
				//if (((num_nodo2 + ncilin) > 0) && ((num_nodo2 + ncilin) <= ncilin)) //relacionado con nodo del cilindro
				//{

				//	k_pis_cil->at(num_nodo1).at(num_nodo2 + ncilin) = _CONDUCTANCIAS_CONDUCT.at(i).k;
				//}
				//if (((num_nodo2 + ncilin - *npis) > 0) && ((num_nodo2 + ncilin - *npis) <= *ncul)) //relacionado con nodos culata
				//{

				//	k_pis_cul->at(num_nodo1).at(num_nodo2 + ncilin - *npis) = _CONDUCTANCIAS_CONDUCT.at(i).k;
				//	k_cul_pis->at(num_nodo2 + ncilin - *npis).at(num_nodo1) = _CONDUCTANCIAS_CONDUCT.at(i).k;
				//}
			}
			else
			{
				num_nodo1 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO1 - (nodos_conv + ngases + ncilin + *npis);
				num_nodo2 = _CONDUCTANCIAS_CONDUCT.at(i).ID_NODO2 - (nodos_conv + ngases + ncilin + *npis);

				if ((num_nodo1 > 0) && (num_nodo1 <= *ncul)) //nodo de la culata
				{
					if ((num_nodo2 > 0) && (num_nodo2 <= *ncul)) //relacionado con nodo culata
					{
						k_cul->at(num_nodo1).at(num_nodo2) = _CONDUCTANCIAS_CONDUCT.at(i).k;
					}
				}
			}

		}

	}


}

void TransmisionCalor::Promedia_h(double *hmean_adm, double *hmean_exh, double *hmean_gas)
{
	//Función que se utiliza para sacar la media de las h
	short filas;
	short cont;
	double hmean;
	double hmean_adm_aux;
	double hmean_exh_aux;
	double hmean_gas_aux;

	filas = _Inicializar.ConfInstrumentacion_NPCL;
	for (cont = 1; cont <= filas; cont++) {

		hmean_adm_aux = hmean_adm_aux + var_inst_nodal.at(cont).h_adm;
		hmean_exh_aux = hmean_exh_aux + var_inst_nodal.at(cont).h_esc;
		hmean_gas_aux = hmean_gas_aux + var_inst_nodal.at(cont).h_gas;
	}

	*hmean_adm = hmean_adm_aux / (filas - 1);
	*hmean_exh = hmean_exh_aux / (filas - 1);
	*hmean_gas = hmean_gas_aux / (filas - 1);
}


void TransmisionCalor::tgas(double *hmean_gas, double *Tmean_gas, double *TGASM)
{

	short filas;
	short cont;
	bool k;

	*Tmean_gas = 0;
	filas = _Inicializar.ConfInstrumentacion_NPCL;

	for (cont = 1; cont <= filas; cont++) {
		*Tmean_gas = *Tmean_gas + var_inst_nodal.at(cont).h_gas * var_inst_nodal.at(cont).T;
	}

	*Tmean_gas = *Tmean_gas / *hmean_gas / (filas - 1);

	*TGASM = *Tmean_gas;
}



double TransmisionCalor::rho_oil(double temp_oil)
{

	return -16572.74497 - 38.71871842 * temp_oil + 0.054641551 * pow(temp_oil, 2) - 3.45327E-05 * pow(temp_oil, 3) + 4398.316892 * log(temp_oil);

}

double TransmisionCalor::mu_oil(double temp_oil)
{
	return 7.48E-05 * exp(1005.2 / (temp_oil - 157.45));
}


double TransmisionCalor::cp_oil(double temp_oil)
{

	return -34040.12013 - 71.24978467 * temp_oil + 0.107361063 * pow(temp_oil, 2) - 6.63966E-05 * pow(temp_oil, 3) + 8670.542844 * log(temp_oil);

}

double TransmisionCalor::k_oil(double temp_oil)
{

	return -13.78897843 - 0.028315003 * temp_oil + 3.80794E-05 * pow(temp_oil, 2) - 2.25321E-08 * pow(temp_oil, 3) + 3.437962885 * log(temp_oil);

}

void TransmisionCalor::calculo_vcil(double Lcil, std::vector< std::vector<double> > * area_inst_nodos_cil)
{
	//area_inst_nodos_cil() dato de salida
	short n;
	short m;
	short cont1;
	short cont2;


	n = _Inicializar.ConfInstrumentacion_NPCL;
	m = _NODOS_CIL.size();
	//UPGRADE_WARNING: El límite inferior de la matriz spiston ha cambiado de 1 a 0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	std::vector<double> spiston(n);
	//UPGRADE_WARNING: El límite inferior de la matriz areacil ha cambiado de 1 a 0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	std::vector<double>areacil(n);

	//ReDim areacil(1 To n, 1 To m) As Double
	//UPGRADE_WARNING: El límite inferior de la matriz area_inst_nodos_cil ha cambiado de 1,1 a 0,0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	// ERROR: Not supported in C#: ReDimStatement

	std::vector< std::vector<double> > area_inst_nodos_cil_aux(n + 1, std::vector<double>(m));


	for (cont1 = 1; cont1 <= n; cont1++) {
		for (cont2 = 1; cont2 < m; cont2++) {
			if (_NODOS_CIL.at(cont2).POS_RAD == 1) {
				if ((Lcil - _NODOS_CIL.at(cont2).L_SUP) > var_inst_nodal.at(cont1).hi) {
					area_inst_nodos_cil_aux.at(cont1).at(cont2) = 0;
				}
				else if ((Lcil - _NODOS_CIL.at(cont2).L_INF) <= var_inst_nodal.at(cont1).hi) {
					area_inst_nodos_cil_aux.at(cont1).at(cont2) = (_NODOS_CIL.at(cont2).L_SUP - _NODOS_CIL.at(cont2).L_INF) * _NODOS_CIL.at(cont2).D_MIN * Constantes::pi * (_NODOS_CIL.at(cont2).ANG_FIN - _NODOS_CIL.at(cont2).ANG_INI) / 360;
				}
				else {
					area_inst_nodos_cil_aux.at(cont1).at(cont2) = (var_inst_nodal.at(cont1).hi - Lcil + _NODOS_CIL.at(cont2).L_SUP) * _NODOS_CIL.at(cont2).D_MIN * Constantes::pi * (_NODOS_CIL.at(cont2).ANG_FIN - _NODOS_CIL.at(cont2).ANG_INI) / 360;
				}

			}
			else {
				area_inst_nodos_cil_aux.at(cont1).at(cont2) = 0;
			}

		}

	}

	*area_inst_nodos_cil = area_inst_nodos_cil_aux;
}

void TransmisionCalor::calculo_vcil_avanzado(double Lcil, std::vector< std::vector<double> > * area_inst_nodos_cil, std::vector< std::vector<std::string> > * matriz_geom_cil, short naxiales, short nradiales, short ncircunf)
{
	//Este método lee la geometía de la matriz que viene en el fichero para modelo nodal avanzado
	//area_inst_nodos_cil() dato de salida
	short n;
	short m;
	short cont1;
	short cont2;
	double L_SUP = 0;
	double L_INF = 0;
	double ANG_INI = 0;
	double ANG_FIN = 0;
	double D_MIN = 0;

	n = _Inicializar.ConfInstrumentacion_NPCL;
	m = naxiales * nradiales * ncircunf + ncircunf;
	std::vector<double> spiston(n);
	std::vector<double>areacil(n);

	std::vector< std::vector<double> > area_inst_nodos_cil_aux(n, std::vector<double>(m));


	for (cont1 = 1; cont1 <= n; cont1++) {
		for (int i = 0; i < matriz_geom_cil->size(); i++)
		{
			cont2 = ObtenerPosicionNodoCilindro(matriz_geom_cil->at(0).at(i).c_str(), naxiales, nradiales, ncircunf);

			L_SUP = atof(matriz_geom_cil->at(2).at(i).c_str());
			L_INF = atof(matriz_geom_cil->at(3).at(i).c_str());
			ANG_INI = atof(matriz_geom_cil->at(4).at(i).c_str());
			ANG_FIN = atof(matriz_geom_cil->at(5).at(i).c_str());
			D_MIN = _Inicializar.ConfMotor_BASE_D;

			if ((Lcil - L_SUP) > var_inst_nodal.at(cont1).hi) {
				area_inst_nodos_cil_aux.at(cont1).at(cont2) = 0;
			}
			else if ((Lcil - L_INF) <= var_inst_nodal.at(cont1).hi) {
				area_inst_nodos_cil_aux.at(cont1).at(cont2) = (L_SUP - L_INF) * D_MIN * Constantes::pi * (ANG_FIN - ANG_INI) / 360;
			}
			else {
				area_inst_nodos_cil_aux.at(cont1).at(cont2) = (var_inst_nodal.at(cont1).hi - Lcil + L_SUP) * D_MIN * Constantes::pi * (ANG_FIN - ANG_INI) / 360;
			}

		}

	}

	*area_inst_nodos_cil = area_inst_nodos_cil_aux;
}

double TransmisionCalor::ha_mean(std::vector<double>* A)
{
	///'''''''''''''''''''''''''''''''''
	//Calcula la conductancia convectiva'
	//(hA) media a partir de h, área y'''
	//temperaturas instantáneas '''''''''

	short n;
	short cont;
	double resultado = 0;

	n = _Inicializar.ConfInstrumentacion_NPCL;
	std::vector<double> h(n + 1);

	for (int k = 1; k <= n; k++) {
		h.at(k) = var_inst_nodal.at(k).h_gas;
	}

	resultado = mult_vector(&h, A) / n;

	return resultado;


}

double TransmisionCalor::mult_vector(std::vector<double>* m1, std::vector<double>* m2)
{

	short cont_f;
	short nnodos;
	short cont_c;
	double A;

	nnodos = m1->size();

	//UPGRADE_WARNING: El límite inferior de la matriz b ha cambiado de 1 a 0. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="0F1C9BE1-AF9D-476E-83B1-17D43BECFF20"'
	std::vector<double> b(nnodos);
	A = 0;
	for (cont_f = 1; cont_f < nnodos; cont_f++) {

		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto m2(cont_f). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto m1(). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		A = A + m1->at(cont_f) * m2->at(cont_f);

	}

	//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto mult_vector. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
	return A;

}



double TransmisionCalor::T_mean(std::vector<double>* A)
{
	///'''''''''''''''''''''''''''''''''
	//Calcula la temperatura media ''''''
	//a partir de h, área y temperaturas'
	//instantáneas ''''''''''''''''''''''

	short n;
	double hA;
	double hAT;
	double result;

	n = _Inicializar.ConfInstrumentacion_NPCL;

	std::vector<double> h(n + 1);
	std::vector<double> T(n + 1);


	for (int k = 1; k <= n; k++) {
		h.at(k) = var_inst_nodal.at(k).h_gas;
		T.at(k) = var_inst_nodal.at(k).T;
	}

	hAT = 0.0;
	hA = 0.0;
	for (int cont = 1; cont <= n; cont++) {
		hAT = h.at(cont) * A->at(cont) * T.at(cont) + hAT;
		hA = hA + h.at(cont) * A->at(cont);
	}


	if (hA != 0.0) {
		result = hAT / hA;
	}
	else {
		result = 0.0;
	}

	return result;
}

double TransmisionCalor::Qrad_mean(std::vector<double>* A_nodo_inst)
{
	//Calcula la Qrad media para un nodo del cilindro ''''''

	short n;
	double qAA = 0;

	n = _Inicializar.ConfInstrumentacion_NPCL;

	std::vector<double> QRAD(n + 1);
	std::vector<double> A_cil(n + 1);

	for (int cont = 1; cont <= n; cont++) {
		QRAD.at(cont) = var_inst_nodal.at(cont).Qrad_cil;
		A_cil.at(cont) = var_inst_nodal.at(cont).area_cil;
	}

	for (int cont = 1; cont <= n; cont++) {
		qAA = QRAD.at(cont) * A_nodo_inst->at(cont) / A_cil.at(cont) + qAA;
	}

	return qAA / n;
}


double TransmisionCalor::hoil_splash(double n, double toil_K)
{

	//Correlacion de Kaplan y Bohac

	const double tref = (273.15 + 90);
	const double href = 350;

	double vref;
	double rho_ref;
	double mu_ref;

	double voil;
	double rho_oil_;
	double mu_oil_;

	rho_ref = rho_oil(tref);
	mu_ref = mu_oil(tref);
	vref = mu_ref / rho_ref;

	rho_oil_ = rho_oil(toil_K);
	mu_oil_ = mu_oil(toil_K);
	voil = mu_oil_ / rho_oil_;

	double pot1 = (vref / voil);
	double pot2 = 0.5;
	double resul_pows = pow(pot1, pot2);
	double result = href * (n / 2000) * resul_pows;

	return result;

}


void TransmisionCalor::diagonal_suma(std::vector< std::vector<double> > *matriz)
{

	double suma;
	short fil;
	short m;
	short n;
	short col;

	m = matriz->size() - 1;
	n = m;
	for (fil = 1; fil <= m; fil++) {
		suma = 0;
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto matriz(fil, fil). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		if (matriz->at(fil).at(fil) != 1)
		{

			for (col = 1; col <= n; col++)
			{
				if ((fil != col))
				{
					//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto matriz(). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
					suma = suma + matriz->at(fil).at(col);
				}
			}
		}
		else
		{
			//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto matriz(). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
			suma = -matriz->at(fil).at(fil);
		}
		//UPGRADE_WARNING: No se puede resolver la propiedad predeterminada del objeto matriz(). Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="6A50421D-15FE-4896-8A1B-2EC21E9037B2"'
		matriz->at(fil).at(fil) = -suma;
	}

}


void TransmisionCalor::invNxN(std::vector< std::vector<double> > *x1, std::vector< std::vector<double> > *matriz_salida)
{
	short m;
	short n;
	short i;
	short j;
	short k;
	double c;
	double valor;

	m = x1->size() - 1;
	n = m;

	std::vector< std::vector<double> > Y(m + 1, std::vector<double>(n + 1));
	std::vector< std::vector<double> > x(m + 1, std::vector<double>(n + 1));
	std::vector< double>  p(n + 1);

	for (i = 1; i <= m; i++) {
		for (j = 1; j <= m; j++) {
			Y.at(i).at(j) = 0;
		}
	}

	complete_triangular_matrix(x1, &x);
	//ImprimirMatriz(x1, "x1");
	//ImprimirMatriz(&x, "x");

	pivot_triangular_matrix(x1, &p);

	for (k = 1; k <= m; k++) {
		c = p.at(k);
		Y.at(k).at(c) = 1;

		for (j = k; j <= n; j++) {
			if (Y.at(j).at(c) != 0) {

				for (i = j + 1; i <= n; i++) {
					valor = Y.at(i).at(c) - Y.at(j).at(c) * x.at(i).at(j);
					Y.at(i).at(c) = valor;
				}
			}
		}

		for (j = n; j >= 1; j += -1) {
			valor = Y.at(j).at(c) / x.at(j).at(j);
			Y.at(j).at(c) = valor;
			if (Y.at(j).at(c) != 0) {
				for (i = 1; i <= j - 1; i++) {
					valor = Y.at(i).at(c) - Y.at(j).at(c) * x.at(i).at(j);
					Y.at(i).at(c) = valor;
				}
			}
		}
	}

	*matriz_salida = Y;
}


void TransmisionCalor::complete_triangular_matrix(std::vector< std::vector<double> >* A, std::vector< std::vector<double> >* matriz_salida)
{
	short m;
	short n;
	short cont;
	short k;
	short i;
	short A0;
	double ptemp;
	double Temp;
	double maxval;
	short j;
	short p;

	n = A->size() - 1;
	m = n;

	std::vector<double> Pivot(m + 1);
	std::vector< std::vector<double> > b(m + 1, std::vector<double>(m + 1));

	b = *A;

	for (cont = 1; cont <= m; cont++) {
		Pivot.at(cont) = cont;
	}

	for (k = 1; k <= m - 1; k++) {
		p = k;
		maxval = abs(b.at(k).at(k));
		for (i = k + 1; i <= m; i++) {
			if ((abs(b.at(i).at(k)) > maxval)) {
				maxval = abs(b.at(i).at(k));
				p = i;
			}
		}
		if ((p != k)) {
			for (A0 = 1; A0 <= n; A0++) {
				Temp = b.at(p).at(A0);
				b.at(p).at(A0) = b.at(k).at(A0);
				b.at(k).at(A0) = Temp;
			}

			ptemp = Pivot.at(k);
			Pivot.at(k) = Pivot.at(p);
			Pivot.at(p) = ptemp;
		}

		if (b.at(k).at(k) != 0) {
			for (i = k + 1; i <= m; i++) {
				b.at(i).at(k) = b.at(i).at(k) / b.at(k).at(k);
			}

			for (j = k + 1; j <= n; j++) {
				for (i = k + 1; i <= m; i++) {
					b.at(i).at(j) = b.at(i).at(j) - b.at(i).at(k) * b.at(k).at(j);
				}
			}
		}
	}

	*matriz_salida = b;

}


void TransmisionCalor::pivot_triangular_matrix(std::vector< std::vector<double> >* A, std::vector<double>* matriz_salida)
{
	short m;
	short n;
	short p;

	short cont;
	short k;
	short i;
	short A0;
	double Temp;

	double ptemp;
	double maxval;
	short j;

	n = A->size() - 1;
	m = n;

	std::vector<double> Pivot(m + 1);
	std::vector< std::vector<double> > b(m + 1, std::vector<double>(m + 1));


	b = *A;

	for (cont = 1; cont <= m; cont++) {
		Pivot.at(cont) = cont;
	}

	for (k = 1; k <= m - 1; k++) {
		p = k;
		maxval = abs(b.at(k).at(k));
		for (i = k + 1; i <= m; i++) {
			if ((abs(b.at(i).at(k)) > maxval)) {
				maxval = abs(b.at(i).at(k));
				p = i;
			}
		}

		if ((p != k)) {
			for (A0 = 1; A0 <= n; A0++) {
				Temp = b.at(p).at(A0);
				b.at(p).at(A0) = b.at(k).at(A0);
				b.at(k).at(A0) = Temp;
			}

			ptemp = Pivot.at(k);
			Pivot.at(k) = Pivot.at(p);
			Pivot.at(p) = ptemp;
		}

		if (b.at(k).at(k) != 0) {
			for (i = k + 1; i <= m; i++) {
				b.at(i).at(k) = b.at(i).at(k) / b.at(k).at(k);
			}

			for (j = k + 1; j <= n; j++) {
				for (i = k + 1; i <= m; i++) {
					b.at(i).at(j) = b.at(i).at(j) - b.at(i).at(k) * b.at(k).at(j);
				}
			}
		}
	}

	*matriz_salida = Pivot;

}


void TransmisionCalor::mult_matriz(std::vector< std::vector<double> >* m1, std::vector<double>* m2, std::vector<double>* matriz_salida)
{
	short cont_f;
	short nnodos;
	short cont_c;
	double A;

	nnodos = m1->size() - 1;

	matriz_salida->resize(nnodos + 1);

	for (cont_f = 1; cont_f <= nnodos; cont_f++) {
		A = 0;
		for (cont_c = 1; cont_c <= nnodos; cont_c++) {
			A = A + m1->at(cont_f).at(cont_c) * m2->at(cont_c);
		}
		matriz_salida->at(cont_f) = A; // _comun->round_doub(A, 2);
	}
}


void TransmisionCalor::calcula_flujos(std::vector< std::vector<double> >* matriz_conductances, std::vector<double>* temperatura_nodos, std::vector< std::vector<double> >* mat_flux)
{

	short nodos_totales;
	short nnodo1;
	short nnodo2;
	bool k;

	nodos_totales = temperatura_nodos->size() - 1;

	for (nnodo1 = 1; nnodo1 <= nodos_totales; nnodo1++) {
		for (nnodo2 = 1; nnodo2 <= nodos_totales; nnodo2++) {
			mat_flux->at(nnodo1).at(nnodo2) = matriz_conductances->at(nnodo1).at(nnodo2) * (temperatura_nodos->at(nnodo1) - temperatura_nodos->at(nnodo2));
		}
	}

}

void TransmisionCalor::guardar_resultado_F_K(short num_cil, std::string tipo_calculo, std::vector< std::vector<double> >* matriz)
{
	//tipo_calculo es F para flujos y K para conductancias
	short i;
	short j;
	short num_nodos;

	num_nodos = matriz->size() - 1;

	if (tipo_calculo == "F") {
		stResultadosNodalF reg;
		for (i = 1; i <= num_nodos; i++) {
			for (j = 1; j <= num_nodos; j++) {
				if (matriz->at(i).at(j) != 0) {
					reg.Z = num_cil;
					reg.ID_NODO1 = i;
					reg.ID_NODO2 = j;
					reg.valor = matriz->at(i).at(j);
					_RESULTADOS_F.push_back(reg);
				}
			}
		}
	}
	else if (tipo_calculo == "K") {
		stResultadosNodalK reg;
		for (i = 1; i <= num_nodos; i++) {
			for (j = 1; j <= num_nodos; j++) {
				if (matriz->at(i).at(j) != 0) {
					reg.Z = num_cil;
					reg.ID_NODO1 = i;
					reg.ID_NODO2 = j;
					reg.valor = matriz->at(i).at(j);
					_RESULTADOS_K.push_back(reg);
				}
			}
		}
	}

}

void TransmisionCalor::guardar_resultado_T(short num_cil, std::vector<double>* matriz)
{
	//tipo_calculo es F para flujos y K para conductancias
	short i;
	short j;
	short num_nodos;

	num_nodos = matriz->size() - 1;

	stResultadosNodalT reg;
	for (i = 1; i <= num_nodos; i++) {
		if (matriz->at(i) != 0) {
			reg.Z = num_cil;
			reg.ID_NODO1 = i;
			reg.valor = matriz->at(i);
			_RESULTADOS_T.push_back(reg);
		}
	}

}

double TransmisionCalor::calcular_Tg_media()
{
	short n;
	short i;
	double Total;
	double hAT;
	double hA;
	n = _Inicializar.ConfInstrumentacion_NPCL;
	hAT = 0;
	hA = 0;
	for (i = 1; i <= n; i++) {
		hAT = hAT + (var_inst_nodal.at(i).h_gas * var_inst_nodal.at(i).area_cil * var_inst_nodal.at(i).T);
		hA = hA + (var_inst_nodal.at(i).h_gas * var_inst_nodal.at(i).area_cil);
	}
	return hAT / hA;
}


double TransmisionCalor::calcular_hA_media()
{
	short n;
	short i;
	double Total;
	n = _Inicializar.ConfInstrumentacion_NPCL;
	Total = 0;
	for (i = 1; i <= n; i++) {
		Total = Total + (var_inst_nodal.at(i).h_gas * var_inst_nodal.at(i).area_cil);
	}

	return Total / n;

}

void TransmisionCalor::concatena_4_matrices(std::vector< std::vector<double> >* SI, std::vector< std::vector<double> >* SD, std::vector< std::vector<double> >* II, std::vector< std::vector<double> >* Id, std::vector< std::vector<double> >*matriz)
{

	short fil1;
	short col1;
	short col2;
	short fil2;
	short fil;

	fil1 = SI->size() - 1;
	col1 = SI->at(0).size() - 1;

	fil2 = II->size() - 1;
	col2 = Id->at(0).size() - 1;
	fil = fil1 + fil2;

	std::vector< std::vector<double> > mleft(fil + 1, std::vector<double>(col1 + 1));
	std::vector< std::vector<double> > mright(fil + 1, std::vector<double>(col2 + 1));

	ImprimirMatriz(&mleft, "mleft");
	concatena_matrices_up_down(SI, II, &mleft);
	ImprimirMatriz(&mleft, "mleft2");
	ImprimirMatriz(&mright, "mright");
	concatena_matrices_up_down(SD, Id, &mright);
	ImprimirMatriz(&mright, "mright2");
	concatena_matrices_left_right(&mleft, &mright, matriz);
	ImprimirMatriz(matriz, "matrizXXXX");
}

void TransmisionCalor::concatena_matrices_up_down(std::vector< std::vector<double> >* m1, std::vector< std::vector<double> >* m2, std::vector< std::vector<double> >*matriz)
{

	short f1;
	short f2;
	short f;
	short c;
	short cont1;
	short cont2;
	short ndim;

	c = m1->at(0).size() - 1;

	/*	c = 1;
		if (ndim != 1)
		c = m1->size() - 1;*/

	f1 = m1->size() - 1;
	f2 = m2->size() - 1;

	f = f1 + f2;

	matriz->resize(f + 1, std::vector<double>(c + 1));

	for (cont1 = 1; cont1 <= f; cont1++) {
		for (cont2 = 1; cont2 <= c; cont2++) {
			if (cont1 <= f1) {
				matriz->at(cont1).at(cont2) = m1->at(cont1).at(cont2);
			}
			else {
				matriz->at(cont1).at(cont2) = m2->at(cont1 - f1).at(cont2);
			}
		}
	}

}


void TransmisionCalor::concatena_matrices_left_right(std::vector< std::vector<double> >* m1, std::vector< std::vector<double> >* m2, std::vector< std::vector<double> >*matriz)
{
	short f;
	short c1;
	short c2;
	short c;
	short cont1;
	short cont2;

	f = m1->size() - 1;
	c1 = m1->at(0).size() - 1;
	c2 = m2->at(0).size() - 1;

	c = c1 + c2;

	matriz->resize(f + 1, std::vector<double>(c + 1));

	for (cont1 = 1; cont1 <= f; cont1++) {
		for (cont2 = 1; cont2 <= c; cont2++) {
			if (cont2 <= c1) {
				matriz->at(cont1).at(cont2) = m1->at(cont1).at(cont2);
			}
			else {
				matriz->at(cont1).at(cont2) = m2->at(cont1).at(cont2 - c1);
			}
		}
	}

}

int TransmisionCalor::ObtenerPosicionNodoCilindro(std::string nodo, int naxiales, int nradiales, int ncircunf)
{
	std::string tipo_nodo;
	short na;
	short nr;
	short nc;
	int pos;

	tipo_nodo = nodo.substr(0, 3);

	if (tipo_nodo == "CIL") // es un nodo de cilindro
	{
		na = atoi(nodo.substr(3, 1).c_str());
		nr = atoi(nodo.substr(4, 1).c_str());
		nc = atoi(nodo.substr(5, 1).c_str());

		if (nr <= nradiales)
		{
			pos = ((na - 1)* nradiales * ncircunf) + ((nr - 1)* ncircunf) + nc;

		}
		else
		{
			pos = ((nr - 1)* naxiales * ncircunf) + nc;
		}
	}
	else
		pos = -1;

	return pos;
}


//******************************************************************
// Cálculo del calor transmitido por Woschni
//******************************************************************
//
//
// Llamado desde CCiclo.Ejecutar_Ciclo, desde CAjusteArrastre.CalculaExpol
// y desde CAjusteCombustion.Ajuste_Presion_Pmi
//
void TransmisionCalor::Calcula_Calor_Woschni(double T, double p)
{
	//double c1;
	double x;
	double A;
	double b;
	double c;
	double factor;
	double cu_prima;
	double k_nodis;
	double ratio_ctms;
	double d_piston;
	double radio;
	double parr;
	double COMB;
	double c1;

	try
	{
		if (_Inicializar.TipoEnsayo == "C") {
			parr = _InicializarAlRCA.PCA * pow((_InicializarAlRCA.VCA / _InicializarPorAngulo.VCIL), 1.36);
			COMB = ((_Inicializar.ConfMotor_BASE_VD * _InicializarAlRCA.TCA) / (_InicializarAlRCA.PCA * _InicializarAlRCA.VCA)) * (p - parr);
			if (COMB < 0)
				COMB = 0;
		}

		if (_InicializarPorAngulo.ciclo_cerrado) {

			if (_Inicializar.ConfMotor_MOTOR_TIPO_ENCENDIDO == "MEC") {
				cu_prima = pow((_Inicializar.ConfMotor_BASE_D / _Inicializar.ConfMotor_PISTON_DB), 2) * (_Inicializar.ConfMotor_PISTON_DB / 2) * (2 * Constantes::pi * _Inicializar.VarMed_N * _Inicializar.ConfMotor_CULATA_CTM);
				k_nodis = exp(-0.200679 * pow(_Inicializar.ConfMotor_CULATA_CTM, 0.431262));
				if (k_nodis > 1)
					k_nodis = 1;
				ratio_ctms = pow((_Inicializar.ConfMotor_PISTON_DB / _Inicializar.ConfMotor_BASE_D), 2) / k_nodis;
				x = _comun->Xp(_InicializarPorAngulo.ang, ratio_ctms);
				_CU = k_nodis * cu_prima * x;
			}
			else {
				_CU = 0;
			}

			A = _Inicializar.Constante_TRANS_CALOR_C_A_WOSCHNI;
			b = _Inicializar.Constante_TRANS_CALOR_C_B_WOSCHNI;
			c = _Inicializar.Constante_TRANS_CALOR_C_C_WOSCHNI;

			factor = A + b * pow(_Inicializar.DatosCalcGeneral_CM, c);
		}
		else {
			if (_Inicializar.ConfMotor_MOTOR_TIPO_ENCENDIDO == "MEC") {
				d_piston = _comun->Dist_Despl_inst(_InicializarPorAngulo.ang, _Inicializar.ConfMotor_BASE_LM, _Inicializar.ConfMotor_BASE_LB, _Inicializar.ConfMotor_BASE_E);
				radio = sqrt(pow(_Inicializar.ConfMotor_PISTON_DB, 2) * (1 - d_piston / _Inicializar.ConfMotor_BASE_S) + ((d_piston / _Inicializar.ConfMotor_BASE_S)) * pow(_Inicializar.ConfMotor_BASE_D, 2)) / 2;
				_CU = _Inicializar.ConfMotor_CULATA_CTM * (2 * Constantes::pi * _Inicializar.VarMed_N) * radio * 0.75 * pow((_Inicializar.ConfMotor_BASE_D / (2 * radio)), 2);
			}
			else {
				_CU = 0;
			}
			factor = 1;
		}

		if (_InicializarPorAngulo.ciclo_cerrado)
		{
			c1 = _InicializarAlRCA.CW1 + _InicializarAlRCA.CW2 * _CU / _Inicializar.DatosCalcGeneral_CM;
		}
		else
			c1 = _Inicializar.ConfMotor_BASE_WR1A + _Inicializar.ConfMotor_BASE_WR1B * _CU / _Inicializar.DatosCalcGeneral_CM;


		if (_Inicializar.Constante_TRANS_CALOR_C_N == 0.8) {
			double pot1 = pow(_Inicializar.ConfMotor_BASE_D, (_Inicializar.Constante_TRANS_CALOR_C_N - 1));
			double pot2 = pow(T, (-0.53));
			double pot3 = pow(p, _Inicializar.Constante_TRANS_CALOR_C_N);
			double vg = (c1 * _Inicializar.DatosCalcGeneral_CM + _Inicializar.ConfMotor_BASE_W2 * COMB);
			double pot4 = pow(vg, _Inicializar.Constante_TRANS_CALOR_C_N);
			//_h = v.C_WH1 * v.d ^ (v.C_N - 1) * T ^ (-0.53) * p ^ v.C_N * vg ^ v.C_N
			_h = _Inicializar.Constante_TRANS_CALOR_C_WH1 * pot1  * pot2 * pot3 * pot4;
		}
		else {
			_h = _Inicializar.Constante_TRANS_CALOR_C_WH1 * pow(_Inicializar.ConfMotor_BASE_D, (_Inicializar.Constante_TRANS_CALOR_C_N - 1)) * pow(T, (0.75 - 1.62 * _Inicializar.Constante_TRANS_CALOR_C_N)) * pow(p, _Inicializar.Constante_TRANS_CALOR_C_N) * pow((c1 * _Inicializar.DatosCalcGeneral_CM + _InicializarAlRCA.CW2 * COMB), _Inicializar.Constante_TRANS_CALOR_C_N);
		}

		//h = v.C_WH1 * v.d ^ (-0.2) * T ^ (-0.53) * p ^ 0.8 * (c1 * CM + C_W2 * COMB) ^ 0.8 'W/m2 K

		_QCUL = _InicializarPorAngulo.deltat * factor * _h * _Inicializar.ConfMotor_CULATA_ACUL * (T - _InicializarAlRCA.TCUL);
		_QPIS = _InicializarPorAngulo.deltat * factor * _h * _Inicializar.ConfMotor_PISTON_AP * (T - _InicializarAlRCA.TPIS);
		_QCIL = _InicializarPorAngulo.deltat * factor * _h * _InicializarPorAngulo.areacil * (T - _InicializarAlRCA.TCIL);
		_QW = _QCUL + _QPIS + _QCIL;

		if (_QW != _QW)
		{
			_QW = _QW;
		}

		_QRAD = 0;

	}
	catch (const std::exception& e)
	{
		_comun->TratarError("Ajustes::Calcula_Calor_Woschni", e.what());
		throw e;
	}
	catch (...)
	{
		_comun->TratarError("Ajustes::Calcula_Calor_Woschni", "Error Desconocido");
		throw;
	}
}

void TransmisionCalor::Calcular_Temp_Media_Nodos(short Z, std::string nombre_nodo, double *temp_nodo)
{
	int id_nodo;
	int cont = 0;
	double T = 0.;

	for (int i = 1; i <= _NODOS_CUL.size(); i++)
	{
		int es_nodo = _NODOS_CUL.at(i).NOMBRE.find(nombre_nodo);
		if (es_nodo >= 0)
		{
			id_nodo = _NODOS_CUL.at(i).ID_NODO;
			break;
		}

		for (int j = 0; j < _RESULTADOS_T.size(); j++)
		{
			if ((id_nodo == _RESULTADOS_T.at(j).ID_NODO1) && (_RESULTADOS_T.at(j).Z == Z))
			{
				cont++;
				T = T + _RESULTADOS_T[j].valor;
				break;
			}
		}
	}

	if (cont > 0)
		*temp_nodo = T / cont;
	else
		*temp_nodo = 0;

}

void TransmisionCalor::ImprimirMatriz(std::vector< std::vector<double> >* m1, std::string nombre)
{
	ofstream fs("c:\\temp\\" + nombre + ".txt");

	// Enviamos una cadena al fichero de salida:

	short f;
	short c1;
	short c2;
	short c;
	short cont1;
	short cont2;

	f = m1->size() - 1;
	c1 = m1->at(0).size() - 1;

	for (cont1 = 1; cont1 <= f; cont1++) {
		for (cont2 = 1; cont2 <= c1; cont2++) {
			fs << m1->at(cont1).at(cont2) << ";";
		}
		fs << endl;
	}
	// Cerrar el fichero, 
	// para luego poder abrirlo para lectura:
	fs.flush();
	fs.close();
}

void TransmisionCalor::ImprimirVector(std::vector<double> * m1, std::string nombre)
{
	ofstream fs("c:\\temp\\" + nombre + ".txt");

	// Enviamos una cadena al fichero de salida:

	short f;
	short c1;
	short c2;
	short c;
	short cont1;
	short cont2;

	f = m1->size() - 1;

	for (cont1 = 1; cont1 <= f; cont1++) {
		fs << m1->at(cont1) << ";";
		fs << endl;
	}
	// Cerrar el fichero, 
	// para luego poder abrirlo para lectura:
	fs.flush();
	fs.close();


}

//}