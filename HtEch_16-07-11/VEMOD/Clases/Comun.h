#include <string>
#include <vector>
//#include "DatosCalculo.h"
//#include "Comunicacion.h"

//namespace CalculoCALMEC
//{
//extern std::string ruta;


#define KMAX 65540 //Doble del número maximo de armónicos para el cálculo de ruido de Tono
#define KMIN (KMAX/2) //Numero maximo de armónicos para el cálculo de ruido de Tono

#ifndef __COMUN_H
#define __COMUN_H


// La clase que almacena funciones comunes
class Comun
{
public:
	struct CInyInicioFinal
	{
		//Public numero As Integer
		double inicio;
		double fin;
		double duracion;
		//se va a asumir que si el inicio está en el ciclo abierto, toda la inyección lo estará
		bool inicio_ciclo_cerrado;
		bool fin_ciclo_cerrado;
	};

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

	struct stResultInyecciones
	{
		double NINY; //clave primaria?
		double AFINY;
		double AINY;
		double DURINY; //Tiempo
		double MFINY; //Antes MINY
		bool inicio_ciclo_cerrado;
		bool fin_ciclo_cerrado;
		double ZONA_INI;
		double ZONA_FIN;
		double T_DELAY; //comentar con Jaime para ver que es ... NORMA
	};

	struct Alfa_Valor
	{
		double eje_x;
		double eje_y;
	};

	class ExcepcionAplicacion : public std::exception
	{
		virtual const char* what() const throw()
		{
			return _texto_error.c_str();
		}
	private:
		std::string _texto_error = "";
	public:
		void ColocarTexto(std::string texto)
		{
			_texto_error = texto;
		}
	} exAplic;

	void TratarError(std::string donde, std::string error);
	void GuardarEnLog(std::string mensaje);
	bool Obtener_Ciclo_medio(std::vector<double>* datos, int NPC, int NCL, std::vector<double>* resultado);
	double modX(double n, int numpuntos);
	int Int(double numero);
	int Busca_Potencia(int numero);
	double Log10(double numero);
	bool desviacion_standard(std::vector<double>* vector, double* med, double* dev, double* cov);
	void FFT_Main(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n);
	void FFT_MainOld(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n);
	void IFFT_Main(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n);
	void IFFT_MainOld(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n);
	void Normaliza(std::vector<double>* sin_norma, std::vector<double>* con_norma, double num_sin, double num_con);
	bool Maximo_Vector_Parabola(std::vector<double>* vector_y, double* Y, double* x, double interv, short puntos_ajuste, int desde, int hasta);
	int ObtenerArmonicoCorte(int cRevCiclo, double rpm, double RegimenCorte0, double RegimenCorte1, int ArmonicoCorte0, int ArmonicoCorte1);
	bool calculo_PMI_PMB(int cGradosCiclo, double vel_ang, std::vector<double>* presionr, double* PMI, double* PMB, double desfase, double vcc, double d, double LM, double LB, double desc, double VD, double inc, double intang, double factp, double facti);
	double Volumen_def(double angulo, double vel_ang, double vcc, double inc, double pre, double factp, double facti, double intang, double D, double LM, double LB, double desc, double VD);
	double Aceleracion(double ang, double vel_ang, double  LM, double LB, double desc);
	double DerivadaVolumen(double angulo, double  LM, double LB, double desc, double D, double factp, double facti, double der_pre, double vel_ang);
	bool Eliminar_Offset_Inyeccion(int NPC, double intang, int cRevCiclo, double cGradosCiclo, std::vector<double>* vector, double desfase_pms);
	bool Calcula_Inicio_Fin_Inyeccion(int NPC, double n, double intang, double cGradosCiclo, double ADETOT, double RCA, double AAE, int NINY, std::vector<double>* ciclo_, std::vector<CInyInicioFinal>& Result);
	//bool Filtra_Señal(std::vector<double>* datos, int num_puntos, int armonico_corte, int ancho_corte, std::vector<double>* resultados, double m1, double m2, double f1, double f2, bool fsel);
	bool Filtra_Señal(std::vector<double>* datos, int num_puntos, int armonico_corte, int ancho_corte, bool derivar, int cRevCiclo, std::vector<double>* resultados, std::vector<double>* resultados_derivada, double m1, double m2, double f1, double f2, bool fsel);
	double modang(double ang, int cGradosCiclo);
	void maximo_vector(std::vector<double>* vector, int *indice);
	bool is_digits(const std::string &str);
	std::string ExtFichero(std::string fichero);
	int UBound(std::vector<double>* vector);
	int UBound(std::vector<Comun::Alfa_Valor>* vector);
	std::wstring s2ws(const std::string& s);
	//void Prepara2(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n1, int n);
	double Log2(double numero);
	//void FftOld(int n, int power, std::vector<double>* r_, std::vector<double>* i_);
	void Fft(int n, int power, std::vector<double>* r_, std::vector<double>* i_);
	void Regresion(std::vector<double>* suma, std::vector<double> *Fc, double *b0, double *b1);
	void SetIdEnsayo(std::string id);
	void SetRutaSalida(std::string ruta);
	std::string getRutaSalida();
	void Interpola(int n, int nf, int grado, std::vector<double> *a1, std::vector<double> *v1, std::vector<double> *a2, std::vector<double> *v2);//De Tono: lo usa para interpolar la presión a 2k puntos (Revisar e incorporar: se ha dejado Normaliza de Calmec en el filtrado)
	void PreparaOld(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n);
	void Prepara(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n);//De Tono: prepara vector de p para calcular fft
	void ffton(int nf, double regimen, std::vector<double>* a1, std::vector<double> *v1, std::vector<double> *frec, std::vector<double> *esp, int cRevCiclo);//De Tono: prepara vector de p y calcula fft para obtener el ruido
	//static void ColocarValor(std::string cadena);
	void EliminarPosicionCeroVector(std::vector<double> *vector);
	void ColocarPosicionCeroVector(std::vector<double> *vector);
	int Obtener_Posicion_en_Vector(double ang, double desfase, int NPC, double intang);
	void Calcula_Modulo(std::vector<double>* vector, std::vector<double>* modulo, int n_puntos);
	double modCGradosCiclo(double ang, int cGradosCiclo);
	double trunc_doub(double val, int digits);
	double round_doub(double value, int digits);
	double Interpolar_valor(double valor_x, std::vector<Comun::Alfa_Valor>* vector, bool interpolar_fuera_intervalo);
	double Interpola(double x, double x1, double x2, double y1, double y2);
	double BuscaCero(double x1, double x2, double y1, double y2);
	void Referencia_en_Y(double intang, int NPC, int cRevCiclo, std::vector<double> *presionf, std::vector<double> *presionr, double valor, double RCA, double ADETOT, double *premax);

	int Busqueda_Dico(double valor_x, double *valor_y, std::vector<Comun::Alfa_Valor>* vector);
	void Busqueda_Dico(double valor_x, double *valor_y, std::vector<double>* vector_x, std::vector<double>* vector_y);

	double ObtenerRcDinamica(double angulo, double vcc, double voldef, double d, double LM, double LB, double e, double VD);
	void Deriva_Señal(std::vector<double> *datos, int num_puntos, int armonico_corte, int ancho_corte, int cRevCiclo, std::vector<double> *resultados);

	//Calculo propiedades
	void C_v(std::string tipo_comb, std::string tipo_comb1, std::string tipo_comb2, double YFUEL1, double YFUEL2, double Xq1_BLEND, double PMC, double PMC1, double T, double* cv_a, double* cv_f, double* cv_q);
	void En_int(std::string tipo_comb, std::string tipo_comb1, std::string tipo_comb2, double YFUEL1, double YFUEL2, double Xq1_BLEND, double PMC, double PMC1, double T, double *u_a, double *u_f, double *u_q);
	double CalculoHCL(double CAHFL, double CBHFL, double TF, double C_TIY);

	void calcula_composicion_RCA_AAE(double MACC, double MAADM, double megr, double MFCC, double Mfa, double Mfadm, double Mainya, double MBB, double mres, double MCC,
		double FE, double *YaRCA, double *YaAAE, double *YqRCA, double *YqAAE, double *YfRCA, double *YfAAE, double *Yaadm, double *Yaesc, double *Yqadm,
		double *Yqesc, double *Yfadm, double *Yfesc);
	double Volumen(double angulo, double vcc, double Diametro, double mani, double Biela, double desc);


	double Xp(double x, double ratio_ctms);
	void Calcula_Blow_By(double pre, double CBB, double temp, double gamma, double r, double *masabb, double pc, double Cfbb, double inc, double d);
	double Obtener_PresionPMI(std::vector<double> *Presion, double desfase, double ang_pmi, double intang);
	bool Impar(short valor);
	void Maximo_Parabola(std::vector<Alfa_Valor> *origen, double intang, double *x, double *Y);
	void Desfasar_Vector(std::vector<double>* vector_origen, std::vector<double>* vector_destino, int desfase, int NPC, double intang);
	double Dist_Despl_inst(double angulo, double mani, double Biela, double desc);
	void Convierte_Inst(std::vector<double>* datos, double primer_angulo, std::vector<double> *result_datos, std::vector<double> *result_angulos, double fact, double delta, double intang);
	double asn(double x);
	double DistPistonCulataDef(double angulo, double vcc, double voldef, double d, double LM, double LB, double e, double h0);
	bool Ajuste_Regresion(std::vector<double> *vector_y, std::vector<double>* vector_x, double* Y, double* x, short n);
	void Calcula_Rendimientos(double Plim, double RCG, double VD, double Pa, double HP, double MFCC, double WIHP, double WILP, double WE, double Qtrans,
		double *rend_comb, double *rend_teo, double *rend_forma, double *rend_bombeo, double *rend_adiab, double *rend_meca);
	double MediaVectorEmpiezaEn1(std::vector<double> *vector);
	double DistanciaEntreAngulos(double ang_ini, double ang_fin, double cGradosCiclo);
	void maximo_vector(std::vector<int>* vector, double* max);
	void minimo_vector(std::vector<int>* vector, double* min);
	double obtener_min(double valor1, double valor2);
	void Puntos_G(double rc, double *g1, double *g2, double lambda);
private:
	//Comunicacion *_com;

	std::string _id_ensayo = ""; //Se apunta aqui para imprimir avisos y errores
	std::string _ruta_salida;
	bool Ordena_Inyecciones(std::vector<CInyInicioFinal>& inyecciones);
	void ifftOld(int n, int power, std::vector<double>* r_, std::vector<double>* i_);
	void ifft(int n, int power, std::vector<double>* r_, std::vector<double>* i_);
	void RecomponeOld(std::vector<double>* vector, std::vector<double>* A, std::vector<double>* b, int n);
	void Recompone(std::vector<double>* vector, std::vector<double>* A, std::vector<double>* b, int n);
	//double modCGradosCiclo(double ang, DatosCalculo* ent);
	std::string ToUpper(std::string str);
	double Polinomio(int ii, int grado, int kk, int nf, std::vector<double> *a1, std::vector<double> *v1, std::vector<double> *a2);
	void Filtra_Deriva(std::vector<double> *vector, int corte, int ancho, int n_puntos, int cRevCiclo);
protected:
	//void Filtra(std::vector<double>* vector, int corte, int ancho, int n_puntos, double m1, double m2, double f1, double f2, bool fsel, short num_cil = 0);
	void Filtra(std::vector<double>* vector, std::vector<double>* vector_derivado, bool derivar, int cRevCiclo, int corte, int ancho, int n_puntos, double m1, double m2, double f1, double f2, bool fsel, short num_cil);
};

#endif
//}