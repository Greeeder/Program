#include "Comun.h"
#include "Constantes.h"
#include <cmath>
#include "windows.h" //Para poder hacer OutputDebugString
#include <ctime>
#include <Eigen/Core>
//#include "global.h"

//namespace CalculoCALMEC
//{

void Comun::TratarError(std::string donde, std::string error)
{
	std::string mensajeCompleto;
	FILE *fp;

	try
	{
		mensajeCompleto = donde + " " + error;

		//mensajeCompleto = g_ruta_comun;

		/*char *path = NULL;
		size_t size;
		path = getcwd(path, size);*/

		// current date/time based on current system
		time_t now = time(0);

		// convert now to string form
		char* dt = ctime(&now);

		std::string diaHora(dt);

		mensajeCompleto = diaHora + ":" + mensajeCompleto;

		if (_ruta_salida != "")
		{
			fp = fopen((_ruta_salida + "\\TRAZA_" + _id_ensayo + ".dat").c_str(), "a");

			fprintf(fp, "%s\n", mensajeCompleto.c_str());

			fclose(fp);
		}

		//Escribir en fichero log los errores y la traza
		/*fich_errores = new Fichero();
		fich_errores->AbrirFichero("F:/LogEjecucionCALMEC.txt");
		fich_errores->EscribirLinea(mensajeCompleto);
		fich_errores->CerrarFichero();*/

		//OutputDebugString(s2ws(mensajeCompleto).c_str());
	}
	catch (const std::exception& e)
	{
		throw;
	}
}

/*void Comun::ColocarValor(std::string cadena)
{
g_ruta_comun = cadena;
}*/

void Comun::GuardarEnLog(std::string mensaje)
{
	FILE *fp;

	try
	{

		/*char *path = NULL;
		size_t size;
		path = getcwd(path, size);*/

		// current date/time based on current system
		time_t now = time(0);

		// convert now to string form
		char* dt = ctime(&now);

		std::string diaHora(dt);

		mensaje = diaHora + ":" + mensaje;

		if (_ruta_salida != "")
		{
			fp = fopen((_ruta_salida + "\\TRAZA_" + _id_ensayo + ".dat").c_str(), "a");

			fprintf(fp, "%s\n", mensaje.c_str());

			fclose(fp);
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::GuardarEnLog", e.what());
		throw;
	}
	catch (...)
	{
		TratarError("Comun::GuardarEnLog", "Error Desconocido");
		throw;
	}
	//OutputDebugString(s2ws(mensajeCompleto).c_str());
}

void Comun::SetIdEnsayo(std::string id)
{
	_id_ensayo = id;
}

void Comun::SetRutaSalida(std::string ruta)
{
	_ruta_salida = ruta;
}

std::string Comun::getRutaSalida()
{
	return _ruta_salida;
}

std::wstring Comun::s2ws(const std::string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}

bool Comun::desviacion_standard(std::vector<double>* vector, double* med, double* dev, double* cov)
{
	int npuntos;
	double suma1;
	double suma2;
	int j;

	try
	{
		npuntos = UBound(vector);

		if (npuntos == 1) {
			*med = vector->at(1);
			*dev = 0;
			return false; //no hay error, pero se sale del método
		}


		*med = 0;

		for (j = 1; j <= npuntos; j++) {
			suma1 = suma1 + pow(vector->at(j), 2);
			suma2 = suma2 + vector->at(j);
			*med = *med + vector->at(j);
		}

		*med = *med / npuntos;
		*dev = sqrt((npuntos * suma1 - pow(suma2, 2)) / (npuntos * (npuntos - 1)));
		*cov = *dev / *med * 100;

		return false;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::desviacion_standard", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::desviacion_standard", "Error Desconocido");
		throw;
	}
}

void Comun::Filtra(std::vector<double>* vector, std::vector<double>* vector_derivado, bool derivar, int cRevCiclo, int corte, int ancho, int n_puntos, double m1, double m2, double f1, double f2, bool fsel, short num_cil)
{
	//Este procedimiento pasa la señal de entrada (vector()) al dominio de la
	//frecuencia, y elimina los armónicos que hay por encima de corte si
	//ancho es cero , si ancho es mayor que cero lo que hace es eliminar los
	//armónicos pero de modo más suave.
	//Este filtro se comporta como un PASO BAJO
	//
	// vector()  ---->  señal de entrada sin filtrar y de salida filtrada
	//vector_derivado() -----> salida derivada
	//derivar -----> true si se quiere derivar el vector enviado
	// corte     ---->  armónico a partir del cual queremos filtrar
	// ancho     ---->  ancho (en armónicos) del filtro "suave"
	// n_puntos  ---->  nº de puntos del vector de entrada y de salida
	//
	// NOTA: es necesario que n_puntos sea potencia de dos

	int cont;
	double i;
	double rel;
	double factor;
	std::vector<double> parte_real;
	std::vector<double> parte_imag;
	double Temp;
	double m;
	double foo;
	short descriptor;

	try
	{

		parte_real.resize(n_puntos / 2 + 1);
		parte_imag.resize(n_puntos / 2 + 1);

		//Devuelve la parte_real y parte_imaginaria
		Comun::FFT_Main(vector, &parte_real, &parte_imag, n_puntos - 1);
		//se pasa al dominio de la frecuencia

		cont = 0;


		if (fsel) {
			for (i = f1; i <= f2; i++) {
				Temp = (Log10(m2) - Log10(m1)) / (Log10(f2) - Log10(f1)) * Log10(i) + (Log10(m1) - Log10(f1) * ((Log10(m2) - Log10(m1)) / (Log10(f2) - Log10(f1))));

				Temp = pow(10, Temp);
				m = sqrt(pow(parte_real.at(i), 2) + pow(parte_imag.at(i), 2));

				parte_real.at(i) = Temp * parte_real.at(i) / m;
				parte_imag.at(i) = Temp * parte_imag.at(i) / m;
			}
		}

		if (ancho != 0) {
			rel = Constantes::pi / ancho;
			for (i = (corte - Int(ancho / 2)); i <= (corte + Int(ancho / 2)); i++) {
				factor = cos(cont * rel) / 2 + 0.5;
				parte_real.at(i) = parte_real.at(i) * factor; //'filtrado "suave"
				parte_imag.at(i) = parte_imag.at(i) * factor;
				cont = cont + 1;
			}
		}

		for (i = (corte + Int(ancho / 2) + 1); i <= (Int(n_puntos / 2)); i++) {
			parte_real.at(i) = 0;
			parte_imag.at(i) = 0;//eliminamos el resto de armónocos
		}
		
		//derivamos
		if (derivar)
		{
			std::vector<double> parte_real_der;
			std::vector<double> parte_imag_der;

			parte_real_der.resize(parte_real.size());
			memcpy(&parte_real_der.at(0), &parte_real.at(0), parte_real.size());

			parte_imag_der.resize(parte_imag.size());
			memcpy(&parte_imag_der.at(0), &parte_imag.at(0), parte_imag.size());


			for (i = 1; i <= n_puntos / 2; i++) {
				Temp = parte_real_der.at(i);
				parte_real_der.at(i) = -1 * (i - 1) * parte_imag_der.at(i);
				parte_imag_der.at(i) = (i - 1) * Temp;
			}

			IFFT_Main(vector_derivado, &parte_real_der, &parte_imag_der, n_puntos);
			//una vez filtrada y derivada pasamos la señal al dominio del tiempo

			for (i = 0; i <= n_puntos - 1; i++) {
				vector_derivado->at(i) = vector_derivado->at(i) / cRevCiclo; //esto se hace porque la derivada sale multiplicada por dos
			}
		}

		IFFT_Main(vector, &parte_real, &parte_imag, n_puntos);
		//pasamos al dominio del tiempo



	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Filtra", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Filtra", "Error Desconocido");
		throw;
	}
}

bool Comun::Maximo_Vector_Parabola(std::vector<double>* vector_y, double* Y, double* x, double interv, short puntos_ajuste, int desde, int hasta)
{
	//Este procedimiento encuentra el máximo en un vector y además hace un
	//ajuste parábolico en torno al máximo y devuelve sus coordenadas
	//
	//    vector_y ---> vector de entrada
	//    y        ---> coordenada y del máximo (ajustado)
	//    x        ---> coordenada x del máximo (ajustado)
	//    interv  ---> entre dos puntos en x qué "distancia" hay ?
	//    puntos_ajuste ---> nº de puntos para el ajuste parabólico
	//    desde ---> empezamos a buscar aquí
	//    hasta ---> acabamos de buscar aquí

	int i;
	double maximo;
	std::vector<double> v_y;
	std::vector<double> v_x;
	int mitad;

	try
	{
		v_y.resize(puntos_ajuste + 1);
		v_x.resize(puntos_ajuste + 1);

		*x = desde;
		maximo = vector_y->at(*x);
		for (i = desde; i <= hasta; i++) {
			if (vector_y->at(i) > maximo) {
				maximo = vector_y->at(i);
				*x = i;
			}
		}

		//aquí falta un control de errores por si nos pasamos

		mitad = Int(puntos_ajuste / 2);

		int PuntoActual;
		for (i = 0; i <= puntos_ajuste; i++) {
			PuntoActual = *x - mitad + i;
			//comprueba que no sobrepasa los limites del vector
			if ((PuntoActual > UBound(vector_y)))
				PuntoActual = PuntoActual - UBound(vector_y);
			if ((PuntoActual < 1))
				PuntoActual = UBound(vector_y) - (-PuntoActual);
			v_y.at(i) = vector_y->at(PuntoActual);
			v_x.at(i) = (*x - mitad + i - 1) * interv;
		}

		Ajuste_Regresion(&v_y, &v_x, Y, x, puntos_ajuste);

		return true;

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Maximo_Vector_Parabola", e.what());
		throw e;
		return false;
	}
	catch (...)
	{
		TratarError("Comun::Maximo_Vector_Parabola", "Error Desconocido");
		throw;
	}

}

bool Comun::Ajuste_Regresion(std::vector<double> *vector_y, std::vector<double>* vector_x, double* Y, double* x, short n)
{
	//Este procedimiento hace un ajuste por mínimos cuadrados
	//de una parábola a un conjunto de pares de puntos,
	//el número de puntos es n y los pares de puntos van en
	//vector_x() y vector_y(), como resultado da las coordenadas
	//del punto máximo de la parábola que más se ajusta a esa nube
	//de puntos, también puede dar los coeficientes del polinomio
	//a, b y c

	double sa3;
	double sa;
	double sa2;
	double sa4;
	double sap;
	double sp;
	double sa2p;
	double s;
	double b;
	double A;
	double c;
	short i;

	try
	{


		sa = 0;
		sa2 = 0;
		sa3 = 0;
		sa4 = 0;


		for (i = 0; i <= n - 1; i++) {
			sa = sa + vector_x->at(i);
			sa2 = sa2 + pow(vector_x->at(i), 2);
			sa3 = sa3 + pow(vector_x->at(i), 3);
			sa4 = sa4 + pow(vector_x->at(i), 4);
		}

		s = (sa4 * sa2 * n) + (sa3 * sa2 * sa) + (sa3 * sa2 * sa) - (sa2 * sa2 * sa2) - (sa3 * sa3 * n) - (sa4 * sa * sa);

		//determinante de la matriz S
		sp = 0;
		sap = 0;
		sa2p = 0;

		for (i = 0; i <= n - 1; i++) {
			sp = sp + vector_y->at(i);
			sap = sap + vector_x->at(i) * vector_y->at(i);
			sa2p = sa2p + pow(vector_x->at(i), 2) * vector_y->at(i);
		}


		//Estos son los coeficientes de la parábola
		A = ((sa2p * sa2 * n) + (sa3 * sa * sp) + (sap * sa * sa2) - (sa2 * sa2 * sp) - (sap * sa3 * n) - (sa * sa * sa2p)) / s;
		b = ((sa4 * sap * n) + (sa2p * sa * sa2) + (sa3 * sp * sa2) - (sa2 * sap * sa2) - (sa3 * sa2p * n) - (sa4 * sp * sa)) / s;
		c = (sp - b * sa - A * sa2) / n;

		//Calculamos la derivada e igualamos a cero para despejar
		//la x donde se ha dado el máximo.
		*x = -b / (2 * A);

		//Sustituimos en la ecuacion de la parabola para sacar la y
		*Y = A * pow(*x, 2) + b * *x + c;

		return true;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Ajuste_Regresion", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Ajuste_Regresion", "Error Desconocido");
		throw;
	}
}


bool Comun::Filtra_Señal(std::vector<double>* datos, int num_puntos, int armonico_corte, int ancho_corte, bool derivar, int cRevCiclo, std::vector<double>* resultados, std::vector<double>* resultados_derivada, double m1, double m2, double f1, double f2, bool fsel)
{
	int l;
	int pot;
	short descriptor;

	try
	{
		resultados->resize(num_puntos + 1);

		if (derivar)
			resultados_derivada->resize(num_puntos + 1);

		pot = Busca_Potencia(num_puntos);
		std::vector<double> normalizado;
		std::vector<double> derivada_normalizado;

		normalizado.resize(pot + 1);

		derivada_normalizado.resize(pot + 1);

		if (armonico_corte > (pot / 2))
		{
			exAplic.ColocarTexto("armonico_corte > pot/2");
			throw exAplic;
		}

		//OJO
		datos->at(0) = datos->at(num_puntos);

		//Devuelve vector normalizado
		Normaliza(datos, &normalizado, num_puntos, pot);

		//OJO paso como ultimo parametro num_cil 0 que es el valor por defecto
		Filtra(&normalizado, &derivada_normalizado, derivar, cRevCiclo, armonico_corte, ancho_corte, pot, m1, m2, f1, f2, fsel, 0);

		//Devuelve los resultados
		Normaliza(&normalizado, resultados, pot, num_puntos);

		if (derivar)
			Normaliza(&derivada_normalizado, resultados_derivada, pot, num_puntos);

		resultados->at(num_puntos) = resultados->at(0);
		
		if (derivar)
			resultados_derivada->at(num_puntos) = resultados_derivada->at(0);


		return false;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Filtra_Señal", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Filtra_Señal", "Error Desconocido");
		throw;
	}
}


double Comun::modX(double n, int numpuntos)
{
	try
	{
		if (n >= 1 & n <= numpuntos) {
			return n;
		}

		if (n > numpuntos) {
			return n - numpuntos;
		}

		if (n < 1) {
			return n + numpuntos;
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::modX", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::modX", "Error Desconocido");
		throw;
	}
}

void Comun::IFFT_Main(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n)
{

	//calcula la transformada inversa rápida de Fourier de una var.
	//compleja proporcionando la var. temporal (vector())

	int nf;
	int n1;
	int power;

	try
	{
		//Elimino la posicion 0 de los vectores que es 0. 
		//y que el codigo de tono espera los vectores desde la posicion 0
		EliminarPosicionCeroVector(p_r);
		EliminarPosicionCeroVector(p_i);

		n1 = n;
		nf = n;
		power = 0;
		while (nf > 1) {
			nf = Int(nf / 2);
			power = power + 1;
		}
		ifft(n1, power, p_r, p_i);

		//ENVIO n - 1 Y NO n es necesario para el codigo de recompone de tono
		n = n - 1;
		Recompone(vector, p_r, p_i, n);

		ColocarPosicionCeroVector(p_r);
		ColocarPosicionCeroVector(p_i);
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::IFFT_Main", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::IFFT_Main", "Error Desconocido");
		throw;
	}
}

void Comun::IFFT_MainOld(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n)
{
	//calcula la transformada inversa rápida de Fourier de una var.
	//compleja proporcionando la var. temporal (vector())

	int nf;
	int n1;
	int power;

	try
	{
		n1 = n;
		nf = n;
		power = 0;
		while (nf > 1) {
			nf = Int(nf / 2);
			power = power + 1;
		}
		ifftOld(n1, power, p_r, p_i);

		RecomponeOld(vector, p_r, p_i, n);

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::IFFT_Main", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::IFFT_Main", "Error Desconocido");
		throw;
	}
}

void Comun::RecomponeOld(std::vector<double>* vector, std::vector<double>* A, std::vector<double>* b, int n)
{
	int i;
	int j;

	try
	{
		j = -1;
		for (i = 1; i <= Int(n / 2); i++) {
			j = j + 1;
			vector->at(j) = A->at(i);
			j = j + 1;
			vector->at(j) = b->at(i);
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Recompone", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Recompone", "Error Desconocido");
		throw;
	}
}

void Comun::Recompone(std::vector<double>* vector, std::vector<double>* A, std::vector<double>* b, int n)
{
	int i;
	int j;

	try
	{
		j = -1;
		for (i = 0; i <= Int(n / 2); i++) {
			j = j + 1;
			vector->at(j) = A->at(i);
			j = j + 1;
			vector->at(j) = b->at(i);
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::RecomponeTono", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::RecomponeTono", "Error Desconocido");
		throw;
	}
}

void Comun::ifft(int n, int pot, std::vector<double>* r_, std::vector<double>* i_)
{
	//Este código se actualizó con el código de Tono 24/11/2015--> Ya se usa para ruido y filtrado
	int i, j, k, r, g;
	double a, b, c, e, p, q, u, v;

	n = n / 2;
	a = Constantes::pi / n;
	p = cos(a);
	q = -1.0*sin(a);
	a = r_->at(0);
	r_->at(0) = a + i_->at(0);
	i_->at(0) = a - i_->at(0);
	c = -1.0;
	e = 0.0;
	for (j = 2; j <= (n / 2); j++)
	{
		a = e*p + c*q;
		c = c*p - e*q;
		e = a;
		k = n - j + 2;
		a = r_->at(j - 1) + r_->at(k - 1);
		b = (i_->at(j - 1) + i_->at(k - 1))*c;
		b -= e*(r_->at(j - 1) - r_->at(k - 1));
		u = i_->at(j - 1) - i_->at(k - 1);
		v = (i_->at(j - 1) + i_->at(k - 1))*e;
		v += c*(r_->at(j - 1) - r_->at(k - 1));
		r_->at(j - 1) = (a + b) / 2.0;
		i_->at(j - 1) = (u - v) / 2.0;
		r_->at(k - 1) = (a - b) / 2.0;
		i_->at(k - 1) = -(u + v) / 2.0;
	}

	i_->at(n / 2) = -1.0*i_->at(n / 2);

	k = 0;
	for (j = 1; j<n; j++)
	{
		i = 2;
		while (k >= (n / i))
		{
			k -= (n / i);
			i += i;
		}
		k += (n / i);
		if (k>j)
		{
			a = r_->at(j);
			r_->at(j) = r_->at(k);
			r_->at(k) = a;
			a = i_->at(j);
			i_->at(j) = i_->at(k);
			i_->at(k) = a;
		}
	}
	g = 1;
	p = 1.0;
	for (i = 1; i < pot; i++)
	{
		c = 1.0;
		e = 0.0;
		q = -1.0*sqrt((1.0 - p) / 2.0);
		if (i == 1)
			p = -1.0*sqrt((1.0 + p) / 2.0);
		else
			p = sqrt((1.0 + p) / 2.0);

		for (r = 1; r <= g; r++)
		{
			j = r;
			while (j <= n)
			{
				k = j + g;
				a = c*r_->at(k - 1) + e*i_->at(k - 1);
				b = e*r_->at(k - 1) - c*i_->at(k - 1);
				r_->at(k - 1) = r_->at(j - 1) - a;
				i_->at(k - 1) = i_->at(j - 1) + b;
				r_->at(j - 1) = r_->at(j - 1) + a;
				i_->at(j - 1) = i_->at(j - 1) - b;
				j += 2 * g;
			}
			a = e*p + c*q;
			c = c*p - e*q;
			e = a;
		}
		g += g;
	}
}



void Comun::ifftOld(int n, int power, std::vector<double>* r_, std::vector<double>* i_)
{
	int i;
	int j;
	int r;
	int g;
	int k;
	double A;
	double b;
	double c;
	double e;
	double p;
	double Q;
	double u;
	double v;

	try
	{

		n = Int(n / 2);
		A = Constantes::pi / n;
		p = cos(A);
		Q = -1 * sin(A);
		A = r_->at(1);
		r_->at(1) = A + i_->at(1);
		i_->at(1) = A - i_->at(1);
		c = -1;
		e = 0;


		for (j = 2; j <= Int(n / 2); j++) {
			A = e * p + c * Q;
			c = c * p - e * Q;
			e = A;
			k = n - j + 2;
			A = r_->at(j) + r_->at(k);
			b = (i_->at(j) + i_->at(k)) * c;
			b = b - (r_->at(j) - r_->at(k)) * e;
			u = i_->at(j) - i_->at(k);
			v = (i_->at(j) + i_->at(k)) * e;
			v = v + (r_->at(j) - r_->at(k)) * c;
			r_->at(j) = (A + b) / 2;
			i_->at(j) = (u - v) / 2;
			r_->at(k) = (A - b) / 2;
			i_->at(k) = -(u + v) / 2;
		}


		i_->at(Int(n / 2) + 1) = -i_->at(Int(n / 2) + 1);


		k = 0;
		for (j = 1; j <= n - 1; j++) {
			i = 2;
			while (k >= Int(n / i)) {
				k = k - Int(n / i);
				i = 2 * i;
			}
			k = k + Int(n / i);
			if (k > j) {
				A = r_->at(j + 1);
				r_->at(j + 1) = r_->at(k + 1);
				r_->at(k + 1) = A;
				A = i_->at(j + 1);
				i_->at(j + 1) = i_->at(k + 1);
				i_->at(k + 1) = A;
			}
		}
		g = 1;
		p = 1;
		for (i = 1; i <= power - 1; i++) {
			c = 1;
			e = 0;
			Q = -1 * sqrt((1 - p) / 2);
			if (i == 1) {
				p = -1 * sqrt((1 + p) / 2);
			}
			else {
				p = sqrt((1 + p) / 2);
			}
			for (r = 1; r <= g; r++) {
				j = r;
				while (j <= n) {
					k = j + g;
					A = c * r_->at(k) + e * i_->at(k);
					b = e * r_->at(k) - c * i_->at(k);
					r_->at(k) = r_->at(j) - A;
					i_->at(k) = i_->at(j) + b;
					r_->at(j) = r_->at(j) + A;
					i_->at(j) = i_->at(j) - b;
					j = j + 2 * g;
				}
				A = e * p + c * Q;
				c = c * p - e * Q;
				e = A;
			}
			g = 2 * g;
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::ifft", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::ifft", "Error Desconocido");
		throw;
	}
}

int Comun::UBound(std::vector<double>* vector)
{
	return vector->size() - 1;
}

int Comun::UBound(std::vector<Comun::Alfa_Valor>* vector)
{
	return vector->size() - 1;
}

bool Comun::Obtener_Ciclo_medio(std::vector<double>* datos, int NPC, int NCL, std::vector<double>* resultado)
{
	//Este procedimiento pasa la señal de entrada (vector()) al dominio de la
	//frecuencia, se queda con los armonicos que son multiplos de NCL y contruye
	//el ciclo medio a partir de estos armonicos
	//
	// vector()  ---->  señal de entrada sin filtrar y de salida filtrada
	// n_puntos  ---->  nº de puntos del vector de entrada y de salida
	// NOTA: es necesario que n_puntos sea potencia de dos

	int num_puntos;
	int pot;
	std::vector<double> normalizado;
	std::vector<double> parte_real;
	std::vector<double> parte_imag;
	std::vector<double> parte_real_ciclounico;
	std::vector<double> parte_imag_ciclounico;
	int cont;
	double referencia;
	int puntosref;
	short posiciones;
	int contArmonicos;

	try

	{
		num_puntos = UBound(datos);

		pot = Busca_Potencia(num_puntos);

		normalizado.resize(pot + 1); //Mas uno para la posicion 0

		datos->at(0) = datos->at(UBound(datos));

		//obtiene una referencia para acercar a cero todos los valores cogiendo el valor
		//medio de los puntos del inicio y del final de la señal para poder realizar correctamente
		//el enventanado de la señal
		referencia = 0;
		puntosref = 5;
		for (cont = 0; cont <= puntosref - 1; cont++) {
			referencia = referencia + datos->at(cont);
		}
		for (cont = num_puntos; cont >= num_puntos - puntosref + 1; cont += -1) {
			referencia = referencia + datos->at(cont);
		}
		referencia = referencia / (puntosref * 2);
		//referencia toda la señal
		for (cont = 0; cont <= num_puntos; cont++) {
			datos->at(cont) = datos->at(cont) - referencia;
		}

		//normaliza el megaciclo
		Normaliza(datos, &normalizado, num_puntos, pot);
		//FILE *fp;
		//fp = fopen("pruebaDebug.txt", "w");
		//for (int i = 0; i <= 32768; i++)
		//{
		//	fprintf(fp, "Iteracion: %d %f %f\n", i, datos->at(i), normalizado.at(i));
		//}
		//fclose(fp); //Cierro el fichero
		//realiza el enventanado de la señal normalizada
		int puntosVentana;
		puntosVentana = 5;
		//enventanado del principio de la señal
		for (cont = 0; cont <= puntosVentana - 1; cont++) {
			normalizado.at(cont) = normalizado.at(cont) * (double)cont / ((double)puntosVentana - 1);
		}
		//enventanado del final de la señal
		for (cont = pot - 1; cont >= (pot - 1) - puntosVentana + 1; cont += -1) {
			normalizado.at(cont) = normalizado.at(cont) * ((double)pot - 1 - cont) / ((double)puntosVentana - 1);
		}

		//realiza la fft de la señal normalizada
		parte_real.resize(pot / 2 + 1); //Mas 1 para dejar sitio al elemento 0
		parte_imag.resize(pot / 2 + 1); //Mas 1 para dejar sitio al elemento 0

		FFT_Main(&normalizado, &parte_real, &parte_imag, pot - 1);
		//se pasa al dominio de la frecuencia

		//recoge los armonicos que son multiplos de NCL porque son los que
		//definen al ciclo unico
		if ((pot / 2) % NCL > 0) {
			posiciones = Int(pot / 2 / NCL) + 1;
		}
		else {
			posiciones = Int(pot / 2 / NCL);
		}

		parte_real_ciclounico.resize(posiciones + 1);
		parte_imag_ciclounico.resize(posiciones + 1);

		contArmonicos = 0;
		for (cont = 1; cont <= pot / 2; cont += NCL) {
			contArmonicos = contArmonicos + 1;
			parte_real_ciclounico.at(contArmonicos) = parte_real.at(cont);
			parte_imag_ciclounico.at(contArmonicos) = parte_imag.at(cont);
		}

		//para realiza la IFFT el numero de puntos de la parte real e imaginaria deben ser
		//multiplos de 2. Se busca la potencia de 2 por debajo del numero de puntos
		pot = Busca_Potencia(UBound(&parte_real_ciclounico));
		if (pot > parte_real_ciclounico.size()) {
			pot = pot / 2;
		}

		//redimensiona los vectores para que tenga pot2 numero de puntos
		parte_real_ciclounico.resize(pot + 1);
		parte_imag_ciclounico.resize(pot + 1);

		//redimensiona el vector normalizado que contendra el resultado
		normalizado.resize(pot * 2 + 1);

		IFFT_Main(&normalizado, &parte_real_ciclounico, &parte_imag_ciclounico, pot * 2);
		//pasamos al dominio del tiempo

		resultado->resize(NPC + 1);
		Normaliza(&normalizado, resultado, UBound(&normalizado), NPC);
		resultado->at(NPC) = resultado->at(0);

		return false;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Obtener_Ciclo_medio", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Obtener_Ciclo_medio", "Error Desconocido");
		throw;
	}
}

void Comun::Normaliza(std::vector<double>* sin_norma, std::vector<double>* con_norma, double num_sin, double num_con)
{
	//Este procedimiento hace la conversión de un vector con un número de
	//puntos (sin_norma) en otro con un número de puntos diferente (con_norma)

	//los parámetros de entrada son:
	//    sin_norma --->   vector de entrada
	//    con_norma --->   vector de salida
	//    num_sin   --->   número de puntos del vector de entrada
	//    num_con   --->   número de puntos del vector de salida

	double intervalo_con;
	double intervalo_sin;
	int f;
	double punto_ant;
	double punto_post;
	double punto_int;
	double int_sin;
	double a_;
	double b;
	double b_;
	double rel;

	try
	{

		int_sin = 1 / num_sin;
		rel = num_sin / num_con;

		for (f = 0; f <= num_con; f++) {
			punto_int = f / num_con;
			punto_ant = Int(f * rel) / num_sin;
			punto_post = punto_ant + int_sin;
			if (round(punto_post * num_sin) < UBound(sin_norma)) {
				a_ = sin_norma->at(round(punto_post * num_sin)) - sin_norma->at(round(punto_ant * num_sin));
			}
			else {
				a_ = sin_norma->at(0) - sin_norma->at(round(punto_ant * num_sin));
			}
			b_ = int_sin;
			b = punto_int - punto_ant;
			con_norma->at(f) = sin_norma->at(round(punto_ant * num_sin)) + a_ * b / b_;
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Normaliza", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Normaliza", "Error Desconocido");
		throw;
	}
}

int Comun::Busca_Potencia(int numero)
{
	int result;
	//esta funcion retorna la potencia de dos mas proxima al
	//numero introducido por arriba
	double x;

	try
	{

		x = Log2((double)numero);
		x = Int(x);
		result = std::pow(2, x);

		if (result < numero) {
			result = std::pow(2, (Int(x) + 1));
		}

		return result;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Busca_Potencia", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Busca_Potencia", "Error Desconocido");
		throw;
	}
}

void Comun::FFT_Main(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n)
{
	//calcula la transformada rápida de Fourier de una var.
	//temporal proporcionando su parte real (p_r)
	//y su parte imaginaria (p_i)

	int nf;
	int power;
	int i;
	int T1;

	try
	{
		//Elimino la posicion 0 de los vectores que es 0. 
		//y que el codigo de tono espera los vectores desde la posicion 0
		EliminarPosicionCeroVector(p_r);
		EliminarPosicionCeroVector(p_i);

		nf = n;
		power = 1;

		while (nf > 1) {
			nf = Int(nf / 2);
			power = power + 1;
		}

		nf = 1;
		for (i = 1; i <= power; i++) {
			nf = nf * 2;
		}

		if ((n + 1 != nf)) {
			//if (permitir_interaccion)
			//{
			//	//MsgBox("Nombre de points incorrect");
			//	//err.Raise((0));

			//	//Enviar mensaje al interfaz grafico
			//	_com->EnviarCadena("MENSAJE");

			//	_com->EnviarEntero(12);

			//	//Espero al respuesta del usuario (interfaz), si es que NO se aborta el cálculo

			//	_com->LeerCadena(&result_mensaje);

			//	if (result_mensaje == "FIN")
			//	{
			//		//Se aborta el calculo

			//		_comun->exAplic.ColocarTexto("Cálculo abortado por condiciones MRESI < 0 o Mret < 0");

			//		throw _comun->exAplic;

			//	}
			//}
		}


		Prepara(vector, p_r, p_i, nf - 1);


		Fft(Int(nf / 2), power, p_r, p_i);

		ColocarPosicionCeroVector(p_r);
		ColocarPosicionCeroVector(p_i);

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::FFT_Main", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::FFT_Main", "Error Desconocido");
		throw;
	}
}

//void Comun::FFT_MainOld(std::vector<double>* vector, std::vector<double>* p_r, std::vector<double>* p_i, int n)
//{
//	//calcula la transformada rápida de Fourier de una var.
//	//temporal proporcionando su parte real (p_r)
//	//y su parte imaginaria (p_i)
//
//	int nf;
//	int power;
//	int i;
//	int T1;
//
//	try
//	{
//		nf = n;
//		power = 1;
//
//		while (nf > 1) {
//			nf = Int(nf / 2);
//			power = power + 1;
//		}
//
//		nf = 1;
//		for (i = 1; i <= power; i++) {
//			nf = nf * 2;
//		}
//
//		if ((n + 1 != nf)) {
//			//if (CMostrarInterfaz)
//			//	MsgBox("Nombre de points incorrect");
//			//err.Raise((0));
//		}
//
//		PreparaOld(vector, p_r, p_i, nf);
//
//		FftOld(Int(nf / 2), power, p_r, p_i);
//
//
//	}
//	catch (const std::exception& e)
//	{
//		TratarError("Comun::FFT_Main", e.what());
//		throw e;
//	}
//	catch (...)
//	{
//		TratarError("Comun::FFT_Main", "Error Desconocido");
//		throw;
//	}
//}

//void Comun::PreparaOld(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n)
//{
//	int i;
//	int j;
//
//	try
//	{
//		j = -1;
//		for (i = 1; i <= Int(n / 2); i++) {
//			j = j + 1;
//			AAA->at(i) = vector->at(j);
//			j = j + 1;
//			bbb->at(i) = vector->at(j);
//		}
//	}
//	catch (const std::exception& e)
//	{
//		TratarError("Comun::Prepara", e.what());
//		throw e;
//	}
//	catch (...)
//	{
//		TratarError("Comun::Prepara", "Error Desconocido");
//		throw;
//	}
//}

void Comun::Prepara(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n)
{
	//Este método nos lo pasó Tono 24/11/2015
	int i;
	int j;

	try
	{
		j = -1;
		for (i = 0; i <= Int(n / 2); i++) {
			j = j + 1;
			AAA->at(i) = vector->at(j);
			j = j + 1;
			bbb->at(i) = vector->at(j);
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::PreparaTono", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::PreparaTono", "Error Desconocido");
		throw;
	}
}

//void Comun::FftOld(int n, int power, std::vector<double>* r_, std::vector<double>* i_)
//{
//
//	int i;
//	int j;
//	int r;
//	int g;
//	int k;
//	double A;
//	double b;
//	double c;
//	double e;
//	double p;
//	double Q;
//	double u;
//	double v;
//
//	try
//	{
//		k = 0;
//		for (j = 1; j <= n - 1; j++) {
//			i = 2;
//			while (k >= Int(n / i)) {
//				k = k - Int(n / i);
//				i = 2 * i;
//			}
//			k = k + Int(n / i);
//			if (k > j) {
//				A = r_->at(j + 1);
//				r_->at(j + 1) = r_->at(k + 1);
//				r_->at(k + 1) = A;
//				A = i_->at(j + 1);
//				i_->at(j + 1) = i_->at(k + 1);
//				i_->at(k + 1) = A;
//			}
//		}
//		g = 1;
//		p = 1;
//
//
//		for (i = 1; i <= power - 1; i++) {
//			c = 1;
//			e = 0;
//			Q = sqrt((1 - p) / 2);
//			if (i == 1) {
//				p = -1 * sqrt((1 + p) / 2);
//			}
//			else {
//				p = sqrt((1 + p) / 2);
//			}
//			for (r = 1; r <= g; r++) {
//				j = r;
//				while (j <= n) {
//					k = j + g;
//					A = c * r_->at(k) + e * i_->at(k);
//					b = e * r_->at(k) - c * i_->at(k);
//
//					r_->at(k) = r_->at(j) - A;
//					i_->at(k) = i_->at(j) + b;
//					r_->at(j) = r_->at(j) + A;
//					i_->at(j) = i_->at(j) - b;
//					//fprintf(fp, "Iteracion: %d %d %d %f %f %f %f %f %f\n", i, j, k, A, b, r_->at(k), r_->at(j), i_->at(k), i_->at(j));
//					j = j + 2 * g;
//				}
//				A = e * p + c * Q;
//				c = c * p - e * Q;
//				e = A;
//			}
//			g = 2 * g;
//		}
//
//
//
//		A = Constantes::pi / n;
//		p = cos(A);
//		Q = sin(A);
//		A = r_->at(1);
//		r_->at(1) = A + i_->at(1);
//		i_->at(1) = A - i_->at(1);
//		r_->at(1) = r_->at(1) / 2;
//		i_->at(1) = i_->at(1) / 2;
//		c = 1;
//		e = 0;
//		for (j = 2; j <= Int(n / 2); j++) {
//			A = e * p + c * Q;
//			c = c * p - e * Q;
//			e = A;
//			k = n - j + 2;
//			A = r_->at(j) + r_->at(k);
//			b = (i_->at(j) + i_->at(k)) * c;
//			b = b - (r_->at(j) - r_->at(k)) * e;
//			u = i_->at(j) - i_->at(k);
//			v = (i_->at(j) + i_->at(k)) * e;
//			v = v + (r_->at(j) - r_->at(k)) * c;
//			r_->at(j) = (A + b) / 2;
//			i_->at(j) = (u - v) / 2;
//			r_->at(k) = (A - b) / 2;
//			i_->at(k) = -(u + v) / 2;
//		}
//
//		i_->at(Int(n / 2) + 1) = -1 * i_->at(Int(n / 2) + 1);
//
//		for (j = 1; j <= n; j++) {
//			r_->at(j) = r_->at(j) / n;
//			i_->at(j) = i_->at(j) / n;
//		}
//	}
//	catch (const std::exception& e)
//	{
//		TratarError("Comun::Fft", e.what());
//		throw e;
//	}
//	catch (...)
//	{
//		TratarError("Comun::Fft", "Error Desconocido");
//		throw;
//	}
//}

void Comun::Fft(int n, int pot, std::vector<double>* r_, std::vector<double>* i_)
{
	//Este código se actualizó con el código de Tono 24/11/2015--> Ya seusa para ruido y filtrado
	int i, j, r, g, k;
	double a, b, c, e, p, q, u, v;
	k = 0;
	for (j = 1; j<n; ++j)
	{
		i = 2;
		while (k >= (n / i))
		{
			k -= (n / i);
			i += i;
		}
		k += (n / i);
		if (k>j)
		{
			a = r_->at(j);
			r_->at(j) = r_->at(k);
			r_->at(k) = a;
			a = i_->at(j);
			i_->at(j) = i_->at(k);
			i_->at(k) = a;
		}
	}
	g = 1;
	p = 1.0;
	for (i = 1; i < pot; i++)
	{
		c = 1.0;
		e = 0.0;
		q = sqrt((1.0 - p) / 2.0);
		if (i == 1)
			p = -1.0*sqrt((1.0 + p) / 2.0);
		else
			p = sqrt((1.0 + p) / 2.0);

		for (r = 1; r < g + 1; r++)
		{
			j = r;
			while (j <= n)
			{
				k = j + g;
				a = c* r_->at(k - 1) + e* i_->at(k - 1);
				b = e* r_->at(k - 1) - c* i_->at(k - 1);
				r_->at(k - 1) = r_->at(j - 1) - a;
				i_->at(k - 1) = i_->at(j - 1) + b;
				r_->at(j - 1) = r_->at(j - 1) + a;
				i_->at(j - 1) = i_->at(j - 1) - b;
				j += 2 * g;
			}
			a = e*p + c*q;
			c = c*p - e*q;
			e = a;
		}
		g += g;
	}


	a = Constantes::pi / n;
	p = cos(a);
	q = sin(a);
	a = r_->at(0);
	r_->at(0) = a + i_->at(0);
	i_->at(0) = a - i_->at(0);
	r_->at(0) = r_->at(0) / 2;
	i_->at(0) = i_->at(0) / 2;
	c = 1.0;
	e = 0.0;
	for (j = 2; j <= n / 2; j++)
	{
		a = e*p + c*q;
		c = c*p - e*q;
		e = a;
		k = n - j + 2;
		a = r_->at(j - 1) + r_->at(k - 1);
		b = c*(i_->at(j - 1) + i_->at(k - 1));
		b -= e*(r_->at(j - 1) - r_->at(k - 1));
		u = i_->at(j - 1) - i_->at(k - 1);
		v = e*(i_->at(j - 1) + i_->at(k - 1));
		v += c*(r_->at(j - 1) - r_->at(k - 1));
		r_->at(j - 1) = (a + b) / 2;
		i_->at(j - 1) = (u - v) / 2;
		r_->at(k - 1) = (a - b) / 2;
		i_->at(k - 1) = -(u + v) / 2;
	}
	i_->at(n / 2) = -1.0* i_->at(n / 2);
	for (j = 1; j <= n; j++)
	{
		r_->at(j - 1) = r_->at(j - 1) / n;
		i_->at(j - 1) = i_->at(j - 1) / n;
	}

}

double Comun::Log2(double numero)
{
	//logaritmo en base 2 de numero
	return log(numero) / log(2.0);
}


double Comun::Log10(double numero)
{
	//logaritmo en base 10 de numero
	return log(numero) / log(10.0);
}

int Comun::Int(double numero)
{
	//Devuelve la parte entera de un numero
	return (int)numero;
}

int Comun::ObtenerArmonicoCorte(int cRevCiclo, double rpm, double RegimenCorte0, double RegimenCorte1, int ArmonicoCorte0, int ArmonicoCorte1)
{
	//Devuelve el armonico de corte interpolando entre los puntos (RegimenCorte0,ArmonicoCorte0)-(RegimenCorte1,ArmonicoCorte1) para el valor Regimen
	int Result;
	double m;
	double b;

	try
	{
		ArmonicoCorte1 = ArmonicoCorte1 / 2 * cRevCiclo;
		ArmonicoCorte0 = ArmonicoCorte0 / 2 * cRevCiclo;

		m = (ArmonicoCorte1 - ArmonicoCorte0) / (RegimenCorte1 - RegimenCorte0);
		b = ArmonicoCorte0 - m * RegimenCorte0;
		Result = (int)((m * rpm) + b);
		if (Result < ArmonicoCorte1)
			Result = ArmonicoCorte1;
		if (Result > ArmonicoCorte0)
			Result = ArmonicoCorte0;

		return Result;

		//ObtenerArmonicoCorte = 30
		//ObtenerArmonicoCorte = 750 / (rpm / 60) 'Para eliminar los picos debido al conducto en motor 2T de Vossloh
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::ObtenerArmonicoCorte", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::ObtenerArmonicoCorte", "Error Desconocido");
		throw;
	}
}



//el ensayo se envía por referencia para eviar sólo el puntero y consumir menos recursos
//double Comun::modCGradosCiclo(DatosCalculo* ent, double ang)
//{
//	double result;
//	//deja el angulo -360 < ang < 360 en 4T y -180 < ang < 180 en 2T

//	try
//	{

//		if (ang >= -(ent->getDatosCalcGeneral().cGradosCiclo / 2) & ang < (ent->getDatosCalcGeneral().cGradosCiclo / 2)) {
//			result = ang;
//		}
//		else {
//			if (ang >= ent->getDatosCalcGeneral().cGradosCiclo / 2) {
//				result = ang - ent->getDatosCalcGeneral().cGradosCiclo;
//			}
//			else {
//				if (ang < -(ent->getDatosCalcGeneral().cGradosCiclo / 2)) {
//					result = ang + ent->getDatosCalcGeneral().cGradosCiclo;
//				}
//			}
//		}

//		return result;
//	}
//	catch (std::exception& e)
//	{
//		Comun::TratarError("Comun::modCGradosCiclo", e.what());
//		throw e;
//	}
//}

bool Comun::calculo_PMI_PMB(int cGradosCiclo, double vel_ang, std::vector<double>* presionr, double* PMI, double* PMB, double desfase, double vcc, double d, double LM, double LB, double e, double VD,
	double inc, double intang, double factp, double facti)
{
	int NPC;
	int i;
	double wtc;
	double wta;
	double ang;
	double presion1;
	double presion2;
	double Volumen1;
	double volumen2;
	double Pm;
	double nada;
	double rel;


	try
	{
		NPC = UBound(presionr);

		wtc = 0;
		wta = 0;
		ang = modCGradosCiclo(desfase, cGradosCiclo);
		for (i = 1; i <= NPC; i++) {
			ang = modCGradosCiclo(ang, cGradosCiclo);
			presion1 = presionr->at(i);
			if (i == NPC)
				presion2 = presionr->at(1);
			else
				presion2 = presionr->at(i + 1);
			Pm = (presion1 + presion2) / 2;
			Volumen1 = Volumen_def(ang, vel_ang, vcc, inc, presion1, factp, facti, intang, d, LM, LB, e, VD);
			volumen2 = Volumen_def(ang + intang, vcc, vel_ang, inc, presion2, factp, facti, intang, d, LM, LB, e, VD);

			if (ang >= -180 & ang < 180) {
				wtc = wtc + Pm * (volumen2 - Volumen1);
			}
			else {
				wta = wta + Pm * (volumen2 - Volumen1);
			}
			ang = ang + intang;
		}
		*PMI = wtc / VD;
		*PMB = wta / VD;

		return false;

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Calculo_PMI_PMB", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Calculo_PMI_PMB", "Error Desconocido");
		throw;
	}
}



double Comun::Volumen(double angulo, double vcc, double Diametro, double mani, double Biela, double desc)
{
	//PMS=0 grados de cigüeñal se corresponden con un volumen igual a Vcc.
	//Devuelve el volumen total sin deformaciones en m3
	try
	{
		return vcc + (Constantes::pi * pow(Diametro, 2) / 4) * Dist_Despl_inst(angulo, mani, Biela, desc);
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Volumen", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Volumen", "Error Desconocido");
		throw;
	}
}


double Comun::Volumen_def(double angulo, double vel_ang, double vcc, double inc, double pre, double factp, double facti, double intang, double D, double LM, double LB, double desc, double VD)
{
	//Este procedimiento calcula la deformacion para una presion determinada y un angulo determinado
	//y devuelve el nuevo volumen de la camara con deformaciones en m3

	double deltah_def;
	double vcc_def;
	double vel;
	double vol;
	double ace;
	double ace_antes;
	double ace_despues;
	double ang; //en radianes


	try
	{

		vol = Volumen(angulo, vcc, D, LM, LB, desc);

		ang = (angulo)* Constantes::pi / 180;

		ace = Aceleracion(ang, vel_ang, LM, LB, desc);

		deltah_def = pre * factp - ace * facti;//Deformacion

		vcc_def = vcc + pow(D, 2) * Constantes::pi * deltah_def / 4;

		return (vol - vcc + vcc_def);

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Volumen_def", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Volumen_def", "Error Desconocido");
		throw;
	}
}


double Comun::Aceleracion(double ang, double vel_ang, double  LM, double LB, double desc)
{
	//Este procedimiento calcula la aceleracion del  pistón en m/s2

	double ace;
	double nominador1;
	double nominador2;
	double nominador3;
	double nominador4;
	double denominador2;
	double denominador;



	try
	{
		nominador1 = pow(LM, 2)* cos(2 * ang);
		nominador2 = desc*LM*sin(ang);
		nominador4 = pow(LM, 4)* (pow(sin(2 * ang), 2) / 4) + pow(desc, 2)* pow(LM, 2) * pow(cos(ang), 2) - 2 * desc*pow(LM, 3)* sin(ang)*pow(cos(ang), 2);
		denominador2 = pow(LB, 2) - pow((LM*sin(ang) - desc), 2);
		nominador3 = nominador4 / denominador2;

		denominador = sqrt(denominador2);

		ace = -pow(vel_ang, 2) * (LM*cos(ang) + (nominador1 + nominador2 + nominador3) / denominador);

		return ace;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Aceleracion", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Aceleracion", "Error Desconocido");
		throw;
	}
}


double Comun::DerivadaVolumen(double angulo, double  LM, double LB, double desc, double D, double factp, double facti, double der_pre, double vel_ang)
{
	//Este procedimiento calcula la aceleracion del  pistón en m/s2
	//angulo viene en grados
	double der_vol;
	double der_ace;
	double delta_alpha;
	double ace_antes;
	double ace_despues;
	double ang;

	try
	{
		ang = angulo* Constantes::pi / 180;
		delta_alpha =0.05* Constantes::pi / 180;

		ace_antes = Aceleracion(ang - delta_alpha, vel_ang, LM, LB, desc);
		ace_despues = Aceleracion(ang + delta_alpha, vel_ang, LM, LB, desc);

		der_ace = (ace_despues - ace_antes) / (2 * delta_alpha);


		der_vol = pow(Constantes::pi * D, 2) / 4 * ((LM * sin(ang) + (LM * sin(ang) - desc)*(LM*cos(ang)) / (sqrt(pow(LB, 2) - pow((LM * sin(ang) - desc), 2)))) + (factp*der_pre) + (facti*der_ace));


		return der_vol;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::DerivadaVolumen", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::DerivadaVolumen", "Error Desconocido");
		throw;
	}
}


bool Comun::Eliminar_Offset_Inyeccion(int NPC, double intang, int cRevCiclo, double cGradosCiclo, std::vector<double>* vector, double desfase_pms)
{
	//El vector puede traer la seña de comando o la tasa de inyeccion
	double acum;
	int cont;
	int pos;
	double A;
	int f;


	try
	{
		A = -desfase_pms;
		acum = 0;
		cont = 0;
		for (f = 1; f <= NPC; f++) {
			if (A > (90 * cRevCiclo) | A < -(90 * cRevCiclo)) {
				acum = acum + vector->at(f);
				cont = cont + 1;
			}
			A = A + intang;
			A = modCGradosCiclo(A, cGradosCiclo);
		}
		acum = acum / cont;
		for (f = 1; f <= NPC; f++) {
			vector->at(f) = vector->at(f) - acum;
		}

		return false;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Eliminar_Offset_Inyeccion", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Eliminar_Offset_Inyeccion", "Error Desconocido");
		throw;
	}
}

void Comun::maximo_vector(std::vector<int>* vector, double* max)
{
	try
	{

		*max = vector->at(0);
		for (int j = 0; j < vector->size(); j++) {
			if (vector->at(j) > *max) {
				*max = vector->at(j);
			}
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::maximo_vector", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::maximo_vector", "Error Desconocido");
		throw;
	}
}

void Comun::minimo_vector(std::vector<int>* vector, double* min)
{
	try
	{
		*min = vector->at(0);
		for (int j = 0; j < vector->size(); j++) {
			if (vector->at(j) < *min) {
				*min = vector->at(j);
			}
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::minimo_vector", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::minimo_vector", "Error Desconocido");
		throw;
	}
}


void Comun::maximo_vector(std::vector<double>* vector, int* indice)
{
	int j;
	double max;

	try
	{
		*indice = 0;
		max = vector->at(0);
		for (j = 0; j < vector->size(); j++) {
			if (vector->at(j) > max) {
				max = vector->at(j);
				*indice = j;
			}
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::maximo_vector", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::maximo_vector", "Error Desconocido");
		throw;
	}
}


bool Comun::Ordena_Inyecciones(std::vector<CInyInicioFinal>& inyecciones)
{
	std::vector<CInyInicioFinal> inyecciones_temp;
	CInyInicioFinal inyeccion_temp;
	CInyInicioFinal inyeccion;
	short i;
	short indice;
	double menor;

	try
	{
		do {

			menor = Constantes::MasInfinito;

			for (i = 0; i < inyecciones.size(); i++) {
				inyeccion = inyecciones.at(i);
				if (inyeccion.inicio < menor) {
					menor = inyeccion.inicio;
					indice = i;
					inyeccion_temp = inyeccion;
				}
			}

			//OJO CON EL INDICE A CONTINUACION !!!!
			inyecciones.erase(inyecciones.begin() + indice);
			inyecciones_temp.push_back(inyeccion_temp);

		} while (!(inyecciones.size() == 0));

		inyecciones = inyecciones_temp;

		return false;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Ordena_Inyecciones", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Ordena_Inyecciones", "Error Desconocido");
		throw;
	}
}

bool Comun::Calcula_Inicio_Fin_Inyeccion(int NPC, double n, double intang, double cGradosCiclo, double ADETOT, double RCA, double AAE, int NINY, std::vector<double>* ciclo_, std::vector<CInyInicioFinal>& Result)
{

	//a este procedimiento se entra con el offset corregido,
	//se dan los angulos de inicio y fin de multiples inyecciones segun
	//algoritmo proporcionado por Jean,
	//los angulos son respecto a la primera medida, de modo que luego habra que tener en cuenta
	//el desfase
	//en el vector inicio() se devuelven los inicios de las inyecciones y en
	//fin() los finales
	//Numero_Iny es el numero total de inyecciones encontradas y
	//com_max es la senal de comando maxima encontrada

	//hay dos problemas:
	//si una inyeccion empieza en el final del ciclo y acaba en el principio esto casca
	//las inyecciones no estan ordenadas, pero luego se ordenan

	CInyInicioFinal injection;
	int p1;
	int j;
	int pos;
	int p2;
	double bajo;
	double alto;
	int punto1;
	int punto2;
	double A;
	double b;
	double valor;
	int indice;
	int cont_busqueda;
	double COM;
	double fin_p;
	double ini_p;
	double tiempo_p;
	std::vector<double> CICLO;
	double COMMAX;
	short numero_iny;
	std::vector<CInyInicioFinal> res;
	double k;
	//ReDim inicio(numero_total_inyecciones) As Double
	//ReDim fin(numero_total_inyecciones) As Double

	try
	{
		CICLO.resize(NPC + 1);

		COMMAX = ciclo_->at(1);
		for (j = 1; j <= NPC; j++) {
			CICLO.at(j) = ciclo_->at(j);
			if (ciclo_->at(j) > COMMAX)
				COMMAX = ciclo_->at(j);
		}
		//pasamos de los puntos negativos, se supone que hemos corregido el offset
		j = 0;
		numero_iny = 0;
		cont_busqueda = 0;

		do {
			//Encontramos el punto maximo del vector en indice
			maximo_vector(&CICLO, &indice);
			COM = ciclo_->at(indice);
			if ((numero_iny >= 1) && (COM < COMMAX * Constantes::CPorcentajeAltoInyeccionValida / 100))
				break;
			bajo = COM * Constantes::CPorcentajeBajoInyeccionValida / 100;
			alto = COM * Constantes::CPorcentajeAltoInyeccionValida / 100;
			//buscamos hacia la izquierda
			j = indice;
			do {
				j = modX(j - 1, NPC);
				valor = CICLO.at(j);
				cont_busqueda = cont_busqueda + 1;
			} while (!((valor < bajo) || (cont_busqueda == NPC)));

			if (cont_busqueda == NPC)
			{
				exAplic.ColocarTexto("Error en Calcula_inicio_fin_inyeccion");
				throw exAplic;
			}


			//return;
			punto1 = j;

			cont_busqueda = 0;
			do {
				j = modX(j + 1, NPC);
				valor = CICLO.at(j);
				cont_busqueda = cont_busqueda + 1;
			} while (!((valor > alto) || (cont_busqueda == NPC)));

			if (cont_busqueda == NPC)
			{
				exAplic.ColocarTexto("Error en Calcula_inicio_fin_inyeccion");
				throw exAplic;
			}

			//goto Trataerror;
			//return;
			punto2 = j;

			if (punto1 > punto2) {
				punto1 = punto1 - NPC;
			}

			A = (CICLO.at(modX(punto2, NPC)) - CICLO.at(modX(punto1, NPC))) / (punto2 - punto1);
			b = -A * punto1 + CICLO.at(modX(punto1, NPC));
			k = modX((-b / A - 1), NPC);
			ini_p = k * intang;

			//buscamos hacia la derecha
			j = indice;
			cont_busqueda = 0;
			do {
				j = modX(j + 1, NPC);
				valor = CICLO.at(j);
				cont_busqueda = cont_busqueda + 1;
			} while (!((valor < bajo) || (cont_busqueda == NPC)));

			if (cont_busqueda == NPC)
			{
				exAplic.ColocarTexto("Error en Calcula_inicio_fin_inyeccion");
				throw exAplic;
			}

			//goto Trataerror;
			//return;
			punto1 = j;

			cont_busqueda = 0;
			do {
				j = modX(j - 1, NPC);
				valor = CICLO.at(j);
				cont_busqueda = cont_busqueda + 1;
			} while (!((valor > alto) || (cont_busqueda == NPC)));

			if (cont_busqueda == NPC)
			{
				exAplic.ColocarTexto("Error en Calcula_inicio_fin_inyeccion");
				throw exAplic;
			}

			//goto Trataerror;
			//return;
			punto2 = j;

			if (punto1 < punto2) {
				punto2 = punto2 - NPC;
			}

			A = (CICLO.at(modX(punto2, NPC)) - CICLO.at(modX(punto1, NPC))) / (punto2 - punto1);
			b = -A * punto1 + CICLO.at(modX(punto1, NPC));
			k = modX((-b / A - 1), NPC);
			fin_p = k * intang;

			if (fin_p > ini_p) {
				tiempo_p = (fin_p - ini_p) / (n * 360);
			}
			else {
				tiempo_p = (cGradosCiclo - (ini_p - fin_p)) / (n * 360);
			}

			//si la inyeccion encontrada dura menos de Tiempo_Min_Iny  es que es falsa
			if (tiempo_p > Constantes::CTiempoMinimoInyeccionValida) {
				numero_iny = numero_iny + 1;
				injection.inicio = modCGradosCiclo(ini_p - ADETOT, cGradosCiclo);
				injection.fin = modCGradosCiclo(fin_p - ADETOT, cGradosCiclo);

				//Esto solamente es necesario hacerlo aquí si la señal de inyección es de tasa (no es el setting).
				//Si es el setting lo rehace al calcular la tasa, y si se está llamando desde el estudio estadístico no se usa.
				if ((injection.inicio < RCA) || (injection.inicio > AAE)) {
					injection.inicio_ciclo_cerrado = false;
				}
				else {
					injection.inicio_ciclo_cerrado = true;
				}

				if ((injection.fin < RCA) || (injection.fin > AAE)) {
					injection.fin_ciclo_cerrado = false;
				}
				else {
					injection.fin_ciclo_cerrado = true;
				}
				res.push_back(injection);

				if (numero_iny == 1)
					COMMAX = COM;
			}
			//eliminamos el trozo, tanto si ha habido exito como si no



			p1 = modX(Int(ini_p / intang) + 1, NPC);
			p2 = modX(Int(fin_p / intang) + 1, NPC);

			if (p2 == 0)
			{
				p2 = 1;
			}

			j = p1;
			while (!(j == p2)) {
				CICLO.at(j) = 0;
				j = modX(j + 1, NPC);
			}
			//hasta acabar con el numero de inyecciones maximo
		} while (!((numero_iny + 1 > NINY)));

		Ordena_Inyecciones(res);

		Result = res;

		//sal_bucle:
		return false;

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Calcula_inicio_fin_inyeccion", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Calcula_inicio_fin_inyeccion", "Error Desconocido");
		throw;
	}
}

double Comun::modang(double ang, int cGradosCiclo) //Esta funcion sirve para que el angulo quede entre -360 y +360
{
	//esta función sirve para que un ángulo de 370 sea igual a 10
	//y un ángulo de 350 sea igual a -10
	try
	{
		if (ang < 180 & ang > -180) {
			return ang;
		}

		ang = ang + 360;

		return modCGradosCiclo(ang, cGradosCiclo);
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::modang", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::modang", "Error Desconocido");
		throw;
	}
}


double Comun::modCGradosCiclo(double ang, int cGradosCiclo)
{
	//deja el angulo -360 < ang < 360 en 4T y -180 < ang < 180 en 2T

	try
	{
		if (ang >= -(cGradosCiclo / 2) & ang < (cGradosCiclo / 2)) {
			return ang;
		}

		if (ang >= cGradosCiclo / 2) {
			return ang - cGradosCiclo;
		}

		if (ang < -(cGradosCiclo / 2)) {
			return ang + cGradosCiclo;
		}

		return 0;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::modCGradosCiclo", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::modCGradosCiclo", "Error Desconocido");
		throw;
	}
}

bool Comun::is_digits(const std::string &str)
{
	return str.find_first_not_of("0123456789.,") == std::string::npos;
}

std::string Comun::ToUpper(std::string str)
{
	for (int pos = 0, sz = str.length(); pos < sz; ++pos)
	{
		if (str[pos] >= 'a' && str[pos] <= 'z') { str[pos] += ('A' - 'a'); }
	}

	return str;
}

//Devuelve la extension de un fichero
std::string Comun::ExtFichero(std::string fichero)
{
	std::size_t pos;
	std::string ext;

	pos = fichero.find_last_of(".");      // obtengo la posicion del ultimo punto

	ext = fichero.substr(pos + 1);     // devuelvo del punto al final

	return ToUpper(ext);
}

//void Comun::Prepara2(std::vector<double>* vector, std::vector<double>* AAA, std::vector<double>* bbb, int n1, int n)
//{
//	std::vector<double> *normalizado;
//	short NO;
//	short j;
//	short i;
//
//	normalizado->resize(n1);
//
//	vector->at(0) = vector->at(n);
//	Normaliza(vector, normalizado, n, n1);
//
//	j = -1;
//	//prepara los valores de entrada para la fft
//	for (i = 1; i <= n1 / 2; i++) {
//		j = j + 1;
//		AAA->at(i) = normalizado->at(j);
//		j = j + 1;
//		bbb->at(i) = normalizado->at(j);
//	}
//
//}

void Comun::Regresion(std::vector<double>* suma, std::vector<double> *Fc, double *b0, double *b1)
{
	short inicio;
	short final;
	double Sxy;
	double Sx;
	double Sx2;
	double Sy;
	short Sn;
	short i;

	i = 0;

	do {
		if (Fc->at(i) > 0.1)
			break; // TODO: might not be correct. Was : Exit Do
		i = i + 1;
	} while (true);
	inicio = i;

	do {
		if (Fc->at(i) > 1)
			break; // TODO: might not be correct. Was : Exit Do
		i = i + 1;
	} while (true);
	final = i;

	Sx = 0;
	Sx2 = 0;
	Sxy = 0;
	Sy = 0;
	Sn = 0;

	std::vector<double> x;
	std::vector<double> Y;

	x.resize(final - inicio + 1);
	Y.resize(final - inicio + 1);


	//asigna valores a las matrices para hacer la regresion
	for (i = 0; i <= final - inicio; i++) {
		x.at(i) = Log10(Fc->at(i + inicio));
		Y.at(i) = 10 * Log10(suma->at(i + inicio) / 4E-10);
	}

	//calcula los coeficientes de minimos cuadrados
	for (i = 0; i <= final - inicio; i++) {
		Sx = Sx + (x.at(i));
		Sx2 = Sx2 + (pow(x.at(i), 2));
		Sxy = Sxy + (x.at(i) * Y.at(i));
		Sy = Sy + (Y.at(i));
		Sn = Sn + 1;
	}

	//calcula la pendiente y el nivel de ruido
	*b0 = ((Sy * Sx2) - (Sx * Sxy)) / ((Sn * Sx2) - (pow((Sx), 2)));
	*b1 = ((Sn * Sxy) - (Sx * Sy)) / ((Sn * Sx2) - (pow((Sx), 2)));

}

double Comun::Interpola(double x, double x1, double x2, double y1, double y2)
{
	//este procedimiento busca la ordenada (Y)  para una abcisa dada (x)
	//mediante interpolacion lineal
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

void Comun::Interpola(int n, int nf, int grado, std::vector<double> *a1, std::vector<double> *v1, std::vector<double> *a2, std::vector<double> *v2)
{
	int i, k;
	double inicio, intervalo;

	try
	{
		inicio = a1->at(0);
		intervalo = (a1->at(nf) - inicio) / (n - 1);
		for (i = 0; i < n; i++)
			a2->at(i) = inicio + intervalo * (double)i;
		k = 0;
		v2->at(0) = v1->at(0);
		v2->at(n - 1) = v1->at(nf);
		for (i = 1; i < n - 1; i++) {
			while ((a1->at(k) < a2->at(i)) && (k < nf))
				k++;
			v2->at(i) = Polinomio(i, grado, k, nf, a1, v1, a2);
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Interpola", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Interpola", "Error Desconocido");
		throw;
	}
}


//Interpola linealmente para encontrar la posición x en la que y = 0
double Comun::BuscaCero(double x1, double x2, double y1, double y2)
{
	double interp = 0;
	try
	{
		interp = x1 - y1 * (x2 - x1) / (y2 - y1);
		return interp;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::BuscaCero", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::BuscaCero", "Error Desconocido");
		throw;
	}
}

double Comun::Polinomio(int ii, int grado, int kk, int nf, std::vector<double> *a1, std::vector<double> *v1, std::vector<double> *a2)
{
	int i, k, n, na, nb;
	std::vector<double> d;
	std::vector<double> x;
	double b;

	d.resize(8);
	x.resize(8);

	/* nb=kk-1;
	na=kk-2; */
	nb = kk;
	na = kk - 1;
	k = 2;
	while (k < grado) {
		if (nb<nf) {
			nb++;
			k++;
		}
		if (na>0) {
			na--;
			k++;
		}
	}
	n = nb - na;
	for (i = na; i <= nb; i++) {
		d.at(i - na) = v1->at(i);
		x.at(i - na) = a1->at(i);
	}
	for (k = 1; k <= n; k++) {
		for (i = 0; i <= (n - k); i++)
			d.at(i) = (d.at(i + 1) - d.at(i)) / (x.at(i + k) - x.at(i));
	}
	b = d.at(0);
	for (i = 1; i <= n; i++)
		b = d.at(i) + (a2->at(ii) - x.at(i)) * b;
	return(b);
}

void Comun::ffton(int nf, double regimen, std::vector<double>* a1, std::vector<double> *v1, std::vector<double> *frec, std::vector<double> *esp, int cRevCiclo)
{
	//Este código de Tono prepara el vector de presion (normaliza a 2k) y calcula la fft para el cálculo de ruido
	bool pasa;
	int n, m, grado;
	int i, pote, test, limite;
	double dumb;
	std::vector<double> a2;
	std::vector<double> v2;
	std::vector<double> gim;
	std::vector<double> gre;

	a2.resize(KMAX);
	v2.resize(KMAX);
	gim.resize(KMIN);//Revisar: que pasa si nf>>KMIN
	gre.resize(KMIN);

	test = nf;
	pote = 1;
	while (test > 1)
	{
		test /= 2;
		pote++;
	}
	test = 1;
	for (i = 1; i<pote + 1; i++)
		test *= 2;
	if (test != nf + 1) {
		pasa = true;
		test = nf;
		m = test;
		n = 2;
		grado = 1;
		while (m>1) {
			m /= 2;
			n *= 2;
		}
		Interpola(n, test, grado, a1, v1, &a2, &v2);//Revisar en le filtrado porque se ha mantenido Nomalizar del calmec viejo 
		Prepara(&v2, &gre, &gim, n);
	}
	else {
		pasa = false;
		Prepara(v1, &gre, &gim, test);
	}

	if (pasa) limite = n;
	else limite = test;

	Fft(limite / 2, pote, &gre, &gim);

	FILE *fp;
	//FILE*fp2;
	fp = fopen((getRutaSalida() + "\\Frecuencias" + _id_ensayo + ".dat").c_str(), "w");

	/*fp2 = fopen("Ciclo.txt", "w");

	for (i = 0; i <= nf; i++)
	{
	fprintf(fp2, "%f %f\n", a1->at(i), v1->at(i));
	}*/

	for (i = 0; i <= limite / 2; i++) {
		frec->at(i) = (double)i * regimen / cRevCiclo;
		/*	*(esp+i) = 20.0*log10(sqrt(*(gre+i)* *(gre+i) + *(gim+i)* *(gim+i))/0.00002); */
		//En el codigo de tono el primer numero era 193.9794=10*log(1/0.0000000002^2) ya que se consideraba la presion en bares
		esp->at(i) = 93.9794 + 10.0*log10(gre.at(i)*gre.at(i) + gim.at(i)* gim.at(i));
		//fprintf(fp, "%f %f\n", frec->at(i), esp->at(i));
		//fre(i) gre(i) gim(i) sqrt(gre(i)^2 + gim(i)^2) log10(fre(i)) log10(sqrt(gre(i)^2 + gim(i)^2))
		//fprintf(fp, "%f %f %f %f %f %f\n", frec->at(i), gre.at(i), gim.at(i), sqrt(pow(gre.at(i), 2) + pow(gim.at(i), 2)), log10(frec->at(i)), log10(sqrt(pow(gre.at(i), 2) + pow(gim.at(i), 2))));
		fprintf(fp, "%f %f\n", frec->at(i), esp->at(i));
	}
	fclose(fp); //Cierro el fichero
	//fclose(fp2);
}

void Comun::EliminarPosicionCeroVector(std::vector<double> *vector)
{
	for (int i = 0; i < vector->size() - 1; i++)
	{
		vector->at(i) = vector->at(i + 1);
	}

	vector->pop_back(); //Quito la ultima posicion
}

void Comun::ColocarPosicionCeroVector(std::vector<double> *vector)
{
	double valor;

	vector->push_back(valor); //Coloco otra posicion

	for (int i = UBound(vector) - 1; i >= 0; i--)
	{
		vector->at(i + 1) = vector->at(i);
	}

	vector->at(0) = 0.;
}

int Comun::Obtener_Posicion_en_Vector(double ang, double desfase, int NPC, double intang)
{
	short pos;
	pos = (desfase + ang) / intang + 1;
	pos = modX(pos, NPC);

	return pos;
}

void Comun::Calcula_Modulo(std::vector<double>* vector, std::vector<double>* modulo, int n_puntos)
{
	//Este procedimiento pasa la señal de entrada (vector()) al dominio de la
	//frecuencia, y calcula el modulo

	// vector()  ---->  señal de entrada
	// n_puntos  ---->  nº de puntos del vector de entrada y de salida
	// modulo ----> vector que contiene los modulos
	//
	// NOTA: es necesario que n_puntos sea potencia de dos


	int cont;

	std::vector<double> parte_real;
	std::vector<double> parte_imag;

	parte_real.resize(n_puntos / 2 + 1);
	parte_imag.resize(n_puntos / 2 + 1);

	FFT_Main(vector, &parte_real, &parte_imag, n_puntos - 1);


	for (cont = 1; cont <= Int(n_puntos / 2); cont++) {
		modulo->at(cont) = sqrt(pow(parte_real.at(cont), 2) + pow(parte_imag.at(cont), 2));
	}

}

double Comun::trunc_doub(double val, int digits)
{
	double temp = 0.0;

	temp = (int)(val * pow(10, digits));

	temp /= pow(10, digits);

	return temp;
}

double Comun::round_doub(double value, int digits)
{
	double factor = pow(10.0, digits - ceil(log10(fabs(value))));
	return round(value * factor) / factor;
}

double Comun::Interpolar_valor(double valor_x, std::vector<Comun::Alfa_Valor>* vector, bool interpolar_fuera_intervalo)
{
	//este procedimiento busca (con orden logaritmico)
	//la abcisa de una senyal dada su ordenada
	//se supone la senyal almacenada en una matriz de dos
	//dimensiones con una ordenada de mayor a menor
	//ademas interpola linealmente el resultado
	//Se indica en interpolar_fuera_intervalo si se realiza una 
	//iterpolacion de una valor de x que no esta entre los 
	//valores del vector

	int der;
	int izq;
	int medio;
	double valor_y;
	double pend; //Pendiente de la recta de los primeros dos puntos o ultimos dos puntos

	try
	{
		//comprueba que no esta fuera de los limites
		if (valor_x <= vector->at(1).eje_x) {
			if (interpolar_fuera_intervalo)
			{
				pend = (vector->at(2).eje_y - vector->at(1).eje_y) / (vector->at(2).eje_x - vector->at(1).eje_x);

				return vector->at(1).eje_y + pend *(valor_x - vector->at(1).eje_x);
			}
			else
			{
				return vector->at(1).eje_y;
			}

		}
		if (valor_x >= vector->at(UBound(vector)).eje_x) {
			if (interpolar_fuera_intervalo)
			{
				pend = (vector->at(UBound(vector)).eje_y - vector->at(UBound(vector) - 1).eje_y) / (vector->at(UBound(vector)).eje_x - vector->at(UBound(vector) - 1).eje_x);

				return vector->at(UBound(vector)).eje_y + pend *(valor_x - vector->at(UBound(vector)).eje_x);
			}
			else
			{
				return vector->at(UBound(vector)).eje_y;
			}
		}

		izq = 1;
		der = UBound(vector);
		medio = Int((izq + der) / 2);

		do {
			if (vector->at(medio).eje_x > valor_x) {
				der = medio;
			}
			else {
				izq = medio;
			}
			medio = Int((izq + der) / 2);
		} while (!((izq == der - 1) | (izq == der)));

		if (izq == der) {
			valor_y = Interpola(valor_x, vector->at(der).eje_x, vector->at(1).eje_x, vector->at(der).eje_y, vector->at(1).eje_y);
		}
		else {
			valor_y = Interpola(valor_x, vector->at(medio).eje_x, vector->at(medio + 1).eje_x, vector->at(medio).eje_y, vector->at(medio + 1).eje_y);
		}

		return valor_y;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Interpolar_valor", e.what());
		throw e;
		return true;
	}
	catch (...)
	{
		TratarError("Comun::Interpolar_valor", "Error Desconocido");
		throw;
	}

}

void Comun::Referencia_en_Y(double intang, int NPC, int cRevCiclo, std::vector<double> *presionf, std::vector<double> *presionr, double valor, double RCA, double ADETOT, double *premax)
{
	//Este procedimiento referencia un vector de presión en Y, pone en el ángulo correspondiente al pmi la presión que se especifica
	//Si RCA es distinto del valor nulo es porque se va ajustar la presión con respecto al ángulo del RCA, sino se impone al PMI
	double sube;
	double prepos;
	int pos;
	short f;

	try
	{
		if (RCA == Constantes::ValorNulo)
			pos = (int)(ADETOT / intang) - (int)(NPC / (2 * cRevCiclo)); // PMI
		else
			pos = (int)(ADETOT / intang) + (int)(RCA / intang); //RCA

		pos = modX(pos, NPC);
		prepos = presionf->at(pos);

		//if (prepos > valor)
		//	sube = prepos - valor;
		//else 
		sube = valor - prepos;

		presionr->resize(NPC + 1);
		*premax = presionf->at(1);

		for (f = 0; f <= NPC; f++)
		{
			presionr->at(f) = presionf->at(f) + sube;

			if (presionr->at(f) > *premax)
			{
				*premax = presionr->at(f);
			}
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Referencia_en_Y", e.what());
		throw e;
		//return true;
	}
	catch (...)
	{
		TratarError("Comun::Referencia_en_Y", "Error Desconocido");
		throw;
	}

}

void Comun::C_v(std::string tipo_comb, std::string tipo_comb1, std::string tipo_comb2, double YFUEL1, double YFUEL2, double Xq1_BLEND, double PMC, double PMC1, double T, double* cv_a, double* cv_f, double* cv_q)
{
	//**********************************************************************
	// Cálculo del Calor específico a volumen constante para aire, combustible
	// y quemados
	// Llamado desde CAjusteCombustión.AjustaPresiónPmi
	// desde CCiclo.Ejecutar_Ciclo y desde comunes.Calcula_Cbb
	//**********************************************************************'
	double cv_f_gasoil;
	double cv_q_gasoil;
	double cv_f_gasolina;
	double cv_q_gasolina;
	double cv_f1;
	double cv_q1;
	double cv_f2;
	double cv_q2;
	double Yq1_BLEND; //fracción másica de los productos estequiométricos del combustible 1 en los productos estequiométricos de la combustión de la mezcla

	//Cálculo del Calor específico a volumen constante para aire, combustible y quemados
	try
	{

		*cv_a = -10.4199 * pow(T, 0.5) + 2522.88 - 67227.1 * pow(T, (-0.5)) + 917124.4 * pow(T, (-1)) - 4174853.6 * pow(T, (-1.5));

		//GASOIL
		cv_f_gasoil = -256.4 + 6.95372 * T - 0.00404715 * pow(T, 2) + 9.10259E-07 * pow(T, 3) + 1458487 * pow(T, (-2));
		cv_q_gasoil = 641.154 + 0.43045 * T - 0.0001125 * pow(T, 2) + 8.979E-09 * pow(T, 3);

		//GASOLINA
		if (T < 1000) {
			cv_f_gasolina = -12291002.2 * pow(T, (-2)) + 227580.041 * pow(T, (-1)) - 1545.51269 + 10.8382363 * T - 0.00837844 * pow(T, 2) + 3.2557E-06 * pow(T, 3) - 4.0429E-10 * pow(T, 4) - 72.78128689;
		}
		else {
			cv_f_gasolina = 984559799 * pow(T, (-2)) - 3394060.95 * pow(T, (-1)) + 5673.52925 + 1.036209 * T - 0.00036926 * pow(T, 2) + 5.2754E-08 * pow(T, 3) - 2.7797E-12 * pow(T, 4) - 72.78128689 + 0.0312498;
		}
		cv_q_gasolina = 655.1103865 + 0.401687993 * T - 6.98091E-05 * pow(T, 2) - 9.31055E-09 * pow(T, 3) + 2.53355E-12 *pow(T, 4);

		if (tipo_comb == "FUEL_BLEND") {
			if (tipo_comb1 == "GASOIL") {
				cv_f1 = cv_f_gasoil;
				cv_q1 = cv_q_gasoil;
			}
			else if (tipo_comb1 == "GASOLINA") {
				cv_f1 = cv_f_gasolina;
				cv_q1 = cv_q_gasolina;
			}
			if (tipo_comb2 == "GASOIL") {
				cv_f2 = cv_f_gasoil;
				cv_q2 = cv_q_gasoil;
			}
			else if (tipo_comb2 == "GASOLINA") {
				cv_f2 = cv_f_gasolina;
				cv_q2 = cv_q_gasolina;
			}
			*cv_f = YFUEL1 * cv_f1 + YFUEL2 * cv_f2;
			Yq1_BLEND = Xq1_BLEND * PMC1 / PMC;
			*cv_q = Xq1_BLEND * cv_q1 + Xq1_BLEND * cv_q2;
		}
		else {
			if (tipo_comb == "GASOIL") {
				*cv_f = cv_f_gasoil;
				*cv_q = cv_q_gasoil;
			}
			else if (tipo_comb == "GASOLINA") {
				*cv_f = cv_f_gasolina;
				*cv_q = cv_q_gasolina;
			}
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::C_v", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::C_v", "Error Desconocido");
		throw;
	}

}

//******************************************************************
// Cálculo de la Energía interna del aire, combustible y quemados
// para una temperatura dada
//******************************************************************
//
// LLamado desde CCiclo.Ejecutar_Ciclo y desde
// CAjusteCombustion.Ajusta_PPMI_y_MRCA
//
void Comun::En_int(std::string tipo_comb, std::string tipo_comb1, std::string tipo_comb2, double YFUEL1, double YFUEL2, double Xq1_BLEND, double PMC, double PMC1, double T, double *u_a, double *u_f, double *u_q)
{
	double u_f_gasoil;
	double u_q_gasoil;

	double u_f_gasolina;
	double u_q_gasolina;

	try
	{
		*u_a = -4193697.9 - 6.94661 * pow(T, (1.5)) + 2522.881 * T - 134454.16 * pow(T, 0.5) + 917124.39 * log(T) + 8349707.14 * pow(T, (-0.5));

		//GASOIL
		u_f_gasoil = -1234157.8 - 256.4 * T + 3.47686 * pow(T, 2) - 0.00134905 * pow(T, 3) + 2.27565E-07 * pow(T, 4) - 1458487 * pow(T, (-1));
		u_q_gasoil = -3251495 + 1028.75 * T - 0.15377 * pow(T, 2) + 6.7895E-05 * pow(T, 3);

		//GASOLINA
		if (T < 1000) {
			u_f_gasolina = (-1 * 3230332.559) + 12291002.16 * pow(T, (-1)) + 227580.0409 * log(T) - 1618.293972 * T + 5.41911816 * pow(T, 2) - 0.002792812 * pow(T, 3) + 8.13916E-07 * pow(T, 4) - 8.08583E-11 * pow(T, 5);
		}
		else {
			u_f_gasolina = (18516600.1 + 72.34589655) - 984559798.9 * pow(T, (-1)) - 3394060.946 * log(T) + 5600.74796 * T + 0.5181045 * pow(T, 2) - 0.00012309 * pow(T, 3) + 1.31188E-08 * pow(T, 4) - 5.5593E-13 * pow(T, 5);
		}
		u_q_gasolina = -3205539.609 + 660.7221242 * T + 0.290764161 * pow(T, 2) - 0.000133742 * pow(T, 3) + 3.25919E-08 * pow(T, 4);

		double u_f1;
		double u_q1;
		double u_f2;
		double u_q2;
		double Yq1_BLEND;
		if (tipo_comb == "FUEL_BLEND") {
			if (tipo_comb1 == "GASOIL") {
				u_f1 = u_f_gasoil;
				u_q1 = u_q_gasoil;
			}
			else if (tipo_comb1 == "GASOLINA") {
				u_f1 = u_f_gasolina;
				u_q1 = u_q_gasolina;
			}
			if (tipo_comb2 == "GASOIL") {
				u_f2 = u_f_gasoil;
				u_q2 = u_q_gasoil;
			}
			else if (tipo_comb2 == "GASOLINA") {
				u_f2 = u_f_gasolina;
				u_q2 = u_q_gasolina;
			}
			*u_f = YFUEL1 * u_f1 + YFUEL2 * u_f2;
			Yq1_BLEND = Xq1_BLEND * PMC1 / PMC;
			//Yq1_BLEND = ens.ConfiguracionFluidos.Xq1_BLEND * ens.ConfiguracionFluidos.PMC / ens.ConfiguracionFluidos.PMC_BLEND
			*u_q = Xq1_BLEND * u_q1 + Xq1_BLEND * u_q2;
		}
		else {
			if (tipo_comb == "GASOIL") {
				*u_f = u_f_gasoil;
				*u_q = u_q_gasoil;
			}
			else if (tipo_comb == "GASOLINA") {
				*u_f = u_f_gasolina;
				*u_q = u_q_gasolina;
			}
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::En_int", e.what());
	}
	catch (...)
	{
		TratarError("Comun::En_int", "Error Desconocido");
		throw;
	}
}

//Calculo de la entalpia de combustible liquido
double Comun::CalculoHCL(double CAHFL, double CBHFL, double TF, double C_TIY)
{
	return CAHFL + CBHFL * (TF + C_TIY);
}

int Comun::Busqueda_Dico(double valor_x, double *valor_y, std::vector<Comun::Alfa_Valor>* vector)
{
	//este procedimiento busca (con orden logaritmico)
	//la abcisa de una senyal dada su ordenada
	//se supone la senyal almacenada en una matriz de dos
	//dimensiones con una ordenada de mayor a menor
	//ademas interpola linealmente el resultado

	int der;
	int izq;
	int medio;

	try
	{
		izq = 1;
		der = UBound(vector);
		medio = Int((izq + der) / 2);

		do {
			if (vector->at(medio).eje_x > (valor_x)) {
				der = medio;
			}
			else {
				izq = medio;
			}
			medio = Int((izq + der) / 2);
		} while (!((izq == der - 1) | (izq == der)));

		if (izq == der) {
			*valor_y = Interpola(valor_x, vector->at(der).eje_x, vector->at(1).eje_x, vector->at(der).eje_y, vector->at(1).eje_y);
		}
		else {
			*valor_y = Interpola(valor_x, vector->at(medio).eje_x, vector->at(medio + 1).eje_x, vector->at(medio).eje_y, vector->at(medio + 1).eje_y);
		}

		return 0;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Busqueda_Dico", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Busqueda_Dico", "Error Desconocido");
		throw;
	}
}

void Comun::Busqueda_Dico(double valor_x, double *valor_y, std::vector<double>* vector_x, std::vector<double>* vector_y)
{
	//este procedimiento busca (con orden logaritmico)
	//la abcisa de una senyal dada su ordenada
	//se supone la senyal almacenada en una matriz de dos
	//dimensiones con una ordenada de mayor a menor
	//ademas interpola linealmente el resultado

	int der;
	int izq;
	int medio;

	try
	{
		izq = 1;
		der = UBound(vector_x);
		medio = Int((izq + der) / 2);

		do {
			if (vector_x->at(medio) > (valor_x)) {
				der = medio;
			}
			else {
				izq = medio;
			}
			medio = Int((izq + der) / 2);
		} while (!((izq == der - 1) | (izq == der)));

		if (izq == der) {
			*valor_y = Interpola(valor_x, vector_x->at(der), vector_x->at(1), vector_y->at(der), vector_y->at(1));
		}
		else {
			*valor_y = Interpola(valor_x, vector_x->at(medio), vector_x->at(medio + 1), vector_y->at(medio), vector_y->at(medio + 1));
		}
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Busqueda_Dico", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Busqueda_Dico", "Error Desconocido");
		throw;
	}
}



double Comun::ObtenerRcDinamica(double angulo, double vcc, double vol_def, double d, double LM, double LB, double e, double VD)
{
	//Esta funcion calcula la rc dinamica en funcion del volumen sin deformar y el volumen deformado

	double vcc_def;
	double vol;

	try
	{
		//calcula el volumen actual sin deformaciones
		vol = Volumen(angulo, vcc, d, LM, LB, e);

		vcc_def = vol_def - vol + vcc;
		return (VD + vcc_def) / vcc_def;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::ObtenerRcDinamica", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::ObtenerRcDinamica", "Error Desconocido");
		throw;
	}
}

void Comun::Deriva_Señal(std::vector<double> *datos, int num_puntos, int armonico_corte, int ancho_corte, int cRevCiclo, std::vector<double> *resultados)
{
	int l;
	int pot;
	std::vector<double> normalizado;

	try
	{
		pot = Busca_Potencia(num_puntos);

		normalizado.resize(pot + 1);

		if (armonico_corte > (pot / 2))
		{
			exAplic.ColocarTexto("armonico_corte > pot/2");
			throw exAplic;
		}

		datos->at(0) = datos->at(num_puntos);


		Normaliza(datos, &normalizado, num_puntos, pot);

		Filtra_Deriva(&normalizado, armonico_corte, ancho_corte, pot, cRevCiclo);

		Normaliza(&normalizado, resultados, pot, num_puntos);


		resultados->at(num_puntos) = resultados->at(0);
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Deriva_Señal", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Deriva_Señal", "Error Desconocido");
		throw;
	}
}

void Comun::Filtra_Deriva(std::vector<double> *vector, int corte, int ancho, int n_puntos, int cRevCiclo)
{
	//Este procedimiento pasa la señal de entrada (vector()) al dominio de la
	//frecuencia, y elimina los armónicos que hay por encima de corte si ancho es cero,
	//si ancho es mayor que cero lo que hace es eliminar los armónicos pero de modo
	//más suave además deriva esta señal resultante
	//el filtro se comporta como un PASO BAJO
	//
	// vector()  ---->  señal de entrada sin filtrar y de salida filtrada y derivada
	// corte     ---->  armónico a partir del cual queremos filtrar
	// ancho     ---->  ancho (en armónicos) del filtro "suave"
	// n_puntos  ---->  nº de puntos del vector de entrada y de salida
	//
	// NOTA: es necesario que n_puntos sea potencia de dos

	int cont;
	double i;
	double rel;
	double factor;
	double Temp;
	std::vector<double> parte_real;
	std::vector<double> parte_imag;


	try
	{
		parte_real.resize(n_puntos / 2 + 1);
		parte_imag.resize(n_puntos / 2 + 1);


		FFT_Main(vector, &parte_real, &parte_imag, n_puntos - 1);
		//se pasa la señal al dominio de la frecuencia



		cont = 0;
		if (ancho != 0) {
			rel = Constantes::pi / ancho;
			for (i = (corte - Int(ancho / 2)); i <= (corte + Int(ancho / 2)); i++) {
				factor = cos(cont * rel) / 2 + 0.5;
				parte_real.at(i) = parte_real.at(i) * factor;
				parte_imag.at(i) = parte_imag.at(i) * factor;
				cont = cont + 1;
			}
		}

		for (i = (corte + Int(ancho / 2) + 1); i <= (Int(n_puntos / 2)); i++) {
			parte_real.at(i) = 0;
			parte_imag.at(i) = 0;
		}

		for (i = 1; i <= n_puntos / 2; i++) {
			Temp = parte_real.at(i);
			parte_real.at(i) = -1 * (i - 1) * parte_imag.at(i);
			parte_imag.at(i) = (i - 1) * Temp;
		}

		IFFT_Main(vector, &parte_real, &parte_imag, n_puntos);
		//una vez filtrada y derivada pasamos la señal al dominio del tiempo


		for (i = 0; i <= n_puntos - 1; i++) {
			vector->at(i) = vector->at(i) / cRevCiclo;
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Filtra_Deriva", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Filtra_Deriva", "Error Desconocido");
		throw;
	}
}

//Calculo de las fracciones masicas al RCA y AAE
void Comun::calcula_composicion_RCA_AAE(double MACC, double MAADM, double megr, double MFCC, double Mfa, double Mfadm, double Mainya, double MBB, double mres, double MCC,
	double FE, double *YaRCA, double *YaAAE, double *YqRCA, double *YqAAE, double *YfRCA, double *YfAAE, double *Yaadm, double *Yaesc, double *Yqadm,
	double *Yqesc, double *Yfadm, double *Yfesc)
{
	double Msum;
	double Mesc;
	double MRCA;
	double MAAE;
	double TEGR;
	double FR;
	//Dim Yaadm As Double 'fracción másica de aire en la admisión
	//Dim Yaesc As Double 'fracción másica de aire en el escape
	//Dim Yfadm As Double 'fracción másica de combustible en la admisión
	//Dim Yfesc As Double 'fracción másica de combustible en el escape
	//Dim Yqadm As Double 'fracción másica de productos quemados en la admisión
	//Dim Yqesc As Double 'fracción másica de productos quemados en el escape
	//Dim MFadm As Double 'masa de combustible inyectada en la admisión. Parte de ella puede cortocircuitarse (para MEP)
	//Dim MFa As Double 'masa de combustible inyectada en la cámara durante el ciclo abierto. Se asume que no se cortocircuita (para MEP y MEC)
	double Mesc_corr;
	//Dim Mainya As Double 'masa de aire inyectada con el sistema neumático en el ciclo abierto

	try
	{
		*Yaadm = 0.;
		*Yaesc = 0.;
		*Yqadm = 0.;
		*Yqesc = 0.;
		*Yfadm = 0.;
		*Yfesc = 0.;
		//Mainya = 0 'por el momento no se considera

		//Todo el planteamiento se hizo asumiendo MCC > 0, si es negativo basta con considerar en el cálculo que el cortocircuito es egr adicional y que no hay cortocircuito como tal.
		if (MCC < 0) {
			megr = megr - MCC;
			MCC = 0.;
		}

		MRCA = MAADM + Mainya + megr + (Mfadm + Mfa) + mres - MCC;
		Msum = MAADM + megr + Mfadm;
		MAAE = MACC + megr + MFCC + mres - MBB - MCC;
		Mesc = MACC + megr + MFCC - MBB - MCC;
		TEGR = megr / (megr + MAADM);
		Mesc_corr = MACC + megr + MFCC - MBB - MCC * TEGR;
		FR = MFCC / MACC / FE;

		*YqAAE = (MFCC + MACC + megr - MCC + MBB / 2 * (mres / MRCA - 1)) - (Mesc / Mesc_corr) * (megr - MCC * TEGR) * (1 - MBB / 2 / MRCA);
		//calculo de las composiciones en función del dosado relativo
		if (FR == 1) {
			*YqAAE = 1;
			*Yqesc = *YqAAE * Mesc / Mesc_corr;
			if (*Yqesc > 1)
				*Yqesc = 1;
			*Yqadm = *Yqesc * TEGR;
			*YqRCA = (*Yqesc * (megr - MCC * TEGR) + mres * *YqAAE) / MRCA;
			*YaAAE = 0.;
			*Yaesc = (Mesc * *YaAAE + MCC * (1 - TEGR)) / Mesc_corr;
			*Yaadm = (MAADM + megr * *Yaesc) / Msum;
			*YaRCA = (MAADM + megr * *Yaesc + mres * *YaAAE - MCC * *Yaadm) / MRCA;
			*YfRCA = 1 - *YqRCA - *YaRCA;
			*YfAAE = 0.;
			*Yfesc = (Mesc * *YfAAE + MCC * Mfadm / Msum) / Mesc_corr;
			*Yfadm = (Mfadm + megr * *Yfesc) / Msum;
		}
		else if (FR < 1) {
			*YfAAE = 0.;
			*Yfesc = (Mesc * *YfAAE + MCC * Mfadm / Msum) / Mesc_corr;
			*Yfadm = (Mfadm + megr * *Yfesc) / Msum;
			*YqAAE = (1 + 1 / FE) * (MFCC - MCC * *Yfadm) / *YqAAE;
			if (*YqAAE > 1) {
				*YqAAE = 1;
			}
			*YaAAE = 1 - *YqAAE - *YfAAE;
			*Yqesc = *YqAAE * Mesc / Mesc_corr;
			if (*Yqesc > 1)
				*Yqesc = 1;
			*YqRCA = (*Yqesc * (megr - MCC * TEGR) + mres * *YqAAE) / MRCA;
			*YfRCA = ((Mfadm + Mfa) + megr * *Yfesc - MCC * *Yfadm) / MRCA;
			*YaRCA = 1 - *YfRCA - *YqRCA;
			*Yqadm = *Yqesc * TEGR;
			*Yaadm = 1 - *Yfadm - *Yqadm;
			*Yaesc = (Mesc * *YaAAE + MCC * (1 - TEGR)) / Mesc_corr;
			*Yaesc = 1 - *Yqesc - *Yfesc;
		}
		else {
			*YaAAE = 0.;
			*Yaesc = (Mesc * *YaAAE + MCC * (1 - TEGR)) / Mesc_corr;
			*Yaadm = (MAADM + megr * *Yaesc) / Msum;
			*YaRCA = (MAADM + megr * *Yaesc + mres * *YaAAE - MCC * *Yaadm) / MRCA;
			*YqAAE = (1 + FE) * (MACC - MCC * *Yaadm - MBB / 2 * *YaRCA) / *YqAAE;
			if (*YqAAE > 1) {
				*YqAAE = 1;
			}
			*YfAAE = 1 - *YqAAE - *YaAAE;
			*Yqesc = *YqAAE * Mesc / Mesc_corr;
			if (*Yqesc > 1)
				*Yqesc = 1;
			*Yqadm = *Yqesc * TEGR;
			*YqRCA = (*Yqesc * (megr - MCC * TEGR) + mres * *YqAAE) / MRCA;
			*YfRCA = 1 - *YaRCA - *YqRCA;
			*Yfesc = (Mesc * *YfAAE + MCC * Mfadm / Msum) / Mesc_corr;
			*Yfadm = (Mfadm + megr * *Yfesc) / Msum;
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::calcula_composicion_RCA_AAE", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::calcula_composicion_RCA_AAE", "Error Desconocido");
		throw;
	}

}

double Comun::Xp(double x, double ratio_ctms)
{
	//curva propuesta por Antonio Gil
	return ratio_ctms + (1 / (pow((cosh(x / 100)), 40) + (ratio_ctms / (1 - ratio_ctms))));

}



//**********************************
//  Calcula las fugas por blow-by
//    pc=presion en el cárter
//    masabb = masa de blow by
//    cbb= coeficiente de blow by
//    cfbb= coeficiente para el calculo de seccion de fugas
//    inc = tiempo que dura un intervalo angular
// Llamado desde CCiclo.Ejecutar_Ciclo, CAjusteArrastre.siciclo,
// comunes.Calcula_Cbb y CajusteArrastre.Calcula_Expol
//***************************************'
void Comun::Calcula_Blow_By(double pre, double CBB, double temp, double gamma, double r, double *masabb, double pc, double Cfbb, double inc, double d)
{
	double ac;

	try
	{
		if ((pc / pre) <= 0.528)
			pc = pre * 0.528;
		if (pc > pre) {
			ac = sqrt((2 * gamma / (gamma - 1)) * (pow((pre / pc), (2 / gamma)) - pow((pre / pc), ((gamma + 1) / gamma)))) * inc;
		}
		else {
			ac = sqrt((2 * gamma / (gamma - 1)) * (pow((pc / pre), (2 / gamma)) - pow((pc / pre), ((gamma + 1) / gamma)))) * inc;
		}
		*masabb = CBB * Cfbb * d * pre * (1 / sqrt(r * temp)) * ac;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Calcula_Blow_By", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Calcula_Blow_By", "Error Desconocido");
		throw;
	}
}

double Comun::Obtener_PresionPMI(std::vector<double> *Presion, double desfase, double ang_pmi, double intang)
{
	int NPC;
	short pos_pmi;
	short pos;

	NPC = UBound(Presion);

	pos_pmi = Obtener_Posicion_en_Vector(ang_pmi, desfase, NPC, intang);

	return Presion->at(pos_pmi);
}

bool Comun::Impar(short valor)
{
	//simplemente impar es true si el numero es impar y si no es false
	if ((valor / 2 != Int(valor / 2)))
		return true;
	else
		return false;

}

void Comun::Maximo_Parabola(std::vector<Alfa_Valor> *origen, double intang, double *x, double *Y)
{
	short indice_ymax;
	short Puntos_defecto;
	short i;
	std::vector<double> vector_y;
	std::vector<double> vector_x;

	//calcula los puntos que se cogeran para el ajuste
	Puntos_defecto = (short)1 / intang * Constantes::CMaximoPresion_VentanaAngular;
	if (!(Impar(Puntos_defecto)))
		Puntos_defecto = Puntos_defecto + 1;

	//obtiene el indice del punto de y maxima
	double ymax;
	ymax = 0;
	indice_ymax = 0;
	for (i = 1; i <= UBound(origen); i++) {
		if ((origen->at(i).eje_y > ymax)) {
			ymax = origen->at(i).eje_y;
			indice_ymax = i;
		}
	}

	//dimensiona los vectores que contienen los puntos para el ajuste
	vector_y.resize(Puntos_defecto - 1);
	vector_x.resize(Puntos_defecto - 1);

	//obtiene los puntos para el ajuste
	short punto;
	for (i = 0; i <= Puntos_defecto - 1; i++) {
		punto = modX(indice_ymax + i - ((Puntos_defecto - 1) / 2), UBound(origen));
		vector_y.at(i) = origen->at(punto).eje_y;
		vector_x.at(i) = origen->at(punto).eje_x;
	}

	//se realiza el ajuste
	Ajuste_Regresion(&vector_y, &vector_x, Y, x, Puntos_defecto);

}

void Comun::Desfasar_Vector(std::vector<double>* vector_origen, std::vector<double>* vector_destino, int desfase, int NPC, double intang)
{
	try
	{
		//rellena el vector angulo con la referencia angular coherente con el desfase y escala con el factor fact el vector
		//delta es el valor con el que se referencia el vector, se usa para la presion de admision, escape y camisa, en las demás vendrá cero
		int n;
		int i;
		short punto;
		int j;

		n = UBound(vector_origen);

		punto = desfase / intang;

		j = 1;

		vector_destino->resize(n);

		if (punto > 0) {

		}
		else {
			punto = n - punto;
		}

		for (i = punto; i <= n; i++) {
			vector_destino->at(j) = vector_origen->at(i);
			j = j + 1;
		}

		for (i = 1; i <= punto - 1; i++) {
			vector_destino->at(j) = vector_origen->at(i);
			j = j + 1;
		}

	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Desfasar_Vector", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Desfasar_Vector", "Error Desconocido");
		throw;
	}

}

double Comun::Dist_Despl_inst(double angulo, double mani, double Biela, double desc)
{
	double ang;
	double result;
	//ang viene en grados y aquí se pasa a radianes.
	//0 grados de cigüeñal
	//se corresponden con un volumen igual a Vcc.
	//se devuelve hi que es el desplazamiento instantáneo (distancia del eje del bulón en el ángulo ang al eje del bulón en el PMS-->Lmax - L en tesis de Octavio)
	try
	{

		ang = (angulo)* Constantes::pi / 180;
		double pot1 = pow((mani + Biela), 2);
		double pot2 = pow(desc, 2);
		double pot3 = pow(Biela, 2);
		double pot4 = pow((mani * sin(ang) - desc), 2);
		result = (sqrt(pot1 - pot2) - (mani * cos(ang) + sqrt(pot3 - pot4)));
		return result;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Dist_Despl_inst", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Dist_Despl_inst", "Error Desconocido");
		throw;
	}
}


void Comun::Convierte_Inst(std::vector<double>* datos, double primer_angulo, std::vector<double> *result_datos, std::vector<double> *result_angulos, double fact, double delta, double intang)
{
	//rellena el vector angulo con la referencia angular coherente con el desfase y escala con el factor fact el vector
	//delta es el valor con el que se referencia el vector, se usa para la presion de admision, escape y camisa, en las demás vendrá cero
	int n;
	int cont;
	int i;
	double ang;

	n = UBound(datos);

	result_datos->resize(n + 1);
	result_angulos->resize(n + 1);

	ang = primer_angulo;
	cont = 0;
	for (i = 1; i <= n; i++) {
		result_angulos->at(i) = ang;
		result_datos->at(i) = datos->at(i) * fact + delta;
		cont = cont + 1;
		ang = primer_angulo + cont * intang;
	}

}

double Comun::asn(double x)
{
	return atan(x / sqrt(1 - pow(x, 2)));
}

double Comun::DistPistonCulataDef(double angulo, double vcc, double voldef, double d, double LM, double LB, double e, double h0)
{
	//Este procedimiento calcula la deformacion para una presion determinada y un angulo determinado
	//da como resultado el nuevo volumen que provoca la presion que hay en la camara
	double h_ahora;
	double deltah_def;
	double vcc_def;
	double vol;

	try
	{
		//calculo la distancia piston culata sin deformaciones
		h_ahora = Dist_Despl_inst(angulo, LM, LB, e);
		//calcula el volumen actual sin deformaciones
		vol = Volumen(angulo, vcc, d, LM, LB, e);
		//calcula la vcc con deformaciones a partir del volumen con deformaciones
		vcc_def = voldef - vol + vcc;
		deltah_def = 4 * (vcc_def - vcc) / (pow(d, 2) * Constantes::pi);
		return (h_ahora + deltah_def + h0);
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::DistPistonCulataDef", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::DistPistonCulataDef", "Error Desconocido");
		throw;
	}
}

void Comun::Calcula_Rendimientos(double Plim, double RCG, double VD, double Pa, double HP, double MFCC, double WIHP, double WILP, double WE, double Qtrans,
	double *rend_comb, double *rend_teo, double *rend_forma, double *rend_bombeo, double *rend_adiab, double *rend_meca)
{
	double VT;
	double Qlib;
	double p2;
	double alfa;
	double beta;
	double WI_teo;
	short descriptor;
	const short Cv = 716;
	const short cp = 716 + 287;
	const double gamma = 1.4;

	try
	{
		//volumen total
		VT = VD * (1 + 1 / (RCG - 1));
		//Calor liberado total
		*rend_comb = 1;
		Qlib = MFCC * HP * *rend_comb;
		p2 = Pa * pow(RCG, 1.4);
		alfa = Plim / p2;
		beta = 1 + (Qlib * 287 / (Pa * VT * pow(RCG, (gamma - 1))) - Cv * (alfa - 1)) * 1 / (cp * alfa);

		if (beta < 1)
			beta = 1;

		*rend_teo = 1 - 1 / pow(RCG, (gamma - 1)) * (alfa * pow(beta, gamma) - 1) / (alfa - 1 + gamma * alfa * (beta - 1));

		WI_teo = Qlib * *rend_teo;
		*rend_forma = WIHP / WI_teo;
		*rend_bombeo = (WIHP + WILP) / WIHP;
		*rend_meca = WE / (WIHP + WILP);
		*rend_adiab = 1 - Qtrans / 100;
	}
	catch (const std::exception& e)
	{
		TratarError("Comun::Calcula_Rendimientos", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Comun::Calcula_Rendimientos", "Error Desconocido");
		throw;
	}
}

double Comun::MediaVectorEmpiezaEn1(std::vector<double> *vector)
{
	double suma = 0.;

	for (int i = 1; i < vector->size(); i++)
	{
		suma = suma + vector->at(i);
	}

	return suma / vector->size() - 1;
}


//esta función te calcula la distancia (positiva) entre ang_ini y ang_fin teniendo en cuenta que pueden estar en diferentes zonas del ciclo
//y por lo tanto sus signos pueden no ser coherentes
double Comun::DistanciaEntreAngulos(double ang_ini, double ang_fin, double cGradosCiclo)
{
	double angulo;
	double result = 0;
	double ang_fin_aux;

	ang_fin_aux = ang_fin;

	if (ang_fin < ang_ini)
	{
		ang_fin_aux = ang_fin + cGradosCiclo;

	}

	result = ang_fin_aux - ang_ini;
	return result;
}


//*************************************************************************
//   Busqueda de los puntos G dada una Relacion de Compresión
//*************************************************************************
void Comun::Puntos_G(double rc, double *g1, double *g2, double lambda)
{
	//esto no es general, para motores con descentramiento habra que construir otro algoritmo
	double p2;
	double p1;
	double p3;
	double v2;
	double v1;
	double v3;
	double res;
	double resolucion;
	bool encontrado;
	double alfa;

	try
	{
		resolucion = Constantes::CPuntosG_Delta;
		encontrado = false;
		alfa = Constantes::CPuntosG_ValorInicial;

		do {
			p1 = (4 * pow(lambda, 2) - pow((sin(alfa)), 2));
			p2 = cos(alfa) * pow(p1, (-1 * 0.5));
			p3 = -((sin(alfa) * (4 * pow(lambda, 2) - 1)) / (pow(p1, 1.5)));
			v1 = 2 / (rc - 1) + 2 * lambda + 1 - cos(alfa) - pow(p1, 0.5);
			v2 = sin(alfa) * (1 + p2);
			v3 = cos(alfa) * (1 + p2) + sin(alfa) * p3;
			res = v3 / v2 - v2 / v1;
			if (res < 0)
				encontrado = true;
			alfa = alfa + resolucion;
		} while (!((alfa > Constantes::pi) | encontrado));
		if (alfa > Constantes::pi)
		{
			TratarError("Ajustes::Puntos_G", "alfa > pi");
			exAplic.ColocarTexto("alfa > pi");
			throw exAplic;
		}

		*g2 = (alfa - resolucion / 2) * 180 / Constantes::pi;
		*g1 = *g2 * -1;
		//g2 = g2 + 360
		//g1 = g1 + 360
	}
	catch (const std::exception& e)
	{
		TratarError("Ajustes::Puntos_G", e.what());
		throw e;
	}
	catch (...)
	{
		TratarError("Ajustes::Puntos_G", "Error Desconocido");
		throw;
	}

}


