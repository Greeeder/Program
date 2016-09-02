#include <string>

//namespace CalculoCALMEC
//{

#ifndef __CONSTANTES_H
#define __CONSTANTES_H

	class Constantes
	{
	public:
		//static const short SW_SHOW = 5;

		static const short VersionBBDDCalmecActual = 19;

		static const double pi;
		static const double e;
		static const double MasInfinito;
		static const double MenosInfinito;
		static const short RU;

		//static const motor2T_2T = 1

		//If motor2T_4T = 0 Then 'motor 4T
		//static const CRevCiclo = 2
		//static const CGradosCiclo = 720
		//ElseIf motor2T_4T = 1 Then ' 2T
		//static const CRevCiclo = 1
		//static const CGradosCiclo = 360
		//End If

		//Tipos de captadores
		static const bool CPiezoElectrico = true;
		static const bool CPiezoResistivo = false;

		//Tipos de medidas
		static const short CPresionCilindro = 1;
		static const short CPresionAdmision = 2;
		static const short CPresionEscape = 3;
		static const short CPresionLinea = 4;
		static const short CPresionCR = 5;
		static const short CLevantamientoAguja = 6;
		static const short CParEfectivo = 7;
		static const short CPresionCamisa = 8;


		//Tipos de resultados instantáneos
		static const short CAngulo = 1;
		static const short CTemperaturaGas = 2;
		static const short CFraccionCalorLiberado = 3;
		static const short CDFraccionCalorLiberado = 4;
		static const short CCalorTransmitido = 5;
		static const short CVolumenCilindro = 6;
		static const short CExponentePolitropico = 7;
		static const short CGamma = 8;
		static const short CCalorTransmitidoRadiacion = 9;


		//Tipos de resultados instantaneos en arrastre
		static const short CLogaritmoPresion = 9;
		static const short CLogaritmoVolumen = 10;
		static const short CCalorExpol = 11;
		static const short CMasaCilArr = 12;

		//Tipos de resultados instantaneos en combustion (distintos de los de arrastre)
		static const short CMasaQuemada = 9;
		static const short CTemperaturadeQuemados = 0;
		static const short CTemperaturasinQuemar = 0;
		static const short CFQLNoAdim = 26;
		static const short CDFQLNoAdim = 27;
		static const short CMasaCilindro = 28;
		static const short CDensidad = 29;
		static const short CMasaBlowBy = 30;
		static const short CRAire = 31;
		static const short CRCombustible = 32;
		static const short CRQuemados = 33;
		static const short CCVAire = 34;
		static const short CCVCombustible = 35;
		static const short CCVQuemados = 36;
		static const short CFraccionMasicaAire = 37;
		static const short CFraccionMasicaCombustible = 38;
		static const short CFraccionMasicaQuemados = 39;
		static const short CMasaInyectada = 40;
		static const short CMasaEvaporada = 41;
		static const short CFlujoCalorTransPiston = 42;
		static const short CFlujoCalorTransCulata = 43;
		static const short CFlujoCalorTransCilindro = 44;
		static const short CRelacionCompresionDinamica = 45;
		static const short CCoefPelicula = 46;
		static const short CQRAD = 47;
		static const short CMasaAireInyectada = 48;
		static const short CMasaFuelDual = 49;


		//resultados cuasiestacionario
		static const short CAnguloCuasi = 47;
		static const short CGastoAdmision = 48;
		static const short CGastoEscape = 49;
		static const short CGastoCortocircuito = 50;
		static const short CCalorCulataCuasi = 51;
		static const short CCalorCilindroCuasi = 52;
		static const short CCalorPistonCuasi = 53;
		static const short CCoefPeliculaCuasi = 54;
		static const short CPresionCilindroCuasi = 55;
		static const short CTempGasCuasi = 56;
		static const short CMasaCilindroCuasi = 57;

		//resultados siciclo
		static const short CAnguloSiciclo = 1;
		static const short CPresionSiciclo = 2;


		//Resultados ruido
		static const short CRuidoCombustion = 1;
		static const short CAtenuacionBloque = 2;
		static const short CRuidoTotal = 3;
		static const short CFrecuenciaNormalizada = 4;

		//Tipos de gráficas
		static const short CGrafPresionCilindroAngulo = 1;
		static const short CGrafTemperaturaAngulo = 2;
		static const short CGrafPresionVolumen = 3;
		static const short CGrafCalorTransmitido = 4;
		static const short CGrafFqlAngulo = 5;
		static const short CGrafLogPLogV = 6;
		static const short CGrafExponentePolitropico = 7;
		static const short CGrafFMQAngulo = 8;
		static const short CGrafRuido = 9;
		static const short CGrafPresionAdmisionAngulo = 10;
		static const short CGrafPresionEscapeAngulo = 11;
		static const short CGrafPresionLineaAngulo = 12;
		static const short CGrafPresionCRAngulo = 13;
		static const short CGrafLevantamientoAngulo = 14;
		static const short CGrafParEfectivoAngulo = 15;
		static const short CGrafGammaAngulo = 16;
		static const short CGrafLevantamientoAnguloModal = 17;
		static const short CGrafSimuladoAnguloModal = 18;
		static const short CGrafMaximoMedioMinimoPresion = 19;
		static const short CGrafMaximoMedioMinimoLevantamiento = 20;
		static const short CGrafTodosLosCiclos = 21;
		static const short CGrafTodosLosCiclosFiltrados = 22;
		static const short CGrafTodosLosCiclosInyeccion = 23;
		static const short CGrafTQ = 25;

		static const short CGrafFQLNoAdimensional = 26;
		static const short CGrafDFQLNoAdimensional = 27;
		static const short CGrafMasaCilindro = 28;
		static const short CGrafDensidad = 29;
		static const short CGrafMasaBlowBy = 30;
		static const short CGrafRAire = 31;
		static const short CGrafRCombustible = 32;
		static const short CGrafRQuemados = 33;
		static const short CGrafCVAire = 34;
		static const short CGrafCVCombustible = 35;
		static const short CGrafCVQuemados = 36;
		static const short CGrafFMAire = 37;
		static const short CGrafFMCombustible = 38;
		static const short CGrafFMQuemados = 39;
		static const short CGrafMasaInyectada = 40;
		static const short CGrafMasaEvaporada = 41;
		static const short CGrafFlujoCalorPiston = 42;
		static const short CGrafFlujoCalorCulata = 43;
		static const short CGrafFlujoCalorCilindro = 44;
		static const short CGrafRCDinamica = 45;
		static const short CGrafCoefPelicula = 46;

		static const short CGrafGastoAdmision = 47;
		static const short CGrafGastoEscape = 48;
		static const short CGrafGastoCortocircuito = 49;

		static const short CGrafGastos = 50;

		//static const CGrafCoefPeliculaCuasi = 54
		//
		//static const CGrafPresionCilAnguloCuasi = 65
		//static const CGrafTemperaturaGasAnguloCuasi = 66
		//static const CGrafMasaCilAnguloCuasi = 67


		static const short CGrafCAnguloCuasi = 70;
		static const short CGrafCGastoAdmision = 71;
		static const short CGrafCGastoEscape = 72;
		static const short CGrafCGastoCortocircuito = 73;
		static const short CGrafCCalorCulataCuasi = 74;
		static const short CGrafCCalorCilindroCuasi = 75;
		static const short CGrafCCalorPistonCuasi = 76;
		static const short CGrafCCoefPeliculaCuasi = 77;
		static const short CGrafCPresionCilindroCuasi = 78;
		static const short CGrafCTempGasCuasi = 79;
		static const short CGrafCMasaCilindroCuasi = 80;


		//static const CGrafFraccionMasicaOxigeno = 55



		//tipos de senyales
		static const short CSPresionCilindroAngulo = 1;
		static const short CSTemperaturaAngulo = 2;
		static const short CSPresionVolumen = 3;
		static const short CSCalorTransmitido = 4;
		static const short CSFqlAngulo = 5;
		static const short CSdFqlAngulo = 6;
		static const short CSLogPLogV = 7;
		static const short CSExponentePolitropico = 8;
		static const short CSFMQAngulo = 9;
		static const short CSRuidoCombustion = 10;
		static const short CSRuidoTotal = 11;
		static const short CSAtenuacion = 12;
		static const short CSPresionAdmisionAngulo = 13;
		static const short CSPresionEscapeAngulo = 14;
		static const short CSPresionLineaAngulo = 15;
		static const short CSPresionCRAngulo = 16;
		static const short CSLevantamientoAngulo = 17;
		static const short CSParEfectivoAngulo = 18;
		static const short CSGammaAngulo = 19;
		static const short CSLevantamientoAnguloFQL = 20;
		static const short CSTemperaturaQuemados = 21;
		static const short CSTemperaturaSinQuemar = 22;
		//presion camisa
		static const short CSPresionCamisaAngulo = 23;

		static const short CSFQLNoAdimensional = 26;
		static const short CSDFQLNoAdimensional = 27;
		static const short CSMasaCilindro = 28;
		static const short CSDensidad = 29;
		static const short CSMasaBlowBy = 30;
		static const short CSRAire = 31;
		static const short CSRCombustible = 32;
		static const short CSRQuemados = 33;
		static const short CSCVAire = 34;
		static const short CSCVCombustible = 35;
		static const short CSCVQuemados = 36;
		static const short CSFMAire = 62;
		static const short CSFMCombustible = 38;
		static const short CSFMQuemados = 39;
		static const short CSMasaInyectada = 40;
		static const short CSMasaEvaporada = 41;
		static const short CSFlujoCalorPiston = 42;
		static const short CSFlujoCalorCulata = 43;
		static const short CSFlujoCalorCilindro = 44;
		static const short CSRCDinamica = 45;

		static const short CSGastoAdmision = 47;
		static const short CSGastoEscape = 48;
		static const short CSGastoCortocircuito = 49;


		static const short CSDerivadaPresionCilindroAngulo = 55;

		static const short CSCalorTransmitidoRadiacion = 56;

		static const short CSFQL_ACT = 57;
		static const short CSFQL_PMX_ACT = 58;

		//static constante para la grafica de Senyal de Comando
		static const short CSSenyaldeComando = 59;

		//static constantes para gráficas de motores MEP
		static const short CSTemperaturaAnguloSQ = 60;
		static const short CSGammaAnguloSQ = 61;

		static const short CSFQLACT_NORM = 63;

		static const short CSCalorExponentePolitropico = 64;





		//Columnas para el fichero PUMA, el valor que tengan será la columna en
		//el fichero puma que contendrá el valor de la variable
		//la static constante debe valer cero si el dato no está disponible en el fichero
		//último cambio sugerido por Santi a fecha 07-04-99 para incluir contaminantes
		static const short ColRegimen = 5;
		static const short ColParEfectivo = 6;
		static const short ColPresionAdmision = 7;
		static const short ColPresionEscape = 9;
		static const short ColTemperaturaAdmision = 8;
		static const short ColTemperaturaEscape = 10;
		static const short ColTemperaturaAguaSalidaMotor = 14;
		static const short ColTemperaturaAceite = 16;
		static const short ColPresionAmbiente = 2;
		static const short ColTemperaturaAmbiente = 3;
		static const short ColHumedadRelativa = 4;
		static const short ColGastoAire = 11;
		static const short ColGastoCombustible = 12;
		static const short ColGastoBlowBy = 18;
		static const short ColGastoEGR = 19;
		//static const ColRegimenTurbo = 0
		//static const ColGastoAgua = 0
		//static const ColPresionEntradaCompresor = 0
		//static const ColPresionSalidaCompresor = 0
		//static const ColPresionEntradaIntercooler = 0
		//static const ColPresionSalidaIntercooler = 0
		//static const ColPresionEntradaTurbina = 0
		//static const ColPresionSalidaTurbina = 0
		static const short ColPresionAceite = 17;
		//static const ColPresionCombustibleSalidaBomba = 0
		static const short ColPresionRail = 28;
		//static const ColTemperaturaEntradaCompresor = 0
		//static const ColTemperaturaSalidaCompresor = 0
		//static const ColTemperaturaEntradaIntercooler = 0
		//static const ColTemperaturaSalidaIntercooler = 0
		//static const ColTemperaturaEntradaTurbina = 0
		//static const ColTemperaturaSalidaTurbina = 0
		static const short ColTemperaturaAguaEntradaMotor = 13;
		static const short ColTemperaturaEntradaBomba = 15;
		static const short ColNOx = 21;
		//static const ColNO = 0
		static const short ColCO = 22;
		static const short ColO2 = 24;
		static const short ColHC = 23;
		static const short ColCO2Admision = 26;
		static const short ColCO2Escape = 25;
		static const short ColHumos = 20;
		//static const ColParticulas = 0
		//static const ColDosado = 0
		static const short ColEGRHoriba = 27;
		static const short ColNumInyecciones = 29;
		static const short ColIniInyeccion1 = 30;
		static const short ColMasInyeccion1 = 31;
		static const short ColIniInyeccion2 = 32;
		static const short ColMasInyeccion2 = 33;
		static const short ColIniInyeccion3 = 34;
		static const short ColMasInyeccion3 = 35;
		static const short ColIniInyeccion4 = 36;
		static const short ColMasInyeccion4 = 37;
		static const short ColIniInyeccion5 = 38;
		static const short ColMasInyeccion5 = 39;
		//Distribucion variable
		static const short ColAAA = 40;
		static const short ColRCA = 41;
		static const short ColLEMA = 42;
		static const short ColAAE = 43;
		static const short ColRCE = 44;
		static const short ColLEME = 45;

		static const short ColGAINY = 46;
		static const short ColPAINY = 47;
		static const short ColGFDUAL = 48;

		static const short ColGOIL = 49;
		static const short ColGastoAgua = 50;
		static const short ColPCOOL = 51;
		static const short ColPresionCombustibleSalidaBomba = 52;
		static const short ColPresionEntradaCompresor = 53;
		static const short ColPresionSalidaCompresor = 54;
		static const short ColPresionEntradaIntercooler = 55;
		static const short ColPresionSalidaIntercooler = 56;
		static const short ColPresionEntradaTurbina = 57;
		static const short ColPresionSalidaTurbina = 58;
		static const short ColRegimenTurbo = 59;
		static const short ColTemperaturaEntradaCompresor = 60;
		static const short ColTemperaturaSalidaCompresor = 61;
		static const short ColTemperaturaEntradaIntercooler = 62;
		static const short ColTemperaturaSalidaIntercooler = 63;
		static const short ColTemperaturaEntradaTurbina = 64;
		static const short ColTemperaturaSalidaTurbina = 65;
		static const short ColDosado = 66;
		static const short ColNO = 67;
		static const short ColParticulas = 68;
		static const short ColOpacidad = 69;



		//Si CGuardaDatos es true, se grabarán ficheros de depuracion
		bool CGuardaDatos;


		//Si CAjusteManual es true, se puede definir que parametros se ajustan en el ajuste en arrastre y combustion
		bool CAjusteManual;
		//Si CLeerArchivoAjuste es true, los parámetros se obtendrán de un fichero .dat
		bool CLeerArchivoAjuste;
		//Si CLeerArchivoTemperaturas es true, las temperaturas de las paredes se obtendrán de un fichero .dat
		bool CLeerArchivoTemperaturas;


		//Si CMenuCMT es true, se activan menus como el de gestion global de ensayos...
		bool CMenuCMT;

		//Si COldCarac es true, el ensayo de arrastre se realiza como antes de esta version, en
		//el que la cw1 y la kdef se ajustaban. En esta version, por defecto, el cw1 y la kdef
		//se dejan fijas y cogen el valor que se indica en la base
		bool COldCarac;

		//Si CPermitirFicherosPUMA es true, se permite la entrada de datos mediante ficheros PUMA
		bool CPermitirFicherosPUMA;

		//CFicherosCSV contiene el directorio donde se han generado los ficheros dat con el e-filemaker
		//para luego eliminar los dat y el proceso de csv a dat sea invisible al usuario
		std::string CFicherosCSV;

		//Si CImportarOsiris es true se importan datos de OSIRIS via OSICOM
		bool CImportarOsiris;
		bool CCalmecFTAuto;

		//Si CMostrarInterfaz es false el programa evita la interaccion con el usuario
		bool CMostrarInterfaz;

		//Si CCalcularACT es false el programa no ejcuta ACT (eso pasa solamente cuando CCalcularPMX y CCalcularRAD son false)
		bool CCalcularACT;

		//Si CCalcularPMX es false el programa no guarda en la base de datos la masa quemada en premezcla
		bool CCalcularPMX;

		//Si CCalcularRAD es false el programa no guarda en la base de datos la radiación
		bool CCalcularRAD;

		//Si CCalcularTQ es false el programa no calcula la temperatura de quemados
		bool CCalcularTQ;

		//Si CGuardarCalorParedes es true se guarda en un fichero dat los calores transmitidos a las paredes
		bool CGuardarCalorParedes;

		//Si CTemperaturaWallManual es true se muestra un form para definir manualmente las temperaturas de pared
		bool CTemperaturaWallManual;

		//Si es True permite mostrar la gráfica de FQL y FQL_PMX  de ACT junto con la FQL no ADIM de CALMEC
		//UPGRADE_NOTE: CACT se actualizó a CACT_Renamed. Haga clic aquí para obtener más información: 'ms-help://MS.VSCC.v80/dv_commoner/local/redirect.htm?keyword="A9E4979A-37FA-4718-9994-97DD76ED70A7"'
		bool CACT_Renamed;

		//Si es True se calculará calmec en todo el ciclo completo
		bool CCicloCompleto;

		//Si es True se buscará el fichero RefPre.txt para referenciar los ficheros de presión de adm, esc y camisa.
		bool CRefPresiones;
		//Si CLeerArchivoRefPresiones es true, ya se ha leido el fichero .dat con las referencias para las presiones
		bool CLeerArchivoRefPresiones;

		bool CPermitirOSIRIS;
		//para que CALMEC trabaje con Osiris solo hace falta cambiar esta variable a true y poner la variable
		//ConOsiris (variable de compilación condicional) a true

		//25/10/2007 variable que indica la selección que ha realizado el usuario en el Formulario FAvisoPMI cuando se muestra interfaz o no.
		//Si es los valores son: A= Opción 1 Continuar /  $valor= Opción 2 Sumar ese valor a PMI+PA+$valor  / C= Opción 3 Eliminar Ensayo
		std::string COpcionPMI;

		//CConfiguracionSiciclo CAnalisisCicloRealIdeal = new CConfiguracionSiciclo();


		//static constantes Varias
		//static constantes para la curva utilizada para el término de velocidad
		static const short Ck1 = 100;
		static const short Ck2 = 2;

		//estas tienen que ver con la localización de inicios y finales de inyección

		static const double CTiempoMinimoInyeccionValida;
		static const short CPorcentajeBajoInyeccionValida = 5;
		static const short CPorcentajeAltoInyeccionValida = 15;

		static const short CMaximoPresion_VentanaAngular = 10;

		static const short CNumeroMaximoGraficasAbiertas = 12;

		//Las siguientes static constantes se utilizan en los procedimientos iterativos
		//que se utilizan para ajustar parámetros
		static const short CAjustaPresionPMICombustion_IncrementoInicial = 10000;
		static const double CAjustaMRCACombustion_IncrementoInicial;
		static const short CAjustaPPMI_MRCA_AnguloInicioInyeccion = -60;
		static const short CAjustaPPMI_MRCA_MaxNumIter = 100;
		static const short CAjustaCbb_MaxNumIter = 100;

		static const short CAjustaPresionPMIArrastre_IncrementoInicial = 5000;
		static const double CAjustaMRCAArrastre_IncrementoInicial;
		static const double CAjustaPresionPMI_MRCAArrastre_Delta;

		static const short CAjustaPresionPMI_MRCAArrastre_MaxNumIter = 100;
		static const short CAjustaDifPresionPMIArrastre_MaxNumIter = 15;



		static const double CTasaDesfasada_Delta;
		static const short CTasaDesfasada_MaxNumIter = 100;



		static const double CAjustaRelacionCompresion_IncrementoInicial;
		static const double CAjustaRelacionCompresion_Delta;
		static const short CAjustaRelacionCompresion_MaxNumIter = 25;


		static const double CAjustaCalorTransmitidoMaximo_IncrementoInicial;
		static const double CAjustaCalorTransmitidoMaximo_Delta;
		static const short CAjustaCalorTransmitidoMaximo_MaxNumIter = 25;




		static const double CPuntosG_ValorInicial;
		static const double CPuntosG_Delta;


		//static constantes que tiene que ver con el dibujado de graficas

		static const short CAnchoCirculoInyeccion = 75;
		static const int CColorRejilla = 0xc0c0c0;
		static const short CModoAngulo = 1;
		static const short CModoTiempo = 2;
		static const short CModoFQL = 3;

		static const short CBordeEjeX = 5;
		static const short CBordeEjeY = 3;

		static const std::string CFuenteTexto;
		static const std::string CFuenteTextoArial;
		static const std::string CFuenteCifra;
		static const short CTamañoTexto = 12;
		static const short CTamañoCifra = 10;
		static const short CTamañoCifraPeque = 8;
		static const short CColorCifra = 0x0;
		static const short CColorTexto = 0x0;
		static const int CColorFondo = 0xffffff;

		//static constantes por defecto si no están en calmec.ini

		//static const CPathBasedeDatos = "\calmecdb\calmet.mdb"
		//static const CPathAyuda = "\aide\AIDE_CALMEC.HLP"


		//static constantes que identifican a las distintas cadenas de medidas de tipo CMT
		static const short CYokoAR1100 = 1;
		static const short CYokoOR1400 = 2;
		static const short CYokoORM1200 = 3;
		static const short CYokoDL708 = 4;
		static const short CYkDL850V = 5;

		//static constantes añadidas para compatibilizar CALMEC con OSICOM, estos son los
		//indices que se utilizarán para hacer referencia a una matriz donde se
		//almacenarán los nombres de las variables OSICOM y los factores de conversión


		static const short CMAX_OSIRIS_VALUES = 231;
		static const short CPorcentajeCalmecNoValido = 10;
		static const double CColaFQlCalmecNoValido;
		static const short CPorcentajePmaxCicloNoValido = 5;
		static const short CPorcentajeCiclosMalosEnsayoNoValido = 20;
		//arrastre
		static const short C_OSI_ARR_N = 1;
		static const short C_OSI_ARR_ME = 2;
		static const short C_OSI_ARR_GA = 3;
		static const short C_OSI_ARR_GF = 4;
		static const short C_OSI_ARR_GBLBY = 5;
		static const short C_OSI_ARR_GEGR = 6;
		static const short C_OSI_ARR_PA = 7;
		static const short C_OSI_ARR_PE = 8;
		static const short C_OSI_ARR_PAMB = 9;
		static const short C_OSI_ARR_TA = 10;
		static const short C_OSI_ARR_TE = 11;
		static const short C_OSI_ARR_TAC = 12;
		static const short C_OSI_ARR_TAMB = 13;
		static const short C_OSI_ARR_HR = 14;
		static const short C_OSI_ARR_TRS = 15;

		//combustion
		static const short C_OSI_COMB_N = 30;
		static const short C_OSI_COMB_NT = 31;
		static const short C_OSI_COMB_MEN = 32;
		static const short C_OSI_COMB_GA = 33;
		static const short C_OSI_COMB_GF = 34;
		static const short C_OSI_COMB_GBLBY = 35;
		static const short C_OSI_COMB_GEGR = 36;
		static const short C_OSI_COMB_GAG = 37;
		static const short C_OSI_COMB_PA = 38;
		static const short C_OSI_COMB_PAEC = 39;
		static const short C_OSI_COMB_PASC = 40;
		static const short C_OSI_COMB_PAEI = 41;
		static const short C_OSI_COMB_PASI = 42;
		static const short C_OSI_COMB_PE = 43;
		static const short C_OSI_COMB_PEET = 44;
		static const short C_OSI_COMB_PEST = 45;
		static const short C_OSI_COMB_PAC = 46;
		static const short C_OSI_COMB_PF = 47;
		static const short C_OSI_COMB_PCR = 48;
		static const short C_OSI_COMB_PAMB = 49;
		static const short C_OSI_COMB_TAMB = 50;
		static const short C_OSI_COMB_TA = 51;
		static const short C_OSI_COMB_TAEC = 52;
		static const short C_OSI_COMB_TASC = 53;
		static const short C_OSI_COMB_TAEI = 54;
		static const short C_OSI_COMB_TASI = 55;
		static const short C_OSI_COMB_TE = 56;
		static const short C_OSI_COMB_TEET = 57;
		static const short C_OSI_COMB_TEST = 58;
		static const short C_OSI_COMB_TF = 59;
		static const short C_OSI_COMB_TRE = 60;
		static const short C_OSI_COMB_TRS = 61;
		static const short C_OSI_COMB_TAC = 62;
		static const short C_OSI_COMB_HR = 63;
		static const short C_OSI_COMB_NOX = 64;
		static const short C_OSI_COMB_NO = 65;
		static const short C_OSI_COMB_CO = 66;
		static const short C_OSI_COMB_CO2A = 67;
		static const short C_OSI_COMB_CO2E = 68;
		static const short C_OSI_COMB_O2 = 69;
		static const short C_OSI_COMB_HC = 70;
		static const short C_OSI_COMB_UB = 71;
		static const short C_OSI_COMB_PART = 72;
		static const short C_OSI_COMB_FH = 73;
		static const short C_OSI_COMB_GEGRH = 74;
		static const short C_OSI_COMB_NINY = 75;
		static const short C_OSI_COMB_ECO = 76;
		static const short C_OSI_COMB_EHC = 77;
		static const short C_OSI_COMB_ENOX = 78;

		//resultados
		static const short C_OSI_RES_PMAX = 200;
		static const short C_OSI_RES_APMAX = 201;
		static const short C_OSI_RES_GRAP = 202;
		static const short C_OSI_RES_AGRAP = 203;
		static const short C_OSI_RES_PMI = 204;
		static const short C_OSI_RES_PMB = 205;
		static const short C_OSI_RES_RTOT = 206;
		static const short C_OSI_RES_TMAX = 207;
		static const short C_OSI_RES_DUCOM = 208;
		static const short C_OSI_RES_AICOM = 209;
		static const short C_OSI_RES_DFQLMAX = 210;
		static const short C_OSI_RES_QT = 211;
		static const short C_OSI_RES_NQPIS = 212;
		static const short C_OSI_RES_NQCUL = 213;
		static const short C_OSI_RES_NQCIL = 214;
		static const short C_OSI_RES_NQPISA = 215;
		static const short C_OSI_RES_NQCULA = 216;
		static const short C_OSI_RES_NQCILA = 217;
		static const short C_OSI_RES_QWPISCF = 218;
		static const short C_OSI_RES_QWCULCF = 219;
		static const short C_OSI_RES_QWCILCF = 220;
		static const short C_OSI_RES_QWPISCO = 221;
		static const short C_OSI_RES_QWCULCO = 222;
		static const short C_OSI_RES_QWCILCO = 223;
		static const short C_OSI_RES_AQ25PMS = 224;
		static const short C_OSI_RES_AQ50PMS = 225;
		static const short C_OSI_RES_AQ75PMS = 226;
		static const short C_OSI_RES_AQ90PMS = 227;
		static const short C_OSI_RES_CSI = 228;
		static const short C_OSI_RES_CSINET = 229;
		static const short C_OSI_RES_RCOM = 230;
		static const short C_OSI_RES_RMECNET = 231;

		//tipo de graficas para señales de creacion de ensayos en transitorio
		static const short CGrafTransGastoAire = 1;
		static const short CGrafTransGastoFuel = 2;
		static const short CGrafTransPresionCilindro = 3;
		static const short CGrafTransRegimen = 4;
		static const short CGrafTransPar = 5;
		static const short CGrafTransPresionAdmision = 6;
		static const short CGrafTransPresionEscape = 7;
		static const short CGrafTransAlfa0 = 8;

		//24/07/07 variable para Corrección de v.C_W2 al llamar a Calcula_Calor_Woschni cuando hay radiación
		//static const C2_ConRadiacion = 0.00054 'primera version SOOT
		//static const C2_ConRadiacion = 0.0005 'segunda version SOOT
		//static const C2_ConRadiacion = 0 'segunda version SOOT Nicolas

		//static constante indicando que es un valor nulo en la base de datos
		//debido a que los double nulos los considera como un 0
		static const short ValorNulo = -9999;


		static const int Azul = 16777088;
		static const int Gris = 14737632;
		static const int Naranja = 0x80c0ff;
		static const int Blanco = 0xffffff;
		static const short Negro = 0x0;
		static const short Rojo = 0xff;
		static const int Verde = 0x80ff80;
		static const int Marron = 0x4080;
		static const int GrisDesactivado = 0xe0e0e0;
		static const int Amarillo = 0x80ffff;


		//Tipos de resultados instantaneos en combustion MEP
		static const short CAngulo_MEP = 201;
		static const short CCVAireQ_MEP = 202;
		static const short CCVAireSQ_MEP = 203;
		static const short CCVCombustibleQ_MEP = 204;
		static const short CCVCombustibleSQ_MEP = 205;
		static const short CCVQuemadosQ_MEP = 206;
		static const short CCVQuemadosSQ_MEP = 207;
		static const short CDensidadQ_MEP = 208;
		static const short CDensidadSQ_MEP = 209;
		static const short CDFraccionCalorLiberado_MEP = 210;
		static const short CDFQLNoAdim_MEP = 211;
		static const short CIncMasaBBQ_MEP = 212;
		static const short CIncMasaBBSQ_MEP = 213;
		static const short CIncMasaLlama_MEP = 214;
		static const short CExponentePolitropico_MEP = 215;
		static const short CFraccionCalorLiberado_MEP = 216;
		static const short CFQLNoAdim_MEP = 217;
		static const short CGammaQ_MEP = 218;
		static const short CGammaSQ_MEP = 219;
		static const short CCoefPeliculaQ_MEP = 220;
		static const short CCoefPeliculaSQ_MEP = 221;
		static const short CMasaQ_MEP = 222;
		static const short CMasaSQ_MEP = 223;
		static const short CCalorTransmitido_MEP = 224;
		static const short CCalorTransmitidoCilQ_MEP = 225;
		static const short CCalorTransmitidoCilSQ_MEP = 226;
		static const short CCalorTransmitidoCulQ_MEP = 227;
		static const short CCalorTransmitidoCulSQ_MEP = 228;
		static const short CCalorTransmitidoPisQ_MEP = 229;
		static const short CCalorTransmitidoPisSQ_MEP = 230;
		static const short CMasaQuemada_MEP = 231;
		static const short CRelacionCompresionDinamica_MEP = 232;
		static const short CMasaTotal_MEP = 233;
		static const short CRadioFrenteLlama_MEP = 234;
		//static const CRQuemados_MEP = 235
		static const short CTemperaturaGasQ_MEP = 236;
		static const short CTemperaturaGasSQ_MEP = 237;
		static const short CVolumenCilindro_MEP = 238;
		static const short CVelocidadCombustion_MEP = 239;
		static const short CVelocidadFrenteLlama_MEP = 240;
		static const short CVolumenQ_MEP = 241;
		static const short CVolumenSQ_MEP = 242;
		static const short CLogaritmoPresion_MEP = 243;
		static const short CLogaritmoVolumen_MEP = 244;
		static const short CCalorTransmitidoQ_MEP = 245;
		static const short CCalorTransmitidoSQ_MEP = 246;


		//CARPETAS PARA CALMEC AUTOMATICO
		//Carpeta para guardar los ensayos que el usuario ha decidido guardar
		static const std::string CarpetaEnsayosGuardadosCA;
		//Carpeta para guardar los ficheros del ensayo que se está calculando actualmente y liberar la carpeta para que
		//el programa de Ali pueda poner nuevos ficheros
		static const std::string CarpetaEnsayoEnCalculoCA;
		static const std::string CarpetaEnsayoCalculadoCA;
		static const std::string CarpetaEnsayoEnPantallaCA;
	};
#endif
//}