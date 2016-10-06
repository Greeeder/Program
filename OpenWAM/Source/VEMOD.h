//---------------------------------------------------------------------------

//#ifndef VEMODH
//#define VEMODH
//---------------------------------------------------------------------------

#ifdef __GNUC__
	#define EXLIB_API
#else
	#define EXLIB_API __declspec(dllexport)
#endif

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

/*!
 * \brief	Loads the model. It must be call at the begining of the simulation.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	18/05/2016
 *
 * \param	icase			The index of the case.
 * \param	FileName		Input data file name
 */

EXTERNC EXLIB_API void LoadModel(const int icase, char *FileName);

/*!
 * \brief	Executes the model during a time step.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	18/05/2016
 *
 * \param	icase			The index of the case.
 * \param	dt			 	The current time [s].
 */

EXTERNC EXLIB_API void RunModel(const int icase, const double t);

/*!
 * \brief	Gets the value of a sensor.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	18/05/2016
 *
 * \param	icase			The index of the case.
 * \param	sensor_id	 	Identifier for the sensor.
 *
 * \return	The sensor value.
 */

EXTERNC EXLIB_API double getSensor(const int icase, const int sensor_id);

/*!
 * \brief	Sets an actuator.
 *
 * \author	F.J. Arnau (farnau@mot.upv.es)
 * \date	18/05/2016
 *
 * \param	icase	   	The index of the case.
 * \param	actuator_id	Identifier for the actuator.
 * \param	value	   	The value.
 */

EXTERNC EXLIB_API void setActuator(const int icase, const int actuator_id, const double value);








//extern "C" __declspec(dllexport)void __stdcall INITCMODELF(void);
//
//extern "C" __declspec(dllexport)void __stdcall LOADCMODELF(char*FileName);
//
//extern "C" __declspec(dllexport)void __stdcall ALLOCATEMODELF();
////
//// extern "C" __declspec(dllexport)void __stdcall UPDATEBOUNDARY(int i, double U0, double U1,
//// double T0, double T1, double P0, double P1, double t);
////
//// extern "C" __declspec(dllexport)void __stdcall INITIATEBOUNDARY(int i, double D0, double D1,
//// double dX);
//
//extern "C" __declspec(dllexport)void __stdcall UPDATEBOUNDARY(int i, double U0, double T0, double P0, double t);
//
//extern "C" __declspec(dllexport)void __stdcall UPDATECONTROLLER(int i, double Val, double t);
//
//extern "C" __declspec(dllexport)void __stdcall INITIATEBOUNDARY(int i, double A);
//
//extern "C" __declspec(dllexport)void __stdcall RUNSTEP(double t);
//
//extern "C" __declspec(dllexport)void __stdcall LOADNEWDATA(int i, double*p, double*T, double*u);
//
//extern "C" __declspec(dllexport)void __stdcall CLOSEMODEL(void);
//
//extern "C" __declspec(dllexport)void __stdcall GETPLOT(int ID, double*val);
//
//extern "C" __declspec(dllexport)void __stdcall UPDATEHTMDATA(int TC, int ID, double val);

//#endif
