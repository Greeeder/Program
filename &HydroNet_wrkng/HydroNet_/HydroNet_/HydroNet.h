#ifndef HydroNet_h
#define HydroNet_h

//#include "stdafx.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>
#include <Math_wam/Math_wam.h>
#include <Includes/Constantes.h>


class HydroNet_model {

public:

	std::string circuit_file;								//!< File contaning the circuit information
	double dt;												//!< Time step
	std::vector<double> flows;								//!< Each branche flow [L / s]
	std::string fluid_type;									//!< Type of fluid
	int max_divisions =20;									//!< Maximum branch division

	/*!
	*\brief Initialice values and execututes the model
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	HydroNet_model();
	~HydroNet_model();
	void HydroNet_Main();



private:

	// DATA STRUCTURES

	struct strObject {
		int ID;												//!< Object's idintificator
		std::string name;									//!< Object's name
		std::string Class;									//!< Object's clas
		double volume;										//!< Objects volume [m ^ 3]
		std::vector<int> adjacent;							//!< Adjacent objects
		int branch;											//!< Branch containing the object
	};
	
	struct strPipe
	{
		int ID;										//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Pipe's name
		int inlet_object;									//!< Pipe's inlet object identifier
		int outlet_object;									//!< Pipe's outlet object identifier
		double length;										//!< Pipe's length [m]
		double diameter;									//!< Pipe's diameter [m]
		double friction_coef;								//!< Pipe's friction coefficient

	};

	struct strValve_fix
	{
		int ID;									//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Valve's name	
		int inlet_object;									//!< Valve's inlet object identifier
		int outlet_object;									//!< Valve's outlet object identifier
		double head_loss;									//!< Head loss [m.c.f. / m ^ 3 / s]
	};

	struct strValve_var
	{
		int ID;									//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Valve's name	
		int inlet_object;									//!< Valve's inlet object identifier
		int outlet_object;									//!< Valve's outlet object identifier
		double opening;										//!< Opening [degree]
		double coef0;										//!< Coefficent order 0 [m.c.f.]
		double coef1;										//!< Coefficent order 1 [m.c.f. / degree]
		double coef2;										//!< Coefficent order 2 [m.c.f. * degree ^ (-2)]
		double coef3;										//!< Coefficent order 3 [m.c.f. * degree ^ (-3)]
	};

	struct strThermostat
	{
		int ID;									//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Thermostat's name	
		int inlet_object;									//!< Thermostat's inlet object identifier
		int outlet_object;									//!< Thermostat's outlet object identifier
		double T_closed;										//!< Start temperaturure for the thermostat change [K]
		double T_open;										//!< Finish temperaturure for the thermostat change [K]
		double Shape_factor;								//!< Slope shape factor
		bool to_open;										//!< 0 to closed; 1 to open'
		double coef_h_0;									//!< Coefficent order 0 [m.c.f.]
		double coef_h_1;									//!< Coefficent order 1 [m.c.f. / m ^ 3 / s]
		double coef_h_2;									//!< Coefficent order 2 [m.c.f. / (m ^ 3 / s) ^ 2]
		double coef_h_3;									//!< Coefficent order 3 [m.c.f. / (m ^ 3 / s) ^ 3]
		double opening;										//!< Opening [degree]
		double Ref_temp;									//!< Temperature of reference to control the thermostat
		int Sensor;											//!< Object in wich the temperature is taken
	};

	struct strPump_turbo
	{
		int ID;										//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Pump's name	
		int inlet_object;									//!< Pups's inlet object identifier
		int outlet_object;									//!< Pump's outlet object identifier
		double coef_N0;										//!< Coefficent order 0 [m.c.f]
		double coef_N1;										//!< Coefficent order 0 N [m.c.f / (rad / s)]
		double coef_N2;										//!< Coefficent order 0 N2 [m.c.f / (rad / s) ^ 2]
		double Q_coef_N0;									//!< Coefficent order 1 [m.c.f / m ^ 3 * s]
		double Q_coef_N1;									//!< Coefficent order 1 N [m.c.f / m ^ 3 * s / (rad/s)]
		double Q_coef_N2;									//!< Coefficent order 1 N2 [m.c.f / m ^ 3 * s / (rad/s)^2]
		double Q2_coef_N0;									//!< Coefficent order 2 [m.c.f /  m ^ 6 * s ^ 2]
		double Q2_coef_N1;									//!< Coefficent order 2 N [m.c.f /  m ^ 6 * s ^ 2 / (rad/s)]
		double Q2_coef_N2;									//!< Coefficent order 2 N2 [m.c.f. / m ^ 6 * s ^ 2 / (rad/s)^2]
		double pump_speed;									//!< Pump speed [rad/s]
		double head_max;									//!< Max head preassure
		double heat;										//!< Heat Energy Exchanged [J]
		double pump_head;									//!< Head preassure [m.c.f.]
		double pump_power;									//!< Pumps power	[W]
	};

	struct strPump_volum
	{
		int ID;										//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Pump's name	
		int inlet_object;									//!< Pups's inlet object identifier
		int outlet_object;									//!< Pump's outlet object identifier
		double coef0;										//!< Coefficent order 0 [m ^ 3 / s]
		double coef1;										//!< Coefficent order 1 [m^3/s / (rad / s)]
		double coef2;										//!< Coefficent order 2 [ m^3/s / (rad / s) ^ 2]
		double coef3;										//!< Coefficent order 3 [m^3/s / (rad / s) ^ 3]
		double pump_speed;									//!< Pump speed [rad/s]
		double head_max;									//!< Max head preassure
		double heat;										//!< Heat Energy Exchanged [J]
		double pump_head;									//!< Head preassure [m.c.f.]
		double pump_power;									//!< Pumps power [W]
	};

	struct strHeat_exch
	{
		int ID;									//!< Identifier
		int obj_index;										//!< Index
		std::string name;									//!< Heat Exchanger name
		int inlet_object;									//!< Heat Exchanger's inlet object identifier
		int outlet_object;									//!< Heat Exchanger's outlet object identifier
		std::string type;									//!< Heat Echanger's type
		double heat;										//!< Heat Energy Exchanged [J]
		double T_out;										//!< Outlet temperaure [K]
		double T_in;										//!< Intlet temperaure [K]
		double hydr_resist;									//!< Hydraulic resistance = head loss / flow^2 [m.c.f. / m ^ 6 * s ^ 2]
	};

	//struct strHeat_exch_Tout
	//{
	//	int ID;									//!< Identifier
	//  int obj_index;										//!< Index
	//	std::string name;									//!< Heat Exchanger name
	//	int inlet_object;									//!< Heat Exchanger's inlet object identifier
	//	int outlet_object;									//!< Heat Exchanger's outlet object identifier
	//	double T_out;										//!< Outlet temperaure [K]
	//	double T_in;										//!< Intlet temperaure [K]
	//	double hydr_resist;									//!< Hydraulic resistance = head loss / flow^2 [m.c.f. / m ^ 6 * s ^ 2]
	//};

	struct strTank
	{
		int ID;												//!< Identifier
		int obj_index;										//!< Index
		std::string name;										//!< Tank's name
		int inlet_object;										//!< Tank's inlet object identifier
		int outlet_object;										//!< Tank's outlet object identifier
		double volume;											//!< Tank's volume [m ^ 3]
		double head_loss;										//!< Head loss [m.c.f. / m3/s]
	};	

	struct strObjects
	{
		int ID;												//!< Object's idintificator
		std::string name;									//!< Object's name
		std::string Class;									//!< Object's clas
		double volume;										//!< Objects volume [m ^ 3]
		int inlet_object;									//!< Tank's inlet object identifier
		int outlet_object;									//!< Tank's outlet object identifier
		std::vector<int> adjacent;							//!< Adjacent objects
	};	

	struct strBranches
	{
		int branch_ind;										//!< branch idintificator
		int branches_cycle;
		//int n_branch;										//!< Number of branches
		std::vector<strObjects> objects;					//!< Objects in Brnach
	};

	struct strSize
	{
		int temperature;									//!< Temperature vector size
		int position;										//!< Position vector size
		std::vector<int> new_temp;							//!< New_temp vector size
		std::vector<int> new_pos;							//!< New_pos vector size
		int overflow_temperature_aux;						//!< Overflow_temperature_aux vector size
		int overflow_temp_volpos_aux;						//!< Overflow_temp_volpos_aux vector size
		std::vector<int> node_inlet_branches;				//!< Node_inlet_branches vector size
		std::vector<int> node_outlet_branches;				//!< Node_outlet_branches vector size
		int branch_ind_order_aux;							//!< Branch_ind_order_aux vector size
		int aux_vec;										//!< Aux_vec vector size
		int to_remove;										//!< To_remove vector size
		int unk_ind_aux;									//!< Unk_ind_aux vector size
		int remain_branch;									//!< Remain_branch vector size
		
	};

	struct strInsertVolum {
		std::vector<double> temperature;					//!< New postion vector size [K]
		std::vector<double> position;						//!< New postion vector size
		int temp_size;										//!< New temperature vector size
		int pos_size;										//!< New postion vector size
	};


	std::vector <strObject> objects;						//!< Vectror with all objects present in the net and their interaction
	std::vector<strBranches> branches;						//!< Vector with branches data
	std::vector<strPump_volum> pump_volum;					//!< Vector with volumetric pumps data
	std::vector<strPump_turbo> pump_turbo;					//!< Vector with turbo pumps data
	std::vector<strValve_var> valve_var;					//!< Vector with variable valves data
	std::vector<strValve_fix>valve_fix;						//!< Vector with fix valves data
	std::vector<strThermostat> thermostat;					//!< Vector with thermostat data
	std::vector<strPipe>pipe;								//!< Vector with pipes data
	//std::vector<strHeat_exch_Tout> heat_exch_Tout;			//!< Vector with heat echangers type temperature oulet data
	std::vector<strHeat_exch> heat_exch;					//!< Vector with heat echangers data
	std::vector<strTank> tank;								//!< Vector with tanks data

	std::vector<std::vector<double>> branch_temp_pos;		//!< Position of the temperaatures in the branches
	std::vector<std::vector<double>> branch_temperature;	//!< Temperaatures in the branches [K]
	std::vector<std::vector<double>> obj_inlet_pos;			//!< Vector of objects inlets
	std::vector<std::vector<double>> obj_outlet_pos;		//!< Vector of objects outlets
	std::vector<std::vector<int>> mesh_branches;			//!< Vector with the maesh branches index ordered for each mesh
	std::vector<std::vector<int>> node_branches;			//!< Vector with the branches conected to each node
	std::vector<std::vector<int>>branches_id;				//!< Vector with the branches objcts identificator ordered for each branch
	std::vector<std::vector<int>>branches_ind;				//!< Vector with the branches objcts index ordered for each branch
//	std::vector<std::vector<int>>branch_htx_Tout;			//!< Vector with the heat exchangers type temperature outlet index in each branch
	std::vector<std::vector<int>>branch_htx;			//!< Vector with the heat exchangers type fix index in each branch
	std::vector<std::vector<int>>branch_pump_volum;			//!< Vector with the volumetric pumps index in each branch
	std::vector<std::vector<int>>branch_pump_turbo;			//!< Vector with the turbo pumps index in each branch
	std::vector<std::vector<int>>branch_valve_var;			//!< Vector with the valve_var index in each branch
	std::vector<std::vector<int>>branch_valve_fix;			//!< Vector with the valve_fix index in each branch
	std::vector<std::vector<int>>branch_pipe;			//!< Vector with the pipe index in each branch
	std::vector<double> branch_volume;						//!< Volume of each branch
	std::vector<double> hydr_resist1;						//!< Branches 1st order hydraulic resistance
	std::vector<double> hydr_resist2;						//!< Branches 2nd order hydraulic resistance
	std::vector<double> positions;							//!< Position of temperatures in branch
	std::vector<double> temperatures;						//!< Temperatures in branch [K]
	std::vector<bool> branch_cycle;							//!< Branche cycle
	std::vector<double> head_loss;							//!< Branches 0 order hydraulic resistance [m.c.f]
	bool flag_pump_speed_change;							//!< Pump's speed changed indicator
	bool flag_thermostat_changes;							//!< Thermostat changed indicator
	std::vector<int> nodes_id;								//!< Objects identificator conected to each node 
	std::vector<int> nodes_ind;								//!< Objects index conected to each node 
	int temp_action;										//!< Type of temperature change for the new volume inserted
	bool flag_valve_change;									//!< Valve changed indicator
	bool flag_first_exec;									//!< Firs execution indicator
	double gravity_acc;										//!< Gravity
	double obj_pos;											//!< Object position
	double start_pos;										//!< Start position
	double end_pos;											//!< End position
	double value;											//!< Key value to determine the new volume temperature
	double volume;											//!< New volume's volum [m ^ 3]
	int n_branch;											//!< Branches number
	int n_nodes;											//!< Nodes number
	int n_mesh;												//!< Meshes number
	int n_obj;												//!< Objects number	
	int n_tanks;											//!< Tanks number
	std::vector<double> bound_flows;						//!< Bound flows [L / s]
	std::vector<int> tanks_ind;								//!< Tanks index
	std::vector<double> head_vol_pump;						//!< Volumetric pump head preassure [m.c.f.]
	double weight_fraction;									//!< weight fraction




	/*!
	*\brief Reads the data file and stores the information im the vector Objects 
	* and aech class vector
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_ReadObj();
		
	/*!
	*\brief Using the data from HydroNet_ReadObj() generates the hydraulic net model data vectors
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void  HydroNet_Components();

	/*!
	*\brief Calls HydroNet_ReadObj() and HydroNet_Componets() and creates the rest of the data variables 
	* used by the other functions.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_Create();

	/*!
	*\brief Using the data from HydroNet_Create() and the initialitation in HydroNet_Maun() solves 
	* the hydralic net.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_Executions();

	/*!
	*\brief Calculates the temperature inside aech branche.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_Temperature();

	/*!
	*\brief Returns the temperature at a branch position.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*\pram obj_pos				Position wich temperature is needed.
	*\pram branch_temp_pos		Vector contanining the position in wich the temperature changes in the branch.
	*\pram branch_temperature	Vector contanining the diferent temperatures in the branch [K].
	*
	*/

	double HydroNet_GetPosTemperature(double obj_pos, std::vector<double> branch_temp_pos, std::vector<double> branch_temperature);

	/*!
	*\brief Returns the temperature at an object.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*\pram obj_index			PObject's index.
	*\pram branch_temp_pos		Vector contanining the position in wich the temperature changes in the branch.
	*\pram branch_temperature	Vector contanining the diferent temperatures in the branch [K].
	*\pram obj_intlet_pos		Vector whith the objects inlet position
	*\pram obj_outlet_pos		Vector whith the objects outlet position
	*
	*/

	double HydroNet_GetObjTemperature(double obj_index, std::vector<double> branch_temp_pos, std::vector<double> branch_temperature);

	/*!
	*\brief Inserts a new volume whth it's temperature
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*\pram positions			Vector contanining the position in wich the temperature changes in the branch.
	*\pram temperatures			Vector contanining the diferent temperatures in the branch [K].
	*\pram start_pos			New volume's start position.
	*\pram end_pos				New volume's end position.
	*\pram temp_action			Way in wich the new volumes's temperature in calculated.
	*\pram fluid_type			Fluid type (e.g. water).
	*\pram value				Key value to calculate the new volumes's temperature.
	*\pram volume				New volume's volume [m ^ 3].
	*\pram max_divisions		Maximum branch divisons.
	*\pram size_position		Vector positions size.
	*\pram size_temperature		Vector temperatures size.
	*
	*/

	strInsertVolum HydroNet_InsertVolume(std::vector<double> positions, std::vector<double> temperatures, double start_pos, double end_pos, int temp_action, std::string fluid_type, double value, double volume, int max_divisions, int size_position, int size_temperature);

	/*!
	*\brief Generates and solves the hydraulic net equations system.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_FlowSolver();

	/*!
	*\brief Calculate the aech object's hydraulic resistance.
	*
	*\author J. Salvador and J. Monfort
	*\date 13/10/2016
	*
	*/

	void HydroNet_HeadLoss();


};

#endif