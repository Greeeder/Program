// HydroNet_Main.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string>		// string use; necessary for std::stoi
#include <vector>


struct strPipe
{
	int pipe_id;					// Identifier
	std::string name;				// Pipe's name
	int inlet_object;				// Pipe's inlet object identifier
	int outlet_object;				// Pipe's outlet object identifier
	double length;					// Pipe's length [meter]
	double diameter;				// Pipe's diameter [meter]
	double friction_coef;			// Pipe's friction coefficient

};

struct strValve_fix
{
	int valve_fix_id;				// Identifier
	std::string name;				// Valve's name	
	int inlet_object;				// Valve's inlet object identifier
	int outlet_oject;				// Valve's outlet object identifier
	double head_loss;				// Head loss [m.c.f. / m3/s]
};

struct strValve_var
{
	int valve_var_id;				// Identifier
	std::string name;				// Valve's name	
	int inlet_object;				// Valve's inlet object identifier
	int outlet_oject;				// Valve's outlet object identifier
	double opening;					// Opening [degree]
	double coef0;					// Coefficent order 0 [m.c.f.]
	double coef1;					// Coefficent order 1 [m.c.f. / degree]
	double coef2;					// Coefficent order 2 [m.c.f. * degree^(-2)]
	double coef3;					// Coefficent order 3 [m.c.f. * degree^(-3)]
};

struct strThermostat
{
	int themostat_id;					// Identifier
	std::string name;					// Thermostat's name	
	int inlet_object;					// Thermostat's inlet object identifier
	int outlet_oject;					// Thermostat's outlet object identifier
	double coef_T_0;					// Coefficent order 0 [-]
	double coef_T_1;					// Coefficent order 1 [K^(-1)]
	double coef_T_2;					// Coefficent order 2 [K^(-2)]
	double coef_T_3;					// Coefficent order 3 [K^(-3)]
	double coef_h_0;					// Coefficent order 0 [m.c.f.]
	double coef_h_1;					// Coefficent order 1 [m.c.f. / m3/s]
	double coef_h_2;					// Coefficent order 2 [m.c.f. / (m3/s)^2]
	double coef_h_3;					// Coefficent order 3 [m.c.f. / (m3/s)^3]
	double opening;						// Thermostat valv's opening
};

struct strPump_turbo
{
	int Pump_turbo;					// Identifier
	std::string name;				// Pump's name	
	int inlet_object;				// Pups's inlet object identifier
	int outlet_oject;				// Pump's outlet object identifier
	double coef_N0;					//coefficient order 0 [m.c.f]
	double coef_N1;					//coefficient order 0 N [m.c.f / (rad/s)]
	double coef_N2;					//coefficient order 0 N2 [m.c.f / (rad/s)^2]
	double Q_coef_N0;				//coefficient order 1 [m.c.f / m3*s]
	double Q_coef_N1;				//coefficient order 1 N [m.c.f / m3*s / (rad/s)]
	double Q_coef_N2;				//coefficient order 1 N2 [m.c.f / m3*s / (rad/s)^2]
	double Q2_coef_N0;				//coefficient order 2 [m.c.f /  m^6*s^2]
	double Q2_coef_N1;				//coefficient order 2 N [m.c.f /  m^6*s^2 / (rad/s)]
	double Q2_coef_N2;				//coefficient order 2 N2 [m.c.f. / m^6*s^2 / (rad/s)^2]
	double pump_speed;				// Pump speed [rad/s]
	double head_max;				// Max head preassure

};

struct strPump_volum
{
	int Pump_volum;					// Identifier
	std::string name;				// Pump's name	
	int inlet_object;				// Pups's inlet object identifier
	int outlet_oject;				// Pump's outlet object identifier
	double coef0;					//coefficient order 0 [m^3/s]
	double coef1;					//coefficient order 1 [m^3/s / (rad/s)]
	double coef2;					//coefficient order 2 [ m^3/s / (rad/s)^2]
	double coef3;					//coefficient order 3 [m^3/s / (rad/s)^3]
	double pump_speed;				// Pump speed [rad/s]
	double head_max;				// Max head preassure

};

struct strHeat_exch_fix
{
	int heat_exch_fix;				// Identifier
	std::string name;				// Heat Exchanger name
	int inlet_object;				// Heat Exchanger's inlet object identifier
	int outlet_object;				// Heat Exchanger's outlet object identifier
	double heat;					// Heat Energy Exchanged [J]
	double T_in;					// Intlet temperaure [K]
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};

struct strHeat_exch_Tout
{
	int heat_exch_Tout;				// Identifier
	std::string name;				// Heat Exchanger name
	int inlet_object;				// Heat Exchanger's inlet object identifier
	int outlet_object;				// Heat Exchanger's outlet object identifier
	double T_out;					// Outlet temperaure [K]
	double T_in;					// Intlet temperaure [K]
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};
/*struct strTank
{
int tank;					// Identifier
std::string name;				// Tank's name
int inlet_object;				// Tank's inlet object identifier
int outlet_oject;				// Tank's outlet object identifier
double head_loss;				// Head loss [m.c.f. / m3/s]
};
*/
struct strObjects
{
	int ID;							// Object's idintificator
	std::string name;				// Object's name
	std::string Class;				// Object's clas
	double volume;					// Objects volume
	int inlet_object;				// Tank's inlet object identifier
	int outlet_oject;				// Tank's outlet object identifier
};
struct strBranches
{
	int branch_ind;				// branch idintificator
	int branches_cycle;
	//int n_branch;				// Number of branches
	std::vector<strObjects> objects;	// Objects in Brnach
};
double HydroNet_GetObjTemperature(double obj_pos, std::vector<double> branch_temp_pos, std::vector<double> branch_temperature) {

	double distance = -1000.0;
	double temperature = -500.0;
	int count = 1;


	//Returns temperature at the specified position

	// Index of the temperature position closest to the specified positionand the difference of position between them

	if (obj_pos == branch_temp_pos.at(1))
		temperature = branch_temperature.at(1);
	else if (obj_pos == branch_temp_pos.at(branch_temp_pos.size() - 1))
		temperature = branch_temperature.at(branch_temperature.size() - 1);
	else {
		while (distance < 0.0) {
			distance = branch_temp_pos.at(count + 1) - obj_pos;		// count+1 because the fert positin is checked before
			count++;
		}
		temperature = branch_temperature.at(count - 1);				// -1 because vector estarts at 0
	}
	return temperature;


}



void HydroNet_Main(std::vector<strPump_volum> pump_volum, std::vector<strPump_turbo> pump_turbo, std::vector<strValve_var> valve_var, std::vector<strThermostat> thermostat, std::vector<strBranches> branches, std::vector<std::vector<double>> obj_inlet_pos, std::vector<std::vector<double>> obj_outlet_pos) {

	std::string circuit_file;
	double gravity_acc;
	std::vector<std::vector<double>> branch_temp_pos,branch_temperature;
	std::vector<double> head_loss, hydr_resist1, hydr_resist2, flows, dt_vec, volum_pump_speed, turbo_pump_speed, valve_opening, thst_sensitivity, thermostat_opening_new;
	double n_branch,dt, obj_pos, thermostat_temp;
	int n_exec;
	bool flag_first_exec, flag_pump_speed_change, flag_valve_change, flag_thermostat_changes;

	// Controls data flow


	//// MODEL CREATION
	// All valves and thermostats are assumed open.

	circuit_file = "Coolant_basic_nodes.txt";
	//circuit_file = 'Circuit_test_1.txt';
	
	///[ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, heat_exch_fix, tank, nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, branch_htx_fix] = HydroNet_Create( circuit_file );



	//// INITIALIZATION

	gravity_acc = 9.81; // gravity acceleration [m/s^2] for pump power

	// Cell to store volumetric position of the points where temperature changes,
	// as a // of the total volume of the branch // temperature between those points is constant
	// branch_temp_pos = cell(n_branch, 1);
	branch_temp_pos.resize(7);
	branch_temp_pos.at(0) = { 0,10,100 };
	branch_temp_pos.at(1) = { 0, 100};
	branch_temp_pos.at(2) = { 0, 100};
	branch_temp_pos.at(3) = { 0, 20, 100};
	branch_temp_pos.at(4) = { 0, 30, 100};
	branch_temp_pos.at(5) = { 0, 100};
	branch_temp_pos.at(6) = { 0, 10, 100};


	// Cell to store temperatures of every volume section in the branch
	// branch_temperature = cell(n_branch, 1);
	branch_temperature.resize(7);
	branch_temperature.at(0) = { 300,250 };
	branch_temperature.at(1) = {350};
	branch_temperature.at(2) = {400};
	branch_temperature.at(3) = {450, 300};
	branch_temperature.at(4) = {500, 400};
	branch_temperature.at(5) = {550};
	branch_temperature.at(6) = {600, 550};

	// Output variables related to branches
	head_loss.assign(n_branch, 0);
	hydr_resist1.assign(n_branch, 0);
	hydr_resist2.assign(n_branch, 0);
	flows.assign(n_branch, 0);



	//// EXECUTIONS
	// Some variables must be recalculated only if other change

	n_exec = 5;

	dt_vec = { 0.1, 0.1, 0.1, 0.1, 0.2 };

	volum_pump_speed = {105, 210, 314, 419, 419};
	turbo_pump_speed = {105, 210, 314, 419, 419};
	valve_opening = {20, 50, 70, 80, 30};

	thst_sensitivity.assign(n_exec, 5); // thermostat sensitivity [%]


	// Execution loop
	flag_first_exec = true;
	for (int count_exec = 0 ; count_exec< n_exec;count_exec++){

		// TIME SPAN
		dt = dt_vec.at(count_exec);


	// PUMP SPEED
	// If the speed of any pump changes, flows must be recalculated

		flag_pump_speed_change = false;
		for (int count=0; count<pump_volum.size; count++)
			if  ( volum_pump_speed.at(count_exec) != pump_volum.at(count).pump_speed )
				flag_pump_speed_change = true;
		
		for (int count = 0; count<pump_turbo.size; count++)
			if (turbo_pump_speed.at(count_exec) != pump_turbo.at(count).pump_speed)
				flag_pump_speed_change = true;

		for (int count = 0; count<pump_volum.size; count++)
		pump_volum.at(count).pump_speed = volum_pump_speed.at(count_exec);
		for (int count = 0; count<pump_turbo.size; count++)
		pump_turbo.at(count).pump_speed = turbo_pump_speed.at(count_exec);


	// VALVES
	// If the valve opening varies, vectors for head losses and thus flows must be recalculated.
	// If the valve is closed, flow in its branch is zero.

		flag_valve_change = false;
		for (int count = 0; count<valve_var.size; count++)
			if (valve_opening.at(count_exec) != valve_var.at(count).opening)
				flag_valve_change = true;
	
		for (int count = 0; count<valve_var.size; count++)
		valve_var.at(count).opening = valve_opening.at(count_exec);


	// THERMOSTATS
	// If the thermostat opening varies because there is a change in the inlet
	// temperature, vectors for head losses and thus flows must be recalculated.
	// Head losses and flows will be recalculated only if the opening variation
	// of any thermostat is larger than the sensitivity of the thermostat[%]
	// If the thermostat is closed, flow in its branch is zero.

		thermostat_opening_new.assign(thermostat.size(),0);
		for (int count_thst = 0, obj_aux, branch_aux; count_thst < thermostat.size();count_thst++) {
			obj_aux = thermostat.at(count_thst).inlet_object;
			for (int count=0; count<branches.size();count++)
				for (int count1 = 0; count<branches.at(count).objects.size(); count++)
					if( branches.at(count).objects.at(count1).ID == obj_aux )
						branch_aux =count ;
			obj_pos = (obj_inlet_pos.at(obj_aux).at(0) + obj_outlet_pos.at(obj_aux).at(0)) / 2;

			thermostat_temp = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos.at(branch_aux ), branch_temperature.at(branch_aux ));

			thermostat_opening_new.at(count_thst) = thermostat.at(count_thst).coef_T_0 + thermostat.at(count_thst).coef_T_1 * thermostat_temp + thermostat.at(count_thst).coef_T_2 * thermostat_temp* thermostat_temp;
		}
	
		flag_thermostat_changes = false;
		for (int count=0; count< thermostat.size();count++)
			if (abs(thermostat_opening_new.at(count) - thermostat.at(count).opening) / thermostat.at(count).opening * 100 > thst_sensitivity.at(count_exec))
				flag_thermostat_changes = true;

		for (int count = 0; count< thermostat.size(); count++)
			thermostat.at(count).opening = thermostat_opening_new.at(count);


	// MAIN EXECUTION FUNCTION CALL
		[flows, branches_ind, branches_id, branch_temp_pos, branch_temperature, branch_htx_fix, branch_htx_Tout, ...
		heat_exch_Tout, heat_exch_fix, obj_inlet_pos, obj_outlet_pos] = HydroNet_Executions(fluid_type, objects, ...
		pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, heat_exch_fix, tank, ...
		nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, ...
		n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, branch_htx_fix, ...
		branch_temp_pos, branch_temperature, gravity_acc, dt, flag_pump_speed_change, flag_valve_change, ...
		flag_thermostat_changes, flag_first_exec, head_loss, hydr_resist1, hydr_resist2, flows);

		flag_first_exec = false;

		display(flows);



	}

	
}

int main(){


    return 0;
}

