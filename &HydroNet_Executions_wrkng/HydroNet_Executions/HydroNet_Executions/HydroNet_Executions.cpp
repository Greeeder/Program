// HydroNet_Executions.cpp : Defines the entry point for the console application.
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
	double pump_head;				// Head preassure
	double pump_power;
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
	double pump_head;				// Head preassure
	double pump_power;

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

template < class T>
T flip(T A) {
	T aux;
	aux.resize(A.size());
	for (int i = 0; i < A.size(); i++)
		aux.at(i) = A.at(A.size() - i - 1);
	return aux;
}


void HydroNet_HeadLoss() {}

void HydroNet_Executions(std::string fluid_type, std::vector<strBranches> branches, std::vector<strPipe>pipe, std::vector<strValve_fix>valve_fix, std::vector<strValve_var>valve_var, std::vector<strThermostat> thermostat, std::vector<strPump_volum> pump_volum, std::vector<strPump_turbo>pump_turbo, std::vector<strHeat_exch_Tout> heat_exch_Tout, std::vector<strHeat_exch_fix> heat_exch_fix, /*std::vector<strTank> tank,*/ std::vector<int> nodes_ind, std::vector<int> nodes_id, int n_nodes, std::vector<std::vector<int>>branches_ind, std::vector<std::vector<int>>branches_id, std::vector<bool> branch_cycle, std::vector<std::vector<int>> mesh_branches, std::vector<std::vector<int>> node_branches, int n_mesh, int n_branch,/*int n_tanks,*/ std::vector<double> branch_volume, std::vector<std::vector<double>> obj_inlet_pos, std::vector<std::vector<double>> obj_outlet_pos, std::vector<std::vector<int>> branch_htx_Tout, std::vector<std::vector<int>> branch_htx_fix, std::vector<std::vector<double>> branch_temp_pos, std::vector<std::vector<double>> branch_temperature, double gravity_acc, double dt, bool flag_pump_speed_change, bool flag_valve_change, bool flag_thermostat_changes, bool flag_first_exec, std::vector<double> head_loss, std::vector<double> hydr_resist1, std::vector<double>hydr_resist2, std::vector<double> flows) {
	
	std::vector<double>bound_flows;
	bool solver_flag=false;
	double head_aux, pos_aux, temperature_aux, obj_inlet_pos_aux;
	std::vector<double> head_vol_pump, hundred;
	hundred.assign(100, n_branch);
	// Main execution of the HydroNet model


	//// HEAD LOSS IN BRANCHES
	// Calculates head loss in every branch

	if (flag_first_exec || flag_valve_change || flag_thermostat_changes)

		HydroNet_HeadLoss(branches, pipe, valve_fix, valve_var, thermostat, pump_turbo, heat_exch_fix, heat_exch_Tout,/* tank,*/ branches_ind, branch_cycle, n_branch, obj_inlet_pos, obj_outlet_pos, branch_temp_pos, branch_temperature);




	//// FLOWS IN BRANCHES and HEAD IN VOLUMETRIC PUMPS

	// Sets boundary conditions: flow = 0
	// - All branches that are not inside a closed cycle
	// - All branches that contain closed valves or thermostats

	bound_flows.assign(NAN, n_branch); // vector filled with NaN to store boundary conditions
	for (int count = 0; count < bound_flows.size(); count++)
		if (!branch_cycle.at(count))
			bound_flows.at(count) = 0;
	//bound_flows(not(branch_cycle)) = 0; // outside a closed cycle
	for (int count = 1, obj_aux, branch_aux;count < valve_var.size();count++) {
		if (valve_var.at(count).opening == 0){
			obj_aux = valve_var.at(count).valve_var_id;
			for (int count1 = 0; count1 < branches_ind.size(); count1++)
				for (int count2 = 0; count2 < branches_ind.at(count1).size(); count2++)
					if (branches_ind.at(count1).at(count2) == obj_aux) {
						branch_aux = count1;
						break;
						break;
					}
			//branch_aux = objects.branch{ obj_aux };
			bound_flows.at(branch_aux) = 0;
		}
	}
	for (int count = 1, obj_aux, branch_aux; count < thermostat.size(); count++) {
		if (thermostat.at(count).opening == 0) {
			obj_aux = thermostat.at(count).themostat_id;
			for (int count1 = 0; count1 < branches_ind.size(); count1++)
				for (int count2 = 0; count2 < branches_ind.at(count1).size(); count2++)
					if (branches_ind.at(count1).at(count2) == obj_aux) {
						branch_aux = count1;
						break;
						break;
					}
			//branch_aux = objects.branch{ obj_aux };
			bound_flows.at(branch_aux) = 0;
		}
	}


	// Loop to solve flows
	if (flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes) {
		solver_flag = true;

		while (solver_flag) {

			solver_flag = false;

			HydroNet_FlowSolver(objects, pump_volum, branches_id, bound_flows, nodes_id, mesh_branches, node_branches, n_mesh, n_branch, head_loss, hydr_resist1, hydr_resist2, n_tanks);

			// Checks whether the obtained head of volumetric pumps are higher than the maximum
			// head of each volumetric pump. In such case, assign flow = 0 as boundary
			// condition and call the solver again
			for (int count = 0, branch_aux; count < pump_volum.size(); count++) {// counter of volumetric pumps
				branch_aux = objects.branch.at(pump_volum.at(count).Pump_volum); // pump branch
				pump_volum.at(count).pump_head = head_vol_pump.at(branch_aux); // store in table
				if (abs(pump_volum.at(count).pump_head) > pump_volum.at(count).head_max) {
					// calculated head is higher than maximum head
					bound_flows.at(branch_aux) = 0; // assign flow = 0 as boundary condition
					solver_flag = true; // re-calculate flows
				}
			}
		}



		// When needed, branch directions are flipped so that flow and branch directions are the same
		for (int count_branch = 0 ;count_branch< flows.size; count_branch++)
			if (flows.at(count_branch) < 0) {
			branches_ind.at(count_branch ) = flip(branches_ind.at(count_branch ));
			branches_id.at(count_branch ) = flip(branches_id.at(count_branch ));
			branch_temp_pos.at(count_branch ) = hundred - flip(branch_temp_pos.at(count_branch ));
			branch_temperature.at(count_branch ) = flip(branch_temperature.at(count_branch ));
			branch_htx_fix.at(count_branch ) = flip(branch_htx_fix.at(count_branch ));
			branch_htx_Tout.at(count_branch )
				= flip(branch_htx_Tout.at(count_branch ));

			for (int count = 0, count_obj; count < branches_ind.at(count_branch).size; count++) {
				std::vector<double> obj_inlet_pos_aux;
				count_obj = branches_ind.at(count_branch).at(count);
				obj_inlet_pos_aux = obj_inlet_pos.at( count_obj ); // saves this vector before overwriting it
				obj_inlet_pos.at(count_obj ) = 100 - obj_outlet_pos.at(count_obj ); // oulets are now inlets
				obj_outlet_pos.at(count_obj ) = 100 - obj_inlet_pos_aux; // inlets are now outlets
			}
		}
		flows = abs(flows); // All flows are positive now

	}



	//// PUMP POWER

	if (flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes){

		// Volumetric pumps
		for (int count = 0, object_aux; count < pump_volum.size(); count++) { // counter of volumetric pumps
		object_aux = pump_volum.obj_index(count); // index of pump in the object's list
		branch_aux = objects.branch.at(object_aux); // branch that contains the pump (only one)
		pos_aux = obj_inlet_pos.at(object_aux) +obj_outlet_pos.at(object_aux) / 2;
		temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos.at(branch_aux), branch_temperature.at(branch_aux));
		// temperature at the pump's middle

		pump_volum.at(count).pump_power = density(fluid_type, temperature_aux) * gravity_acc * head_vol_pump.at(branch_aux) * flows(branch_aux) / efficiency(head_vol_pump(branch_aux), flows(branch_aux));
		}

	// Turbopumps
		for (int count = 0, object_aux, branch_aux; count < pump_turbo.size(); count++) { // counter of turbopumps
		object_aux = pump_turbo.at(count).Pump_turbo; // index of pump in the object's list
		branch_aux = objects.branch.at( object_aux ); // branch that contains the pump (only one)
		pos_aux = obj_inlet_pos.at( object_aux ) +obj_outlet_pos.at( object_aux ) / 2;
		temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos.at( branch_aux ), branch_temperature.at(branch_aux ));
		// temperature at the pump's middle

		// Head curve of pump: head(speed, flow)
		head_aux = (pump_turbo.at(count).coef_N0 + pump_turbo.at(count).coef_N1 * pump_turbo.at(count).pump_speed +pump_turbo.at(count).coef_N2 * pump_turbo.at(count).pump_speed ^ 2) + flows.at(branch_aux) * (pump_turbo.at(count).Q_coef_N0 + pump_turbo.at(count).Q_coef_N1 * pump_turbo.at(count).pump_speed + pump_turbo.at(count).Q_coef_N2 * pump_turbo.at(count).pump_speed ^ 2) + flows.at(branch_aux) ^ 2 * (pump_turbo.at(count).Q2_coef_N0 + pump_turbo.at(count).Q2_coef_N1 * pump_turbo.at(count).pump_speed + pump_turbo.at(count).Q2_coef_N2 * pump_turbo.at(count).pump_speed ^ 2);

		pump_turbo.at(count).pump_power = density(fluid_type, temperature_aux) * gravity_acc * head_aux * flows.at(branch_aux) / efficiency(head_vol_pump.at(branch_aux), flows.at(branch_aux));
	}
	}



	//// TEMPERATURES

	HydroNet_Temperature(fluid_type, objects, obj_inlet_pos, obj_outlet_pos, heat_exch_Tout, heat_exch_fix, nodes_ind, branches_ind, node_branches, n_nodes, n_branch, branch_volume, branch_temp_pos, branch_temperature, branch_htx_Tout, branch_htx_fix, flows, dt);


}

int main()
{
    return 0;
}

