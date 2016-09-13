// HydroNet_Temperature.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

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
	double heat;					// Heat Exchanged
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

void HydroNet_InsertVolume() {
}

double HydroNet_GetObjTemperature() {

}

double vec_abs_sum(std::vector<double> in) {
	double sum = 0;
	for (int i = 0; i < in.size(); i++)
		sum = sum + abs(in.at(i));
	return sum;
}
double vec_sum_elements(std::vector<double> in,std::vector<int> elem) {
	double sum = 0;
	for (int i = 0; i < elem.size(); i++)
		sum = sum + in.at(elem.at(i));
	return sum;
}

double density(/*std::string fluid_type, double temperature*/)
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 1000;	// for tests: water [Kg/m3]
}




void HydroNet_Temperature(std::string fluid_type, std::vector<strBranches> branches, std::vector<double> obj_inlet_pos, std::vector<double> obj_outlet_pos, std::vector<strHeat_exch_Tout> heat_exch_Tout, std::vector<strHeat_exch_fix>heat_exch_fix, std::vector<int> nodes_id, std::vector<std::vector<int>>branches_id, std::vector<std::vector<int>> node_branches, int n_nodes, int n_branch, std::vector<double> branch_volume, std::vector<std::vector<double>>branch_temp_pos, std::vector<std::vector<double>>branch_temperature, std::vector<std::vector<int>>branch_htx_Tout, std::vector<std::vector<int>>branch_htx_fix, std::vector<double> flows, double dt) {

	// Finds positions where temperature changes and temperatures in those ranges.
	// Moves positions and temperatures according to flows and, simultaneously,
	// enters heat from heat exchangers. Returns refreshed ranges and temperatures
	// in every branch as well as temperatures for the inlet of heat exchangers.
	// Merges ranges (volumes) and its temperatures if there are too many of them in a branch.

	int max_divisions, position_count, position, node_count, end_node_ind, start_node_ind, obj_pos, node_aux, pos_aux;
	std::vector<std::vector<double>>new_pos, new_temp;
	std::vector<std::vector<int>>  node_inlet_branches, node_outlet_branches;
	std::vector<double> overflow, overflow_temperature, overflow_temp_volpos_aux, overflow_temperature_aux;
	double delta_pos, start_pos, end_pos, value, temp_action, sum_aux, temp_aux, flow_sum, extra_div, sum_mass, sum_temp_mass;
	bool flows_flag = false;
	std::vector<bool>flag_inlet, to_remove;
	std::vector<double> node_volume_in, node_mass, node_temperature, volume_share, new_pos_aux, aux_vec;
	std::vector<int>branch_ind_order_aux;




	// Maximum divisions inside each branch; if there are more, some are merged
	max_divisions = 20;



	/// VOLUME THAT STAYS INSIDE THE SAME BRANCH AND OVERFLOWING VOLUME

	new_pos.resize(n_branch); // temporary cell to store positions in the current instant
	new_temp.resize(n_branch); // temporary cell to store temperatures in the current instant
	overflow.assign(n_branch, 0); // vector to store volumes that move out from its original branch
	for (int count = 0; count < branch_volume.size(); count++)
		if (branch_volume.at(count) > 0)
			overflow.at(count) = flows.at(count) * dt;

	overflow_temperature.assign(n_branch, 0);//ones(n_branch, 1) * (-999999); // vector to store temperatures of overflowed volumes


	// Branches without volume
	for (int count_branch = 0; count_branch < branch_volume.size(); count_branch++) {//count_branch = transpose(find(branch_volume == 0))
		if (branch_volume.at(count_branch) == 0) {
			new_pos.at(count_branch).insert(new_pos.at(count_branch).begin(), 0);
			new_temp.at(count_branch) = branch_temperature.at(count_branch);
		}
	}


	// Branch loop excluding branches without volume
	for (int count_branch = 0; count_branch < branch_volume.size(); count_branch++) {//count_branch = transpose(find(branch_volume == 0))
		if (branch_volume.at(count_branch) > 0) {
			delta_pos = flows.at(count_branch) * dt / branch_volume.at(count_branch) * 100;
			// \% of branch volume that flow moved

		// NO FLOW MOVEMENT: refreshes temperatures inside the branch
			if (delta_pos == 0) {// flow is zero, no movement

				new_pos.at(count_branch) = branch_temp_pos.at(count_branch);
				new_temp.at(count_branch) = branch_temperature.at(count_branch);

				// If there are heat exchangers, inserts volume for it and refreshes temperature
				// Heat exchanger of type T_out
				for (int count_obj_htx = 0, ind_aux; count_obj_htx < branch_htx_Tout.at(count_branch).size(); count_obj_htx++) {//count_obj_htx = branch_htx_Tout{ count_branch }

					ind_aux = heat_exch_Tout.at(count_obj_htx).heat_exch_Tout;
					// index of the heat exchanger in the table objects

					start_pos = obj_inlet_pos.at(ind_aux); // object's inlet position
					end_pos = obj_outlet_pos.at(ind_aux); // object's outlet position
					temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = heat_exch_Tout.at(count_obj_htx).T_out; // final temperature

					HydroNet_InsertVolume(/*new_pos{ count_branch }, new_temp{ count_branch }, start_pos, end_pos, temp_action, fluid_type, value*/);
				}
				// Heat exchanger of type fixed heat
				for (int count_obj_htx = 0, ind_aux; count_obj_htx < branch_htx_fix.at(count_branch).size(); count_obj_htx++) {//count_obj_htx = branch_htx_fix{ count_branch }

					ind_aux = heat_exch_fix.at(count_obj_htx).heat_exch_fix;
					// index of the heat exchanger in the table objects

					start_pos = obj_inlet_pos.at(ind_aux); // object's inlet position
					end_pos = obj_outlet_pos.at(ind_aux); // object's outlet position
					temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = heat_exch_fix.at(count_obj_htx).heat; // heat to add

					HydroNet_InsertVolume(/*new_pos{ count_branch }, new_temp{ count_branch }, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch)*/);
				}
			}
		}



		////////////////////////////////////
		// THERE IS FLOW MOVEMENT INSIDE THE BRANCH
		else

			position_count = 2; // position counter

			// REFRESHES TEMPERATURES INSIDE THE BRANCH WITHOUT HEAT EXCHANGE
		while (((branch_temp_pos.at(count_branch).at(position_count) + delta_pos) <= 100) && ((position_count) < (branch_temp_pos.at(count_branch).size()))) {//and ((branch_temp_pos{ count_branch }(position_count)+delta_pos) <= 10(position_count) < numel(branch_temp_pos{ count_branch })) // old position exists inside the vector
		// Positions are moved along the branch according to flowing flow

			position = branch_temp_pos.at(count_branch).at(position_count);

			// Flow must be positive (same direction as branch)
			new_pos.at(count_branch).push_back(position + delta_pos); // add new position at the end of the vector
			new_temp.at(count_branch).push_back(branch_temperature.at(count_branch).at(position_count - 1));

			position_count = position_count + 1;
		}
		// New position is in another branch OR there are no more positions in the branch
		// Closes branch
		new_pos.at(count_branch).insert(new_pos.at(count_branch).begin(), 0); // add vector start 
		new_pos.at(count_branch).push_back(100); // add vector end
		new_temp.at(count_branch).push_back(branch_temperature.at(count_branch).at(position_count - 1));
		// add new temperature at the end of the vector


	////////////////////////////////////
	// VOLUME OVERFLOW: makes volume buffer without taking into account heat exchangers
	// First step to find the temperature of the buffer volume
		if (position_count <= branch_temp_pos.at(count_branch).size()) {
			// old position exists inside the vector, thus after the previous in-branch
			// while loop the new position has to be in another branch: OVERFLOW

			// Makes a new pair of vectors like branch_temp_pos and branch_temperature
			// for the volume that goes out of the branch
			overflow_temp_volpos_aux = { 100 }; // initialization: branch end
				// position specified as \// of the branch volume
			overflow_temperature_aux;

			for (int count = 0; count < branch_temp_pos.size(); count++)
				if (position_count == branch_temp_pos.at(count).size()) {

					overflow_temp_volpos_aux.push_back(branch_temp_pos.at(count).at(position_count) + delta_pos);
					// position specified as \// of the branch volume (> 100%)
					overflow_temperature_aux.push_back(branch_temperature.at(count).at(position_count - 1));
				}
		}


		////////////////////////////////////
		// HEAT EXCHANGERS
		// If there are heat exchangers, checks whether the volume out of them remains inside the branch or it
		// goes out, inserts volume for it, refreshes temperature inside the branch and stores overflowed volume.
		// There can only be heat exchangers if branch volume is > 0.
		// Second step to find the temperature of the buffer volume

		// Heat exchanger of type T_out
		for (int count = 0, count_obj_htx, ind_aux; count < branch_htx_Tout.at(count_branch).size(); count++) {// loop starting from the closest exchanger to the branch's end
			count_obj_htx = branch_htx_Tout.at(count_branch).at(count);

			ind_aux = heat_exch_Tout.at(count_obj_htx).heat_exch_Tout;
			// index of the heat exchanger in the table objects

			temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			value = heat_exch_Tout.at(count_obj_htx).T_out; // new temperature

			start_pos = obj_outlet_pos.at(ind_aux); // object's outlet position is the start position
			end_pos = start_pos + delta_pos;

			if (end_pos > 100) {
				start_pos = 100;

				HydroNet_InsertVolume(/*overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value*/);

				// For the volume inside the branch
				start_pos = obj_outlet_pos.at(ind_aux); // object's outlet position is the start position
				end_pos = 100;
			}

			HydroNet_InsertVolume(/*new_pos{ count_branch }, new_temp{ count_branch }, start_pos, end_pos, temp_action, fluid_type, value*/);
		}

		// Heat exchanger of type fixed heat
		for (int count = 0, count_obj_htx, ind_aux; count < branch_htx_fix.at(count_branch).size(); count++) { // loop starting from the closest exchanger to the branch's end
			count_obj_htx = branch_htx_fix.at(count_branch).at(count);
			ind_aux = heat_exch_fix.at(count_obj_htx).heat_exch_fix;
			// index of the heat exchanger in the table objects

			temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			value = heat_exch_fix.at(count_obj_htx).heat;

			start_pos = obj_outlet_pos.at(ind_aux); // object's outlet position is the start position
			end_pos = start_pos + delta_pos;

			if (end_pos > 100) {

				value = heat_exch_fix.at(count_obj_htx).heat * (end_pos - 100) / (end_pos - start_pos);
				// heat to add averaged for the volume that goes out of the branch

				start_pos = 100;

				HydroNet_InsertVolume(/*overflow_temp_volpos_aux,	overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch)*/);

				// For the volume inside the branch
				start_pos = obj_outlet_pos.at(ind_aux); // object's outlet position is the start position
				end_pos = 100;
				value = heat_exch_fix.at(count_obj_htx).heat * (100 - start_pos) / (end_pos - start_pos);
				// heat to add averaged for the volume that stays inside the branch
			}

			HydroNet_InsertVolume(/*new_pos{ count_branch }, new_temp{ count_branch }, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch)*/);
		}


		////////////////////////////////////
		// TEMPERATURE OF THE BUFFER VOLUME (OVERFLOW)
		// Final averaging of temperatures in the buffer volume

		start_pos = 100;
		end_pos = overflow_temp_volpos_aux.at(overflow_temp_volpos_aux.size());
		temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature

		HydroNet_InsertVolume(/*overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type*/);
	}




	//// SEPARES BRANCHES AT THE INLET AND OUTLET OF EVERY NODE
	// Separates branches that are inlets to nodes and those that are outlets to nodes from node_branches
	// All branches should start and end at a node!!
	for (int count = 0; count < flows.size(); count++)
		if (flows.at(count) > 0)
			flows_flag = true;

	if (flows_flag == true) { // if there is no movement in the branch, there is no need for this

		node_inlet_branches = node_branches;
		node_outlet_branches = node_branches;
		//node_volume_out = zeros(n_nodes, 1);
			// stores the sum of all volumes leaving every node to check if outlet branches are overflowed
			// no used because, later, volume_share is more precise
		node_count = 0;
		for (int count; node_count < n_nodes; node_count++) { // node loop
			flag_inlet.assign(false, node_branches.at(node_count).size()); // marks inlet branches to this node
			count = 0; // counter of branches in every node
			for (int count1 = 0, count_branch; count1 < node_branches.at(node_count).size(); count1++) {// alternative: branch_count = 1 : numel(node_branches{node_count})

				count_branch = node_branches.at(node_count).at(count1);
				count = count + 1;
				if (branches_id.at(count_branch).at(branches_id.at(count_branch).size()) == nodes_id.at(node_count))
					// flow goes into the node from the branch: it is an inlet branch to the node
					flag_inlet.at(count) = true; // mark branch as inlet to node
		//         else // flow goes into the branch from the node: it is an outlet branch from the node
		//             node_volume_out(node_count) = node_volume_out(node_count) + branch_volume(branch_count);
			}

			// Separates inlet and outlet branches to each node
			for (int aux = 0; aux < node_inlet_branches.size(); aux++)
				if (flag_inlet.at(aux) == true)
					node_inlet_branches.at(node_count).at(aux) = node_inlet_branches.at(node_count).at(aux);
				else
					node_inlet_branches.at(node_count).erase(node_inlet_branches.at(node_count).begin() + aux);
			for (int aux = 0; aux < node_outlet_branches.size(); aux++)
				if (flag_inlet.at(aux) == false)
					node_outlet_branches.at(node_count).at(aux) = node_outlet_branches.at(node_count).at(aux);
				else
					node_outlet_branches.at(node_count).erase(node_outlet_branches.at(node_count).begin() + aux);
		}
	}


	//// MANAGES OVERFLOWED VOLUMES
	// All branches should start and end at a node!!
	// The number of branches without volume must be minimized!!

	while (vec_abs_sum(overflow) != 0) { // while there are overflows

		// Makes an enthalpy balance of all volumes mixing in every node
		node_volume_in.assign(n_nodes, 0); // stores the sum of all volumes mixing in every node
		node_mass.assign(n_nodes, 0); // temporary vector to store mass in every node
		node_temperature.assign(n_nodes, 0); // stores temperatures in every node after mixing

		// Volume and mass in nodes
		for (int count = 0, count_branch; count < overflow.size(); count++) { // 1 : n_branch // branch loop
			if (overflow.at(count) > 0)
				count_branch = count;
			// do not consider branches that do not have overflows (first iteration: do not considerer nodes without volume)
			// branch_count = 1 : n_branch would be valid too (overflow = 0 does not add volume or mass)
			for (int aux = 0; aux < branches_id.at(count_branch).size(); aux++)
				if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size()))
					end_node_ind = aux; // index of the node at the end of the branch

			node_volume_in.at(end_node_ind) = node_volume_in.at(end_node_ind) + overflow.at(count_branch);
			node_mass.at(end_node_ind) = node_mass.at(end_node_ind) + density(/*fluid_type, overflow_temperature.at(count_branch)*/) * overflow.at(count_branch);
		}

		// Temperature in nodes, provisional
		// Processes nodes that have overflow buffer volumes > 0

		for (int count = 0, count_node; count < node_mass.size(); count++) { // to avoid dividing by zero
			if (node_mass.at(count) > 0)
				count_node = count;
			for (int count1 = 0, count_branch; count1 < node_inlet_branches.at(count_node).size(); count1++) {
				count_branch = node_inlet_branches.at(count_node).at(count1);
				for (int aux = 0; aux < branches_id.at(count_branch).size(); aux++)
					if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size()))
						end_node_ind = aux; // index of the node at the end of the branch
				node_temperature.at(end_node_ind) = node_temperature.at(end_node_ind) + overflow_temperature.at(count_branch) * density(/*fluid_type, overflow_temperature.at(count_branch)*/) *overflow.at(count_branch) / node_mass.at(end_node_ind);
				// T branch to node * branch mass to node / total mass to node
			}
		}

		// Control goes to buffer volumes in the nodes node_volume_in
		for (int count = 0; count < overflow.size(); count++) {
			overflow.at(count) = 0;
			overflow_temperature.at(count) = 0;
		}


		// Passes buffer volumes that are at the inlet node of branches with
		// volume = 0 to the end node

		sum_aux = 0; // sum of volumes at the start of branches without volume
		for (int count = 0, count_branch; count < flows.size(); count++) {
			if ((branch_volume.at(count) == 0) && (flows.at(count) > 0)) {
				count_branch = count;
				// loop of branches without volume but with flow movement

				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = nodes_id.at(aux);
				sum_aux = sum_aux + node_volume_in.at(start_node_ind);
			}
		}

		while (sum_aux > 0) {
			// there are buffer volumes at the start of branches without volume; they must be past through

			// Priority: the node at the start of the branch has no inlet branches
			// without volume (they also would get passed a volume from behind)

			for (int count = 0, count_branch; count < flows.size(); count++)
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0))
					branch_ind_order_aux.push_back(count); // branches without volume but with flow movement

			for (int count = 0, count_branch; count < flows.size(); count++) {
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0))
					count_branch = count;
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = nodes_id.at(aux);
				for (int aux = 0; aux < node_inlet_branches.at(start_node_ind).size(); aux++)
					if (branch_volume.at(node_inlet_branches.at(start_node_ind).at(aux)) == 0) {
						// any of the inlet branches to the start node has no volume
						for (int count1 = 0; count < branch_ind_order_aux.size(); count1++) {
							if (branch_ind_order_aux.at(count1) == count_branch) {
								branch_ind_order_aux.erase(branch_ind_order_aux.begin() + count1); // remove the branch
								branch_ind_order_aux.push_back(count_branch); // put the branch at the bottom of the list
							}
						}
					}
			}

			// Moves volumes according to priority
			for (int count = 0, count_branch; count < branch_ind_order_aux.size(); count++) {
				count_branch = branch_ind_order_aux.at(count);
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = nodes_id.at(aux);// index of the node at the start of the branch
				for (int aux = 0; aux < branches_id.at(count_branch).size(); aux++)
					if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size()))
						end_node_ind = aux; // index of the node at the end of the branch

				node_volume_in.at(end_node_ind) = node_volume_in.at(end_node_ind) + node_volume_in.at(start_node_ind);
				temp_aux = (node_temperature.at(end_node_ind) * node_mass.at(end_node_ind) + node_temperature.at(start_node_ind) * node_mass.at(start_node_ind)) / (node_mass.at(end_node_ind) + node_mass.at(start_node_ind));
				// balance between volume already present and incoming volume
				node_temperature.at(end_node_ind) = temp_aux;
				node_mass.at(end_node_ind) = node_mass.at(end_node_ind) + node_mass.at(start_node_ind);
				node_volume_in.at(start_node_ind) = 0;
				node_mass.at(start_node_ind) = 0;
				node_temperature.at(start_node_ind) = 0;
				new_temp.at(count_branch).push_back(temp_aux); // assigns temperature to the branch
			}

			sum_aux = 0; // sum of volumes at the start of branches without volume
			for (int count = 0, count_branch; count < flows.size(); count++) {
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0))
					count_branch = count;
				// loop of branches without volume but with flow movement
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = nodes_id.at(aux); // index of the node at the start of the branch
				sum_aux = sum_aux + node_volume_in.at(start_node_ind);
			}
		}



		// Distributes volume among branches
		for (int count = 0; count < node_volume_in.size(); count++)
			if (node_volume_in.at(count) > 0) { // loop of nodes that have a buffer volume
				node_count = count;

				flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count)); // sum of flows leaving the node
				for (int aux = 0; aux < node_outlet_branches.at(node_count).size(); aux++)
					volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
				// volume that goes to every branch that is an outlet to the node (proportional to flow)
				for (int count_branch = 0, branch_aux; count_branch < node_outlet_branches.at(node_count).size(); count_branch++) { // loop of branches that are outlets to the node
					branch_aux = node_outlet_branches.at(node_count).at(count_branch); // branch index
					if (volume_share.at(count_branch) > 0) {// processes only branches where volume actually goes in
						if (volume_share.at(count_branch) <= branch_volume.at(branch_aux)) {
							// volume in the branch is bigger than incoming volume
							// also avoids dividing by zero

							// Adds position and temperature at the beginning
							new_pos.at(branch_aux).insert(new_pos.at(branch_aux).begin(), (volume_share.at(count_branch) / branch_volume.at(branch_aux) * 100));
							new_temp.at(branch_aux).insert(new_pos.at(branch_aux).begin(), node_temperature.at(node_count));

							// Removes positions lower than or equal to the new one
							to_remove.assign(false, new_pos.at(branch_aux).size());
							for (int pos_aux = 1; pos_aux < new_pos.at(branch_aux).size(); pos_aux++) {
								if (new_pos.at(branch_aux).at(pos_aux) <= new_pos.at(branch_aux).at(0)) {
									to_remove.at(pos_aux) = true;
								}
							}
							for (int aux = 0; aux < to_remove.size(); aux++)
								if (to_remove.at(aux) == false)
									new_pos.at(branch_aux) = new_pos.at(branch_aux);

							// Does not remove temperature that is after the last position to remove
							for (int count1 = 0; count1 < to_remove.size(); count1++) /*: -1 : 1*/ {
								if (to_remove.at(count) == true) {
									to_remove.at(count) = false;
									break;
								}
							}
							to_remove.erase(to_remove.begin() + to_remove.size());
							for (int aux = 0; aux < to_remove.size(); aux++)
								if (to_remove.at(aux) == false)
									new_temp.at(branch_aux) = new_temp.at(branch_aux);

							// Adds position 0 at the beginning
							new_pos.at(branch_aux).insert(new_pos.at(branch_aux).begin(), 0);
						}

						else { // volume in the branch is smaller than incoming volume

							// Fills branch
							new_pos.at(branch_aux) = { 0, 100 };
							new_temp.at(branch_aux).at(0) = node_temperature.at(node_count);

							// Stores the excess volume
							overflow.at(branch_aux) = volume_share.at(count_branch) - branch_volume.at(branch_aux);
							overflow_temperature.at(branch_aux) = node_temperature.at(node_count);
						}
					}
				}
			}

		for (int count = 0; count < node_volume_in.size(); count++)
			node_volume_in.at(count) = 0; // just in case it is used after the while loop
				// control passes to overflow in branches
	}




	//// MERGES VOLUMES INSIDE A BRANCH
	// If there are more volumes than the value specified in max_divisions

	for (int count_branch = 0; count_branch < n_branch; count_branch++) {
		extra_div = (new_pos.at(count_branch).size() - 1) - max_divisions; // number of divisions over max_divisions
		while (extra_div > 0) {
			new_pos_aux = new_pos.at(count_branch);

			// Finds smaller combination of adjacent volumes
			for (int aux = 0; aux < new_pos.size() - 3; aux++)
				aux_vec.push_back(new_pos_aux.at(3 + aux) - new_pos_aux.at(aux));
			pos_aux = std::min_element(aux_vec.at(0), aux_vec.at(aux_vec.size()));

			// Obtains temperature by means of an enthalpy balance
			temp_aux = (new_temp.at(count_branch).at(pos_aux) * density(/*fluid_type, new_temp(pos_aux)*/) * (new_pos.at(count_branch).at(pos_aux + 1) - new_pos.at(count_branch).at(pos_aux)) +new_temp.at(count_branch).at(pos_aux + 1) * density(/*fluid_type, new_temp(pos_aux + 1)*/) * (new_pos.at(count_branch).at(pos_aux + 2) - new_pos.at(count_branch).at(pos_aux + 1))) / (density(/*fluid_type, new_temp(pos_aux)*/) * (new_pos_aux.at(pos_aux + 1) - new_pos_aux.at(pos_aux)) + density(/*fluid_type, new_temp(pos_aux + 1)*/) * (new_pos_aux.at(pos_aux + 2) - new_pos_aux.at(pos_aux + 1)));
			// Merges volumes
			new_pos.at(count_branch).erase(new_pos.at(count_branch).begin() + (pos_aux + 1)); // remove position in the middle
			new_temp.at(count_branch).erase(new_temp.at(count_branch).begin() + (pos_aux + 1)); // remove second temperature
			new_temp.at(count_branch).at(pos_aux) = temp_aux; // enter temperature from balance

			extra_div = extra_div - 1; // next division to merge
		}
	}



	//// CALCULATES INLET TEMPERATURE TO HEAT EXCHANGERS
	// Asumption: the flow that will pass through the heat exchanger in the next
	// instant is equal to the branch flow in the current instant.
	// Exception: if flow = 0 in the current instant, inserts a volume for the
	// heat exchanger and obtains its temperature.
	// Backwards direction from the outlet of the heat exchanger.

	// Heat exchangers of type Tout
	for (int count_htx = 0, obj_aux, branch_aux; count_htx < heat_exch_Tout.size(); count_htx++) {

		obj_aux = heat_exch_Tout.at(count_htx).heat_exch_Tout; // index of the heat exchanger in objects table
		//branch_aux = objects.branch.at(obj_aux); // branch that contains the heat exchanger

		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}
		if (flows.at(branch_aux) == 0) { // no movement in the branch

			// Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux);
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux }, start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the heat exchanger
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_Tout.at(count_htx).T_in = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux }, start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_Tout.at(count_htx).T_in = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
		}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			// First, obtains temperature of the volume between the start of the
			// branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux }, start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
			heat_exch_Tout.at(count_htx).T_in = temp_aux;

			// Initializes enthalpy balance
			sum_mass = density(/*fluid_type, temp_aux*/) * (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

			// Stores volume that will come from previous branches
			overflow.at(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (vec_abs_sum(overflow) != 0) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow.at(count) > 0)
						count_branch = count;

					node_aux = branches_id.at(count_branch).at(1); // node at the start of the branch
					flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count)); // sum of flows going into the branch

					for (int aux = 0; aux < node_outlet_branches.at(node_count).size(); aux++)
						volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
					// volume that comes from every inlet branch (proportional to flow)
					for (int count_branch1 = 0, branch_aux1; count_branch1 < node_inlet_branches.at(node_aux).size(); count_branch1++) { // loop of inlet branches
						branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
						if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
							// volume in the branch is bigger than outcoming volume
							// also avoids dividing by zero

							// Inserts a volume for the affected area and obtains its temperature
							start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
							end_pos = 100;
							temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
							HydroNet_InsertVolume(/*new_pos{ branch_aux1 }, new_temp{ branch_aux1 }, start_pos, end_pos, temp_action, fluid_type*/);

							// Reads the temperature in the middle of the affected volume
							obj_pos = (start_pos + end_pos) / 2;
							temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);

							// Enters temperature and mass in the enthalpy balance
							sum_mass = sum_mass + density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);
							sum_temp_mass = sum_temp_mass + temp_aux * density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);
						}

						else {// volume in the branch is smaller than outcoming volume

							// Inserts a volume for the entire branch and obtains its temperature
							start_pos = 0;
							end_pos = 100;
							temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
							HydroNet_InsertVolume(/*new_pos{ branch_aux1 }, new_temp{ branch_aux1 }, start_pos, end_pos, temp_action, fluid_type*/);

							// Reads the temperature in the middle of the branch
							obj_pos = (start_pos + end_pos) / 2;
							temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);

							// Enters temperature and mass in the enthalpy balance
							sum_mass = sum_mass + density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);
							sum_temp_mass = sum_temp_mass + temp_aux * density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);

							// Stores the excess volume
							overflow.at(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
						}
					}
				}
			}
			for (int count = 0; count < overflow.size(); count++)
				overflow.at(count) = 0;// leaves everything as found

			// Stores resulting inlet temperature
			heat_exch_Tout.at(count_htx).T_in = sum_temp_mass / sum_mass;
		}
	}


	// Heat exhangers of type fixed heat
	for (int count_htx, branch_aux, obj_aux = 0; count_htx < heat_exch_fix.size(); count_htx++) {


		obj_aux = heat_exch_fix.at(count_htx).heat_exch_fix; // index of the heat exchanger in objects table
		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}// branch that contains the heat exchanger

		if (flows.at(branch_aux) == 0) { // no movement in the branch

			// Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux);
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux }, start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the heat exchanger
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_fix.at(count_htx).T_in = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux }, start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_fix.at(count_htx).T_in = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
		}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			// First, obtains temperature of the volume between the start of the
			// branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			HydroNet_InsertVolume(/*new_pos{ branch_aux }, new_temp{ branch_aux },start_pos, end_pos, temp_action, fluid_type*/);

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);
			heat_exch_fix.at(count_htx).T_in = temp_aux;

			// Initializes enthalpy balance
			sum_mass = density(/*fluid_type, temp_aux*/) * (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

			// Stores volume that will come from previous branches
			overflow.at(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (vec_abs_sum(overflow) != 0) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow.at(count) > 0)
						count_branch = count;

					node_aux = branches_id.at(count_branch).at(1); // node at the start of the branch
					flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count)); // sum of flows going into the branch

					for (int aux = 0; aux < node_outlet_branches.at(node_count).size(); aux++)
						volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
					// volume that comes from every inlet branch (proportional to flow)
					for (int count_branch1 = 0, branch_aux1; count_branch1 < node_inlet_branches.at(node_aux).size(); count_branch1++) { // loop of inlet branches
						branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
						if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
							// volume in the branch is bigger than outcoming volume
							// also avoids dividing by zero

							// Inserts a volume for the affected area and obtains its temperature
							start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
							end_pos = 100;
							temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
							HydroNet_InsertVolume(/*new_pos{ branch_aux1 }, new_temp{ branch_aux1 }, start_pos, end_pos, temp_action, fluid_type*/);

							// Reads the temperature in the middle of the affected volume
							obj_pos = (start_pos + end_pos) / 2;
							temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);

							// Enters temperature and mass in the enthalpy balance
							sum_mass = sum_mass + density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);
							sum_temp_mass = sum_temp_mass + temp_aux * density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);
						}

						else {// volume in the branch is smaller than outcoming volume

							// Inserts a volume for the entire branch and obtains its temperature
							start_pos = 0;
							end_pos = 100;
							temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
							HydroNet_InsertVolume(/*new_pos{ branch_aux1 }, new_temp{ branch_aux1 }, start_pos, end_pos, temp_action, fluid_type*/);

							// Reads the temperature in the middle of the branch
							obj_pos = (start_pos + end_pos) / 2;
							temp_aux = HydroNet_GetObjTemperature(/*obj_pos, branch_temp_pos_aux, branch_temp_aux*/);

							// Enters temperature and mass in the enthalpy balance
							sum_mass = sum_mass + density(/*fluid_type, temp_aux) * volume_share(count_branch1*/);
							sum_temp_mass = sum_temp_mass + temp_aux * density(/*fluid_type, temp_aux*/) * volume_share.at(count_branch1);

							// Stores the excess volume
							overflow.at(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
						}
					}
				}
			}

			for (int count = 0; count < overflow.size(); count++)
				overflow.at(count) = 0;// leaves everything as found

			// Stores resulting inlet temperature
			heat_exch_fix.at(count_htx).T_in = sum_temp_mass / sum_mass;

		}
	}
}


int main()
{
	std::vector<strBranches> branches;
	std::vector<strPump_volum> pump_volum;
	std::vector<std::vector<int>>branches_id;
	std::vector<double> bound_flows;
	std::vector<int> nodes_id;
	std::vector<std::vector<int>> mesh_branches;
	std::vector<std::vector<int>> node_branches;
	int n_mesh;
	int n_nodes;
	int n_branch;
	std::vector<double> head_loss;
	std::vector<double> hydr_resist1;
	std::vector<double> hydr_resist2;
	int n_tanks;
	std::string fluid_type;
	std::vector<double> flows;
	std::vector<double> branch_volume;
	std::vector<std::vector<int>>branch_htx_Tout;
	std::vector<std::vector<int>>branch_htx_fix;
	double dt;
	std::vector<double> obj_inlet_pos;
	std::vector<double> obj_outlet_pos;
	std::vector<strHeat_exch_Tout> heat_exch_Tout;
	std::vector<strHeat_exch_fix>heat_exch_fix;
	std::vector<std::vector<double>> branch_temp_pos;
	std::vector<std::vector<double>> branch_temperature;

	fluid_type = "water";

	branches_id.resize(6);

	branches_id.at(0) = { 12,13,1,2,3 };
	branches_id.at(1) = { 12,11,10,8,7 };
	branches_id.at(2) = { 7,6 };
	branches_id.at(3) = { 6,4,3 };
	branches_id.at(4) = { 6,14,5,3 };
	branches_id.at(5) = { 12,5,9,7 };

	bound_flows.assign(6, NAN);
	bound_flows.at(1) = 0;

	nodes_id = { 3,6,7,12 };
	n_nodes = 4;

	mesh_branches.resize(6);

	mesh_branches.at(0) = { 0,1,2,3 };
	mesh_branches.at(1) = { 0,1,2,4 };
	mesh_branches.at(2) = { 1,5 };
	mesh_branches.at(3) = { 0,5,2,3 };
	mesh_branches.at(4) = { 0,5,2,4 };
	mesh_branches.at(5) = { 3,4 };

	node_branches.resize(4);

	node_branches.at(0) = { 0,3,4 };
	node_branches.at(1) = { 2,3,4 };
	node_branches.at(2) = { 1,2,5 };
	node_branches.at(3) = { 0,1,5 };

	n_mesh = 6;

	n_branch = 6;

	head_loss = { 0,0,0,0,0,1 };

	hydr_resist1 = { 0,0,0,0,0,0 };

	hydr_resist2 = { 2.120926015593668E06,1.016417857504270E06,0,50000,100000,5.1642E5 };

	n_tanks = 0;

	pump_volum.resize(6);

	pump_volum.at(0).Pump_volum = 1;
	pump_volum.at(0).coef0 = 1;
	pump_volum.at(0).coef1 = 0;
	pump_volum.at(0).coef2 = 0;
	pump_volum.at(0).coef3 = 0;
	pump_volum.at(0).head_max = 1000000;
	pump_volum.at(0).inlet_object = 1;
	pump_volum.at(0).outlet_oject = 2;
	pump_volum.at(0).pump_speed = 105;

	branches.resize(6);


	branches.at(0).branch_ind = 1;
	branches.at(1).branch_ind = 2;
	branches.at(2).branch_ind = 3;
	branches.at(4).branch_ind = 4;
	branches.at(4).branch_ind = 5;
	branches.at(5).branch_ind = 6;


	branches.at(0).objects.resize(3);
	branches.at(1).objects.resize(3);
	branches.at(2).objects.resize(1);
	branches.at(3).objects.resize(1);
	branches.at(4).objects.resize(2);
	branches.at(5).objects.resize(2);

	branches.at(0).objects.at(0).Class = "pump_volum";
	branches.at(0).objects.at(0).ID = 1;
	branches.at(0).objects.at(0).name = "Main_Pump";
	branches.at(0).objects.at(0).volume = 1E-3;
	branches.at(0).objects.at(0).inlet_object = 13;
	branches.at(0).objects.at(0).outlet_oject = 2;

	branches.at(0).objects.at(1).Class = "pipe";
	branches.at(0).objects.at(1).ID = 2;
	branches.at(0).objects.at(1).name = "Enigine_Inlet";
	branches.at(0).objects.at(1).volume = 8.8357E-5;
	branches.at(0).objects.at(1).inlet_object = 1;
	branches.at(0).objects.at(1).outlet_oject = 3;

	branches.at(0).objects.at(2).Class = "pipe";
	branches.at(0).objects.at(2).ID = 13;
	branches.at(0).objects.at(2).name = "Return_Pipe";
	branches.at(0).objects.at(2).volume = 6.2832E-4;
	branches.at(0).objects.at(2).inlet_object = 12;
	branches.at(0).objects.at(2).outlet_oject = 1;

	branches.at(1).objects.at(0).Class = "thermostat";
	branches.at(1).objects.at(0).ID = 8;
	branches.at(1).objects.at(0).name = "Thermostat_rad";
	branches.at(1).objects.at(0).volume = 0;
	branches.at(1).objects.at(0).inlet_object = 7;
	branches.at(1).objects.at(0).outlet_oject = 10;

	branches.at(1).objects.at(1).Class = "heat_exch_fix";
	branches.at(1).objects.at(1).ID = 10;
	branches.at(1).objects.at(1).name = "Radiator";
	branches.at(1).objects.at(1).volume = 5E-4;
	branches.at(1).objects.at(1).inlet_object = 8;
	branches.at(1).objects.at(1).outlet_oject = 11;

	branches.at(1).objects.at(2).Class = "pipe";
	branches.at(1).objects.at(2).ID = 11;
	branches.at(1).objects.at(2).name = "Radiator_outlet";
	branches.at(1).objects.at(2).volume = 3.1416E-4;
	branches.at(1).objects.at(2).inlet_object = 10;
	branches.at(1).objects.at(2).outlet_oject = 12;

	branches.at(2).objects.at(0).Class = "node";




	branches.at(3).objects.at(0).Class = "heat_exch_fix";
	branches.at(3).objects.at(0).ID = 4;
	branches.at(3).objects.at(0).name = "Block_cyl1";
	branches.at(3).objects.at(0).volume = 5E-4;
	branches.at(3).objects.at(0).inlet_object = 3;
	branches.at(3).objects.at(0).outlet_oject = 6;

	branches.at(4).objects.at(0).Class = "heat_exch_fix";
	branches.at(4).objects.at(0).ID = 5;
	branches.at(4).objects.at(0).name = "Head_cyl1";
	branches.at(4).objects.at(0).volume = 5E-4;
	branches.at(4).objects.at(0).inlet_object = 3;
	branches.at(4).objects.at(0).outlet_oject = 14;

	branches.at(4).objects.at(1).Class = "heat_exch_fix";
	branches.at(4).objects.at(1).ID = 14;
	branches.at(4).objects.at(1).name = "Head_cyl2";
	branches.at(4).objects.at(1).volume = 5E-4;
	branches.at(4).objects.at(1).inlet_object = 5;
	branches.at(4).objects.at(1).outlet_oject = 6;

	branches.at(5).objects.at(0).Class = "thermostat";
	branches.at(5).objects.at(0).ID = 9;
	branches.at(5).objects.at(0).name = "Thermostat_byp";
	branches.at(5).objects.at(0).volume = 0;
	branches.at(5).objects.at(0).inlet_object = 7;
	branches.at(5).objects.at(0).outlet_oject = 15;

	branches.at(5).objects.at(1).Class = "pipe";
	branches.at(5).objects.at(1).ID = 5;
	branches.at(5).objects.at(1).name = "Bypass_Pipe";
	branches.at(5).objects.at(1).volume = 3.1416E-4;
	branches.at(5).objects.at(1).inlet_object = 9;
	branches.at(5).objects.at(1).outlet_oject = 12;

	HydroNet_Temperature(fluid_type, branches, obj_inlet_pos, obj_outlet_pos, heat_exch_Tout, heat_exch_fix, nodes_id, branches_id, node_branches, n_nodes, n_branch, branch_volume, branch_temp_pos, branch_temperature, branch_htx_Tout, branch_htx_fix, flows, dt);

	return 0;
}