// HydroNet_Create.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <string>

void HydroNet_ReadObj(std::string circuit_file);

void HydroNet_Create(std:: string circuit_file) {
	// Creates a model for the specified circuit
	int n_obj;

	//// ANALYSIS OF CIRCUIT OBJECTS
	// Reads circuit configuration and specifications of objects and stores data in tables

	HydroNet_ReadObj(circuit_file);

	n_obj = objects.size();



	//// NODES
	// List of all nodes

	nodes_ind = find(strcmp(objects.class, 'node'));
	nodes_id = objects.id(nodes_ind);
	n_nodes = size(nodes_ind, 1);

	// Using a loop
	// nodes_ind = zeros(n_obj, 1); // stores position in object list
	// nodes_id = zeros(n_obj, 1); // stores id
	//     // maximum size: all objects are nodes
	//
	// node_count = 0;
	// for count = 1 : n_obj // objects loop
	//     if strcmp(objects.class{count}, 'node')
	//         node_count = node_count + 1;
	//         nodes_ind(node_count) = count;
	//         nodes_id(node_count) = objects.id(count);
	//     end
	// end
	// nodes_ind = nodes_ind(1 : node_count); // adjust size
	// nodes_id = nodes_id(1 : node_count);



	//// TANKS
	// List of all tanks

	tanks_ind = zeros(n_obj, 1); // stores tank ids
		// maximum size: all objects are tanks

	tank_count = 0;
	for count = 1 : n_obj // objects loop
		if strcmp(objects.class { count }, 'tank') {
			tank_count = tank_count + 1;
			tanks_ind(tank_count) = count;
		}

	tanks_ind = tanks_ind(1 : tank_count); // adjust size

	n_tanks = length(tanks_ind);



	//// ADD OBJECT'S POSITION INTO CERTAIN TABLES
	// Index of the object in the object's table

	valve_var.obj_index = find(strcmp(objects.class, 'valve_var'));
	thermostat.obj_index = find(strcmp(objects.class, 'thermostat'));
	pump_volum.obj_index = find(strcmp(objects.class, 'pump_volum'));
	pump_turbo.obj_index = find(strcmp(objects.class, 'pump_turbo'));
	heat_exch_fix.obj_index = find(strcmp(objects.class, 'heat_exch_fix'));
	heat_exch_Tout.obj_index = find(strcmp(objects.class, 'heat_exch_Tout'));
	tank.obj_index = find(strcmp(objects.class, 'tank'));



	//// CIRCUIT COMPONENTS
	// Obtains nodes, branches and meshes in the circuit.
	// Determines branches connected to nodes and branches that are part of a mesh.
	// All valves and thermostats are assumed open.

	HydroNet_Components ( objects, n_obj, nodes_ind, pump_volum, pump_turbo, heat_exch_Tout, tank );



	//// BRANCH VOLUME

	branch_volume.assign(n_branch, 0);

	for (int count = 0; count < n_branch; count++) {// branch loop
		sum_aux = 0;
		for count1 = branches_ind{ count } // loop of objects in the branch
			sum_aux = sum_aux + objects.volume(count1);

		branch_volume(count) = sum_aux;
	}



	//// INLET AND OUTLET OF OBJECTS
	// Stores volumetric position of the inlet and outlet of objects, as a \// of the total branch volume

	// Inlet positions
	obj_inlet_pos_cell = cell(n_branch, 1); // branch / objects format
	obj_inlet_pos = cell(n_obj, 1); // "vector" of objects format
	for (int count_branch = 0; count_branch < n_branch; count_branch++) {// branch loop
		inlet_aux.assign(branches_ind.at(count_branch).size(), 0); // auxiliary vector to store inlets in branch
		for (int count = 0; count < branches_ind.at(count_branch).size(); count++) // initialization loop
			count_obj = branches_ind.at(count_branch).at(count);
		obj_inlet_pos.at( count_obj ) = [obj_inlet_pos.at( count_obj ) 0]; // all objects of the branch have the inlet at 0//
		
		if (branch_volume(count_branch) != 0) { // to avoid dividing by zero if branch volume is 0
			for (int count_obj = 1; count_obj < branches_ind.at(count_branch).size(); count_obj++) // loop of objects in the branch
				inlet_aux(count_obj) = inlet_aux(count_obj - 1) + objects.volume(branches_ind.at(count_branch).at(count_obj - 1)) / branch_volume.at(count_branch) * 100;
			// \// of the previous object inlet + volume previous object / branch volume * 100
			obj_inlet_pos.at(branches_ind.at(count_branch).at(count_obj)) = [obj_inlet_pos.at( branches_ind.at(count_branch).at(count_obj) )(1:(end - 1)) inlet_aux(count_obj)];
		}
			
		obj_inlet_pos_cell.at( count_branch ) = inlet_aux; // stores inlets of the current branch
	}


	// Outlet positions
	obj_outlet_pos = cell(n_obj, 1);
	for (int count = 0; count < n_obj;count++) {
		if (objects.volume.at(count) == 0)
			obj_outlet_pos.at( count ) = obj_inlet_pos.at( count);
		else // volume > 0
			obj_outlet_pos.at( count ) = obj_inlet_pos.at( count ) +objects.volume.at(count) / branch_volume.at(objects.branch.at( count )) * 100;
	}


	//// INDEXES OF HEAT EXCHANGERS IN EVERY BRANCH

	branch_htx_fix = cell(n_branch, 1);
	branch_htx_Tout = cell(n_branch, 1);

	for (int count_branch = 0; count_branch < n_branch; count_branch++) {
		branch_htx_fix.at(count_branch) = transpose(find(cell2mat(objects.branch(heat_exch_fix.obj_index)) == count_branch));
		// index in the object's table of heat exchangers of type fix heat that are present in the branch
		// Note: a heat exchanger cannot be in two branches but a branch can contain several heat exchangers
		if (isempty(branch_htx_fix.at(count_branch)))
			branch_htx_fix.at(count_branch) = {}; // to avoid problems with [0x1 double]

		branch_htx_Tout.at( count_branch ) = transpose(find(cell2mat(objects.branch(heat_exch_Tout.obj_index)) == count_branch));
		// index in the object's table of heat exchangers of type T_out that are present in the branch
		if (isempty(branch_htx_Tout(count_branch)))
			branch_htx_Tout.at(count_branch) = {}; // to avoid problems with [0x1 double]

		// Sort according to closeness to the end of the branch
		[~, ind_aux] = sort(obj_inlet_pos_cell.at(count_branch).at(branch_htx_fix.at(count_branch)), 'descend');
		// sorted indexes of heat exchangers in obj_inlet_pos
		branch_htx_fix.at(count_branch) = branch_htx_fix.at(count_branch).at(ind_aux);
		// reorder heat exchangers
		// Repeat for branch_htx_Tout
		[~, ind_aux] = sort(obj_inlet_pos_cell.at(count_branch).at(branch_htx_Tout.at(count_branch)), 'descend');
		branch_htx_Tout.at(count_branch) = branch_htx_Tout.at(count_branch).at(ind_aux);
	}





}

int main()
{
    return 0;
}

