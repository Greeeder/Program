// HydroNet_Components.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <string>
//#include <Eigen/Dense>

struct strObject {
	int ID;							// Object's idintificator
	std::string name;				// Object's name
	std::string Class;				// Object's clas
	double volume;					// Objects volume
	std::vector<int> adjacent;		// Adjacent objects
};

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
	int outlet_object;					// Thermostat's outlet object identifier
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
	int outlet_object;				// Pump's outlet object identifier
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
	int outlet_object;				// Pump's outlet object identifier
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
int outlet_object;				// Tank's outlet object identifier
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
	std::vector<int> adjacent;		// Adjacent objects
};
struct strBranches
{
	int branch_ind;				// branch idintificator
	int branches_cycle;
	//int n_branch;				// Number of branches
	std::vector<strObjects> objects;	// Objects in Brnach
};
//std::vector<std::vector<bool>> mirror (std::vector<std::vector<bool>> A){
//	int n = A.size();
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < n; j++)
//			A.at(j).at(i) = A.at(i).at(j);
//	return A;
//
//}

bool all(std::vector<bool> A) {
	bool all = true;
	for (int i = 0; i < A.size(); i++)
		if (A.at(i) == false) {
			all = false;
			break;
		}
	return all;
}

bool any(std::vector<bool> A) {
	bool any = false;
	for (int i = 0; i < A.size(); i++)
		if (A.at(i) == true) {
			any = true;
			break;
		}
	return any;
}
bool in(int a, std::vector<int> A) {
	bool in = false;

	for (int i = 0; i < A.size(); i++)
		if (a == A.at(i)) {
			in = true;
			break;
		}
	return in;
}

std::vector<int> flip(std::vector<int> A) {
	std::vector<int> aux ;

	for (int i = 0; i < A.size(); i++)
		aux.push_back(A.at(A.size()-1 - i));

	return aux;


}

void  HydroNet_Components(/*std::vector<strBranches> branches,*/ std:: vector <strObject> objects,int n_obj, std::vector<int> nodes_id, std::vector<strPump_volum> pump_volum, std::vector<strPump_turbo> pump_turbo, std::vector<strHeat_exch_Tout> heat_exch_Tout/*, std::vector<strTank> tank*/) {

	int edge_count, start_obj, path_count, obj_count, row_aux, column, prev_row, n_paths, branch_count, n_mesh, n_new_paths;
	std::vector<std::vector<bool>> edges, mat_aux;
	std::vector<int> path, branch_aux, vec_aux1,vec_vec_aux,path_aux;
	std::vector<std::vector<int>>paths_ind, mesh_ind, branches_ind, branch_mesh, node_branches, mesh_branches, objects_branches, branches_id, ind_outlet_aux;
	bool flag_path,first;
	std::vector<bool>path_cycle, to_remove, vec_aux, branch_cycle;
	//Eigen::RowVectorXi vec_aux;
	std::string node = "node", tank = "tank";
	std::vector<strBranches> branches;

	/*std::vector<int> tanks_id = {};
	for (int count=0;count<tank.size();count++)
		tanks_id.push_back(tank.at(count).tank)*/


	// Obtains all information about a network and its components - meshes,
	// nodes and branches.

	// Stores all nodes in the network(index and id)
	// Finds all paths using a Depth - first search algorithm.
	// Checks whether paths are closed cycles(mesh).
	// Stores all objects in each mesh(index and id)
	// Finds all branches(paths without bifurcations crossed by a single flow)
	// Stores all objects in each branch(index and id)
	// Checks whether branches form part of closed cycles(meshes).
	// Finds branches that form part of each mesh
	// Finds branches connected to each node

	// Inputs are the table of objects, tables for classes valve and thermostat
	// (for checking opening) and temperature at the inlet of thermostats.

	// Outputs are :
	//   Nodes in the network(object indices)
	// Nodes in the network(object id)
	// Objects in every mesh(object indices)
	// Objects in every mesh(object id)
	// Objects in every branch(object indices)
	// Objects in every branch(object id)
	// Logical vector that marks whether branches form part of closed cycles
	//   Branches that form part of each mesh
	// 	Branches connected to each node
	//   Number of meshes
	//   Number of branches


	// Initialization for testing
	// circuit_file = 'Coolant_basic_nodes.txt';
	// //circuit_file = 'Test_nodes.txt';
	//[fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, ...
	//     heat_exch_fix, heat_exch_Tout, tank] = HydroNet_ReadObj(circuit_file);
	// temperature = 0;



	//// EDGES
	// Looks for all connections among objects(indirected).

	edge_count = 0; // counter of edges
	edges.resize(n_obj);
	for (int count = 0; count < n_obj; count++)
		edges.at(count).assign(n_obj, false); // matrix that contains false if there is no connection between

	// objects and true if they are connected(rows and columns are the indices of objects)

	//for (int count = 0; count < n_obj; count++) {// objects loop(no need to check last object)
	//// counter is equal to index
	//	for (int brch = 0; brch < branches.size(); brch++)
	//		for (int count1 = 0, id_aux, ind_aux; branches.at(brch).objects.at(count).adjacent.size(); count1++) { // connections loop
	//			id_aux = branches.at(brch).objects.at(count).adjacent.at(count1); // id of adjacent object
	//			for (int count2 = 0; count2 < branches.at(brch).objects.size(); count2++)
	//				if (branches.at(brch).objects.at(count2).ID == id_aux)
	//					ind_aux = branches.at(brch).objects.at(count2).ID;
	//			if (ind_aux > count) {// the adjacent object has not been visited
	//			// the edge has not been stored
	//				edge_count = edge_count + 1; // next edge
	//				edges.at(count).at(ind_aux) = true; // stores edge
	//			}
	//		}
	//} // edges is a triangular matrix

	//edges = mirror(edges); // full symmetrical matrix(false OR true = true)

	
	for (int obj = 0; obj < n_obj; obj++)
		for (int adj = 0; adj <objects.at(obj).adjacent.size(); adj++) 
			edges.at(objects.at(obj).ID).at(objects.at(obj).adjacent.at(adj)) = true;
				
	for (int row = 0; row < edges.size(); row++)
		for (int col = 0; col < edges.at(row).size(); col++)
			if (edges.at(row).at(col)==true)
				edge_count++;

	edge_count = edge_count / 2;
			

	// Alternative: Store list in a cell
	// edges = cell((obj_num * (obj_num - 1)) / 2, 1); // stores all edges
	// maximum size : all objects are nodes and are connected with each other
	// obj_id = zeros(obj_num, 1);
	// obj_count = 0; // counter of visited objects
	// for count = 1 : obj_num - 1 // objects loop(no need to check last object)
	// id_aux = objects.id(count); // id of current object
	//     obj_count = obj_count + 1; // next visited object
	//     obj_id(obj_count) = id_aux; // store object as visited
	//     for count1 = 1 : numel(objects.adjacent{ count }) // connections loop
	//         id_aux1 = objects.adjacent{ count }(count1); // id of adjacent object
	//         if not(any((id_aux1 == obj_id)) // the adjacent object is not in the list of visited objects
		//             // the edge has not been stored
		//             edge_count = edge_count + 1; // next edge
		//             edges{ edge_count } = [id_aux id_aux1]; // stores edge
		//         end
		//     end
		// end
		// edges = edges(1 : edge_count); // adjust size



		//// PATHS

	// Selects starting object : first node of the list
	if (nodes_id.size() == 0) { // no nodes found
		start_obj = 1; // index of first object in the list
		//paths_ind.push_back({}); // stores all paths(indices) // only one path
	}
	else {
		start_obj = nodes_id.at(0); // index of first node
		//paths_ind.resize(edge_count); // stores all paths(indices)
			// maximum size : same as number of edges
	}


	// Finds all paths using a Depth - first search algorithm.

	path_count =0; // counter of paths
	mat_aux = edges; // auxiliary matrix to mark visited objects
	//path.assign(n_obj + 1, 0); // stores indices of objects in a path
	// +1 because in a closed cycle the starting object is also the end object
	//path.at(0) = start_obj; // initialize path
	path = { start_obj };
	obj_count = 1; // counter of objects along the path
	row_aux = start_obj; // initialize row(starting object)
	column = -1;


	/*while ( ! ((row_aux == start_obj) && (column ==(n_obj-1)) ) ) {
		// while there are unexplored paths from the starting object(did not reach the end of the starting row)
		
		column = column + 1; // next column

		if ( mat_aux.at(row_aux).at(column) == true) { // there is a connection

			mat_aux.at(row_aux).at(column) = false; // mark edge as explored
			prev_row = row_aux;
			row_aux = column; // move to connection's row
			flag_path = false; // path not finished
			obj_count = obj_count + 1; // next object along the path
			path.at(obj_count-1)= row_aux; // store index of object in the path
			mat_aux.at(row_aux).at(prev_row) = false; // going back the same way is not allowed

			//for (int count = 0; count < mat_aux.at(row_aux).size(); count++) {
				if (!any(mat_aux.at(row_aux)) ) // there are no more ways to go
					flag_path = true; // end path
				else
					for (int count1 = 0;  count1 < obj_count - 2 ; count1++) {
						if (path.at(count1) == row_aux) {// the new object was already in the path - the cycle is closed
							flag_path = true;
							break;// end path
						}
							//         elseif strcmp(objects.class(row_aux), 'valve_var') // is a valve and may be closed
							//             if  valve_var.opening(valve_var.id == objects.id(row_aux)) == 0 // it is closed
							//                 flag_path = true; // end path
							//             end
							//         elseif strcmp(objects.class(row_aux), 'thermostat') // is a thermostat and may be closed
							//             pos_aux = find(thermostat.id == objects.id(row_aux)); // position in thermostat table
							//             if thermostat.T_coef0(pos_aux) + thermostat.T_coef1(pos_aux) * temperature ...
							//                     +thermostat.T_coef2(pos_aux) * temperature ^ 2 == 0 // it is closed
							//                 flag_path = true; // end path
							//             end
					}
			//}

			if (flag_path) { // pat_h has ended
				path.resize(obj_count);
				path.shrink_to_fit();// adjust size 
				path_count = path_count + 1; // next path
				paths_ind.at(path_count-1)= path; // save path
				path.pop_back(); // remove last object from path - step back
				path.resize(n_obj+1);
				obj_count = obj_count - 1; // previous object along the path
				if (row_aux == start_obj)// the starting row never gets fully re - initialized
					mat_aux.at(row_aux).at(prev_row) = true; // just allow the same way back for future paths
				else
					mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object

				row_aux = path.at(obj_count); // return to upper object
			}
			else // path goes further
				column = -1; // stay in the new row and start from the beginning

		}
		while ((column >= (n_obj-1)) && (row_aux != start_obj)) {// loop to step back in a path
		// end of the row and it is not the starting object // skip fully explored objects
			mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object
			column = path.at(obj_count - 1); // continue the upper object where it was left
			path.pop_back(); // remove last object from path - step back
			path.resize(n_obj);
			obj_count = obj_count - 1; // previous object along the path
			row_aux = path.at(obj_count+1); // return to upper object
			if (path.size() > 1)
				mat_aux.at(row_aux).at(path.at(obj_count)) = false; // the previous object along the path is not allowed
		}
	}*/
	while (!(row_aux == start_obj && column == n_obj)) {
		// while there are unexplored paths from the starting object(did not reach the end of the starting row)



		column = column + 1; // next column

		if (column<n_obj && mat_aux.at(row_aux).at(column)) { // there is a connection

			mat_aux.at(row_aux).at(column) = false; // mark edge as explored
			prev_row = row_aux;
			row_aux = column; // move to connection's row
			flag_path = false; // path not finished
			obj_count = obj_count + 1; // next object along the path
			path.push_back(row_aux); // store index of object in the path
			mat_aux.at(row_aux).at(prev_row) = false; // going back the same way is not allowed


			path_aux = path;
			path_aux.pop_back();


			if (!any(mat_aux.at(row_aux))) // there are no more ways to go
				flag_path = true; // end path



			else if (in(row_aux, path_aux)) // the new object was already in the path - the cycle is closed
				flag_path = true; // end path
								  //         elseif strcmp(objects.class(row_aux), 'valve_var') // is a valve and may be closed
								  //             if  valve_var.opening(valve_var.id == objects.id(row_aux)) == 0 // it is closed
								  //                 flag_path = true; // end path
								  //             end
								  //         elseif strcmp(objects.class(row_aux), 'thermostat') // is a thermostat and may be closed
								  //             pos_aux = find(thermostat.id == objects.id(row_aux)); // position in thermostat table
								  //             if thermostat.T_coef0(pos_aux) + thermostat.T_coef1(pos_aux) * temperature ...
								  //                     +thermostat.T_coef2(pos_aux) * temperature ^ 2 == 0 // it is closed
								  //                 flag_path = true; // end path
								  //             end


			if (flag_path) {// path has ended
				path_count = path_count + 1; // next path
				paths_ind.push_back(path); // save path
				path.pop_back(); // remove last object from path - step back
				obj_count = obj_count - 1; // previous object along the path
				if (row_aux == start_obj) // the starting row never gets fully re - initialized
					mat_aux.at(row_aux).at(prev_row) = true; // just allow the same way back for future paths
				else
					mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object

				row_aux = path.at(path.size() - 1); // return to upper object

				column = -1;
			}
			else // path goes further
				column = -1; // stay in the new row and start from the beginning
		}

		while ((column >= n_obj - 1) && (row_aux != start_obj)) { // loop to step back in a path
																  // end of the row and it is not the starting object // skip fully explored objects
			mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object
			column = path.at(path.size() - 1); // continue the upper object where it was left
			path.pop_back(); // remove last object from path - step back
			obj_count = obj_count - 1; // previous object along the path
			row_aux = path.at(path.size() - 1); // return to upper object
			if (path.size() > 1)
				mat_aux.at(row_aux).at(path.at(path.size() - 2)) = false; // the previous object along the path is not allowed

		}
	}

	n_paths = path_count; // number of paths
	paths_ind.resize(n_paths);
	path.shrink_to_fit();// adjust size
	path_aux.clear();
	path_aux.shrink_to_fit();
	path_aux = {};



	// Find and isolate closed paths(cycles)
	path_cycle.assign(n_paths, false);
	// vector that indicates whether the path forms a mesh(closed cycle)
	for (int count_path = 0; count_path < n_paths; count_path++) {// path loop
	// Criterium to know if the paths are cycles : the last object and one of the other objects in the path are the same
	// Two cases : pure cycle(first object and last object are the same)
	// branched cycle(first objects are not a part of the cycle)
		if (paths_ind.at(count_path).at(0) == paths_ind.at(count_path).at(paths_ind.at(count_path).size()-1)) // pure cycle
			path_cycle.at(count_path) = true;
		else {//any(paths_ind{ count }(2 : end - 2) == paths_ind{ count }(end)) // branched cycle
			for (int count_obj = 1; count_obj < paths_ind.at(count_path).size() - 3; count_obj++) {// minus first, last and next - to - last
				if (paths_ind.at(count_path).at(count_obj) == paths_ind.at(count_path).at(paths_ind.at(count_path).size()-1)) { // cycle starts here
					path_cycle.at(count_path) = true;
					// Create new path for the branch that is out of the cycle
					first = true;
					for (int aux_count_obj = 0; aux_count_obj < count_obj+1; aux_count_obj++) {
						if (first) {
							first = false;
							paths_ind.push_back({});
						}

						paths_ind.at(paths_ind.size() - 1 ).push_back(paths_ind.at(count_path).at(aux_count_obj));
					}
					// Convert branched cycle in a pure cycle
					for (int aux_count_obj = count_obj; aux_count_obj < paths_ind.at(count_path).size(); aux_count_obj++)
						path_aux.push_back(paths_ind.at(count_path).at(aux_count_obj));
					paths_ind.at(count_path) = path_aux;
					path_aux.clear();
					path_aux.shrink_to_fit();
					path_aux = {};
					break; // path processed; go to next
				}
			}
		}
	}
	n_paths = paths_ind.size(); // refresh number of paths
	n_new_paths = n_paths - path_cycle.size();
	for (int count = 0; count <n_new_paths; count++)
		path_cycle.push_back(false);


	// Removes repeated paths and cycles(all cycles are duplicated once - two directions)
	to_remove.assign(n_paths, false); // stores paths to be removed
	for (int count = 0; count < n_paths; count++) { // path loop
		if (to_remove.at(count) == false) { // the path is not marked to be deleted
			for (int count1 = count + 1; count1 < n_paths; count1++) {// compare with next paths
				if (paths_ind.at(count1).size() == paths_ind.at(count).size()) {
					// to compare for equality, vectors must have the same length
					//                 if or (paths_ind{ count1 } == paths_ind{ count }, paths_ind{ count1 } == flip(paths_ind{ count }))
					// // paths are equal OR path is equal to the inverse path
					//                     to_remove(count1) = true; // mark to be removed
					//                 end
					// all objects of one branch are present in the other one
					// (even if start / end node is not the sa)me)
					vec_aux.assign(paths_ind.at(count).size(), false);
					for (int count_obj = 0; count_obj < paths_ind.at(count).size(); count_obj++) {
						if (in(paths_ind.at(count).at(count_obj),paths_ind.at(count1))) {
							vec_aux.at(count_obj) = true;
						}
					}
					if (all(vec_aux))// all objects are present
						to_remove.at(count1) = true; // mark to be removed
				}
			}
		}
	}
	
	for (int count = 0,aux = 0; count < to_remove.size(); count++)
		if (to_remove.at(count) == true) {
			
			paths_ind.erase(paths_ind.end() - (paths_ind.size()-count + aux) ); // remove repeated paths and cycles
			paths_ind.shrink_to_fit();
			path_cycle.erase(path_cycle.end() - (path_cycle.size() - count + aux) );
			path_cycle.shrink_to_fit();
			aux++;
		}
	n_paths = paths_ind.size(); // refresh number of paths


	// // Convert paths to change indices by id's
	// paths_id = cell(path_count, 1); // store all paths(id)
	// for count = 1 : path_count // path loop
	//     for count1 = 1 : numel(paths_ind{ count }) // objects in a path
	//         paths_id {count}(count1) = objects.id(paths_ind{ count }(count1)); // assign id
	//     end
	// end



	//// MESHES(cycles)
	for (int count = 0; count < paths_ind.size(); count++) {
		if (path_cycle.at(count) == true)
			mesh_ind.push_back(paths_ind.at(count)); // paths that are cycles(indices)
			// mesh_id = paths_id(path_cycle); // (id)
	}

	// Convert meshes to change indices by id's
	// mesh_id = cell(size(mesh_ind, 1), 1); // store all meshes(id)
	// for count = 1 : size(mesh_id, 1) // mesh loop
	//     for count1 = 1 : numel(mesh_ind{ count }) // objects in a mesh
	//         mesh_id {	count}(count1) = objects.id(mesh_ind{ count }(count1)); // assign id
	//     end
	// end



	//// BRANCHES
	// Paths without bifurcations crossed by a single flow.
	// Branches that are not included in a mesh will be obtained too because they may be
	// traversed by a flow if there is a boundary condition at one end(flow or head).

	branches_ind.resize(edge_count * n_paths); // stores unique branches of all paths(indices)
	// maximum size : number of edges * number of paths
	branch_cycle.assign(edge_count * n_paths, false); // stores whether the branch is inside a mesh
	branch_mesh.resize(edge_count * n_paths); // stores the mesh where the branch is contained(indices)
	// is a cell because it will store all meshes that contain that branch
	// if the branch is not inside a mesh, the path is not stored
	// it is used to make a cell later that stores branches in each mesh


	// Find branches in paths by splitting them at nodes and tanks
	branch_count = 0; // branch counter
	branch_aux = {};
	for (int count = 0, count1; count < n_paths; count++) { // path loop
		start_obj = 0; // starting index for the first branch
		for (count1 = 1; count1 < (paths_ind.at(count).size() -1); count1++) { // counter of objects in a path excluding first and last
			if (/*(objects.at(paths_ind.at(count).at(count1)).Class == node)--->*/ in(paths_ind.at(count).at(count1),nodes_id)/* || in(paths_ind.at(count).at(count1), tanks_id) /*<-----(objects.at(paths_ind.at(count).at(count1)).Class == tank)*/) {
				branch_count = branch_count + 1; // next branch

				for (int aux = start_obj; aux <= count1; aux++)
					branch_aux.push_back(paths_ind.at(count).at(aux)); // save branch

				branches_ind.at(branch_count-1) = branch_aux;
				branch_aux.clear();
				branch_aux.shrink_to_fit();
				branch_aux = {};
				start_obj = count1; // last object of the present branch is the first object of the next one
				branch_cycle.at(branch_count-1) = path_cycle.at(count); // save whether the branch is inside a mesh
				if (path_cycle.at(count) == true) // path is a mesh
					for (int aux = 0; aux < path_cycle.size(); aux++)
						if (path_cycle.at(aux) == true)
							if (aux == count) {
								branch_mesh.at(branch_count - 1) = { count }; // save origin mesh(index)
								break;
							}

			}
		}
		if (paths_ind.at(count).size() == 2) {// branches with only 2 objects
			count1 = -1; // needed because the loop above was avoided
		}
		count1 = count1 + 1; // next object // count1 = end of path
		branch_count = branch_count + 1; // next branch
		for (int aux = start_obj; aux < count1; aux++)
			branch_aux.push_back(paths_ind.at(count).at(aux)); // save branch
		branches_ind.at(branch_count-1) = branch_aux;
		branch_aux.clear();
		branch_aux.shrink_to_fit();
		branch_aux = {};
		branch_cycle.at(branch_count-1) = path_cycle.at(count); // save whether the branch is inside a mesh
		if (path_cycle.at(count) == true) {// path is a mesh
			for (int aux = 0; aux < path_cycle.size(); aux++)
				if (path_cycle.at(aux) == true)
					if (aux == count) {
						branch_mesh.at(branch_count - 1) = { count }; // save origin mesh(index)
						break;
					}
		}
	}
	branches_ind.resize(branch_count); // adjust size
	branch_cycle.resize(branch_count);
	branch_mesh.resize(branch_count);
	branches_ind.shrink_to_fit();
	branch_cycle.shrink_to_fit();
	branch_mesh.shrink_to_fit();

	// Removes repeated branches
	to_remove.clear();
	to_remove.shrink_to_fit();
	to_remove.assign(branch_count,false); // stores branches to be removed
	for (int count = 0; count < branch_count;count++) { // branch loop
		if (to_remove.at(count)== false) { // the branch is not marked to be deleted
			for (int count1 = count + 1; count1 < branch_count;count1++) {// compare with next branches
				if (branches_ind.at(count1).size() == branches_ind.at( count ).size()) {
					// to compare for equality, vectors must have the same length
					if (branches_ind.at( count1 ) == branches_ind.at(count ) || branches_ind.at( count1) ==flip(branches_ind.at(count ))) {
						// branches are equal OR branch is equal to the inverse branch
						to_remove.at(count1) = true; // mark to be removed
						if (branch_cycle.at(count)==false && branch_cycle.at(count1)==true) {
							// first branch is not in a cycle but second(and same) branch is
							branch_cycle.at(count) = true; // branch is inside a cycle
						}
						for (int aux=0; aux<branch_mesh.at(count1).size();aux++)
						branch_mesh.at(count ).push_back( branch_mesh.at(count1 ).at(aux)); // add mesh
					}
				}
			}
		}
	}
	for (int count = 0,aux=0; count < to_remove.size(); count++)
		if (to_remove.at(count) == true) {
			branches_ind.erase(branches_ind.end() - (branches_ind.size() - count + aux)); // remove repeated branches
			branches_ind.shrink_to_fit();
			branch_cycle.erase(branch_cycle.end() - (branch_cycle.size() - count + aux));
			branch_cycle.shrink_to_fit();
			branch_mesh.erase(branch_mesh.end() - (branch_mesh.size() - count + aux));
			branch_mesh.shrink_to_fit();
			aux++;
		}
		branch_count =branches_ind.size(); // refresh number of branches


	// Finds and merges split branches.
	// If there are no nodes, the starting object may be inside a mesh.In that case,
	// it is in the middle of a branch but the algorithms have splitted it in two parts.
		if (nodes_id.size() == 0) { // no nodes found(first condition)
			for (int count = 0; count < branch_count; count++) {// branch loop
				if (branches_ind.at(count).at(0) == 1 || branches_ind.at(count).at(branches_ind.at(count).size()-1) == 1) {
					// starting object is in the branch(its index is 1)
					if (branch_cycle.at(count) == true) { // is in a mesh(second condition)
						for (int count1 = count; count1 < branch_count; count1++) {// looks for the other side
							if (branches_ind.at(count1).at(0) == 1 || branches_ind.at(count1).at(branches_ind.at(count1).size()-1) == 1) {
								// Join branches
								if (branches_ind.at(count).at(branches_ind.at(count).size()-1) == 1 && branches_ind.at(count).at(0) == 1)
									// starting object is the last of first branch and the first of the second brandlch
									for (int aux = 1; aux < branches_ind.at(count1).size(); aux++)
										branches_ind.at(count).push_back(branches_ind.at(count1).at(aux)); // join directly // do not repeat starting object
								else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size()-1) == 1)
									// starting object is the last of first branch and the last of the second branch
									for (int aux = 1; aux < branches_ind.at(count1).size(); aux++)
										branches_ind.at(count).push_back(branches_ind.at(count1).at(branches_ind.at(count1).size() - aux)); // flip second branch
								else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size()-1) == 1) {
									// starting object is the first of first branch and the first of the second branch
									for (int aux = 0; aux < branches_ind.at(count1).size(); aux++) {
										branch_aux.push_back(branches_ind.at(count).at(branches_ind.at(count).size()-1 - aux));
									}
									for (int aux = 0; aux < branches_ind.at(count1).size(); aux++) {
										branch_aux.push_back(branches_ind.at(count1).at(aux));
									}
									branches_ind.at(count) = branch_aux; // flip first branch

								}

								else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size()-1) == 1)
									// starting object is the first of first branch and the last of the second branch
									branch_aux = branches_ind.at(count1);
								for (int aux = 1; aux < branches_ind.at(count1).size(); aux++) {
									branch_aux.push_back(branches_ind.at(count).at(aux));
									branches_ind.at(count) = branch_aux; // join second branch and then first

									branches_ind.erase(branches_ind.begin() + count1); // remove second branch
									branch_cycle.erase(branch_cycle.cbegin() + count1);
									branch_mesh.erase(branch_mesh.begin() + count1);
									branch_count = branch_count - 1; // refresh number of branches
									break; // all done; exit algorithm
								}

							}
						}
						break;// exit algorithm
					}
				}
			}
		}

	// Stores branches connected to each node
	node_branches.resize(nodes_id.size());
	for (int count = 0; count< nodes_id.size(); count++) { // counter of nodes
			for ( int count1 = 0 ; count1<  branch_count;count1++){ // branch loop
				if (nodes_id.at(count) == branches_ind.at( count1).at(0) || nodes_id.at(count) == branches_ind.at( count1 ).at(branches_ind.at(count1).size() - 1)) {
					// branch starts or ends at the node
					node_branches.at(count).push_back(count1); // add index of branch
				}
			}
		}
		


	// Stores branches that form each mesh
	n_mesh=mesh_ind.size();
	mesh_branches.resize(n_mesh);
	for (int count = 0; count < n_mesh; count++)
		mesh_branches.at(count) = {};

	for (int count = 0; count < n_mesh;count++)// mesh loop
		for (int count1 = 0 ; count1 < branch_count ;count1++)// branch loop
			if (in(count , branch_mesh.at( count1 )))// the branch belongs to the current mesh
				mesh_branches.at(count).push_back( count1); // add index of accordant branch
				

	// Re - order branches
	for (int count_mesh = 1; count_mesh < n_mesh; count_mesh++) { // mesh loop
		if (mesh_branches.at(count_mesh).size() > 2) {// if there are only one or two branches, the current order is accepted
			vec_aux1.assign(mesh_branches.at(count_mesh).size(), 0); // stores the ordered branches in the mesh
			vec_aux1.at(1) = mesh_branches.at(count_mesh).at(1); // the first branch it finds is chosen as the first of the mesh
			// Looks for the following branches
			for (int count_branch1 = 0,ind_aux; count_branch1 < vec_aux1.size() - 1; count_branch1++){// after last branch comes first branch again(cycle)
				ind_aux = vec_aux1.at(count_branch1); // current branch index
				for (int count_branch2 = 1, ind_aux1; count_branch2 < mesh_branches.at(count_mesh ).size(); count_branch2++){ // branches to compare with // first branch cannot be
					ind_aux1 = mesh_branches.at(count_mesh ).at(count_branch2); // index of branch to compare
					if (ind_aux != ind_aux1){// does not ckeck itself
						if (branches_ind.at(ind_aux1).at(0) == branches_ind.at(ind_aux).at(branches_ind.at(ind_aux).size()-1) ){
							// second branch starts with the last object of first branch
							vec_aux1.at(count_branch1 + 1) = ind_aux1; // is the next element
							break; // next branch
						}
						else if (branches_ind.at( ind_aux1 ).at(branches_ind.at(ind_aux1).size()-1) == branches_ind.at(ind_aux ).at(branches_ind.at(ind_aux).size()-1)) {
							// second branch ends with the last object of first branch
							branch_aux.clear();
							for (int aux = 0; aux < branches_ind.at(ind_aux1).size(); aux++)
								branch_aux.push_back(branches_ind.at(ind_aux1).at(branches_ind.at(ind_aux1).size()-1 - aux));

							branches_ind.at(ind_aux1) = branch_aux; // inverts branch
							vec_aux1.at(count_branch1 + 1) = ind_aux1; // is the next element
							break; // next branch
						}
					}
				}
			}

			mesh_branches.at(count_mesh ) = vec_aux1; // ordered branches in the mesh
		}
	}


	// Adds every object's branch to the objects table
	objects_branches.resize(n_obj);
	for (int count = 0; count < n_obj; count++)
		objects_branches.at(count) = {};

	for (int count = 0; count < branch_count;count++) // branch loop
		for (int count1 = 0, aux; count1 < branches_ind.at(count).size(); count1++) { // objects in a branch loop
			aux = branches_ind.at(count).at(count1);
			objects_branches.at( aux).push_back(count);
		}


	// Converts branches to change indices by id's
	branches_id.resize(branch_count); // store all branches(id)
	for (int count = 0; count < branch_count; count++)
		branches_id.at(count) = {};
	for (int count = 0; count < branch_count;count++) // branch loop
		for (int count1 = 0; count1 < branches_ind.at(count).size();count1++) // objects in a branch
			branches_id.at(count).push_back(objects.at(branches_ind.at( count ).at(count1)).ID); // assign id
	//branches_id = branches_ind;

	// Inverts branches that are in the wrong direction
	// Flow direction is determined in some objects(and their branches) :
	// pumps and heat exchangers of type 'T out'
	for (int aux = 0,i; aux < pump_volum.size() + pump_turbo.size()-1 + heat_exch_Tout.size()-1/*+tank.size()-1*/; aux++) {
		if (aux < pump_volum.size() - 1) {
			i = aux;
			ind_outlet_aux.push_back({ pump_volum.at(i).Pump_volum,pump_volum.at(i).outlet_object });
		}
		else if (pump_volum.size() - 1 < aux < pump_turbo.size() - 1) {
			i = aux- (pump_volum.size()- 1);
			ind_outlet_aux.push_back({ pump_turbo.at(i).Pump_turbo,pump_turbo.at(i).outlet_object });
		}
		else if (pump_volum.size() - 1+pump_turbo.size() - 1 < aux <heat_exch_Tout.size()-1 ) {
			i = aux - (pump_volum.size() - 1 + pump_turbo.size() - 1);
			ind_outlet_aux.push_back({ heat_exch_Tout.at(i).heat_exch_Tout,heat_exch_Tout.at(i).outlet_object });
		}
		/*else if (pump_volum.size() - 1 + pump_turbo.size() - 1 +heat_exch_Tout.size() - 1 < aux < tank.size()) {
			i = aux - (pump_volum.size() - 1 + pump_turbo.size() - 1);
			ind_outlet_aux.push_back({ tank.at(i).tank,tank.at(i).outlet_object });
		}*/
	}
	//ind_outlet_aux.at(0) = { pump_volum.at(0).name, pump_volum.at(0).outlet_oject; pump_turbo.at(0).name ,  pump_turbo.at(0).outlet_oject; heat_exch_Tout.at(0).name , heat_exch_Tout.at(0).outlet_object/*; tank.obj_index tank.at(0).outlet*/ };
	// column 1: object index of all volumetric pumps, turbopumps, tanks and heat exchangers of type 'T out'
	// column 2 : id ot the outlets of those objects
	for (int count = 0, ind_aux, ind_aux1; count< ind_outlet_aux.size() ;count++){ // counter of objects that have a direction
		branch_aux = objects_branches.at( ind_outlet_aux.at(count).at(0) ); // branch that contains the object
		for( int aux = 0;aux<branches_ind.at(branch_aux.at(0)).size();aux++)
			if (branches_ind.at(branch_aux.at(0)).at(aux) == ind_outlet_aux.at(count).at(0)) {
				ind_aux = aux; // position of the object in the branch
				break;
			}
		for (int aux = 0; aux<branches_ind.at(branch_aux.at(0)).size(); aux++)
			if (branches_ind.at(branch_aux.at(0)).at(aux) == ind_outlet_aux.at(count).at(1)) {
				ind_aux1 = aux; // position of the outlet object in the branch
				break;
			}																	// Compare direction of branch and object
		if (ind_aux1 < ind_aux){ // the outlet is located earlier in the branch // opposite direction
							  // this works for tanks in the return branch because(empty < x) == false
			for (int aux = 0; aux < branches_id.at(branch_aux.at(0)).size(); aux++)
				vec_vec_aux.push_back(branches_id.at(branch_aux.at(0)).at(branches_id.at(branch_aux.at(0)).size()-1 - aux));
			branches_id.at( branch_aux.at(0) ) = vec_vec_aux; // invert branch
			vec_vec_aux.clear();
			vec_vec_aux.shrink_to_fit();
			vec_vec_aux = {};
			for (int aux = 0; aux < branches_ind.at(branch_aux.at(0)).size(); aux++)
				vec_vec_aux.push_back(branches_ind.at(branch_aux.at(0)).at(branches_ind.at(branch_aux.at(0)).size()-1 - aux));
			branches_ind.at(branch_aux.at(0) )= vec_vec_aux;
			vec_vec_aux.clear();
			vec_vec_aux.shrink_to_fit();
			vec_vec_aux = {};
		}
	}
	
	// branches vector contruction

	branches.resize(branches_id.size());

	for (int count = 0; count < branches.size(); count++) {
		branches.at(count).branch_ind = count;
		branches.at(count).objects.resize(branches_id.at(count).size());
		for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++) {
			branches.at(count).objects.at(count1).ID = branches_id.at(count).at(count1);
			for (int aux = 0; aux < objects.size(); aux++) {
				if (objects.at(aux).ID == branches.at(count).objects.at(count1).ID) {
					branches.at(count).objects.at(count1).Class = objects.at(aux).Class;
					branches.at(count).objects.at(count1).volume = objects.at(aux).volume;
					branches.at(count).objects.at(count1).name = objects.at(aux).name;
					branches.at(count).objects.at(count1).adjacent = objects.at(aux).adjacent;
					
				}
			}
		}
	}




}


int main()
{
	std::vector<strBranches> branches;
	int n_obj;
	std::vector<int> nodes_id;
	std::vector<strPump_volum> pump_volum;
	std::vector<strPump_turbo> pump_turbo;
	std::vector<strHeat_exch_Tout> heat_exch_Tout;
	/* std::vector<strTank> tank;*/
	std::vector <strObject> objects;

	n_obj = 15;

	nodes_id = { 2,5,6,11 };
	
	pump_volum.resize(6);

	pump_volum.at(0).Pump_volum = 1;
	pump_volum.at(0).coef0 = 1;
	pump_volum.at(0).coef1 = 0;
	pump_volum.at(0).coef2 = 0;
	pump_volum.at(0).coef3 = 0;
	pump_volum.at(0).head_max = 1000000;
	pump_volum.at(0).inlet_object = 1;
	pump_volum.at(0).outlet_object = 2;
	pump_volum.at(0).pump_speed = 105;
	branches.resize(6);


	branches.at(0).branch_ind = 0;
	branches.at(1).branch_ind = 1;
	branches.at(2).branch_ind = 2;
	branches.at(4).branch_ind = 3;
	branches.at(4).branch_ind = 4;
	branches.at(5).branch_ind = 5;


	branches.at(0).objects.resize(5);
	branches.at(1).objects.resize(5);
	branches.at(2).objects.resize(2);
	branches.at(3).objects.resize(3);
	branches.at(4).objects.resize(4);
	branches.at(5).objects.resize(4);

	branches.at(0).objects.at(0).Class = "pump_volum";
	branches.at(0).objects.at(0).ID = 0;
	branches.at(0).objects.at(0).name = "Main_Pump";
	branches.at(0).objects.at(0).volume = 1E-3;
	branches.at(0).objects.at(0).inlet_object = 12;
	branches.at(0).objects.at(0).outlet_oject = 1;
	branches.at(0).objects.at(0).adjacent = { 12,1 };

	branches.at(0).objects.at(1).Class = "pipe";
	branches.at(0).objects.at(1).ID = 1;
	branches.at(0).objects.at(1).name = "Enigine_Inlet";
	branches.at(0).objects.at(1).volume = 8.8357E-5;
	branches.at(0).objects.at(1).inlet_object = 0;
	branches.at(0).objects.at(1).outlet_oject = 2;
	branches.at(0).objects.at(1).adjacent = { 0,2 };

	branches.at(0).objects.at(2).Class = "pipe";
	branches.at(0).objects.at(2).ID = 12;
	branches.at(0).objects.at(2).name = "Return_Pipe";
	branches.at(0).objects.at(2).volume = 6.2832E-4;
	branches.at(0).objects.at(2).inlet_object = 11;
	branches.at(0).objects.at(2).outlet_oject = 0;
	branches.at(0).objects.at(2).adjacent = { 11,0 };

	branches.at(0).objects.at(3).Class = "node";
	branches.at(0).objects.at(3).ID = 11;
	branches.at(0).objects.at(3).adjacent = { 10,12,14 };

	branches.at(0).objects.at(3).Class = "node";
	branches.at(0).objects.at(3).ID = 2;
	branches.at(0).objects.at(3).adjacent = { 1,3,4 };

	branches.at(1).objects.at(0).Class = "thermostat";
	branches.at(1).objects.at(0).ID = 7;
	branches.at(1).objects.at(0).name = "Thermostat_rad";
	branches.at(1).objects.at(0).volume = 0;
	branches.at(1).objects.at(0).inlet_object = 6;
	branches.at(1).objects.at(0).outlet_oject = 9;
	branches.at(1).objects.at(0).adjacent = { 6,9 };

	branches.at(1).objects.at(1).Class = "heat_exch_fix";
	branches.at(1).objects.at(1).ID = 9;
	branches.at(1).objects.at(1).name = "Radiator";
	branches.at(1).objects.at(1).volume = 5E-4;
	branches.at(1).objects.at(1).inlet_object = 7;
	branches.at(1).objects.at(1).outlet_oject = 10;
	branches.at(1).objects.at(1).adjacent = { 7,10 };

	branches.at(1).objects.at(2).Class = "pipe";
	branches.at(1).objects.at(2).ID = 10;
	branches.at(1).objects.at(2).name = "Radiator_outlet";
	branches.at(1).objects.at(2).volume = 3.1416E-4;
	branches.at(1).objects.at(2).inlet_object = 9;
	branches.at(1).objects.at(2).outlet_oject = 11;
	branches.at(1).objects.at(2).adjacent = { 9,11 };

	branches.at(1).objects.at(3).Class = "node";
	branches.at(1).objects.at(3).ID = 11;
	branches.at(1).objects.at(3).adjacent = { 10,12,14 };

	branches.at(1).objects.at(4).Class = "node";
	branches.at(1).objects.at(4).ID = 6;
	branches.at(1).objects.at(4).adjacent = { 5,7,8 };

	branches.at(2).objects.at(0).Class = "node";
	branches.at(2).objects.at(0).ID = 5;
	branches.at(2).objects.at(0).adjacent = { 3,6,13 };

	branches.at(2).objects.at(1).Class = "node";
	branches.at(2).objects.at(1).ID = 6;
	branches.at(2).objects.at(1).adjacent = { 5,7,8 };




	branches.at(3).objects.at(0).Class = "heat_exch_fix";
	branches.at(3).objects.at(0).ID = 3;
	branches.at(3).objects.at(0).name = "Block_cyl1";
	branches.at(3).objects.at(0).volume = 5E-4;
	branches.at(3).objects.at(0).inlet_object = 2;
	branches.at(3).objects.at(0).outlet_oject = 5;
	branches.at(3).objects.at(0).adjacent = { 2,5 };

	branches.at(3).objects.at(1).Class = "node";
	branches.at(3).objects.at(1).ID = 2;
	branches.at(3).objects.at(1).adjacent = { 1,3,4 };

	branches.at(3).objects.at(2).Class = "node";
	branches.at(3).objects.at(2).ID = 5;
	branches.at(3).objects.at(2).adjacent = { 3,6,13 };

	branches.at(4).objects.at(0).Class = "heat_exch_fix";
	branches.at(4).objects.at(0).ID = 4;
	branches.at(4).objects.at(0).name = "Head_cyl1";
	branches.at(4).objects.at(0).volume = 5E-4;
	branches.at(4).objects.at(0).inlet_object = 2;
	branches.at(4).objects.at(0).outlet_oject = 13;
	branches.at(4).objects.at(0).adjacent = { 2,13 };

	branches.at(4).objects.at(1).Class = "heat_exch_fix";
	branches.at(4).objects.at(1).ID = 13;
	branches.at(4).objects.at(1).name = "Head_cyl2";
	branches.at(4).objects.at(1).volume = 5E-4;
	branches.at(4).objects.at(1).inlet_object = 4;
	branches.at(4).objects.at(1).outlet_oject = 5;
	branches.at(4).objects.at(1).adjacent = { 4,5 };

	branches.at(4).objects.at(2).Class = "node";
	branches.at(4).objects.at(2).ID = 2;
	branches.at(4).objects.at(2).adjacent = { 1,3,4 };

	branches.at(4).objects.at(3).Class = "node";
	branches.at(4).objects.at(3).ID = 5;
	branches.at(4).objects.at(3).adjacent = { 3,6,13 };

	branches.at(5).objects.at(0).Class = "thermostat";
	branches.at(5).objects.at(0).ID = 8;
	branches.at(5).objects.at(0).name = "Thermostat_byp";
	branches.at(5).objects.at(0).volume = 0;
	branches.at(5).objects.at(0).inlet_object = 6;
	branches.at(5).objects.at(0).outlet_oject = 14;
	branches.at(5).objects.at(0).adjacent = { 6,14 };


	branches.at(5).objects.at(1).Class = "pipe";
	branches.at(5).objects.at(1).ID = 14;
	branches.at(5).objects.at(1).name = "Bypass_Pipe";
	branches.at(5).objects.at(1).volume = 3.1416E-4;
	branches.at(5).objects.at(1).inlet_object = 8;
	branches.at(5).objects.at(1).outlet_oject = 11;
	branches.at(5).objects.at(1).adjacent = { 8,11 };


	branches.at(5).objects.at(2).Class = "node";
	branches.at(5).objects.at(2).ID = 11;
	branches.at(5).objects.at(2).adjacent = { 10,12,14 };

	branches.at(5).objects.at(2).Class = "node";
	branches.at(5).objects.at(2).ID = 6;
	branches.at(5).objects.at(2).adjacent = { 5,7,8 };

	objects.resize(15);

	objects.at(0).ID = 0;
	objects.at(0).name = "Main_Pump";
	objects.at(0).Class = "pump_volum";
	objects.at(0).volume = 10E-3;
	objects.at(0).adjacent = { 12,1 };

	objects.at(1).ID = 1;
	objects.at(1).name = "Engine_Inlet";
	objects.at(1).Class = "pipe";
	objects.at(1).volume = 8.8235E-5;
	objects.at(1).adjacent = { 0,2 };

	objects.at(2).ID = 2;
	objects.at(2).name = "Splitter_Engine";
	objects.at(2).Class = "node";
	objects.at(2).volume = 0;
	objects.at(2).adjacent = { 1,3,4 };

	objects.at(3).ID = 3;
	objects.at(3).name = "Block_cyl1";
	objects.at(3).Class = "heat_exch_fix";
	objects.at(3).volume = 5E-4;
	objects.at(3).adjacent = { 2,5 };

	objects.at(4).ID = 4;
	objects.at(4).name = "Head_cyl1";
	objects.at(4).Class = "heat_exch_fix";
	objects.at(4).volume = 5E-4;
	objects.at(4).adjacent = { 2,13 };

	objects.at(5).ID = 5;
	objects.at(5).name = "Mixer_Engine";
	objects.at(5).Class = "node";
	objects.at(5).volume = 0;
	objects.at(5).adjacent = { 3,13,6 };

	objects.at(6).ID = 6;
	objects.at(6).name = "Splitter_Thermo";
	objects.at(6).Class = "node";
	objects.at(6).volume = 0;
	objects.at(6).adjacent = { 5,7,8 };

	objects.at(7).ID = 7;
	objects.at(7).name = "Thermostat_rad";
	objects.at(7).Class = "thermostat";
	objects.at(7).volume = 0;
	objects.at(7).adjacent = { 6,9 };

	objects.at(8).ID = 8;
	objects.at(8).name = "Thermostat_byp";
	objects.at(8).Class = "thermostat";
	objects.at(8).volume = 0;
	objects.at(8).adjacent = { 6,14 };

	objects.at(9).ID = 9;
	objects.at(9).name = "Radiator";
	objects.at(9).Class = "heat_exch_fix";
	objects.at(9).volume = 5E-4;
	objects.at(9).adjacent = { 7,10 };

	objects.at(10).ID = 10;
	objects.at(10).name = "Radiator_Outlet";
	objects.at(10).Class = "pump_volum";
	objects.at(10).volume = 3.1416E-4;
	objects.at(10).adjacent = { 9,11 };

	objects.at(11).ID = 11;
	objects.at(11).name = "Mixer_Thermostat";
	objects.at(11).Class = "node";
	objects.at(11).volume = 0;
	objects.at(11).adjacent = { 14,10,12 };

	objects.at(12).ID = 12;
	objects.at(12).name = "Return_pipe";
	objects.at(12).Class = "pipe";
	objects.at(12).volume = 6.2834E-4;
	objects.at(12).adjacent = { 11,0 };

	objects.at(13).ID = 13;
	objects.at(13).name = "Head_cyl2";
	objects.at(13).Class = "heat_exch_fix";
	objects.at(13).volume = 5E-4;
	objects.at(13).adjacent = { 4,5 };

	objects.at(14).ID = 14;
	objects.at(14).name = "Bypass_Pipe";
	objects.at(14).Class = "pipe";
	objects.at(14).volume = 3.1416E-4;
	objects.at(14).adjacent = { 8,11 };





	HydroNet_Components(/*branches,*/objects, n_obj,  nodes_id,  pump_volum,  pump_turbo, heat_exch_Tout/*,tank*/);
    return 0;
}

