// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <iostream>

bool any(std::vector<bool> A) {
	bool any = false;
	for (int i = 0; i < A.size(); i++)
		if (A.at(i) == true)
			any = true;
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

int main()
{
	std::vector<std::vector<bool>> edges, mat_aux;
	int obj_count,start_obj,row_aux,column,path_count, n_obj, prev_row;
	std::vector<int> path, nodes_id, path_aux;
	std::vector<std::vector<int>>paths_ind;
	bool flag_path,flag_aux;

	nodes_id = { 2,5,6,11 };

	edges.resize(15);

	for (int i = 0; i < 15; i++)
		edges.at(i).assign(15,false);

	edges.at(0).at(1) = true;
	edges.at(0).at(12) = true;
	edges.at(1).at(0) = true;
	edges.at(1).at(2) = true;
	edges.at(2).at(1) = true;
	edges.at(2).at(3) = true;
	edges.at(2).at(4) = true;
	edges.at(3).at(2) = true;
	edges.at(3).at(5) = true;
	edges.at(4).at(2) = true;
	edges.at(4).at(13) = true;
	edges.at(5).at(3) = true;
	edges.at(5).at(6) = true;
	edges.at(5).at(13) = true;
	edges.at(6).at(5) = true;
	edges.at(6).at(7) = true;
	edges.at(6).at(8) = true;
	edges.at(7).at(6) = true;
	edges.at(7).at(9) = true;
	edges.at(8).at(6) = true;
	edges.at(8).at(14) = true;
	edges.at(9).at(7) = true;
	edges.at(9).at(10) = true;
	edges.at(10).at(9) = true;
	edges.at(10).at(11) = true;
	edges.at(11).at(10) = true;
	edges.at(11).at(12) = true;
	edges.at(11).at(14) = true;
	edges.at(12).at(0) = true;
	edges.at(12).at(11) = true;
	edges.at(13).at(4) = true;
	edges.at(13).at(5) = true;
	edges.at(14).at(8) = true;
	edges.at(14).at(11) = true;



														std::cout << "\t" << "\t";
														for (int i = 0; i < 15; i++) {
															if (i<10)
																std::cout << i << " " << " ";
															else
																std::cout << i << " " ;
														}
														std::cout << "\n" << "\n";
	
														for (int j = 0; j < 15; j++) {
															std::cout << "\t" << j << "\t";
															for (int k = 0; k < 15; k++) 
																std::cout << (edges.at(j).at(k)) << " " << " ";
															
															std::cout << "\n";
														}
														std::cout << "\n" << "\n" << "\n";
	
	n_obj = 15;
	obj_count = 0;
	start_obj = 2;
	row_aux = 2;
	column = -1;
	mat_aux = edges;
	path = { start_obj };
	path_count=0;
	

	
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

		while ((column >= n_obj-1) && (row_aux != start_obj)) { // loop to step back in a path
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



	//while (!((column == 15) && (start_object==row)))
	//{
	//	if (column!=15 && mat_aux.at(row).at(column)) {
	//		obj_count++;
	//		path.push_back(column);
	//		mat_aux.at(row).at(column) = false;
	//		mat_aux.at(column).at(row) = false;
	//		row = column;
	//		column = 0;
	//		flag = true;
	//	}
	//	else {
	//		column++;
	//	}
	//	if ((start_object==row && flag) || column >= 15 ) {
	//		path_count++;
	//		path_ind.push_back(path);

	//		for (int i = 0; i < path.size(); i++) {
				std::cout << "\t" << "\t";
				for (int ii = 0; ii < 15; ii++) {
					if (ii<10)
						std::cout << ii << " " << " ";
					else
						std::cout << ii << " ";
				}
				std::cout << "\n" << "\n";

				for (int j = 0; j < 15; j++) {
					std::cout << "\t" << j << "\t";
					for (int k = 0; k < 15; k++)
						std::cout << (mat_aux.at(j).at(k)) << " " << " ";

					std::cout << "\n";
				}
				std::cout << "\n" << "\n" << "\n";
	//			if (!in(i,nodes_id)) 
	//				mat_aux.at(path.at( i)) = edges.at(path.at( i));
	//				
	//			else {
	//				if (!any(mat_aux.at(path.at(i)))) {
	//					mat_aux.at(path.at(i)) = edges.at(path.at(i));
	//				}
	//				else
	//					mat_aux.at(path.at(i)).at(path.at(i - 1)) = true;

	//			}



	//			//if (!any(mat_aux.at(path.at( i)))) 
	//			//	mat_aux.at(path.at( i)) = edges.at(path.at( i));
	//			//	
	//			//else if (path.at( i) != start_object && flag_aux){
	//			//	mat_aux.at(path.at( i)).at(path.at( i+1)) = true;
	//			//	//flag_aux = false;
	//			//}
	//			//else if (i<path.size()-1) 
	//			//	mat_aux.at(path.at(i)).at(path.at(i + 1)) = true;
	//		/*	else if (i<path.size() -1)
	//				mat_aux.at(path.at(i)).at(path.at(i + 1)) = true;*/
	//			/*else if (i!=1)

	//				mat_aux.at(path.at(path.size() - i)).at(path.at(path.size() - i +1)) = true;*/
	//			
	//			//if (/*any(mat_aux.at(path.at(i))) &&*/ (path.at(i) != start_object) && flag_aux) {

	//			//	mat_aux.at(path.at(i)).at(path.at(i + 1)) = true;
	//			//	flag_aux = false;
	//			//}

	//			// else //if (!any(mat_aux.at(path.at(i))))
	//			//	mat_aux.at(path.at(i)) = edges.at(path.at(i));


	//		}
	//		path.clear();
	//		path = { start_object };
	//		row = start_object;
	//		column = 0;
	//		flag = false;
	//		flag_aux = true;
	//	}
	//}

return 0;
}

