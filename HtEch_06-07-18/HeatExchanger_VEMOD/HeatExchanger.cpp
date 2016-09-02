#include <cmath>
//#include <iostream>
#include <fstream>
#include <vector>
#include <string>		// std::stoi
#include <algorithm>	// to avoid 'min/max not a member of std' error
#include "HeatExchanger.h"


void HeatExchanger_model::create_HeatExchangers(int n_htxs)
{
	int line;
	int count;
	//int N;

	std::vector<std::string> HeatExchangers_data_archives(n_htxs);
	std::vector<std::string> HeatExchangers_table_archives(n_htxs);

	HeatExchangers_data_archives.at(0) = "HeatExchanger_inputs.txt";

	HeatExchangers_table_archives.at(0) = "HeatExchanger_efficiency_table.txt";				// Cambio para prieva de "HeatExchanger_efficiency_table.txt" a "HtEchngr_Eff_Tbl.txt"
	
	HeatExchangers_data_archives.at(1) = "HeatExchanger_inputs_2.txt";

	HeatExchangers_table_archives.at(1) = "HeatExchanger_efficiency_table_2.txt";
	
	HeatExchangers.resize(n_htxs);															// Crate the same number of entrys in the vector as number of Heat Exchangers
	inputs.resize(n_htxs);																	// Crate the same number of entrys in the vector as number of Heat Exchangers
	results.resize(n_htxs);																	// Crate the same number of entrys in the vector as number of Heat Exchangers
		
	// for (every file / heat exchanger) {
	// Strings, integers, doubles and booleans
	// Read lines of data
	//HeatExchanger_file.open("HeatExchanger_inputs.txt", std::ifstream::out );
	
	for (count = 0; count < n_htxs; count++) {


		std::ifstream HeatExchanger_file(HeatExchangers_data_archives.at(count), std::ifstream::out);
																				
		std::vector<std::string> HeatExchanger_data;
		std::string data;
				

		while (HeatExchanger_file) {
			std::getline(HeatExchanger_file, data);			//Reads an element of HeatExchangers_file

			HeatExchanger_data.push_back(data);				// Save a new element at the end of HeatExchanger_data
			
		}
		HeatExchanger_file.close();
		
	

		// String data goes directly to structure
		HeatExchangers.at(count).name = HeatExchanger_data[2];
		HeatExchangers.at(count).fluid_type_1 = HeatExchanger_data[3];
		HeatExchangers.at(count).fluid_type_2 = HeatExchanger_data[4];

		// Integer data is converted and stored
		HeatExchangers.at(count).ID_htx = std::stoi(HeatExchanger_data[0]);
		HeatExchangers.at(count).ID_htx_HydroNet = std::stoi(HeatExchanger_data[1]);
		HeatExchangers.at(count).exchanger_type = std::stoi(HeatExchanger_data[8]);
		HeatExchangers.at(count).n_shell_passes = std::stoi(HeatExchanger_data[9]);

		// Double data is converted and stored
		HeatExchangers.at(count).hydraulic_resistance = std::stod(HeatExchanger_data[5]);
		HeatExchangers.at(count).head_loss = std::stod(HeatExchanger_data[6]);
		HeatExchangers.at(count).UA = std::stod(HeatExchanger_data[7]);

		// Boolean data is converted and stored
		if (HeatExchanger_data[10] == "1") {
			HeatExchangers.at(count).flag_efficiency_table = true;
		}
		else if (HeatExchanger_data[10] == "0") {
			HeatExchangers.at(count).flag_efficiency_table = false;
		}

		if (HeatExchangers.at(count).flag_efficiency_table == true) {
			// Efficiency table: arrays
			// Read all data in the table
			int lines;

			std::string line_scan;

			HeatExchanger_file.open(HeatExchangers_table_archives.at(count), std::ifstream::out);


			lines = 0;

			while (!HeatExchanger_file.eof()) {
				getline(HeatExchanger_file, line_scan);
				lines++;	// lines counter
				}

			HeatExchanger_file.close();

			std::vector<std::string> efficiency_data;
			std::vector<std::string> NTU_data;
			std::vector<std::string> Cr_data;

			HeatExchanger_file.open(HeatExchangers_table_archives.at(count), std::ifstream::out);
			


			while (HeatExchanger_file) {

				std::getline(HeatExchanger_file, data, '\t');				//Reads an element of HeatExchangers_file, stops at TAB

				NTU_data.push_back(data);									// Save a new element at the end of NTU_data

				std::getline(HeatExchanger_file, data, '\t');				//Reads an element of HeatExchangers_file, stops at TAB
				
				Cr_data.push_back(data);									// Save a new element at the end of Cr_data

				std::getline(HeatExchanger_file, data, '\n');				//Reads an element of HeatExchangers_file, stops at END OF LINE

				efficiency_data.push_back(data);							// Save a new element at the end of efficiency_data

			}
			HeatExchanger_file.close();




			// Array data is converted and stored
			HeatExchangers.at(count).NTU.resize(lines);						// Crate the same number of entrys in the vector as number of lines at the table
			HeatExchangers.at(count).Cr.resize(lines);						// Crate the same number of entrys in the vector as number of lines at the table
			HeatExchangers.at(count).efficiency.resize(lines);				// Crate the same number of entrys in the vector as number of lines at the table


			line = 0;
			while (line < lines) {
				HeatExchangers.at(count).NTU[line] = std::stod(NTU_data[line]);

				HeatExchangers.at(count).Cr[line] = std::stod(Cr_data[line]);

				HeatExchangers.at(count).efficiency[line] = std::stod(efficiency_data[line]);

				line++;
			}

			// Rows in the table: NTU, Cr & efficicency
			HeatExchangers.at(count).n_rows_table = line;
		}
		results.at(count).HeatExchanger_ID = HeatExchangers.at(count).ID_htx;
	}
	//}
	
}


void HeatExchanger_model::calculate_heat_exchange()
{
	int count;
	int index_aux;
	int n_htxch;						//Heat exchangers number
	int ID_htxch_index;					// Heat Exchanger index
	double mCp_1;
	double mCp_2;
	double mCp_min;
	double mCp_max;
	double NTU;
	double Cr;
	double efficiency;
	efficiency = 0;
	ID_htxch_index = 0;

	
	n_htxch = HeatExchangers.size();
	
	for (count = 0;count<n_htxch;count++){
	
	for (index_aux = 0;index_aux<n_htxch;index_aux++){
		
		if (inputs.at(count).HeatExchanger_ID == HeatExchangers.at(index_aux).ID_htx) {
			ID_htxch_index = index_aux;

			break;
		}
		
		
	}
	
	

	mCp_1 = inputs.at(count).mass_flow_1 * cp(HeatExchangers.at(ID_htxch_index).fluid_type_1, inputs.at(count).Temperature_inlet_1);
	mCp_2 = inputs.at(count).mass_flow_2 * cp(HeatExchangers.at(ID_htxch_index).fluid_type_2, inputs.at(count).Temperature_inlet_2);
	mCp_min = std::min(mCp_1, mCp_2);
	mCp_max = std::max(mCp_1, mCp_2);

	NTU = HeatExchangers.at(ID_htxch_index).UA / mCp_min;
	Cr = mCp_min / mCp_max;

	if (HeatExchangers.at(ID_htxch_index).flag_efficiency_table) {	// look up efficiency in a table
		efficiency = lookup_table_2dim(HeatExchangers.at(ID_htxch_index).n_rows_table, NTU, Cr, &HeatExchangers.at(ID_htxch_index).NTU, &HeatExchangers.at(ID_htxch_index).Cr, &HeatExchangers.at(ID_htxch_index).efficiency);
	}
	else {	// no table; obtain efficiency analiticaly considering the heat exchanger characteristics
		switch (HeatExchangers.at(ID_htxch_index).exchanger_type) {
		case 1:	// Simple co-current
			efficiency = (1 - exp(-NTU * (1 + Cr))) / (1 + Cr);
			break;
		case 2:	// Simple counter-current
			if (Cr == 1) {
				efficiency = 1 / (1 + 1 / NTU);
			}
			else {
				efficiency = (1 - exp(-NTU * (1 - Cr))) / (1 - Cr * exp(-NTU * (1 - Cr)));
			}
			break;
		case 3:	// Shell and tube
			efficiency = 1 / (2 * (1 + Cr + sqrt(1 + Cr * Cr) * (1 + exp(-NTU * sqrt(1 + Cr * Cr)) / (1 - exp(-NTU * sqrt(1 + Cr * Cr))))));
			if (HeatExchangers.at(ID_htxch_index).n_shell_passes > 1) {	// multiple passes through shell (for 1, the previous expression is obtained)
				efficiency = (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchangers.at(ID_htxch_index).n_shell_passes) - 1)
					/ (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchangers.at(ID_htxch_index).n_shell_passes) - Cr);
			}
			break;
		case 4:   // Crossed flows (simple pass, fluids not mixed)
			efficiency = 1 - exp(1 / Cr * pow(NTU, 0.22) * (exp(-Cr * pow(NTU, 0.78)) - 1));
			break;
		}
	}

	results.at(count).exchanged_heat_1 = efficiency * mCp_min * (inputs.at(count).Temperature_inlet_2 - inputs.at(count).Temperature_inlet_1); // (+) heat gained from 2 to 1; (-) heat lost from 1 to 2
	results.at(count).exchanged_heat_2 = -results.at(count).exchanged_heat_1; // (+) heat gained from 1 to 2; (-) heat lost from 2 to 1

	results.at(count).Temperature_outlet_1 = inputs.at(count).Temperature_inlet_1 + results.at(count).exchanged_heat_1 / mCp_1;
	results.at(count).Temperature_outlet_2 = inputs.at(count).Temperature_inlet_2 + results.at(count).exchanged_heat_2 / mCp_2;
	}

}


double HeatExchanger_model::lookup_table_2dim(int n_rows, double x, double y, std::vector<double>* x_column, std::vector<double> *y_column, std::vector<double>* z_column)
{

	int count;
	int comb1;
	int comb2;
	int comb3;
	int comb4;
	int smallest;
	std::vector<double> dist_x;
	std::vector<double> dist_y;
	std::vector<double> dist_xy;
	double z;
	z = 0;
	dist_x.resize(n_rows);						// Assigns size to vector
	dist_y.resize(n_rows);						// Assigns size to vector
	dist_xy.resize(n_rows);						// Assigns size to vector

	std::vector<double> X_column;
	std::vector<double> Y_column;
	std::vector<double> Z_column;
	
	X_column=*x_column;
	Y_column=*y_column;
	Z_column=*z_column;
	
	// If x and y are outside the table, bring them inside the limits
	if (x > *std::max_element(X_column.begin(),X_column.end()) ){
		x = *std::max_element(X_column.begin(),X_column.end());
	}
	else if (x < *std::min_element(X_column.begin(),X_column.end())) {
		x = *std::min_element(X_column.begin(),X_column.end());
	}
	if (y > *std::max_element(Y_column.begin(),Y_column.end())) {
		y = *std::max_element(Y_column.begin(),Y_column.end());
	}
	else if (y < *std::min_element(Y_column.begin(),Y_column.end())) {
		y = *std::min_element(Y_column.begin(),Y_column.end());
	}

	// Obtain distances to x and y
	for (count = 0; count < n_rows; count++) {
		dist_x[count] = abs(x - X_column[count]);
	}
	for (count = 0; count < n_rows; count++) {
		dist_y[count] = abs(y - Y_column[count]);
	}

	// Calculate distances to x y combinations (rows) from the queried values
	// Hypotenuse from cathetus
	for (count = 0; count < n_rows; count++) {
		dist_xy[count] = sqrt(dist_x[count] * dist_x[count] + dist_y[count] * dist_y[count]);
	}

	// Find the 4 closest points of combination x y to the queried combination
	if (n_rows > 3) {
		smallest = 0;
		for (count = 1; count < n_rows; count++) {
			if (dist_xy[count] < dist_xy[smallest]) {
				smallest = count;
			}
		}
		comb1 = smallest;	// index of the closest point

		if (comb1 == 0) {
			smallest = 1;
		}
		else {
			smallest = 0;
		}
		for (count = smallest + 1; count < n_rows; count++) {
			if ((dist_xy[count] < dist_xy[smallest]) && (count != comb1)) {
			smallest = count;
			}
		}
		comb2 = smallest;	// index of the second closest point

		if (((comb1 == 0) || (comb2 == 0)) && ((comb1 == 1) || (comb2 == 1))) {
			smallest = 2;
		}
		else if ((comb1 == 0) || (comb2 == 0)) {
			smallest = 1;
		}
		else {
			smallest = 0;
		}
		for (count = smallest + 1; count < n_rows; count++) {
			if ((dist_xy[count] < dist_xy[smallest]) && (count != comb1) && (count != comb2)) {
				smallest = count;
			}
		}
		comb3 = smallest;	// index of the third closest point

		if (((comb1 == 0) || (comb2 == 0) || (comb3 == 0)) && ((comb1 == 1) || (comb2 == 1) || (comb3 == 1)) && ((comb1 == 2) || (comb2 == 2) || (comb3 == 2))) {
			smallest = 3;
		}
		else if (((comb1 == 0) || (comb2 == 0) || (comb3 == 0)) && ((comb1 == 1) || (comb2 == 1) || (comb3 == 1))) {
			smallest = 2;
		}
		else if (((comb1 == 0) || (comb2 == 0) || (comb3 == 0))) {
			smallest = 1;
		}
		else {
			smallest = 0;
		}
		for (count = smallest + 1; count < n_rows; count++) {
			if ((dist_xy[count] < dist_xy[smallest]) && (count != comb1) && (count != comb2) && (count != comb3)) {
				smallest = count;
			}
		}
		comb4 = smallest;	// index of the fourth closest point


		// Interpolate with distances to obtain the queried value
		z = (Z_column[comb1] / dist_xy[comb1] + Z_column[comb2] / dist_xy[comb2] + Z_column[comb3] / dist_xy[comb3] + Z_column[comb4] / dist_xy[comb4])
			/ (1 / dist_xy[comb1] + 1 / dist_xy[comb2] + 1 / dist_xy[comb3] + 1 / dist_xy[comb4]);
	}
	else if (n_rows == 3) {
		z = (Z_column[0] / dist_xy[0] + Z_column[1] / dist_xy[1] + Z_column[2] / dist_xy[2]) / (1 / dist_xy[0] + 1 / dist_xy[1] + 1 / dist_xy[2]);
	}
	else if (n_rows == 2) {
		z = (Z_column[0] / dist_xy[0] + Z_column[1] / dist_xy[1]) / (1 / dist_xy[0] + 1 / dist_xy[1]);
	}
	else if  (n_rows == 1) {
		z = Z_column[0];
	}

	return z;
}



double HeatExchanger_model::cp(std::string fluid_type, double Temperature_inlet)
{
	return 4181.3;	// for tests: water [J/kg/K]
}

int main() {
	int n_htxs;
	int count;


	n_htxs = 2;


	HeatExchanger_model objeto;
	

	objeto.create_HeatExchangers(n_htxs);

	for (count = 0; count < n_htxs; count++ ){

		objeto.inputs.at(count).mass_flow_1 = 0.01;				// Test - Eliminate when working
		objeto.inputs.at(count).mass_flow_2 = 0.5;					// Test - Eliminate when working
		objeto.inputs.at(count).Temperature_inlet_1 = 100;			// Test - Eliminate when working
		objeto.inputs.at(count).Temperature_inlet_2 = 20;			// Test - Eliminate when working

		
	objeto.calculate_heat_exchange();
}
	return 0;
}
