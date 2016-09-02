#include <cmath>		// for sqrt
#include <fstream>		// to read file
#include <vector>
#include <string>		// string use; necessary for std::stoi
#include <algorithm>	// to avoid 'min/max not a member of std' error
#include "HeatExchanger.h"


void HeatExchanger_model::create_HeatExchangers(int n_htxs)
{
	int count;	// counter of heat exchangers
	int lines;	// number of lines in the efficiency table
	int line;	// counter for efficiency table data
	std::string line_scan;	// content of a line as a string

	// Names of data files
	std::vector<std::string> HeatExchangers_data_archives(n_htxs);	// properties
	std::vector<std::string> HeatExchangers_table_archives(n_htxs);	// efficiency table


	// For testing; only one file/heat exchanger
	HeatExchangers_data_archives.at(0) = "HeatExchanger_inputs.txt";				
	HeatExchangers_table_archives.at(0) = "HeatExchanger_efficiency_table.txt";
	

	// Vector dimensions: number of heat exchangers
	HeatExchangers.resize(n_htxs);
	inputs.resize(n_htxs);
	results.resize(n_htxs);
	
	
	// Heat exchanger loop
	for (count = 0; count < n_htxs; count++) {

		std::vector<std::string> HeatExchanger_data;	// vector to store file data as strings
		std::ifstream HeatExchanger_file(HeatExchangers_data_archives.at(count), std::ifstream::out);


		// Read data
		while (!HeatExchanger_file.eof()) {
			std::getline(HeatExchanger_file, line_scan);	// read a line of HeatExchangers_file
			HeatExchanger_data.push_back(line_scan);		// store the new line at the end of HeatExchanger_data
		}
		HeatExchanger_file.close();
		
		
		// String data goes directly to structure
		HeatExchangers.at(count).name = HeatExchanger_data[1];
		HeatExchangers.at(count).fluid_type_1 = HeatExchanger_data[2];
		HeatExchangers.at(count).fluid_type_2 = HeatExchanger_data[3];

		// Integer data is converted and stored
		HeatExchangers.at(count).ID_htx = std::stoi(HeatExchanger_data[0]);
		HeatExchangers.at(count).exchanger_type = std::stoi(HeatExchanger_data[5]);
		HeatExchangers.at(count).n_shell_passes = std::stoi(HeatExchanger_data[6]);

		// Double data is converted and stored
		HeatExchangers.at(count).UA = std::stod(HeatExchanger_data[4]);

		// Boolean data is converted and stored
		if (HeatExchanger_data[7] == "1") {
			HeatExchangers.at(count).flag_efficiency_table = true;
		}
		else if (HeatExchanger_data[7] == "0") {
			HeatExchangers.at(count).flag_efficiency_table = false;
		}


		// Efficiency table data
		if (HeatExchangers.at(count).flag_efficiency_table == true){	// efficiency must be obtained from a table

			std::vector<std::string> efficiency_data;
			std::vector<std::string> NTU_data;
			std::vector<std::string> Cr_data;


			// Read data
			HeatExchanger_file.open(HeatExchangers_table_archives.at(count), std::ifstream::out);
			while (!HeatExchanger_file.eof()) {

				std::getline(HeatExchanger_file, line_scan, '\t');	// read an element of HeatExchangers_file, stop at TAB
				NTU_data.push_back(line_scan);						// store the new element at the end of NTU_data

				std::getline(HeatExchanger_file, line_scan, '\t');	// read an element of HeatExchangers_file, stop at TAB
				Cr_data.push_back(line_scan);						// store the new element at the end of Cr_data

				std::getline(HeatExchanger_file, line_scan, '\n');	// read an element of HeatExchangers_file, stop at END OF LINE
				efficiency_data.push_back(line_scan);				// store the new element at the end of efficiency_data
			}
			HeatExchanger_file.close();


			// Vector dimensions
			lines = efficiency_data.size();
			HeatExchangers.at(count).NTU.resize(lines);
			HeatExchangers.at(count).Cr.resize(lines);
			HeatExchangers.at(count).efficiency.resize(lines);


			// Table data is converted and stored
			line = 0;
			while (line < lines) {
				HeatExchangers.at(count).NTU[line] = std::stod(NTU_data[line]);
				HeatExchangers.at(count).Cr[line] = std::stod(Cr_data[line]);
				HeatExchangers.at(count).efficiency[line] = std::stod(efficiency_data[line]);
				line++;
			}
			HeatExchangers.at(count).n_rows_table = lines;
			
		}

		// Store heat exchanger ID in results vector
		results.at(count).HeatExchanger_ID = HeatExchangers.at(count).ID_htx;
	}
	
}


void HeatExchanger_model::calculate_heat_exchange()
{
	int count;				// counter of heat exchangers
	int n_htxch;			// number of heat exchangers
	int prop_htxch_index;	// index of a heat exchanger in the properties vector
	double mCp_1;			// Mass flow times heat capacity of the first fluid
	double mCp_2;			// Mass flow times heat capacity of the second fluid
	double mCp_min;			// Minimum mass flow times heat capacity of both fluids
	double mCp_max;			// Maximum mass flow times heat capacity of both fluids
	double NTU;				// Parameter NTU of heat ezchanger
	double Cr;				// Ratio of m*Cp
	double efficiency = 0;	// Heat exchanger efficiency
	
	
	// Heat exchanger loop
	n_htxch = HeatExchangers.size();
	for (count = 0; count < n_htxch; count++){
		
		// Check if indexes of inputs and properties vectors are the same, using ID
		if (inputs.at(count).HeatExchanger_ID == HeatExchangers.at(count).ID_htx) {	// same ID
			prop_htxch_index = count;
		}
		else {	// different ID, find the right index for the properties vector
			prop_htxch_index = 0;
			while ((prop_htxch_index < n_htxch) && (inputs.at(count).HeatExchanger_ID != HeatExchangers.at(prop_htxch_index).ID_htx)) {
				// while indexes are different and there are heat exchangers to check
				prop_htxch_index++;
			}
		}
		

		// Calculation of m*Cp and heat exchange constants
		mCp_1 = inputs.at(count).mass_flow_1 * cp(HeatExchangers.at(prop_htxch_index).fluid_type_1, inputs.at(count).Temperature_inlet_1);
		mCp_2 = inputs.at(count).mass_flow_2 * cp(HeatExchangers.at(prop_htxch_index).fluid_type_2, inputs.at(count).Temperature_inlet_2);
		mCp_min = std::min(mCp_1, mCp_2);
		mCp_max = std::max(mCp_1, mCp_2);

		NTU = HeatExchangers.at(prop_htxch_index).UA / mCp_min;
		Cr = mCp_min / mCp_max;


		// Efficiency obtaining
		if (HeatExchangers.at(prop_htxch_index).flag_efficiency_table) {	// look up efficiency in a table
			efficiency = lookup_table_2dim(HeatExchangers.at(prop_htxch_index).n_rows_table, NTU, Cr, &HeatExchangers.at(prop_htxch_index).NTU, &HeatExchangers.at(prop_htxch_index).Cr, &HeatExchangers.at(prop_htxch_index).efficiency);
		}
		else {	// no table; obtain efficiency analiticaly considering the heat exchange characteristics
			
			switch (HeatExchangers.at(prop_htxch_index).exchanger_type) {
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
				if (HeatExchangers.at(prop_htxch_index).n_shell_passes > 1) {	// multiple passes through shell (for 1, the previous expression is obtained)
					efficiency = (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchangers.at(prop_htxch_index).n_shell_passes) - 1)
						/ (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchangers.at(prop_htxch_index).n_shell_passes) - Cr);
				}
				break;
			case 4:   // Crossed flows (simple pass, fluids not mixed)
				efficiency = 1 - exp(1 / Cr * pow(NTU, 0.22) * (exp(-Cr * pow(NTU, 0.78)) - 1));
				break;
			case 5:   // Crossed flows (simple pass, fluid 1 mixed, fluid 2 not mixed)
				if (mCp_1 == mCp_max) {
					efficiency = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))));
				}
				else {
					efficiency = 1 - exp(-1 / Cr * (1 - exp(-Cr * NTU)));
				}
				break;
			case 6:   // Crossed flows (simple pass, fluid 1 not mixed, fluid 2 mixed)
				if (mCp_2 == mCp_max) {
						  efficiency = 1 / Cr * (1 - exp(-Cr * (1 - exp(-NTU))));
				}
				else {
					efficiency = 1 - exp(-1 / Cr * (1 - exp(-Cr * NTU)));
				}
				break;
			case 7:   // Crossed flows (simple pass, both fluids mixed)
				efficiency = 1 / (1 / (1 - exp(-NTU)) + Cr / (1 - exp(-Cr * NTU)) - 1 / NTU);
				break;
			}
		}
		
		
		// Calculation of exhanged heat
		results.at(count).exchanged_heat_1 = efficiency * mCp_min * (inputs.at(count).Temperature_inlet_2 - inputs.at(count).Temperature_inlet_1); // (+) heat gained from 2 to 1; (-) heat lost from 1 to 2
		results.at(count).exchanged_heat_2 = -results.at(count).exchanged_heat_1; // (+) heat gained from 1 to 2; (-) heat lost from 2 to 1

		// Calculation of outlet temperatures
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
	double z = 0;

	// Data vectors
	std::vector<double> X_column;
	std::vector<double> Y_column;
	std::vector<double> Z_column;
	X_column=*x_column;
	Y_column=*y_column;
	Z_column=*z_column;

	// Distance vectors
	std::vector<double> dist_x;
	std::vector<double> dist_y;
	std::vector<double> dist_xy;
	dist_x.resize(n_rows);
	dist_y.resize(n_rows);
	dist_xy.resize(n_rows);


	
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


double HeatExchanger_model::cp(std::string fluid_type, double temperature)
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 4181.3;	// for tests: water [J/kg/K]
}


int main() {

	int n_htxs;	// number of heat exchangers
	int count;	// counter of executions

	n_htxs = 1;	// for test

	// Simulates model initialization
	HeatExchanger_model objeto;
	objeto.create_HeatExchangers(n_htxs);

	// Simulates calls of VEMOD to the model
	for (count = 0; count < n_htxs; count++ ){

		// Test inputs - remove when working
		objeto.inputs.at(count).HeatExchanger_ID = 1;
		objeto.inputs.at(count).mass_flow_1 = 0.01;
		objeto.inputs.at(count).mass_flow_2 = 0.5;
		objeto.inputs.at(count).Temperature_inlet_1 = 100;
		objeto.inputs.at(count).Temperature_inlet_2 = 20;

		objeto.calculate_heat_exchange();	// calculation call
}
	return 0;
}
