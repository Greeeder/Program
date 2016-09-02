#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>		// std::stoi
#include <algorithm>	// to avoid 'min/max not a member of std' error
#include "HeatExchanger.h"


void HeatExchanger_model::create_HeatExchangers()
{
	int line;
	int count;
	std::ifstream HeatExchanger_file;

	// for (every file / heat exchanger) {

	// Strings, integers, doubles and booleans
	// Read lines of data
	HeatExchanger_file.open("HeatExchanger_inputs.txt", std::ios::out);
	std::string HeatExchanger_data[11];
	line = 0;
	for (count = 0; count <= 11; count++) {
		std::getline(HeatExchanger_file, HeatExchanger_data[line]);	// Save line in HeatExchanger_data
		line++;
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
	else if(HeatExchanger_data[10] == "0") {
		HeatExchangers.at(count).flag_efficiency_table = false;
	}


	// Efficiency table: arrays
	// Read all data in the table
	int block;
	int lines;

	std::string line_scan;
	
	HeatExchanger_file.open("HeatExchanger_efficiency_table.txt", std::ios::out);
	
	line = 0;
	
	while (!HeatExchanger_file.eof()) {
		getline(HeatExchanger_file,line_scan);
		lines++;	// lines counter
	}
	
	std::string efficiency_data[line];
	std::string NTU_data[line];
	std::string Cr_data[lines];
	
	block = 0;
	for (count = 0; count < lines; count++) {
		getline(HeatExchanger_file, NTU_data[block],'\t');
		
		getline(HeatExchanger_file, Cr_data[block],'\t');
		
		getline(HeatExchanger_file, efficiency_data[block],'\n');
	
		block++;
	}
	HeatExchanger_file.close();
	
	
	

	// Array data is converted and stored
	HeatExchangers.at(count).NTU.resize(lines);
	HeatExchangers.at(count).Cr.resize(lines);
	HeatExchangers.at(count).efficiency.resize(lines);

	
	line = 0;
	while (line<lines) {
		HeatExchangers.at(count).NTU[line] = std::stod(NTU_data[line]);
	
		HeatExchangers.at(count).Cr[line] = std::stod(Cr_data[line]);
	
		HeatExchangers.at(count).efficiency[line] = std::stod(efficiency_data[line]);
	
		line++;
	}

	// Rows in the table: NTU, Cr & efficicency
	HeatExchangers.at(count).n_rows_table = line;

	results.at(count).HeatExchanger_ID = HeatExchangers.at(count).ID_htx;

	//}
}


void HeatExchange_refresh (int HeatExchanger_ID, std::vector< stHeatExchanger > & HeatExchangers){

	int count;

	count = 0;

	while (HeatExchanger_ID != HeatExchangers.at(count).ID_htx){
	count++;	
	}
	
}


void HeatExchanger_model::calculate_heat_exchange(int HeatExchanger_ID,
	double Temperature_inlet_1, double Temperature_inlet_2,	double mass_flow_1, double mass_flow_2)
{
	int count;
	int index_aux;
	double mCp_1;
	double mCp_2;
	double mCp_min;
	double mCp_max;
	double NTU;
	double Cr;
	double efficiency;

	for (count)
	index_aux = 

	mCp_1 = mass_flow_1 * cp(HeatExchanger.fluid_type_1, Temperature_inlet_1);
	mCp_2 = mass_flow_2 * cp(HeatExchanger.fluid_type_2, Temperature_inlet_2);
	mCp_min = std::min(mCp_1, mCp_2);
	mCp_max = std::max(mCp_1, mCp_2);

	NTU = HeatExchanger.UA / mCp_min;
	Cr = mCp_min / mCp_max;

	if (HeatExchanger.flag_efficiency_table) {	// look up efficiency in a table
		efficiency = lookup_table_2dim(HeatExchanger.n_rows_table, NTU, Cr, HeatExchanger.NTU, HeatExchanger.Cr, HeatExchanger.efficiency);
	}
	else {	// no table; obtain efficiency analiticaly considering the heat exchanger characteristics
		switch (HeatExchanger.exchanger_type) {
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
		case3:	// Shell and tube
			efficiency = 1 / (2 * (1 + Cr + sqrt(1 + Cr * Cr) * (1 + exp(-NTU * sqrt(1 + Cr * Cr)) / (1 - exp(-NTU * sqrt(1 + Cr * Cr))))));
			if (HeatExchanger.n_shell_passes > 1) {	// multiple passes through shell (for 1, the previous expression is obtained)
				efficiency = (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchanger.n_shell_passes) - 1)
					/ (pow(((1 - efficiency * Cr) / (1 - efficiency)), HeatExchanger.n_shell_passes) - Cr);
			}
			break;
		case4: // Crossed flows (simple pass, fluids not mixed)
			efficiency = 1 - exp(1 / Cr * pow(NTU, 0.22) * (exp(-Cr * pow(NTU, 0.78)) - 1));
			break;
		}
	}

	results.at(index_aux).exchanged_heat_1 = efficiency * mCp_min * (Temperature_inlet_2 - Temperature_inlet_1); // (+) heat gained from 2 to 1; (-) heat lost from 1 to 2
	results.at(index_aux).exchanged_heat_2 = -exchanged_heat1; // (+) heat gained from 1 to 2; (-) heat lost from 2 to 1

	results.at(index_aux).Temperature_outlet_1 = Temperature_inlet_1 + exchanged_heat1 / mCp_1;
	results.at(index_aux).Temperature_outlet_2 = Temperature_inlet_2 + exchanged_heat2 / mCp_2;

}

double lookup_table_2dim(int n_rows, double x, double y, double x_column[], double y_column[], double z_column[])
{

	int count;
	int comb1;
	int comb2;
	int comb3;
	int comb4;
	int smallest;
	double dist_x[n_rows];
	double dist_y[n_rows];
	double dist_xy[n_rows];
	double z;

	// If x and y are outside the table, bring them inside the limits
	if (x > std::max_element(x_column)) {
		x = std::max_element(x_column);
	}
	else if (x < std::min_element(x_column)) {
		x = std::min_element(x_column);
	}
	if (y > std::max_element(y_column)) {
		y = std::max_element(y_column);
	}
	else if (y < std::min_element(y_column)) {
		y = std::min_element(y_column);
	}

	// Obtain distances to x and y
	for (count = 0; count <= n_rows; count++) {
		dist_x[count] = abs(x - x_column[count]);
	}
	for (count = 0; count <= n_rows; count++) {
		dist_y[count] = abs(y - y_column[count]);
	}

	// Calculate distances to x y combinations (rows) from the queried values
	// Hypotenuse from cathetus
	for (count = 0; count <= n_rows; count++) {
		dist_xy[count] = sqrt(dist_x[count] * dist_x[count] + dist_y[count] * dist_y[count]);
	}

	// Find the 4 closest points of combination x y to the queried combination
	if (n_rows > 3) {
		smallest = 0;
		for (count = 1; count <= n_rows; count++) {
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
		for (count = smallest + 1; count <= n_rows; count++) {
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
		for (count = smallest + 1; count <= n_rows; count++) {
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
		for (count = smallest + 1; count <= n_rows; count++) {
			if ((dist_xy[count] < dist_xy[smallest]) && (count != comb1) && (count != comb2) && (count != comb3)) {
				smallest = count;
			}
		}
		comb4 = smallest;	// index of the fourth closest point


		// Interpolate with distances to obtain the queried value
		z = (z_column[comb1] / dist_xy[comb1] + z_column[comb2] / dist_xy[comb2] + z_column[comb3] / dist_xy[comb3] + z_column[comb4] / dist_xy[comb4])
			/ (1 / dist_xy[comb1] + 1 / dist_xy[comb2] + 1 / dist_xy[comb3] + 1 / dist_xy[comb4]);
	}
	else if (n_rows == 3) {
		z = (z_column[0] / dist_xy[0] + z_column[1] / dist_xy[1] + z_column[2] / dist_xy[2]) / (1 / dist_xy[0] + 1 / dist_xy[1] + 1 / dist_xy[2]);
	}
	else if (n_rows == 2) {
		z = (z_column[0] / dist_xy[0] + z_column[1] / dist_xy[1]) / (1 / dist_xy[0] + 1 / dist_xy[1]);
	}
	else if  (n_rows == 1) {
		z = z_column[0];
	}

	return z;
}

double cp(std::string fluid_type, double Temperature_inlet)
{
	return 4181.3;	// for tests: water [J/kg/K]
}
