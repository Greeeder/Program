#include <cmath>		// for sqrt
#include <fstream>		// to read file
#include <vector>
#include <string>		// string use; necessary for std::stoi
#include <algorithm>	// to avoid 'min/max not a member of std' error
#include "HeatExchanger.h"


void HeatExchanger_model::create_HeatExchanger()
{
	int count;	// counter of heat exchangers
	int lines;	// number of lines in the efficiency table
	int line;	// counter for efficiency table data
	std::string line_scan;	// content of a line as a string

	// Names of data files
	std::string HeatExchanger_data_archives;	// properties
	std::string HeatExchanger_table_archives;	// efficiency table


	// For testing; only one file/heat exchanger
	HeatExchanger_data_archives = "HeatExchanger_inputs.txt";				
	HeatExchanger_table_archives = "HeatExchanger_efficiency_table.txt";
	

		
	// Heat exchanger loop
	

		std::vector<std::string> HeatExchanger_data;	// vector to store file data as strings
		std::ifstream HeatExchanger_file(HeatExchanger_data_archives, std::ifstream::out);


		// Read data
		while (!HeatExchanger_file.eof()) {
			std::getline(HeatExchanger_file, line_scan);	// read a line of HeatExchangers_file
			HeatExchanger_data.push_back(line_scan);		// store the new line at the end of HeatExchanger_data
		}
		HeatExchanger_file.close();
		
		
		// String data goes directly to structure
		name = HeatExchanger_data[2];
		fluid_type_1 = HeatExchanger_data[3];
		fluid_type_2 = HeatExchanger_data[4];

		// Integer data is converted and stored
		ID_htx = std::stoi(HeatExchanger_data[0]);
		ID_htx_HydroNet = std::stoi(HeatExchanger_data[1]);
		exchanger_type = std::stoi(HeatExchanger_data[8]);
		n_shell_passes = std::stoi(HeatExchanger_data[9]);

		// Double data is converted and stored
		hydraulic_resistance = std::stod(HeatExchanger_data[5]);
		head_loss = std::stod(HeatExchanger_data[6]);
		UA = std::stod(HeatExchanger_data[7]);

		// Boolean data is converted and stored
		if (HeatExchanger_data[10] == "1") {
			flag_efficiency_table = true;
		}
		else if (HeatExchanger_data[10] == "0") {
			flag_efficiency_table = false;
		}


		// Efficiency table data
		if (flag_efficiency_table == true){	// efficiency must be obtained from a table

			std::vector<std::string> efficiency_data;
			std::vector<std::string> NTU_data;
			std::vector<std::string> Cr_data;


			// Read data
			HeatExchanger_file.open(HeatExchanger_table_archives, std::ifstream::out);
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
			NTU.resize(lines);
			Cr.resize(lines);
			efficiency.resize(lines);


			// Table data is converted and stored
			line = 0;
			while (line < lines) {
				NTU[line] = std::stod(NTU_data[line]);
				Cr[line] = std::stod(Cr_data[line]);
				efficiency[line] = std::stod(efficiency_data[line]);
				line++;
			}
			n_rows_table = lines;
			
	

		// Store heat exchanger ID in results vector
		HeatExchanger_ID =ID_htx;
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
	double _NTU;				// Parameter NTU of heat ezchanger
	double _Cr;				// Ratio of m*Cp
	double _efficiency = 0;	// Heat exchanger efficiency
	
	
	// Heat exchanger loop
	

		
		
		

		// Calculation of m*Cp and heat exchange constants
		mCp_1 = mass_flow_1 * cp(fluid_type_1, Temperature_inlet_1);
		mCp_2 = mass_flow_2 * cp(fluid_type_2, Temperature_inlet_2);
		mCp_min = std::min(mCp_1, mCp_2);
		mCp_max = std::max(mCp_1, mCp_2);

		_NTU = UA / mCp_min;
		_Cr = mCp_min / mCp_max;


		// Efficiency obtaining
		if (flag_efficiency_table) {	// look up efficiency in a table
			_efficiency = lookup_table_2dim(n_rows_table, _NTU, _Cr, &NTU, &Cr, &efficiency);
		}
		else {	// no table; obtain efficiency analiticaly considering the heat exchange characteristics
			
			switch (exchanger_type) {
			case 1:	// Simple co-current
				_efficiency = (1 - exp(-_NTU * (1 + _Cr))) / (1 + _Cr);
				break;
			case 2:	// Simple counter-current
				if (_Cr == 1) {
					_efficiency = 1 / (1 + 1 / _NTU);
				}
				else {
					_efficiency = (1 - exp(-_NTU * (1 - _Cr))) / (1 - _Cr * exp(-_NTU * (1 - _Cr)));
				}
				break;
			case 3:	// Shell and tube
				_efficiency = 2 /  (1 + _Cr + sqrt(1 + _Cr * _Cr) * (1 + exp(-_NTU * sqrt(1 + _Cr * _Cr)) / (1 - exp(-_NTU * sqrt(1 + _Cr * _Cr)))));
				if (n_shell_passes > 1) {	// multiple passes through shell (for 1, the previous expression is obtained)
					_efficiency = (pow(((1 - _efficiency * _Cr) / (1 - _efficiency)), n_shell_passes) - 1)
						/ (pow(((1 - _efficiency * _Cr) / (1 - _efficiency)), n_shell_passes) - _Cr);
				}
				break;
			case 4:   // Crossed flows (simple pass, fluids not mixed)
				_efficiency = 1 - exp(1 / _Cr * pow(_NTU, 0.22) * (exp(-_Cr * pow(_NTU, 0.78)) - 1));
				break;
			case 5:   // Crossed flows (simple pass, fluid 1 mixed, fluid 2 not mixed)
				if (mCp_1 == mCp_max) {
					_efficiency = 1 / _Cr * (1 - exp(-_Cr * (1 - exp(-_NTU))));
				}
				else {
					_efficiency = 1 - exp(-1 / _Cr * (1 - exp(-_Cr * _NTU)));
				}
				break;
			case 6:   // Crossed flows (simple pass, fluid 1 not mixed, fluid 2 mixed)
				if (mCp_2 == mCp_max) {
						  _efficiency = 1 / _Cr * (1 - exp(-_Cr * (1 - exp(-_NTU))));
				}
				else {
					_efficiency = 1 - exp(-1 / _Cr * (1 - exp(-_Cr * _NTU)));
				}
				break;
			case 7:   // Crossed flows (simple pass, both fluids mixed)
				_efficiency = 1 / (1 / (1 - exp(-_NTU)) + _Cr / (1 - exp(-_Cr * _NTU)) - 1 / _NTU);
				break;
			}
		}
		
		
		// Calculation of exhanged heat
		exchanged_heat_1 = _efficiency * mCp_min * (Temperature_inlet_2 - Temperature_inlet_1); // (+) heat gained from 2 to 1; (-) heat lost from 1 to 2
		exchanged_heat_2 = -exchanged_heat_1; // (+) heat gained from 1 to 2; (-) heat lost from 2 to 1

		// Calculation of outlet temperatures
		Temperature_outlet_1 = Temperature_inlet_1 + exchanged_heat_1 / mCp_1;
		Temperature_outlet_2 = Temperature_inlet_2 + exchanged_heat_2 / mCp_2;
	

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

	
	int count;	// counter of executions

	

	// Simulates model initialization
	HeatExchanger_model objeto;
	objeto.create_HeatExchanger();

	// Simulates calls of VEMOD to the model
	

		// Test inputs - remove when working
		objeto.HeatExchanger_ID = 1;
		objeto.mass_flow_1 = 0.01;
		objeto.mass_flow_2 = 0.5;
		objeto.Temperature_inlet_1 = 100;
		objeto.Temperature_inlet_2 = 20;

		objeto.calculate_heat_exchange();	// calculation call

	return 0;
}
