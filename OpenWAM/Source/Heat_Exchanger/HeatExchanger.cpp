#include <cmath>		// for sqrt
#include <fstream>		// to read file
#include <vector>
#include <string>		// string use; necessary for std::stod
#include <algorithm>	// to avoid 'min/max not a member of std' error
#include "HeatExchanger.h"


HeatExchanger_model::HeatExchanger_model() {

}


void HeatExchanger_model::calculate_heat_exchange()
{

	double mCp_1;			// Mass flow times heat capacity of the first fluid
	double mCp_2;			// Mass flow times heat capacity of the second fluid
	double mCp_min;			// Minimum mass flow times heat capacity of both fluids
	double mCp_max;			// Maximum mass flow times heat capacity of both fluids
	double NTU;				// Parameter NTU of heat exchanger
	double Cr;				// Ratio of m*Cp
	
	
	if ((mass_flow_1 == 0.0) || (mass_flow_2 == 0.0)) {	// any flow is zero
		// Efficiency method makes no sense

		// Calculation of exhanged heat
		exchanged_heat_1 = 0.0;
		exchanged_heat_2 = 0.0;

		// Calculation of outlet temperatures
		Temperature_outlet_1 = Temperature_inlet_1;
		Temperature_outlet_2 = Temperature_inlet_2;
	}

	else {

		// Calculation of m*Cp and heat exchange constants

		mCp_1 = mass_flow_1 * fluid_type_1->FunCp(Temperature_inlet_1);
		mCp_2 = mass_flow_2 * fluid_type_2->FunCp(Temperature_inlet_2);
		mCp_min = std::min(mCp_1, mCp_2);
		mCp_max = std::max(mCp_1, mCp_2);

		NTU = UA / mCp_min;
		Cr = mCp_min / mCp_max;


		// Obtain efficiency
		if (flag_efficiency_table) {	// look up efficiency in a table
			efficiency = lookup_table_2dim(n_rows_table, NTU, Cr, &NTU_col, &Cr_col, &efficiency_col);
		}
		else {	// no table; obtain efficiency analytically considering the heat exchange characteristics

			switch (exchanger_type) {
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
				efficiency = 2 / (1 + Cr + sqrt(1 + Cr * Cr) * (1 + exp(-NTU * sqrt(1 + Cr * Cr)) / (1 - exp(-NTU * sqrt(1 + Cr * Cr)))));
				if (n_shell_passes > 1) {	// multiple passes through shell (for 1, the previous expression is obtained)
					efficiency = (pow(((1 - efficiency * Cr) / (1 - efficiency)), n_shell_passes) - 1)
						/ (pow(((1 - efficiency * Cr) / (1 - efficiency)), n_shell_passes) - Cr);
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
		exchanged_heat_1 = efficiency * mCp_min * (Temperature_inlet_2 - Temperature_inlet_1); // (+) heat gained from 2 to 1; (-) heat lost from 1 to 2
		exchanged_heat_2 = -exchanged_heat_1; // (+) heat gained from 1 to 2; (-) heat lost from 2 to 1

		heat_1 = exchanged_heat_1 * dt;
		heat_2 = exchanged_heat_2 * dt;

		// Calculation of outlet temperatures
		Temperature_outlet_1 = Temperature_inlet_1 + exchanged_heat_1 / mCp_1;
		Temperature_outlet_2 = Temperature_inlet_2 + exchanged_heat_2 / mCp_2;
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


//double HeatExchanger_model::cp(std::string fluid_type, double temperature)
//{
//	// dummy fuction to be ignored when the model is integrated in VEMOD
//
//	return 4181.3;	// for tests: water [J/kg/K]
//}
