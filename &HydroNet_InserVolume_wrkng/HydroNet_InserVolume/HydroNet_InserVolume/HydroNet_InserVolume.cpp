// HydroNet_InserVolume.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>

double cp(std::string fluid_type, double temperature)
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 4181.3;	// for tests: water [J/kg/K]
}

double density(std::string fluid_type, double temperature)
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 1000;	// for tests: water [Kg/m3]
}




/*function [ positions, temperatures ] = HydroNet_InsertVolume( positions, temperatures, ...
    start_pos, end_pos, temp_action, fluid_type, value, volume )
% Obtains range of positions where a volume is to be placed inside a branch and replaces
% current volumes by the new volume. In addition, a temperature for the volume is calculated.
% Fluid_type is only mandatory for temp_action = 0 and 1 (must calculate outlet temperature).
% Value is only mandatory for temp_action = 1 and 2 (not just averaging).
% Volume is only mandatory for temp_action = 1 (impose heat).

*/

void HydroNet_InsertVolume(std::vector<double> positions, std::vector<double> temperatures, double start_pos, double end_pos, int temp_action, std::string fluid_type, double value, double volume) {


	int start_pos_ind;
	int end_pos_ind;
	int last_element;
	int count = 0;
	int count_pos;
	double temp_to_insert;
	double mass_start;
	double mass_end;
	double temp_mass_end;
	double mass_sum;
	double temp_mass_sum;
	double temp_mass_start;
	bool flag_start_exists = false;
	bool flag_end_exists = false;
	bool flag_inserted = false;
	std::vector<double> _temperatures;
	std::vector<double> _positions;

	_temperatures.resize(1);
	_positions.resize(1);

	

	// Obtains range of positions where a volume is to be placed inside a branch and replaces current volumes by the new volume.In addition, a temperature for the volume is calculated.
	//Fluid_type is only mandatory for temp_action = 0 and 1 (must calculate outlet temperature).
	// Value is only mandatory for temp_action = 1 and 2 (not just averaging).
	// Volume is only mandatory for temp_action = 1 (impose heat).

	if (start_pos == end_pos)
		return;					// to avoid dividing by zero


	// Finds start position inside positions vector

	while (count < positions.size() && start_pos < positions.at(count)) {

		count++;
	}
	if (start_pos == positions.at(count)) {      // Start position allready exixt in possitions
		flag_start_exists = true;
		start_pos_ind = count;
	}
	else
		start_pos_ind = count - 1;				// Position of the highest value in the vector position lower than start_position

	// Finds end position inside positions vector

	while (count < positions.size() && end_pos < positions.at(count)) {

		count++;
	}
	if (end_pos == positions.at(count)) {		// End position allready exixt in possitions
		flag_end_exists = true;
		end_pos_ind = count;
	}
	else
		end_pos_ind = count - 1;					// Position of the shortest value in the vector position lower than end_position

	//Obtains volume temperature

	if (temp_action == 2)						// impose temperature given by 'value'

		temp_to_insert = value;

	// 0 or 1: calculates mean temperature with an enthalpy balance
	// Branch volume will not be used now because it is not mandatory since
	// it is a constant that multiplies and divides everything

	// Temperatures at both ends of the volume, multiplied by density
	// and by the percentage range that they occupy inside the volume
	// Those volumes are taken incomplete because they are cut by the
	// new bounds start_pos and end_pos
	else {
		if (flag_start_exists == true) {					// start_pos already exists in the vector
			mass_start = 0;
			temp_mass_start = 0;
		}
		else if (start_pos_ind == end_pos_ind) {			// both are inside the same volume initially
			mass_start = density(fluid_type, temperatures.at(start_pos_ind - 1)) * (end_pos - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind - 1) * mass_start;
		}
		else {
			mass_start = density(fluid_type, temperatures.at(start_pos_ind - 1)) * (positions.at(start_pos_ind) - start_pos);
			// Not a real mass.Volume is not included because it appears both in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind - 1) * mass_start;
		}

		if (flag_end_exists == true) {							 // end_pos already exists in the vector
			mass_end = 0;
			temp_mass_end = 0;
		}

		else if (start_pos_ind == end_pos_ind) {				//both are inside the same volume initially
			mass_end = 0;
			temp_mass_end = 0;
		}
		else {
			mass_end = density(fluid_type, temperatures.at(end_pos_ind - 1)) * (end_pos - positions.at(end_pos_ind - 1));
			temp_mass_end = temperatures.at(end_pos_ind - 1) * mass_end;
		}

		mass_sum = mass_start + mass_end;						// denominator
		temp_mass_sum = temp_mass_start + temp_mass_end;	// numerator
	}

	// Volumes completely inside new volume

	if (flag_end_exists == true)
		last_element = (end_pos_ind - 1);

	else
		last_element = (end_pos_ind - 2);

	for (count_pos = start_pos_ind; count_pos <= last_element; count_pos++) {
		mass_sum = mass_sum + density(fluid_type, temperatures.at(count_pos)) * (positions.at(count_pos + 1) - positions.at(count_pos));
		temp_mass_sum = temp_mass_sum + temperatures.at(count_pos) * density(fluid_type, temperatures.at(count_pos)) * (positions.at(count_pos + 1) - positions.at(count_pos));
	}


	// Enter averaged temperature
	temp_to_insert = temp_mass_sum / mass_sum;

	if (temp_action == 1)										// calculate temperature imposing exchanged heat

		mass_sum = mass_sum * volume;								// volume needed to calculate real mass
	temp_to_insert = temp_to_insert + value / mass_sum / cp(fluid_type, temp_to_insert);




	// Obtains volume temperature

	if (temp_action == 2)												// impose temperature given by 'value'

		temp_to_insert = value;

	else {// 0 or 1: calculates mean temperature with an enthalpy balance
		 // Branch volume will not be used now because it is not mandatory since
		 // it is a constant that multiplies and divides everything

		 // Temperatures at both ends of the volume, multiplied by density
		 // and by the percentage range that they occupy inside the volume
		 // Those volumes are taken incomplete because they are cut by the
		 // new bounds start_pos and end_pos

		if (flag_start_exists == true) {								// start_pos already exists in the vector
			mass_start = 0;
			temp_mass_start = 0;
		}
		else if (start_pos_ind == end_pos_ind) {						 // both are inside the same volume initially
			mass_start = density(fluid_type, temperatures.at(start_pos_ind - 1)) * (end_pos - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind - 1) * mass_start;
		}

		else {
			mass_start = density(fluid_type, temperatures.at(start_pos_ind - 1)) * (positions.at(start_pos_ind) - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind - 1) * mass_start;
		}


		if (flag_end_exists == true) {													// end_pos already exists in the vector
			mass_end = 0;
			temp_mass_end = 0;
		}
		else if (start_pos_ind == end_pos_ind) {										 // both are inside the same volume initially
			mass_end = 0;
			temp_mass_end = 0;
		}
		else {
			mass_end = density(fluid_type, temperatures.at(end_pos_ind - 1)) * (end_pos - positions.at(end_pos_ind - 1));
			temp_mass_end = temperatures.at(end_pos_ind - 1) * mass_end;
		}


		mass_sum = mass_start + mass_end;													// denominator
		temp_mass_sum = temp_mass_start + temp_mass_end;										// numerator


		 // Volumes completely inside new volume

		if (flag_end_exists == true)
			last_element = (end_pos_ind - 1);
		else
			last_element = (end_pos_ind - 2);


		for (count_pos = start_pos_ind; count_pos <= last_element; count_pos++) {
			mass_sum = mass_sum + density(fluid_type, temperatures.at(count_pos)) * (positions.at(count_pos + 1) - positions.at(count_pos));
			temp_mass_sum = temp_mass_sum + temperatures.at(count_pos) * density(fluid_type, temperatures.at(count_pos)) * (positions.at(count_pos + 1) - positions.at(count_pos));
		}


		// Enter averaged temperature
		temp_to_insert = temp_mass_sum / mass_sum;

		if (temp_action == 1) {	       								// calculate temperature imposing exchanged heat

			mass_sum = mass_sum * volume;							// volume needed to calculate real mass
			temp_to_insert = temp_to_insert + value / mass_sum / cp(fluid_type, temp_to_insert);
		}
	}
	// Inserts calculated temperature and removes intermediate temperatures in the branch
	for (int count = 1, aux = 1; flag_start_exists == false && flag_end_exists == false && count < temperatures.size() + 2
		|| flag_start_exists == true && flag_end_exists == false && count < temperatures.size() + 1
		|| flag_start_exists == false && flag_end_exists == true && count < temperatures.size() + 1
		|| flag_start_exists == true && flag_end_exists == true && count < temperatures.size(); count++) {

		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count--;
			aux--;
		}
		if (positions.at(count) < start_pos) {										// Saves values of temperatures until the new temperature
			_temperatures.push_back(temperatures.at(aux));
			aux++;
		}

		else if (flag_inserted == true ) {											// Saves values of temperatures after the new temperature


			_temperatures.push_back(temperatures.at(aux));
			aux++;
		}

		else {																//Saves value of new temperature betwin two previusly known
			_temperatures.push_back(temp_to_insert);
			flag_inserted = true;
			if (flag_start_exists == true && flag_end_exists == true)
				aux++;

			if (flag_start_exists == false && flag_end_exists == false)
				_temperatures.push_back(temperatures.at(aux - 1));;
		}
		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count++;
			aux++;
		}

	}
	temperatures.clear();

	temperatures = _temperatures;


	//Inserts specified positions and removes intermediate positions in the branch

	flag_inserted = false;

	for (int count = 1, aux = 1; flag_start_exists == false && flag_end_exists == false && count < positions.size() + 2
		|| flag_start_exists == true && flag_end_exists == false && count < positions.size()
		|| flag_start_exists == false && flag_end_exists == true && count < positions.size()
		|| flag_start_exists == true && flag_end_exists == true && count < positions.size(); count++) {

		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count--;

		}

		if (flag_inserted != true && positions.at(count) < start_pos) {										// Saves values of positions until the new temperature
			_positions.push_back(positions.at(aux));
			aux++;
		}

		else if (flag_inserted == true) {											// Saves values of positions after the new temperature


			_positions.push_back(positions.at(aux));
			aux++;
		}

		else {																//Saves value of new temperature betwin two previusly known


			_positions.push_back(start_pos);
			_positions.push_back(end_pos);
			flag_inserted = true;
			if (flag_start_exists == true)
				aux++;

			if (flag_end_exists == true)
				aux++;

			if (flag_start_exists == true && flag_end_exists == true)
				count++;
		}
		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count++;

		}

	}
	positions.clear();

	positions = _positions;

	return void(temperatures);
}

int main()
{
	std::vector<double> positions;
	std::vector<double> temperatures;
	double start_pos;
	double end_pos;
	int temp_action;
	std::string fluid_type;
	double value;
	double volume;


	HydroNet_InsertVolume(positions, temperatures, start_pos, end_pos, temp_action, fluid_type, value, volume);
	temperatures;

    return 0;
}

