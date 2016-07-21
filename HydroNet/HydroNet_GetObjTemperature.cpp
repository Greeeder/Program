#include "stdafx.h"
#include <vector>

using namespace std;

double HydroNet_GetObjTemperature(double obj_pos,std::vector<double> branch_temp_pos,std::vector<double> branch_temperature) {
	
	double distance = -1000.0;
	double temperature=-500.0;
	int count=0;


	//Returns temperature at the specified position

	// Index of the temperature position closest to the specified positionand the difference of position between them

	if (obj_pos == 0.0)
		temperature = branch_temperature.at(0);
	else if (obj_pos == 100.0)
		temperature = branch_temperature.at(branch_temperature.size()-1);
	else {
		while (distance < 0.0) {
			distance = branch_temp_pos.at(count+1) - obj_pos;		// count+1 because the fert positin is checked before
			count++;
		}
		temperature = branch_temperature.at(count - 1);				// -1 because vector estarts at 0
	}
	return temperature;

	
}

int main()
{
	double obj_pos;
	double temperature;
	std::vector<double> branch_temp_pos;
	std::vector<double> branch_temperature;

	temperature = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos, branch_temperature);
	
    return 0;
}

