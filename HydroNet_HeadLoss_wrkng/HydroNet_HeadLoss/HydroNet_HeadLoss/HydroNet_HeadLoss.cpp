#include <vector>
#include <string>		
#define PI 3.14159265358979323846  


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
	int outlet_oject;					// Thermostat's outlet object identifier
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
	int outlet_oject;				// Pump's outlet object identifier
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
	int outlet_oject;				// Pump's outlet object identifier
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
	double heat;					// Heat Exchanged
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};

struct strHeat_exch_Tout
{
	int heat_exch_Tout;				// Identifier
	std::string name;				// Heat Exchanger name
	int inlet_object;				// Heat Exchanger's inlet object identifier
	int outlet_object;				// Heat Exchanger's outlet object identifier
	double T_out;					// Outlet temperaure [K]
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};
/*struct strTank
{
	int tank;					// Identifier
	std::string name;				// Tank's name	
	int inlet_object;				// Tank's inlet object identifier
	int outlet_oject;				// Tank's outlet object identifier
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
};
struct strBranches
{
	int branch_ind;				// branch idintificator
	int branches_cycle;
	//int n_branch;				// Number of branches
	std::vector<strObjects> objects;	// Objects in Brnach
};

double HydroNet_GetObjTemperature(double obj_pos, double branch_temp_pos, double branch_temperature) {
	double z;
	return z;
}

void HydroNet_HeadLoss(int n_branch, std::vector<strPipe> pipe, std::vector<strValve_fix> valve_fix, std::vector<strValve_var> valve_var, std::vector<strThermostat> thermostat, std::vector<strPump_turbo> pump_turbo, std::vector<strPump_volum> pump_volum, std::vector<strHeat_exch_fix> heat_exch_fix, std::vector<strHeat_exch_Tout> heat_exch_Tout, std::vector<strBranches> branches, int temperature/* std::vector<strTank> tank,*/, std::vector<double> obj_inlet_pos, std::vector<double> obj_outlet_pos, std::vector<double> branch_temp_pos, std::vector<double> branch_temperature) {

	std::string _class;
	int count;
	int aux = 0;
	double opening;
	double temperature_aux;
	double obj_pos;
	int _ID;
	int index;
	std::vector <double> head_loss;
	std::vector <double> hydr_resist1;
	std::vector <double> hydr_resist2;

	head_loss.resize(n_branch + 1);							// head loss(m.c.f.) - independent of flow

	hydr_resist1.resize(n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 3 * s)
																			// proportional to flow - for turbopumps
	hydr_resist2.resize(n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 6 * s ^ 2)
																			//proportional to flow ^ 2


	// Ads the values of head_loss, hydr_resist1 and hydr_resist2 for each branch
for (count=1; count<= n_branch; count++){

	head_loss.at(count) = 0.0;
	hydr_resist1.at(count) = 0.0;
	hydr_resist2.at(count) = 0.0;
	opening = 0.0;
	aux = 1;

	while (aux<= branches.at(count).objects.size())
	{
		_class = branches.at(count).objects.at(aux).Class;
		_ID = branches.at(count).objects.at(aux).ID;

	}

	if (_class == "pipe") {

		while (index<pipe.size() && _ID!=pipe.at(index).pipe_id){

			index++;
		}


		hydr_resist2.at(count) = hydr_resist2.at(count) + 8 * pipe.at(index).friction_coef * pipe.at(index).length / (PI*PI / 9.81 / pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter);
	}

	else if (_class == "valve_fix") {

		while (index<valve_fix.size() && _ID != valve_fix.at(index).valve_fix_id) {

			index++;
		}

		head_loss.at(count) = head_loss.at(count) + valve_fix.at(index).head_loss;
	}

	else if (_class == "valve_var") {

		while (index < valve_var.size() && _ID != valve_var.at(index).valve_var_id) {

			index++;
		}

		head_loss.at(count) = head_loss.at(count) + valve_var.at(index).coef0 + valve_var.at(index).coef1 * valve_var.at(index).opening + valve_var.at(index).coef2 * valve_var.at(index).opening*valve_var.at(index).opening + valve_var.at(index).coef3 * valve_var.at(index).opening*valve_var.at(index).opening*valve_var.at(index).opening;
	}
	else if (_class == "thermostat") {

		while (index<thermostat.size() && _ID != thermostat.at(index).themostat_id) {

			index++;;
		}
		obj_pos = (obj_inlet_pos.at( index) +obj_outlet_pos.at( index )) / 2;																	//position inside the branch as a % of the branch volume
		
		temperature_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos.at( count ),branch_temperature.at(count ));						// temperature in the thermostat

		opening = opening + thermostat.at(index).coef_T_0 + thermostat.at(index).coef_T_1 * temperature_aux + thermostat.at(index).coef_T_2 * temperature_aux*temperature_aux + thermostat.at(index).coef_T_3 * temperature_aux*temperature_aux*temperature_aux;

		head_loss.at(count) = head_loss.at(count) + thermostat.at(index).coef_h_0 + thermostat.at(index).coef_h_1 * opening + thermostat.at(index).coef_h_2 * opening*opening + thermostat.at(index).coef_h_3 * opening*opening*opening;

	}

	else if (_class == "pump_turbo") {

		while (index<pump_turbo.size() && _ID != pump_turbo.at(index).Pump_turbo) {

			index++;
		}


		head_loss.at(count) = head_loss.at(count) -( pump_turbo.at(index).coef_N0 + pump_turbo.at(index).coef_N1 * pump_turbo.at(index).pump_speed + pump_turbo.at(index).coef_N2 * pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

		hydr_resist1.at(count) = hydr_resist1.at(count) -( pump_turbo.at(index).Q_coef_N0 + pump_turbo.at(index).Q_coef_N1* pump_turbo.at(index).pump_speed + pump_turbo.at(index).Q_coef_N2* pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

		hydr_resist2.at(count) = hydr_resist2.at(count) -( pump_turbo.at(index).Q2_coef_N0 + pump_turbo.at(index).Q2_coef_N1* pump_turbo.at(index).pump_speed + pump_turbo.at(index).Q2_coef_N2* pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

	}

	else if (_class == "pump_volum") {

		while (index < pump_volum.size() && _ID != pump_volum.at(index).Pump_volum) {

			index++;
		}
		hydr_resist1.at(count) = hydr_resist1.at(count) - (pump_volum.at(index).coef0 + pump_volum.at(index).coef1* pump_volum.at(index).pump_speed + pump_volum.at(index).coef2* pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed + pump_volum.at(index).coef3* pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed);
	}

	else if (_class == "heat_exch_fix") {

		while (index < heat_exch_fix.size() && _ID != heat_exch_fix.at(index).heat_exch_fix) {

			index++;
		}
		hydr_resist2.at(count) = hydr_resist2.at(count) + heat_exch_fix.at(index).hydr_resist;
	}

	else if (_class == "heat_exch_Tout") {

		while (index < heat_exch_Tout.size() && _ID != heat_exch_Tout.at(index).heat_exch_Tout) {

			index++;
		}
		hydr_resist2.at(count) = hydr_resist2.at(count) + heat_exch_Tout.at(index).hydr_resist;
	}

	/*
	else if (_class == "tank"){

	while (index < valve_var.size() && _ID != valve_var.at(index).pipe_id) {

			index++
		}

		head_loss.at(count) = head_loss.at(count) + tank.at(index).head_loss;
		}
	*/
	aux++;
	}
}


int main()
{




	//std::vector<strObjects> objects;      //Objects in each brnach
	std::vector<strPipe> pipe;
	std::vector<strValve_fix> valve_fix;
	std::vector<strValve_var> valve_var;
	std::vector<strThermostat> thermostat;
	std::vector<strPump_turbo> pump_turbo;
	std::vector<strPump_volum> pump_volum;
	std::vector<strHeat_exch_fix> heat_exch_fix;
	std::vector<strHeat_exch_Tout> heat_exch_Tout;
	//std::vector<strTank> tank;
	std::vector<double> obj_inlet_pos;
	std::vector<double> obj_outlet_pos;
	std::vector<double> branch_temp_pos;
	std::vector<double> branch_temperature;
	
	std::vector<strBranches> branches;			// Branches vector
	int temperature;
	int n_branch;								// Number of branches

	HydroNet_HeadLoss(n_branch, pipe, valve_fix, valve_var, thermostat, pump_turbo, pump_volum, heat_exch_fix, heat_exch_Tout,branches, temperature/* tank,*/, obj_inlet_pos, obj_outlet_pos, branch_temp_pos, branch_temperature);

	return 0;
}

