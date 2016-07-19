#include <cmath>
//#include <iostream>
//#include <fstream>
#include <vector>
#include <string>		
#include <algorithm>	
//#include "HeatExchanger.h"

using namespace std;

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
	int Pump_turbo;					// Identifier
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

struct strTank
{
	int tank;					// Identifier
	std::string name;				// Tank's name	
	int inlet_object;				// Tank's inlet object identifier
	int outlet_oject;				// Tank's outlet object identifier
	double head_loss;				// Head loss [m.c.f. / m3/s]
};

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
	int n_branch;				// Number of branches
};

void HydroNet_HEadLoss(int id, std::vector<strObjects> objects, std::vector<strPipe> pipe, std::vector<strValve_fix> valve_fix, std::vector<strValve_var> valve_var, std::vector<strThermostat> thermostat, std::vector<strPump_turbo> pump_turbo, std::vector<strPump_volum> pump_volum, std::vector<strHeat_exch_fix> heat_exch_fix, std::vector<strHeat_exch_Tout> heat_exch_Tout, std::vector<strTank> tank, std::vector<strBranches> branches, int temperature) {

	std::string _class;
	int count;
	int aux;
	double opening;
	std::vector <double> head_loss;
	std::vector <double> hydr_resist1;
	std::vector <double> hydr_resist2;

	head_loss.resize(branches.at(0).n_branch + 1);							// head loss(m.c.f.) - independent of flow

	hydr_resist1.resize(branches.at(0).n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 3 * s)
																			// proportional to flow - for turbopumps
	hydr_resist2.resize(branches.at(0).n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 6 * s ^ 2)
																			//proportional to flow ^ 2

for (count=1; count<= n_branch; count++){
	while (aux< objects.size() && id != objects.at(aux).ID)
	{
		_class = objects.at(aux).Class;
	}

	if (_class == "pipe") 	
		hydr_resist2.at(count) = hydr_resist2.at(count) + 8 * pipe.at(aux).friction_coef * pipe.at(aux).length / (pi ^ 2) / 9.81 / (pipe.at(aux).diameter ^ 5);

	else if (_class == "valve_fix")		
		head_loss.at(count) = head_loss.at(count) + valve_fix.at(aux).head_loss;
		
	else if (_class == "valve_var")
		head_loss.at(count) = head_loss.at(count) + valve_var.at(aux).coef0 + valve_var.at(aux).coef1 * valve_var.at(aux).opening + valve_var.at(aux).coef2 * valve_var.at(aux).opening ^ 2 + valve_var.at(aux).coef3 * valve_var.at(aux).opening^3;
	
	else if (_class == "thermostat"){
	
	
		opening = opening + thermostat.at(aux).coef_T_0 + thermostat.at(aux).coef_T_1 * temperature + thermostat.at(aux).coef_T_2 * temperature^2 + thermostat.at(aux).coef_T_3 * temperature^3;
		
		head_loss.at(count) = head_loss.at(count) + thermostat.at(aux).coef_h_0 + thermostat.at(aux).coef_h_1 * opening + thermostat.at(aux).coef_h_2 * opening^2 + thermostat.at(aux).coef_h_3 * opening^3;
		
	}
	
	else if (_class == "pump_turbo"){
	
	
		head_loss.at(count) = head_loss.at(count) + pump_turbo.at(aux).coef_N0 + pump_turbo.at(aux).coef_N1 * pump_turbo.at(aux).pump_speed + pump_turbo.at(aux).coef_N2 * pump_turbo.at(aux).pump_speed^2 ;
		
		hydr_resist1.at(count) = hydr_resist1.at(count) +pump_turbo.at(aux).Q_coef_N0 + pump_turbo.at(aux).Q_coef_N1* pump_turbo.at(aux).pump_speed + pump_turbo.at(aux).Q_coef_N2* pump_turbo.at(aux).pump_speed^2;
		
		hydr_resist2.at(count) = hydr_resist2.at(count) +pump_turbo.at(aux).Q2_coef_N0 + pump_turbo.at(aux).Q2_coef_N1* pump_turbo.at(aux).pump_speed + pump_turbo.at(aux).Q2_coef_N2* pump_turbo.at(aux).pump_speed^2;

	}
	
	else if (_class == "pump_volum")
		hydr_resist1.at(count) = hydr_resist1.at(count) + pump_volum.at(aux).coef0 + pump_volum.at(aux).coef1* pump_volum.at(aux).pump_speed + pump_volum.at(aux).coef2* pump_volum.at(aux).pump_speed ^2 + pump_volum.at(aux).coef3* pump_volum.at(aux).pump_speed ^ 3;
	}
}

/*
void HydroNet_HeadLoss(std::vector<strObjects>* objects, std::vector<strPipe>* pipe, std::vector<strValve_fix>* valve_fix, std::vector<strValve_fix>* valve_var, std::vector<strThermostat>* thermostat, std::vector<strPump_turbo>* pump_turbo, int heat_exch_fix, int heat_exch_Tout, int tank, std::vector<int>* branches_ind, int branch_cycle, int  n_branch, int temperature) {
	int count;
	int count1;
	int id_aux;
	double opening;
	std::vector<double> hydr_resist2;

	// Calculates head loss in every branch of an hydraulic circuit
	std::vector <double> head_loss;
	std::vector <double> hydr_resist1;
	std::vector <double> hydr_resist2;
	std::string class_aux;
	head_loss.resize(n_branch + 1);							// head loss(m.c.f.) - independent of flow

	hydr_resist1.resize(n_branch + 1);						 // hydraulic resistance(m.c.f. / m ^ 3 * s)

															 // proportional to flow - for turbopumps

	hydr_resist2.resize(n_branch + 1);					// hydraulic resistance(m.c.f. / m ^ 6 * s ^ 2)
														//proportional to flow ^ 2


														// Branch loop
	for (count = 0; count <= n_branch; count++) {// goes only to branches that are part of a cycle

		if (branch_cycle == branches_ind.at(count) {

			//for count = n_branch // analyzes all branches
			for (count1 = branches_ind.at(count) {// objects in branch(indices)
												  //Mirar Concordancia id_aux
				id_aux = objects.id(count1); // id of current object
					class_aux = objects.at(count1); // class of current object

													// I can't use a switch because of the heat exchangers
					if (class_aux == "pipe") { // pipe : hydraulic resistance depends on characteristics
											   //Mirar Concordancia id_aux
						hydr_resist2.at(count) = hydr_resist2.at(count) + 8 * pipe.friction_coef.at(id_aux) * pipe.length.at(id_aux) / (pi ^ 2) / 9.81 / (pipe.diameter.at(id_aux) ^ 5);
					}

					else if (class_aux == "valve_fix") { // valve with fixed head loss
														 //Mirar Concordancia id_aux
						head_loss.at(count) = head_loss.at(count) + valve_fix.head_loss.at(id_aux);
					}

					else if (class_aux == "valve_var") { // valve with head loss dependant on opening
														 //ind_aux = find(valve_var.id == id_aux); // index of current object in class table
						head_loss.at(count) = head_loss.at(count) + valve_var.coef0.at(ind_aux) + valve_var.coef1.at(ind_aux) * valve_var.opening.at(ind_aux) + valve_var.coef2.at(ind_aux) * valve_var.opening.at(ind_aux) ^ 2 + valve_var.coef3.at(ind_aux) * valve_var.opening.at(ind_aux) ^ 3;

					}

					else if (class_aux == "thermostat") {
						// valve with head loss dependant on opening dependant on temperature
						//Mirar Concordancia id_aux
						//ind_aux = find(thermostat.id == id_aux);
						opening = thermostat.T_coef0.at(ind_aux) + thermostat.T_coef1.at(ind_aux) * temperature + thermostat.T_coef2.at(ind_aux) * temperature ^ 2 + thermostat.T_coef3.at(ind_aux) * temperature ^ 3;
						head_loss.at(count) = head_loss.at(count) + thermostat.h_coef0(ind_aux) + thermostat.h_coef1(ind_aux) * opening + thermostat.h_coef2(ind_aux) * opening ^ 2 + thermostat.h_coef3(ind_aux) * opening ^ 3;
					}

					else if (class_aux == "pump_turbo") {
						// turbopump: positive head->negative loss or resistance
						//Mirar Concordancia id_aux
						//ind_aux = find(pump_turbo.id == id_aux);
						head_loss.at(count) = head_loss.at(count) - (pump_turbo.coef_N0.at(ind_aux) + pump_turbo.coef_N1.at(ind_aux) * 	pump_turbo.pump_speed.at(idx_aux) + pump_turbo.coef_N2.at(ind_aux) * pump_turbo.pump_speed.at(ind_aux) ^ 2);
						hydr_resist1.at(count) = hydr_resist1.at(count) - (pump_turbo.Q_coef_N0.at(ind_aux) + pump_turbo.Q_coef_N1.at(ind_aux) * pump_turbo.pump_speed.at(ind_aux) + pump_turbo.Q_coef_N2.at(ind_aux) * pump_turbo.pump_speed.at(ind_aux) ^ 2);
						hydr_resist2.at(count) = hydr_resist2.at(count) - (pump_turbo.Q2_coef_N0.at(ind_aux) + pump_turbo.Q2_coef_N1.at(ind_aux) *pump_turbo.pump_speed.at(ind_aux) + pump_turbo.Q2_coef_N2.at(ind_aux) * pump_turbo.pump_speed.at(ind_aux) ^ 2);
					}

					else if (class_aux == "heat_exch") {// class starts with 'heat_exch'
														// heat exchanger : hydraulic resistance given

														//Mirar Concordancia id_aux
														//ind_aux = find(eval(strcat(class_aux, '.id')) == id_aux); // index of current object in class table
						hydr_resist2.at(count) = hydr_resist2.at(count) + eval(strcat(class_aux, '.hydr_resist(', num2str(ind_aux), ')'));

						//         elseif strcmp(class_aux, 'tank')
						//% tank: %%%%%%%%%%%%%%%%%% TO DO!
					}
			}
		}
	}

}
*/

int main()
{




	std::vector<strObjects> objects;
	std::vector<strPipe> pipe;
	std::vector<strValve_fix> valve_fix;
	std::vector<strValve_var> valve_var;
	std::vector<strThermostat> thermostat;
	std::vector<strPump_turbo> pump_turbo;
	std::vector<strPump_volum> pump_volum;
	std::vector<strHeat_exch_fix> heat_exch_fix;
	std::vector<strHeat_exch_Tout> heat_exch_Tout;
	std::vector<strTank> tank;
	std::vector<strBranches> branches;
	int temperature;

	HydroNet_HeadLoss(&objects, &pipe, &valve_fix, &valve_var, &thermostat, &pump_turbo, heat_exch_fix, heat_exch_Tout, tank, &branches_ind, branch_cycle, n_branch, temperature);

	return 0;
}

