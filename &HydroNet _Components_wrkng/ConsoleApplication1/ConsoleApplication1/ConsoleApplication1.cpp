// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <string>

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
	int outlet_object;					// Thermostat's outlet object identifier
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
	int outlet_object;				// Pump's outlet object identifier
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
	int outlet_object;				// Pump's outlet object identifier
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
	double T_in;					// Intlet temperaure [K]
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};

struct strHeat_exch_Tout
{
	int heat_exch_Tout;				// Identifier
	std::string name;				// Heat Exchanger name
	int inlet_object;				// Heat Exchanger's inlet object identifier
	int outlet_object;				// Heat Exchanger's outlet object identifier
	double T_out;					// Outlet temperaure [K]
	double T_in;					// Intlet temperaure [K]
	double hydr_resist;				// Hydraulic resistance = head loss / flow^2 [m.c.f. / m^6*s^2]
};
/*struct strTank
{
int tank;					// Identifier
std::string name;				// Tank's name
int inlet_object;				// Tank's inlet object identifier
int outlet_object;				// Tank's outlet object identifier
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
	std::vector<int> adjacent;		// Adjacent objects
};
struct strBranches
{
	int branch_ind;				// branch idintificator
	int branches_cycle;
	//int n_branch;				// Number of branches
	std::vector<strObjects> objects;	// Objects in Brnach
};

int main()
{
	std::vector<strBranches> branches;
	std::vector<std::vector<int>> edges;
	int n_obj;

	n_obj = 15;

	edges.resize(n_obj);
	for (int i = 0; i < n_obj; i++)
		edges.at(i).assign(n_obj, 0);




	branches.resize(6);


	branches.at(0).branch_ind = 0;
	branches.at(1).branch_ind = 1;
	branches.at(2).branch_ind = 2;
	branches.at(4).branch_ind = 3;
	branches.at(4).branch_ind = 4;
	branches.at(5).branch_ind = 5;


	branches.at(0).objects.resize(5);
	branches.at(1).objects.resize(5);
	branches.at(2).objects.resize(2);
	branches.at(3).objects.resize(3);
	branches.at(4).objects.resize(4);
	branches.at(5).objects.resize(4);

	branches.at(0).objects.at(0).Class = "pump_volum";
	branches.at(0).objects.at(0).ID = 0;
	branches.at(0).objects.at(0).name = "Main_Pump";
	branches.at(0).objects.at(0).volume = 1E-3;
	branches.at(0).objects.at(0).inlet_object = 12;
	branches.at(0).objects.at(0).outlet_oject = 1;
	branches.at(0).objects.at(0).adjacent = {12,1};

	branches.at(0).objects.at(1).Class = "pipe";
	branches.at(0).objects.at(1).ID = 1;
	branches.at(0).objects.at(1).name = "Enigine_Inlet";
	branches.at(0).objects.at(1).volume = 8.8357E-5;
	branches.at(0).objects.at(1).inlet_object = 0;
	branches.at(0).objects.at(1).outlet_oject = 2;
	branches.at(0).objects.at(1).adjacent = { 0,2 };

	branches.at(0).objects.at(2).Class = "pipe";
	branches.at(0).objects.at(2).ID = 12;
	branches.at(0).objects.at(2).name = "Return_Pipe";
	branches.at(0).objects.at(2).volume = 6.2832E-4;
	branches.at(0).objects.at(2).inlet_object = 11;
	branches.at(0).objects.at(2).outlet_oject = 0;
	branches.at(0).objects.at(2).adjacent = { 11,0 };

	branches.at(0).objects.at(3).Class = "node";
	branches.at(0).objects.at(3).ID = 11;
	branches.at(0).objects.at(3).adjacent = { 10,12,14 };

	branches.at(0).objects.at(3).Class = "node";
	branches.at(0).objects.at(3).ID = 2;
	branches.at(0).objects.at(3).adjacent = { 1,3,4 };

	branches.at(1).objects.at(0).Class = "thermostat";
	branches.at(1).objects.at(0).ID = 7;
	branches.at(1).objects.at(0).name = "Thermostat_rad";
	branches.at(1).objects.at(0).volume = 0;
	branches.at(1).objects.at(0).inlet_object = 6;
	branches.at(1).objects.at(0).outlet_oject = 9;
	branches.at(1).objects.at(0).adjacent = { 6,9 };

	branches.at(1).objects.at(1).Class = "heat_exch_fix";
	branches.at(1).objects.at(1).ID = 9;
	branches.at(1).objects.at(1).name = "Radiator";
	branches.at(1).objects.at(1).volume = 5E-4;
	branches.at(1).objects.at(1).inlet_object = 7;
	branches.at(1).objects.at(1).outlet_oject = 10;
	branches.at(1).objects.at(1).adjacent = { 7,10 };

	branches.at(1).objects.at(2).Class = "pipe";
	branches.at(1).objects.at(2).ID = 10;
	branches.at(1).objects.at(2).name = "Radiator_outlet";
	branches.at(1).objects.at(2).volume = 3.1416E-4;
	branches.at(1).objects.at(2).inlet_object = 9;
	branches.at(1).objects.at(2).outlet_oject = 11;
	branches.at(1).objects.at(2).adjacent = { 9,11 };

	branches.at(1).objects.at(3).Class = "node";
	branches.at(1).objects.at(3).ID = 11;
	branches.at(1).objects.at(3).adjacent = { 10,12,14 };

	branches.at(1).objects.at(4).Class = "node";
	branches.at(1).objects.at(4).ID = 6;
	branches.at(1).objects.at(4).adjacent = { 5,7,8 };

	branches.at(2).objects.at(0).Class = "node";
	branches.at(2).objects.at(0).ID = 5;
	branches.at(2).objects.at(0).adjacent = { 3,6,13};

	branches.at(2).objects.at(1).Class = "node";
	branches.at(2).objects.at(1).ID = 6;
	branches.at(2).objects.at(1).adjacent = { 5,7,8 };




	branches.at(3).objects.at(0).Class = "heat_exch_fix";
	branches.at(3).objects.at(0).ID = 3;
	branches.at(3).objects.at(0).name = "Block_cyl1";
	branches.at(3).objects.at(0).volume = 5E-4;
	branches.at(3).objects.at(0).inlet_object = 2;
	branches.at(3).objects.at(0).outlet_oject = 5;
	branches.at(3).objects.at(0).adjacent = { 2,5 };

	branches.at(3).objects.at(1).Class = "node";
	branches.at(3).objects.at(1).ID = 2;
	branches.at(3).objects.at(1).adjacent = {1,3,4 };

	branches.at(3).objects.at(2).Class = "node";
	branches.at(3).objects.at(2).ID = 5;
	branches.at(3).objects.at(2).adjacent = { 3,6,13 };

	branches.at(4).objects.at(0).Class = "heat_exch_fix";
	branches.at(4).objects.at(0).ID = 4;
	branches.at(4).objects.at(0).name = "Head_cyl1";
	branches.at(4).objects.at(0).volume = 5E-4;
	branches.at(4).objects.at(0).inlet_object = 2;
	branches.at(4).objects.at(0).outlet_oject = 13;
	branches.at(4).objects.at(0).adjacent = { 2,13 };

	branches.at(4).objects.at(1).Class = "heat_exch_fix";
	branches.at(4).objects.at(1).ID = 13;
	branches.at(4).objects.at(1).name = "Head_cyl2";
	branches.at(4).objects.at(1).volume = 5E-4;
	branches.at(4).objects.at(1).inlet_object = 4;
	branches.at(4).objects.at(1).outlet_oject = 5;
	branches.at(4).objects.at(1).adjacent = { 4,5 };

	branches.at(4).objects.at(2).Class = "node";
	branches.at(4).objects.at(2).ID = 2;
	branches.at(4).objects.at(2).adjacent = { 1,3,4 };

	branches.at(4).objects.at(3).Class = "node";
	branches.at(4).objects.at(3).ID = 5;
	branches.at(4).objects.at(3).adjacent = { 3,6,13 };

	branches.at(5).objects.at(0).Class = "thermostat";
	branches.at(5).objects.at(0).ID = 8;
	branches.at(5).objects.at(0).name = "Thermostat_byp";
	branches.at(5).objects.at(0).volume = 0;
	branches.at(5).objects.at(0).inlet_object = 6;
	branches.at(5).objects.at(0).outlet_oject = 14;
	branches.at(5).objects.at(0).adjacent = { 6,14 };


	branches.at(5).objects.at(1).Class = "pipe";
	branches.at(5).objects.at(1).ID = 14;
	branches.at(5).objects.at(1).name = "Bypass_Pipe";
	branches.at(5).objects.at(1).volume = 3.1416E-4;
	branches.at(5).objects.at(1).inlet_object = 8;
	branches.at(5).objects.at(1).outlet_oject = 11;
	branches.at(5).objects.at(1).adjacent = { 8,11 };


	branches.at(5).objects.at(2).Class = "node";
	branches.at(5).objects.at(2).ID = 11;
	branches.at(5).objects.at(2).adjacent = { 10,12,14 };

	branches.at(5).objects.at(2).Class = "node";
	branches.at(5).objects.at(2).ID = 6;
	branches.at(5).objects.at(2).adjacent = { 5,7,8 };


	for (int brnch = 0; brnch < branches.size(); brnch++)
		for (int obj = 0; obj < branches.at(brnch).objects.size(); obj++)
			for (int adj = 0; adj < branches.at(brnch).objects.at(obj).adjacent.size(); adj++)
				edges.at(branches.at(brnch).objects.at(obj).ID).at(branches.at(brnch).objects.at(obj).adjacent.at(adj)) = 1;











    return 0;
}

