// HydroNet_.cpp : Defines the entry point for the application.
//

#include "stdafx.h"
#include "HydroNet.h"


bool all(std::vector<bool> A) {				// Returns true if every component in the vector is true, if not returns false
	bool all = true;
	for (int i = 0; i < A.size(); i++)
		if (A.at(i) == false) {
			all = false;
			break;
		}
	return all;
}

bool any(std::vector<bool> A) {					// Returns true if there is any componet in the vector true
	bool any = false;
	for (int i = 0; i < A.size(); i++)
		if (A.at(i) == true) {
			any = true;
			break;
		}
	return any;
}

bool in(int a, std::vector<int> A) {			// a is part of the vector A
	bool in = false;

	for (int i = 0; i < A.size(); i++)
		if (a == A.at(i)) {
			in = true;
			break;
		}
	return in;
}				

std::vector<double> roots(double a, double b, double c) {			//Roots a*x^2+b*x+c=0
	double r1, r2;
	std::vector<double> r(2);
	r1 = abs((-b + sqrt(b*b - 4 * a*c)) / (2 * a));
	r2 = abs((-b - sqrt(b*b - 4 * a*c)) / (2 * a));

	r.at(0) = Math_wamH::Max(r1, r2);
	r.at(1) = Math_wamH::Min(r1, r2);

	return r;
}

double vec_sum_elements(std::vector<double> in, std::vector<int> elem, int elem_size) {		// Adds the vecctor in elements presents at vector elem
	double sum = 0;
	for (int i = 0; i < elem_size; i++)
		sum = sum + in.at(elem.at(i));
	return sum;
}

double density(std::string fluid_type, double temperature)				
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 1000;	// for tests: water [Kg/m3]
}

double cp(std::string fluid_type, double temperature)
{
	// dummy fuction to be ignored when the model is integrated in VEMOD

	return 4181.3;	// for tests: water [J/kg/K]
}

template < class T>
T flip(T A) {					// Flips vector A
	T aux;
	aux.resize(A.size());
	for (int i = 0; i < A.size(); i++)
		aux.at(i) = A.at(A.size() - i - 1);
	return aux;
}

class Solver {
	// Functor
private:
	std::vector<double> A, B, C, D;
	int n;



public:
	Solver(Eigen::VectorXd x0, Eigen::VectorXd data) {

		//Matrix recostruction

		int count = 0;

		n = x0.size();
		A.resize(n);
		B.resize(n*n);
		C.resize(n*n);
		D.resize(n*n);


		for (count; count < n; count++)
			A[count] = data[count];
		for (count; count < n*n + n; count++)
			B[count - n] = data[count];
		for (count; count < 2 * n*n + n; count++)
			C[count - n*n - n] = data[count];
		for (count; count < 3 * n*n + n; count++)
			D[count - 2 * n*n - n] = data[count];

	};


	int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const;
	int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const;
	int inputs() const;// inputs is the dimension of x.
	int values() const; // "values" is the dimension of F 

};

int Solver::operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {

	double temp, temp1, temp2;
	assert(fvec.size() == n);
	//Implementation of Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional sistems
	for (int k = 0; k < n; k++) {

		temp = 0;
		temp = A[k];
		for (int j = 0; j<n; j++)
			temp = temp + B[j + k * n] * x[j] / abs(x[j]);

		temp1 = 0;
		for (int j = 0; j<n; j++)
			temp1 = temp1 + C[j + k * n] * x[j];


		temp2 = 0;
		for (int j = 0; j<n; j++)
			temp2 = temp2 + D[j + k * n] * abs(x[j])*x[j];

		fvec[k] = temp + temp1 + temp2;
	}
	return 0;
}

int Solver::df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const {

	assert(fjac.rows() == n);
	assert(fjac.cols() == n);
	//Implementation of Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional sistems Jacobian
	for (int k = 0; k < n; k++)
		for (int j = 0; j < n; j++) {

			fjac(k, j) = 2 * D[n*k + j] * x[j] + C[n*k + j];
		}
	return (0);

}
int Solver::inputs() const { return n; } // inputs is the dimension of x.
int Solver::values() const { return n; } // "values" is the number of f_i 

std::vector<double> NonLinearSolver(Eigen::VectorXd _sol, std::vector<std::vector<double>> coef1_mat, std::vector<std::vector<double>> coef2_mat, std::vector<std::vector<double>> const_mat, std::vector<double> const_vec)
{
	/*
	A is an independent constant vector
	B is a vector dependent on the solution sing
	C ia a vector containing the lineal coefficients
	D ia a vector containing the quadratic coefficients
	Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional
	*/

	int n = _sol.size();
	std::vector <double> sol;
	int info;
	Eigen::VectorXd data(10 * n);
	//int count = 0;
	std::vector<double> _data;
	_data.resize(10 * n);


	//The matrix an vector's coeficients are stored in a single Eigen vector
	for (int count = 0; count < n; count++) {
		data[count] = const_vec[count];
		for (int count1 = 0; count1 < n; count1++) {
			data[(count + 1)*n + count1] = const_mat.at(count)[count1];
			data[(count + 4)*n + count1] = coef1_mat.at(count).at(count1);
			data[(count + 7)*n + count1] = coef2_mat.at(count).at(count1);
		}


	}

	// To avoid inestabilities when the starting point is close to 0
	for (int j = 0; j < _sol.size(); j++)
		if (abs(_sol[j]) < 1E-9)
			_sol[j] = _sol[j] / abs(_sol[j] * 1E-9);

	//Computation
	Solver functor(_sol, data);
	Eigen::HybridNonLinearSolver<Solver> solver(functor);
	solver.diag.setConstant(n, 1.);
	solver.useExternalScaling = true;
	info = solver.hybrj1(_sol);


	//The results are stored in the vector _results
	sol.resize(_sol.size());
	for (int count = 0; count < _sol.size(); count++)
		sol.at(count) = _sol(count);

	return sol;
}


double efficiency(double a, double b) {						// Fake efficiency calculation
	return 1;
}


void HydroNet_model::HydroNet_HeadLoss() {						// Calculates each object's hydraulic resistance

	std::string _class;
	int count;
	int aux = 0;
	double opening;
	double temperature_aux;
	//double obj_pos;
	int _ID;
	int index=0;
	//std::vector <double> head_loss;
	//std::vector <double> hydr_resist1;
	//std::vector <double> hydr_resist2;

	head_loss.resize(n_branch + 1);							// head loss(m.c.f.) - independent of flow

	hydr_resist1.resize(n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 3 * s)
															// proportional to flow - for turbopumps
	hydr_resist2.resize(n_branch + 1);						// hydraulic resistance(m.c.f. / m ^ 6 * s ^ 2)
															//proportional to flow ^ 2


															// Ads the values of head_loss, hydr_resist1 and hydr_resist2 for each branch
	for (count = 1; count < n_branch; count++) {

		head_loss.at(count) = 0.0;
		hydr_resist1.at(count) = 0.0;
		hydr_resist2.at(count) = 0.0;
		opening = 0.0;
		aux = 1;

		while (aux < branches.at(count).objects.size())
		{
			_class = branches.at(count).objects.at(aux).Class;
			_ID = branches.at(count).objects.at(aux).ID;



			if (_class == "pipe") {

				while (index < pipe.size() && _ID != pipe.at(index).pipe_id) {				// Looks for the object

					index++;
				}
				// Hydro resistance calculation 
				hydr_resist2.at(count) = hydr_resist2.at(count) + 8 * pipe.at(index).friction_coef * pipe.at(index).length / (__cons::Pi*__cons::Pi / 9.81 / pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter * pipe.at(index).diameter);
			}

			else if (_class == "valve_fix") {

				while (index < valve_fix.size() && _ID != valve_fix.at(index).valve_fix_id) {				// Looks for the object

					index++;
				}
				// Hydro resistance calculation 	
				head_loss.at(count) = head_loss.at(count) + valve_fix.at(index).head_loss;
			}

			else if (_class == "valve_var") {

				while (index < valve_var.size() && _ID != valve_var.at(index).valve_var_id) {							// Looks for the object

					index++;
				}
				// Hydro resistance calculation 	
				head_loss.at(count) = head_loss.at(count) + valve_var.at(index).coef0 + valve_var.at(index).coef1 * valve_var.at(index).opening + valve_var.at(index).coef2 * valve_var.at(index).opening*valve_var.at(index).opening + valve_var.at(index).coef3 * valve_var.at(index).opening*valve_var.at(index).opening*valve_var.at(index).opening;
			}
			else if (_class == "thermostat") {

				while (index < thermostat.size() && _ID != thermostat.at(index).themostat_id) {						// Looks for the object

					index++;;
				}
				// Hydro resistance calculation 	
				obj_pos = (obj_inlet_pos.at(count).at(index) + obj_outlet_pos.at(count).at(index)) / 2;																	//position inside the branch as a % of the branch volume

				temperature_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos.at(count), branch_temperature.at(count));						// temperature in the thermostat

				opening = thermostat.at(index).opening;

				head_loss.at(count) = head_loss.at(count) + thermostat.at(index).coef_h_0 + thermostat.at(index).coef_h_1 * opening + thermostat.at(index).coef_h_2 * opening*opening + thermostat.at(index).coef_h_3 * opening*opening*opening;

			}

			else if (_class == "pump_turbo") {

				while (index < pump_turbo.size() && _ID != pump_turbo.at(index).Pump_turbo) {							// Looks for the object

					index++;
				}

				// Hydro resistance calculation 	
				head_loss.at(count) = head_loss.at(count) - (pump_turbo.at(index).coef_N0 + pump_turbo.at(index).coef_N1 * pump_turbo.at(index).pump_speed + pump_turbo.at(index).coef_N2 * pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

				hydr_resist1.at(count) = hydr_resist1.at(count) - (pump_turbo.at(index).Q_coef_N0 + pump_turbo.at(index).Q_coef_N1* pump_turbo.at(index).pump_speed + pump_turbo.at(index).Q_coef_N2* pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

				hydr_resist2.at(count) = hydr_resist2.at(count) - (pump_turbo.at(index).Q2_coef_N0 + pump_turbo.at(index).Q2_coef_N1* pump_turbo.at(index).pump_speed + pump_turbo.at(index).Q2_coef_N2* pump_turbo.at(index).pump_speed*pump_turbo.at(index).pump_speed);

			}

			else if (_class == "pump_volum") {

				while (index < pump_volum.size() && _ID != pump_volum.at(index).Pump_volum) {					// Looks for the object

					index++;
				}
				// Hydro resistance calculation 	
				hydr_resist1.at(count) = hydr_resist1.at(count) - (pump_volum.at(index).coef0 + pump_volum.at(index).coef1* pump_volum.at(index).pump_speed + pump_volum.at(index).coef2* pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed + pump_volum.at(index).coef3* pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed*pump_volum.at(index).pump_speed);
			}

			else if (_class == "heat_exch") {

				while (index < heat_exch.size() && _ID != heat_exch.at(index).heat_exch) {				// Looks for the object

					index++;
				}
				// Hydro resistance calculation 	
				hydr_resist2.at(count) = hydr_resist2.at(count) + heat_exch.at(index).hydr_resist;
			}

			//else if (_class == "heat_exch_Tout") {

			//	while (index < heat_exch_Tout.size() && _ID != heat_exch_Tout.at(index).heat_exch_Tout) {				// Looks for the object

			//		index++;
			//	}
			//	// Hydro resistance calculation 	
			//	hydr_resist2.at(count) = hydr_resist2.at(count) + heat_exch_Tout.at(index).hydr_resist;
			//}

			/*
			else if (_class == "tank"){

			while (index < valve_var.size() && _ID != valve_var.at(index).pipe_id) {				// Looks for the object

			index++
			}
			// Hydro resistance calculation 	
			head_loss.at(count) = head_loss.at(count) + tank.at(index).head_loss;
			}
			*/
			aux++;
			index = 0;
		}
	}
}

void HydroNet_model::HydroNet_Main() {

	std::vector<double>  dt_vec, volum_pump_speed, turbo_pump_speed, valve_opening, thst_sensitivity, thermostat_opening_new;
	double  thermostat_temp;
	int n_exec;
	double Tmin, Tmax, y_end, x0, k ;
	

	// Controls data flow


	//// MODEL CREATION
	// All valves and thermostats are assumed open.

	//circuit_file = "data//test.txt";
	circuit_file = "data//Coolant_basic_nodes - copia.txt";

	HydroNet_Create( );



	//// INITIALIZATION

	gravity_acc = 9.81; // gravity acceleration [m/s^2] for pump power

						// Cell to store volumetric position of the points where temperature changes,
						// as a // of the total volume of the branch // temperature between those points is constant
						// branch_temp_pos = cell(n_branch, 1);
	branch_temp_pos.resize(7);
	branch_temp_pos.at(0) = { 0,10,100 };
	branch_temp_pos.at(1) = { 0, 100 };
	branch_temp_pos.at(2) = { 0, 100 };
	branch_temp_pos.at(3) = { 0, 20, 100 };
	branch_temp_pos.at(4) = { 0, 30, 100 };
	branch_temp_pos.at(5) = { 0, 100 };
	branch_temp_pos.at(6) = { 0, 10, 100 };


	// Cell to store temperatures of every volume section in the branch
	// branch_temperature = cell(n_branch, 1);
	branch_temperature.resize(7);
	branch_temperature.at(0) = { 290,290 };
	branch_temperature.at(1) = { 290 };
	branch_temperature.at(2) = { 290 };
	branch_temperature.at(3) = { 290, 290 };
	branch_temperature.at(4) = { 290, 290 };
	branch_temperature.at(5) = { 290 };
	branch_temperature.at(6) = { 290, 290 };

	// Output variables related to branches
	head_loss.assign(n_branch, 0);
	hydr_resist1.assign(n_branch, 0);
	hydr_resist2.assign(n_branch, 0);
	flows.assign(n_branch, 0);



	//// EXECUTIONS
	// Some variables must be recalculated only if other change

	n_exec = 5;

	dt_vec = { 0.1, 0.1, 0.1, 0.1, 0.2 };

	volum_pump_speed = { 105, 210, 314, 419, 419 };
	turbo_pump_speed = { 105, 210, 314, 419, 419 };
	valve_opening = { 20, 50, 70, 80, 30 };

	thst_sensitivity.assign(n_exec, 5); // thermostat sensitivity [%]


										// Execution loop
	flag_first_exec = true;
	for (int count_exec = 0; count_exec< n_exec; count_exec++) {

		// TIME SPAN
		dt = dt_vec.at(count_exec);


		// PUMP SPEED
		// If the speed of any pump changes, flows must be recalculated

		flag_pump_speed_change = false;
		for (int count = 0; count<pump_volum.size(); count++)
			if (volum_pump_speed.at(count_exec) != pump_volum.at(count).pump_speed)
				flag_pump_speed_change = true;

		for (int count = 0; count<pump_turbo.size(); count++)
			if (turbo_pump_speed.at(count_exec) != pump_turbo.at(count).pump_speed)
				flag_pump_speed_change = true;

		for (int count = 0; count<pump_volum.size(); count++)
			pump_volum.at(count).pump_speed = volum_pump_speed.at(count_exec);
		for (int count = 0; count<pump_turbo.size(); count++)
			pump_turbo.at(count).pump_speed = turbo_pump_speed.at(count_exec);


		// VALVES
		// If the valve opening varies, vectors for head losses and thus flows must be recalculated.
		// If the valve is closed, flow in its branch is zero.

		flag_valve_change = false;
		for (int count = 0; count<valve_var.size(); count++)
			if (valve_opening.at(count_exec) != valve_var.at(count).opening)
				flag_valve_change = true;

		for (int count = 0; count<valve_var.size(); count++)
			valve_var.at(count).opening = valve_opening.at(count_exec);


		// THERMOSTATS
		// If the thermostat opening varies because there is a change in the inlet
		// temperature, vectors for head losses and thus flows must be recalculated.
		// Head losses and flows will be recalculated only if the opening variation
		// of any thermostat is larger than the sensitivity of the thermostat[%]
		// If the thermostat is closed, flow in its branch is zero.

		thermostat_opening_new.assign(thermostat.size(), 0);
		for (int count_thst = 0, obj_aux, branch_aux; count_thst < thermostat.size(); count_thst++) {
			obj_aux = thermostat.at(count_thst).inlet_object;
			for (int count = 0; count<branches.size(); count++)
				for (int count1 = 0; count1<branches.at(count).objects.size(); count1++)
					if (branches.at(count).objects.at(count1).ID == obj_aux)
						branch_aux = count;
			Tmin = thermostat.at(count_thst).T_min;
			Tmax = thermostat.at(count_thst).T_max;
			y_end = thermostat.at(count_thst).Shape_factor;

			x0 = (Tmax + Tmin) / 2;
			k = -log(y_end) / (Tmin - x0);
			obj_pos = (obj_inlet_pos.at(obj_aux).at(0) + obj_outlet_pos.at(obj_aux).at(0)) / 2;

			thermostat_temp = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos.at(branch_aux), branch_temperature.at(branch_aux));
			if (thermostat.at(count_thst).to_open)
				thermostat_opening_new.at(count_thst) = 1. / (1 + exp(-k * (thermostat_temp - x0)));
			else
				thermostat_opening_new.at(count_thst) = 1 - 1. / (1 + exp(-k * (thermostat_temp - x0)));
		}

		flag_thermostat_changes = false;
		for (int count = 0; count< thermostat.size(); count++)
			if (abs(thermostat_opening_new.at(count) - thermostat.at(count).opening) / thermostat.at(count).opening * 100 > thst_sensitivity.at(count_exec))
				flag_thermostat_changes = true;

		for (int count = 0; count< thermostat.size(); count++)
			thermostat.at(count).opening = thermostat_opening_new.at(count);


		// MAIN EXECUTION FUNCTION CALL
		HydroNet_Executions();

		flag_first_exec = false;

		
	}


}

void HydroNet_model::HydroNet_ReadObj() {
	// Reads circuit from configuration file and prepares the data in arrays.
	//   Input is the name of the file that contains the information of the
	//   circuit. Outputs are an object's list with id, name, class and joints
	//   (table  'objects') and one table for all objects of each class (except
	//   nodes). In addition an analysis is performed to check if all objects
	//   are connected and identifiers are right.


	
	std::string _aux;
	int adj_size, pump_volum_count = 0, pump_turbo_count = 0,/* heat_exch_Tout_count = 0,*/ heat_exch_count = 0, pipe_count = 0, valve_fix_count = 0, valve_var_count = 0, thermostat_count = 0, tank_count = 0;
	std::ifstream circuit_text;

	n_obj = 0;


	circuit_text.open(circuit_file, std::ifstream::in);					// File opening

	std::getline(circuit_text, fluid_type, '\t');

	while (!circuit_text.eof()) {									// Reads and stores the object data present in the file

		n_obj++;

		objects.resize(n_obj);

		circuit_text.ignore(1000,'#');

		std::getline(circuit_text, _aux, '\t');

		objects.at(n_obj - 1).ID = std::stoi(_aux)-1;

		circuit_text.ignore(1000, '\n');

		std::getline(circuit_text, objects.at(n_obj - 1).name, '\t');

		circuit_text.ignore(1000, '\n');

		std::getline(circuit_text, objects.at(n_obj - 1).Class, '\t');

		// Knowing the object's class, reads and stores the data folowing each objects structuture

		if (objects.at(n_obj - 1).Class == "node") {

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			adj_size = std::stoi(_aux);

			objects.at(n_obj - 1).adjacent.resize(adj_size);

			for (int count = 0; count < adj_size; count++) {

				circuit_text.ignore(1000, '\n');

				std::getline(circuit_text, _aux, '\t');

				objects.at(n_obj - 1).adjacent.at(count) = std::stoi(_aux)-1;

			}
		}
		else if (objects.at(n_obj - 1).Class == "pump_volum") {

			pump_volum_count++;

			pump_volum.resize(pump_volum_count);

			pump_volum.at(pump_volum_count-1).name = objects.at(n_obj - 1).name;
			pump_volum.at(pump_volum_count-1).Pump_volum = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			pump_volum.at(pump_volum_count-1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			pump_volum.at(pump_volum_count-1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).volume = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).pump_speed = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).head_max = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).coef0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).coef1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).coef2 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).coef3 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).pump_head = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_volum.at(pump_volum_count-1).pump_power = std::stod(_aux);


		}
		else if (objects.at(n_obj - 1).Class == "pump_turbo") {

			pump_turbo_count++;

			pump_turbo.resize(pump_volum_count);

			pump_turbo.at(pump_turbo_count-1).name = objects.at(n_obj - 1).name;
			pump_turbo.at(pump_turbo_count - 1).Pump_turbo = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			pump_turbo.at(pump_turbo_count - 1).inlet_object = std::stoi(_aux) - 1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			pump_turbo.at(pump_turbo_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).volume = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).pump_speed = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).coef_N0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).coef_N1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).coef_N2 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q_coef_N0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q_coef_N1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q_coef_N2 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q2_coef_N0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q2_coef_N1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).Q2_coef_N2 = std::stod(_aux);


			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).pump_head = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pump_turbo.at(pump_turbo_count - 1).pump_power = std::stod(_aux);


		}
		else if (objects.at(n_obj - 1).Class == "pipe") {

			pipe_count++;

			pipe.resize(pipe_count);

			pipe.at(pipe_count-1).name = objects.at(n_obj - 1).name;
			pipe.at(pipe_count - 1).pipe_id = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			pipe.at(pipe_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			pipe.at(pipe_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pipe.at(pipe_count - 1).length = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pipe.at(pipe_count - 1).diameter = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			pipe.at(pipe_count - 1).friction_coef = std::stod(_aux);

			objects.at(n_obj - 1).volume = __cons::Pi*pipe.at(pipe_count - 1).diameter*pipe.at(pipe_count - 1).length;


		}
		else if (objects.at(n_obj - 1).Class == "heat_exch") {

			heat_exch_count++;

			heat_exch.resize(heat_exch_count);

			heat_exch.at(heat_exch_count-1).name = objects.at(n_obj - 1).name;
			heat_exch.at(heat_exch_count - 1).heat_exch = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			heat_exch.at(heat_exch_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			heat_exch.at(heat_exch_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).volume = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch.at(heat_exch_count - 1).type = (_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch.at(heat_exch_count - 1).heat = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch.at(heat_exch_count - 1).T_out = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch.at(heat_exch_count - 1).hydr_resist = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch.at(heat_exch_count - 1).T_in = std::stod(_aux);


		}
		/*else if (objects.at(n_obj - 1).Class == "heat_exch_Tout") {

			heat_exch_Tout_count++;

			heat_exch_Tout.resize(heat_exch_Tout_count);

			heat_exch_Tout.at(heat_exch_Tout_count-1).name = objects.at(n_obj - 1).name;
			heat_exch_Tout.at(heat_exch_Tout_count - 1).heat_exch_Tout = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			heat_exch_Tout.at(heat_exch_Tout_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			heat_exch_Tout.at(heat_exch_Tout_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).volume = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch_Tout.at(heat_exch_Tout_count - 1).T_in = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch_Tout.at(heat_exch_Tout_count - 1).T_out = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			heat_exch_Tout.at(heat_exch_Tout_count - 1).hydr_resist = std::stod(_aux);


		}*/
		else if (objects.at(n_obj - 1).Class == "thermostat") {

			thermostat_count++;

			thermostat.resize(thermostat_count);

			thermostat.at(thermostat_count-1).name = objects.at(n_obj - 1).name;
			thermostat.at(thermostat_count - 1).themostat_id = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			thermostat.at(thermostat_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			thermostat.at(thermostat_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).T_min = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).T_max = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).Shape_factor = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).to_open = std::stoi(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).opening = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).coef_h_0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).coef_h_1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).coef_h_2 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			thermostat.at(thermostat_count - 1).coef_h_3 = std::stod(_aux);

		}
		else if (objects.at(n_obj - 1).Class == "valve_var") {

			valve_var_count++;

			valve_var.resize(valve_var_count);

			valve_var.at(valve_var_count-1).name = objects.at(n_obj - 1).name;
			valve_var.at(valve_var_count - 1).valve_var_id = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			valve_var.at(valve_var_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			valve_var.at(valve_var_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_var.at(valve_var_count - 1).opening = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_var.at(valve_var_count - 1).coef0 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_var.at(valve_var_count - 1).coef1 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_var.at(valve_var_count - 1).coef2 = std::stod(_aux);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_var.at(valve_var_count - 1).coef3 = std::stod(_aux);


		}
		else if (objects.at(n_obj - 1).Class == "valve_fix") {

			valve_fix_count++;

			valve_fix.resize(valve_fix_count);

			valve_fix.at(valve_fix_count-1).name = objects.at(n_obj - 1).name;
			valve_fix.at(valve_fix_count - 1).valve_fix_id = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			valve_fix.at(valve_fix_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			valve_fix.at(valve_fix_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			valve_fix.at(valve_fix_count - 1).head_loss = std::stod(_aux);

		}
		else if (objects.at(n_obj - 1).Class == "tank") {

			tank_count++;

			tank.resize(tank_count);

			tank.at(tank_count-1).name = objects.at(n_obj - 1).name;
			tank.at(tank_count - 1).tank = objects.at(n_obj - 1).ID;


			objects.at(n_obj - 1).adjacent.resize(2);

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(0) = std::stoi(_aux)-1;
			tank.at(tank_count - 1).inlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			objects.at(n_obj - 1).adjacent.at(1) = std::stoi(_aux)-1;
			tank.at(tank_count - 1).outlet_object = std::stoi(_aux)-1;

			circuit_text.ignore(1000, '\n');

			std::getline(circuit_text, _aux, '\t');

			tank.at(tank_count - 1).head_loss = std::stod(_aux);

		}
	}

	n_obj--;

	objects.resize(n_obj);
	objects.shrink_to_fit();

	}

void HydroNet_model::HydroNet_Create() {
	// Creates a model for the specified circuit
	int tank_count;
	double sum_aux;
	std::vector<int> inlet_aux;
	std::vector<std::vector<int>> obj_inlet_pos_cell;

	//// ANALYSIS OF CIRCUIT OBJECTS
	// Reads circuit configuration and specifications of objects and stores data in tables

	HydroNet_ReadObj();

	n_nodes = 0;

	//// NODES
	// List of all nodes
	for (int count = 0; count<n_obj; count++)
		if (objects.at(count).Class == "node") {
			nodes_ind.resize(n_nodes + 1);
			nodes_id.resize(n_nodes + 1);
			nodes_ind.at(n_nodes) = count;
			nodes_id.at(n_nodes) = objects.at(count).ID;
			n_nodes++;
		}

	
	//// TANKS
	// List of all tanks

	tanks_ind.assign(n_obj, 0); // stores tank ids
								// maximum size: all objects are tanks

	tank_count = 0;
	for (int count = 0; count < n_obj; count++) // objects loop
		if (objects.at(count).Class == "tank") {
			tank_count = tank_count + 1;
			tanks_ind.at(tank_count) = count;
		}

	tanks_ind.resize(tank_count); // adjust size

	n_tanks = tanks_ind.size();

	
	//// CIRCUIT COMPONENTS
	// Obtains nodes, branches and meshes in the circuit.
	// Determines branches connected to nodes and branches that are part of a mesh.
	// All valves and thermostats are assumed open.

	HydroNet_Components();



	//// BRANCH VOLUME

	branch_volume.assign(n_branch, 0);

	for (int count = 0; count < n_branch; count++) {													// branch loop
		sum_aux = 0;
		for (int count1 = 0, aux; count1 < branches_ind.at(count).size(); count1++) {				 // loop of objects in the branch
			aux = branches_ind.at(count).at(count1);
			sum_aux = sum_aux + objects.at(count1).volume;
		}
		branch_volume.at(count) = sum_aux;
	}



	//// INLET AND OUTLET OF OBJECTS
	// Stores volumetric position of the inlet and outlet of objects, as a \// of the total branch volume

	// Inlet positions
	obj_inlet_pos_cell.resize(n_branch); // branch / objects format
	obj_inlet_pos.resize(n_obj); // "vector" of objects format
	for (int count_branch = 0; count_branch < n_branch; count_branch++) {// branch loop
		inlet_aux.assign(branches_ind.at(count_branch).size(), 0); // auxiliary vector to store inlets in branch
		for (int count = 0, count_obj; count < branches_ind.at(count_branch).size(); count++) { // initialization loop
			count_obj = branches_ind.at(count_branch).at(count);
			obj_inlet_pos.at(count_obj).push_back(objects.at(count_obj).adjacent.at(0)); // all objects of the branch have the inlet at 0//
		}
		if (branch_volume.at(count_branch) != 0) { // to avoid dividing by zero if branch volume is 0
			for (int count_obj = 1; count_obj < branches_ind.at(count_branch).size(); count_obj++) { // loop of objects in the branch
				inlet_aux.at(count_obj) = inlet_aux.at(count_obj - 1) + objects.at(branches_ind.at(count_branch).at(count_obj - 1)).volume / branch_volume.at(count_branch) * 100;
				// \// of the previous object inlet + volume previous object / branch volume * 100
				for (int count1 = 1; count1 <objects.at(count_obj).adjacent.size(); count1++)
					obj_inlet_pos.at(branches_ind.at(count_branch).at(count_obj)).push_back(objects.at(count_obj).adjacent.at(count1));
				obj_inlet_pos.at(branches_ind.at(count_branch).at(count_obj)).push_back(inlet_aux.at(count_obj));
			}
		}

		obj_inlet_pos_cell.at(count_branch) = inlet_aux; // stores inlets of the current branch
	}


	// Outlet positions
	obj_outlet_pos.resize(n_obj);
	for (int count = 0; count < n_obj; count++) {
		if (objects.at(count).volume == 0)
			obj_outlet_pos.at(count) = obj_inlet_pos.at(count);
		else // volume > 0
			for (int count1 = 0; count1 <obj_inlet_pos.at(count).size(); count1++)
				obj_outlet_pos.at(count).push_back(obj_inlet_pos.at(count).at(count1) + objects.at(count).volume / branch_volume.at(objects.at(count).branch) * 100);
	}


	//// INDEXES OF HEAT EXCHANGERS & PUMPS IN EVERY BRANCH

	branch_htx.resize(n_branch);
	//branch_htx_Tout.resize(n_branch);
	branch_pump_volum.resize(n_branch);
	branch_pump_turbo.resize(n_branch);
	//Initialisation of vectors
	for (int count_branch = 0; count_branch < n_branch; count_branch++) {
		branch_htx.at(count_branch) = {};
		//branch_htx_Tout.at(count_branch) = {};
		branch_pump_volum.at(count_branch) = {};
		branch_pump_turbo.at(count_branch) = {};
	}

	for (int count_branch = 0; count_branch < n_branch; count_branch++) {
		for (int count_obj = 0; count_obj < branches.at(count_branch).objects.size(); count_obj++) {
			if (branches.at(count_branch).objects.at(count_obj).Class == "heat_exch")								// Asigns the heat exchanger identificator when the object is a Heat Echanger type fix
				branch_htx.at(count_branch).push_back(branches.at(count_branch).objects.at(count_obj).ID);
			//if (branches.at(count_branch).objects.at(count_obj).Class == "heat_exch_Tout")								// Asigns the heat exchanger identificator when the object is a Heat Echanger type Tout
				//branch_htx_Tout.at(count_branch).push_back(branches.at(count_branch).objects.at(count_obj).ID);
			if (branches.at(count_branch).objects.at(count_obj).Class == "pump_volum")									// Asigns the pump identificator when the object is a Volumetric Pump
				branch_pump_volum.at(count_branch).push_back(branches.at(count_branch).objects.at(count_obj).ID);
			if (branches.at(count_branch).objects.at(count_obj).Class == "pump_turbo")									// Asigns the pump identificator when the object is a Turbo Pump
				branch_pump_turbo.at(count_branch).push_back(branches.at(count_branch).objects.at(count_obj).ID);
		}
	}

}

void HydroNet_model::HydroNet_Components() {

	int edge_count, start_obj, path_count, obj_count, row_aux, column, prev_row, n_paths, branch_count, n_new_paths;
	std::vector<std::vector<bool>> edges, mat_aux;
	std::vector<int> path, branch_aux, vec_aux1, vec_vec_aux, path_aux;
	std::vector<std::vector<int>>paths_ind, mesh_ind,  branch_mesh, objects_branches,  ind_outlet_aux;
	bool flag_path, first;
	std::vector<bool>path_cycle, to_remove, vec_aux;
	std::string node = "node", tank = "tank";
	/*std::vector<int> tanks_id = {};
	for (int count=0;count<tank.size();count++)
	tanks_id.push_back(tank.at(count).tank)*/


	// Obtains all information about a network and its components - meshes,
	// nodes and branches.

	// Stores all nodes in the network(index and id)
	// Finds all paths using a Depth - first search algorithm.
	// Checks whether paths are closed cycles(mesh).
	// Stores all objects in each mesh(index and id)
	// Finds all branches(paths without bifurcations crossed by a single flow)
	// Stores all objects in each branch(index and id)
	// Checks whether branches form part of closed cycles(meshes).
	// Finds branches that form part of each mesh
	// Finds branches connected to each node

	// Inputs are the table of objects, tables for classes valve and thermostat
	// (for checking opening) and temperature at the inlet of thermostats.

	// Outputs are :
	//   Nodes in the network(object indices)
	// Nodes in the network(object id)
	// Objects in every mesh(object indices)
	// Objects in every mesh(object id)
	// Objects in every branch(object indices)
	// Objects in every branch(object id)
	// Logical vector that marks whether branches form part of closed cycles
	//   Branches that form part of each mesh
	// 	Branches connected to each node
	//   Number of meshes
	//   Number of branches


	
	//// EDGES
	// Looks for all connections among objects(indirected).

	edge_count = 0; // counter of edges
	edges.resize(n_obj);
	for (int count = 0; count < n_obj; count++)
		edges.at(count).assign(n_obj, false); // matrix that contains false if there is no connection between

											  // objects and true if they are connected(rows and columns are the indices of objects)

											
	for (int obj = 0; obj < n_obj; obj++)				// Creates a matrix n_obj * n_obj indicating if the objects are adjacent
		for (int adj = 0; adj <objects.at(obj).adjacent.size(); adj++)
			edges.at(objects.at(obj).ID).at(objects.at(obj).adjacent.at(adj)) = true;


	//edges count
	for (int row = 0; row < edges.size(); row++)
		for (int col = 0; col < edges.at(row).size(); col++)
			if (edges.at(row).at(col) == true)
				edge_count++;

	edge_count = edge_count / 2;

	
	//// PATHS

	// Selects starting object : first node of the list
	if (nodes_id.size() == 0) { // no nodes found
		start_obj = 1; // index of first object in the list
	}
	else {
		start_obj = nodes_id.at(0); // index of first node
	}


	// Finds all paths using a Depth - first search algorithm.

	path_count = 0; // counter of paths
	mat_aux = edges; // auxiliary matrix to mark visited objects
	path = { start_obj };// Path initialitacio
	obj_count = 1; // counter of objects along the path
	row_aux = start_obj; // initialize row(starting object)
	column = -1;


	
	while (!(row_aux == start_obj && column == n_obj)) {
		// while there are unexplored paths from the starting object(did not reach the end of the starting row)
		
		column = column + 1; // next column

		if (column<n_obj && mat_aux.at(row_aux).at(column)) { // there is a connection

			mat_aux.at(row_aux).at(column) = false; // mark edge as explored
			prev_row = row_aux;
			row_aux = column; // move to connection's row
			flag_path = false; // path not finished
			obj_count = obj_count + 1; // next object along the path
			path.push_back(row_aux); // store index of object in the path
			mat_aux.at(row_aux).at(prev_row) = false; // going back the same way is not allowed


			path_aux = path;
			path_aux.pop_back();


			if (!any(mat_aux.at(row_aux))) // there are no more ways to go
				flag_path = true; // end path
			

			else if (in(row_aux, path_aux)) // the new object was already in the path - the cycle is closed
				flag_path = true; // end path
								  

			if (flag_path) {// path has ended
				path_count = path_count + 1; // next path
				paths_ind.push_back(path); // save path
				path.pop_back(); // remove last object from path - step back
				obj_count = obj_count - 1; // previous object along the path
				if (row_aux == start_obj) // the starting row never gets fully re - initialized
					mat_aux.at(row_aux).at(prev_row) = true; // just allow the same way back for future paths
				else
					mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object

				row_aux = path.at(path.size() - 1); // return to upper object

				column = -1;
			}
			else // path goes further
				column = -1; // stay in the new row and start from the beginning
		}

		while ((column >= n_obj - 1) && (row_aux != start_obj)) { // loop to step back in a path
																  // end of the row and it is not the starting object // skip fully explored objects
			mat_aux.at(row_aux) = edges.at(row_aux); // re - initialize connections of object
			column = path.at(path.size() - 1); // continue the upper object where it was left
			path.pop_back(); // remove last object from path - step back
			obj_count = obj_count - 1; // previous object along the path
			row_aux = path.at(path.size() - 1); // return to upper object
			if (path.size() > 1)
				mat_aux.at(row_aux).at(path.at(path.size() - 2)) = false; // the previous object along the path is not allowed

		}
	}

	n_paths = path_count; // number of paths
	paths_ind.resize(n_paths);
	path.shrink_to_fit();// adjust size
	path_aux.clear();
	path_aux.shrink_to_fit();
	path_aux = {};



	// Find and isolate closed paths(cycles)
	path_cycle.assign(n_paths, false);
	// vector that indicates whether the path forms a mesh(closed cycle)
	for (int count_path = 0; count_path < n_paths; count_path++) {// path loop
																  // Criterium to know if the paths are cycles : the last object and one of the other objects in the path are the same
																  // Two cases : pure cycle(first object and last object are the same)
																  // branched cycle(first objects are not a part of the cycle)
		if (paths_ind.at(count_path).at(0) == paths_ind.at(count_path).at(paths_ind.at(count_path).size() - 1)) // pure cycle
			path_cycle.at(count_path) = true;
		else {
			for (int count_obj = 1; count_obj < paths_ind.at(count_path).size() - 3; count_obj++) {// minus first, last and next - to - last
				if (paths_ind.at(count_path).at(count_obj) == paths_ind.at(count_path).at(paths_ind.at(count_path).size() - 1)) { // cycle starts here
					path_cycle.at(count_path) = true;
					// Create new path for the branch that is out of the cycle
					first = true;
					for (int aux_count_obj = 0; aux_count_obj < count_obj + 1; aux_count_obj++) {
						if (first) {
							first = false;
							paths_ind.push_back({});
						}

						paths_ind.at(paths_ind.size() - 1).push_back(paths_ind.at(count_path).at(aux_count_obj));
					}
					// Convert branched cycle in a pure cycle
					for (int aux_count_obj = count_obj; aux_count_obj < paths_ind.at(count_path).size(); aux_count_obj++)
						path_aux.push_back(paths_ind.at(count_path).at(aux_count_obj));
					paths_ind.at(count_path) = path_aux;
					path_aux.clear();
					path_aux.shrink_to_fit();
					path_aux = {};
					break; // path processed; go to next
				}
			}
		}
	}
	n_paths = paths_ind.size(); // refresh number of paths
	n_new_paths = n_paths - path_cycle.size();
	for (int count = 0; count <n_new_paths; count++)
		path_cycle.push_back(false);


	// Removes repeated paths and cycles(all cycles are duplicated once - two directions)
	to_remove.assign(n_paths, false); // stores paths to be removed
	for (int count = 0; count < n_paths; count++) { // path loop
		if (to_remove.at(count) == false) { // the path is not marked to be deleted
			for (int count1 = count + 1; count1 < n_paths; count1++) {// compare with next paths
				if (paths_ind.at(count1).size() == paths_ind.at(count).size()) {
					vec_aux.assign(paths_ind.at(count).size(), false);
					for (int count_obj = 0; count_obj < paths_ind.at(count).size(); count_obj++) {
						if (in(paths_ind.at(count).at(count_obj), paths_ind.at(count1))) {
							vec_aux.at(count_obj) = true;
						}
					}
					if (all(vec_aux))// all objects are present
						to_remove.at(count1) = true; // mark to be removed
				}
			}
		}
	}

	for (int count = 0, aux = 0; count < to_remove.size(); count++)
		if (to_remove.at(count) == true) {

			paths_ind.erase(paths_ind.end() - (paths_ind.size() - count + aux)); // remove repeated paths and cycles
			paths_ind.shrink_to_fit();
			path_cycle.erase(path_cycle.end() - (path_cycle.size() - count + aux));
			path_cycle.shrink_to_fit();
			aux++;
		}
	n_paths = paths_ind.size(); // refresh number of paths

	if (n_paths == 1)
		paths_ind.at(0).pop_back();

	
								//// MESHES(cycles)
	for (int count = 0; count < paths_ind.size(); count++) {
		if (path_cycle.at(count) == true)
			mesh_ind.push_back(paths_ind.at(count)); // paths that are cycles(indices)
	}

	
	//// BRANCHES
	// Paths without bifurcations crossed by a single flow.
	// Branches that are not included in a mesh will be obtained too because they may be
	// traversed by a flow if there is a boundary condition at one end(flow or head).

	branches_ind.resize(edge_count * n_paths); // stores unique branches of all paths(indices)
											   // maximum size : number of edges * number of paths
	branch_cycle.assign(edge_count * n_paths, false); // stores whether the branch is inside a mesh
	branch_mesh.resize(edge_count * n_paths); // stores the mesh where the branch is contained(indices)
											  // is a cell because it will store all meshes that contain that branch
											  // if the branch is not inside a mesh, the path is not stored
											  // it is used to make a cell later that stores branches in each mesh


											  // Find branches in paths by splitting them at nodes and tanks
	branch_count = 0; // branch counter
	branch_aux = {};
	for (int count = 0, count1; count < n_paths; count++) { // path loop
		start_obj = 0; // starting index for the first branch
		for (count1 = 1; count1 < (paths_ind.at(count).size() - 1); count1++) { // counter of objects in a path excluding first and last
			if ( in(paths_ind.at(count).at(count1), nodes_id)/* || in(paths_ind.at(count).at(count1), tanks_id) */) {
				branch_count = branch_count + 1; // next branch

				for (int aux = start_obj; aux <= count1; aux++)
					branch_aux.push_back(paths_ind.at(count).at(aux)); // save branch

				branches_ind.at(branch_count - 1) = branch_aux;
				branch_aux.clear();
				branch_aux.shrink_to_fit();
				branch_aux = {};
				start_obj = count1; // last object of the present branch is the first object of the next one
				branch_cycle.at(branch_count - 1) = path_cycle.at(count); // save whether the branch is inside a mesh
				if (path_cycle.at(count) == true) // path is a mesh
					for (int aux = 0; aux < path_cycle.size(); aux++)
						if (path_cycle.at(aux) == true)
							if (aux == count) {
								branch_mesh.at(branch_count - 1) = { count }; // save origin mesh(index)
								break;
							}

			}
		}
		if (paths_ind.at(count).size() == 2) {// branches with only 2 objects
			count1 = -1; // needed because the loop above was avoided
		}
		count1 = count1 + 1; // next object // count1 = end of path
		branch_count = branch_count + 1; // next branch
		for (int aux = start_obj; aux < count1; aux++)
			branch_aux.push_back(paths_ind.at(count).at(aux)); // save branch
		branches_ind.at(branch_count - 1) = branch_aux;
		branch_aux.clear();
		branch_aux.shrink_to_fit();
		branch_aux = {};
		branch_cycle.at(branch_count - 1) = path_cycle.at(count); // save whether the branch is inside a mesh
		if (path_cycle.at(count) == true) {// path is a mesh
			for (int aux = 0; aux < path_cycle.size(); aux++)
				if (path_cycle.at(aux) == true)
					if (aux == count) {
						branch_mesh.at(branch_count - 1) = { count }; // save origin mesh(index)
						break;
					}
		}
	}
	branches_ind.resize(branch_count); // adjust size
	branch_cycle.resize(branch_count);
	branch_mesh.resize(branch_count);
	branches_ind.shrink_to_fit();
	branch_cycle.shrink_to_fit();
	branch_mesh.shrink_to_fit();

	// Removes repeated branches
	to_remove.clear();
	to_remove.shrink_to_fit();
	to_remove.assign(branch_count, false); // stores branches to be removed
	for (int count = 0; count < branch_count; count++) { // branch loop
		if (to_remove.at(count) == false) { // the branch is not marked to be deleted
			for (int count1 = count + 1; count1 < branch_count; count1++) {// compare with next branches
				if (branches_ind.at(count1).size() == branches_ind.at(count).size()) {
					// to compare for equality, vectors must have the same length
					if (branches_ind.at(count1) == branches_ind.at(count) || branches_ind.at(count1) == flip(branches_ind.at(count))) {
						// branches are equal OR branch is equal to the inverse branch
						to_remove.at(count1) = true; // mark to be removed
						if (branch_cycle.at(count) == false && branch_cycle.at(count1) == true) {
							// first branch is not in a cycle but second(and same) branch is
							branch_cycle.at(count) = true; // branch is inside a cycle
						}
						for (int aux = 0; aux<branch_mesh.at(count1).size(); aux++)
							branch_mesh.at(count).push_back(branch_mesh.at(count1).at(aux)); // add mesh
					}
				}
			}
		}
	}
	for (int count = 0, aux = 0; count < to_remove.size(); count++)
		if (to_remove.at(count) == true) {
			branches_ind.erase(branches_ind.end() - (branches_ind.size() - count + aux)); // remove repeated branches
			branches_ind.shrink_to_fit();
			branch_cycle.erase(branch_cycle.end() - (branch_cycle.size() - count + aux));
			branch_cycle.shrink_to_fit();
			branch_mesh.erase(branch_mesh.end() - (branch_mesh.size() - count + aux));
			branch_mesh.shrink_to_fit();
			aux++;
		}
	branch_count = branches_ind.size(); // refresh number of branches


										// Finds and merges split branches.
										// If there are no nodes, the starting object may be inside a mesh.In that case,
										// it is in the middle of a branch but the algorithms have splitted it in two parts.
	if (nodes_id.size() == 0) { // no nodes found(first condition)
		for (int count = 0; count < branch_count; count++) {// branch loop
			if (branches_ind.at(count).at(0) == 1 || branches_ind.at(count).at(branches_ind.at(count).size() - 1) == 1) {
				// starting object is in the branch(its index is 1)
				if (branch_cycle.at(count) == true) { // is in a mesh(second condition)
					for (int count1 = count+1, aux_size; count1 < branch_count; count1++) {// looks for the other side
						if (branches_ind.at(count1).at(0) == 1 || branches_ind.at(count1).at(branches_ind.at(count1).size() - 1) == 1) {
							// Join branches
							if (branches_ind.at(count).at(branches_ind.at(count).size() - 1) == 1 && branches_ind.at(count).at(0) == 1) {
								// starting object is the last of first branch and the first of the second brandlch
								aux_size = branches_ind.at(count1).size();
								for (int aux = 1; aux < aux_size; aux++)
									branches_ind.at(count).push_back(branches_ind.at(count1).at(aux)); // join directly // do not repeat starting object
							}
							else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size() - 1) == 1) {
								// starting object is the last of first branch and the last of the second branch
								aux_size = branches_ind.at(count1).size();
								for (int aux = 1; aux < aux_size; aux++)
									branches_ind.at(count).push_back(branches_ind.at(count1).at(branches_ind.at(count1).size() - aux)); // flip second branch
							}
							else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size() - 1) == 1) {
								// starting object is the first of first branch and the first of the second branch
								aux_size = branches_ind.at(count1).size();
								for (int aux = 0; aux < aux_size; aux++) {
									branch_aux.push_back(branches_ind.at(count).at(branches_ind.at(count).size() - 1 - aux));
								}
								for (int aux = 0; aux < aux_size; aux++) {
									branch_aux.push_back(branches_ind.at(count1).at(aux));
								}
								branches_ind.at(count) = branch_aux; // flip first branch

							}

							else if (branches_ind.at(count).at(1) == 1 && branches_ind.at(count).at(branches_ind.at(count).size() - 1) == 1) {
								// starting object is the first of first branch and the last of the second branch
								branch_aux = branches_ind.at(count1);
								aux_size = branches_ind.at(count1).size();
								for (int aux = 1; aux < aux_size; aux++) {
									branch_aux.push_back(branches_ind.at(count).at(aux));
									branches_ind.at(count) = branch_aux; // join second branch and then first

									branches_ind.erase(branches_ind.begin() + count1); // remove second branch
									branch_cycle.erase(branch_cycle.cbegin() + count1);
									branch_mesh.erase(branch_mesh.begin() + count1);
									branch_count = branch_count - 1; // refresh number of branches
									break; // all done; exit algorithm
								}
							}
						}
					}
					break;// exit algorithm
				}
			}
		}
	}

	// Stores branches connected to each node
	node_branches.resize(nodes_id.size());
	for (int count = 0; count< nodes_id.size(); count++) { // counter of nodes
		for (int count1 = 0; count1< branch_count; count1++) { // branch loop
			if (nodes_id.at(count) == branches_ind.at(count1).at(0) || nodes_id.at(count) == branches_ind.at(count1).at(branches_ind.at(count1).size() - 1)) {
				// branch starts or ends at the node
				node_branches.at(count).push_back(count1); // add index of branch
			}
		}
	}



	// Stores branches that form each mesh
	n_mesh = mesh_ind.size();
	mesh_branches.resize(n_mesh);
	for (int count = 0; count < n_mesh; count++)
		mesh_branches.at(count) = {};

	for (int count = 0; count < n_mesh; count++)// mesh loop
		for (int count1 = 0; count1 < branch_count; count1++)// branch loop
			if (in(count, branch_mesh.at(count1)))// the branch belongs to the current mesh
				mesh_branches.at(count).push_back(count1); // add index of accordant branch


														   // Re - order branches
	for (int count_mesh = 1; count_mesh < n_mesh; count_mesh++) { // mesh loop
		if (mesh_branches.at(count_mesh).size() > 2) {// if there are only one or two branches, the current order is accepted
			vec_aux1.assign(mesh_branches.at(count_mesh).size(), 0); // stores the ordered branches in the mesh
			vec_aux1.at(1) = mesh_branches.at(count_mesh).at(1); // the first branch it finds is chosen as the first of the mesh
																 // Looks for the following branches
			for (int count_branch1 = 0, ind_aux; count_branch1 < vec_aux1.size() - 1; count_branch1++) {// after last branch comes first branch again(cycle)
				ind_aux = vec_aux1.at(count_branch1); // current branch index
				for (int count_branch2 = 1, ind_aux1; count_branch2 < mesh_branches.at(count_mesh).size(); count_branch2++) { // branches to compare with // first branch cannot be
					ind_aux1 = mesh_branches.at(count_mesh).at(count_branch2); // index of branch to compare
					if (ind_aux != ind_aux1) {// does not ckeck itself
						if (branches_ind.at(ind_aux1).at(0) == branches_ind.at(ind_aux).at(branches_ind.at(ind_aux).size() - 1)) {
							// second branch starts with the last object of first branch
							vec_aux1.at(count_branch1 + 1) = ind_aux1; // is the next element
							break; // next branch
						}
						else if (branches_ind.at(ind_aux1).at(branches_ind.at(ind_aux1).size() - 1) == branches_ind.at(ind_aux).at(branches_ind.at(ind_aux).size() - 1)) {
							// second branch ends with the last object of first branch
							branch_aux.clear();
							for (int aux = 0; aux < branches_ind.at(ind_aux1).size(); aux++)
								branch_aux.push_back(branches_ind.at(ind_aux1).at(branches_ind.at(ind_aux1).size() - 1 - aux));

							branches_ind.at(ind_aux1) = branch_aux; // inverts branch
							vec_aux1.at(count_branch1 + 1) = ind_aux1; // is the next element
							break; // next branch
						}
					}
				}
			}

			mesh_branches.at(count_mesh) = vec_aux1; // ordered branches in the mesh
		}
	}


	// Adds every object's branch to the objects table
	objects_branches.resize(n_obj);
	for (int count = 0; count < n_obj; count++)
		objects_branches.at(count) = {};

	for (int count = 0; count < branch_count; count++) // branch loop
		for (int count1 = 0, aux; count1 < branches_ind.at(count).size(); count1++) { // objects in a branch loop
			aux = branches_ind.at(count).at(count1);
			objects_branches.at(aux).push_back(count);
		}


	// Converts branches to change indices by id's
	branches_id.resize(branch_count); // store all branches(id)
	for (int count = 0; count < branch_count; count++)
		branches_id.at(count) = {};
	for (int count = 0; count < branch_count; count++) // branch loop
		for (int count1 = 0; count1 < branches_ind.at(count).size(); count1++) // objects in a branch
			branches_id.at(count).push_back(objects.at(branches_ind.at(count).at(count1)).ID); // assign id
																							   //branches_id = branches_ind;

																							   // Inverts branches that are in the wrong direction
																							   // Flow direction is determined in some objects(and their branches) :
																							   // pumps and heat exchangers of type 'T out'
	for (int aux = 0, i; aux < pump_volum.size() + pump_turbo.size()  + heat_exch.size() /*+tank.size()*/; aux++) {
		if (aux < pump_volum.size()) {
			i = aux;
			ind_outlet_aux.push_back({ pump_volum.at(i).Pump_volum,pump_volum.at(i).outlet_object });
		}
		else if (pump_turbo.size()!=0 && pump_volum.size() < aux < pump_turbo.size() ) {
			i = aux - (pump_volum.size() - 1);
			ind_outlet_aux.push_back({ pump_turbo.at(i).Pump_turbo,pump_turbo.at(i).outlet_object });
		}
		else if (heat_exch.size()!=0 &&pump_volum.size()  + pump_turbo.size() < aux <heat_exch.size() ) {
			if (heat_exch.at(i).type == "T_out") {
				i = aux - (pump_volum.size() - 1 + pump_turbo.size() - 1);
				ind_outlet_aux.push_back({ heat_exch.at(i).heat_exch,heat_exch.at(i).outlet_object });
			}
		}
		/*else if (pump_volum.size() - 1 + pump_turbo.size() - 1 +heat_exch_Tout.size() - 1 < aux < tank.size()) {
		i = aux - (pump_volum.size() - 1 + pump_turbo.size() - 1);
		ind_outlet_aux.push_back({ tank.at(i).tank,tank.at(i).outlet_object });
		}*/
	}
	
	for (int count = 0, ind_aux, ind_aux1; count< ind_outlet_aux.size(); count++) { // counter of objects that have a direction
		branch_aux = objects_branches.at(ind_outlet_aux.at(count).at(0)); // branch that contains the object
		for (int aux = 0; aux<branches_ind.at(branch_aux.at(0)).size(); aux++)
			if (branches_ind.at(branch_aux.at(0)).at(aux) == ind_outlet_aux.at(count).at(0)) {
				ind_aux = aux; // position of the object in the branch
				break;
			}
		for (int aux = 0; aux<branches_ind.at(branch_aux.at(0)).size(); aux++)
			if (branches_ind.at(branch_aux.at(0)).at(aux) == ind_outlet_aux.at(count).at(1)) {
				ind_aux1 = aux; // position of the outlet object in the branch
				break;
			}																	// Compare direction of branch and object
		if (ind_aux1 < ind_aux) { // the outlet is located earlier in the branch // opposite direction
								  // this works for tanks in the return branch because(empty < x) == false
			for (int aux = 0; aux < branches_id.at(branch_aux.at(0)).size(); aux++)
				vec_vec_aux.push_back(branches_id.at(branch_aux.at(0)).at(branches_id.at(branch_aux.at(0)).size() - 1 - aux));
			branches_id.at(branch_aux.at(0)) = vec_vec_aux; // invert branch
			vec_vec_aux.clear();
			vec_vec_aux.shrink_to_fit();
			vec_vec_aux = {};
			for (int aux = 0; aux < branches_ind.at(branch_aux.at(0)).size(); aux++)
				vec_vec_aux.push_back(branches_ind.at(branch_aux.at(0)).at(branches_ind.at(branch_aux.at(0)).size() - 1 - aux));
			branches_ind.at(branch_aux.at(0)) = vec_vec_aux;
			vec_vec_aux.clear();
			vec_vec_aux.shrink_to_fit();
			vec_vec_aux = {};
		}
	}

	// branches vector contruction

	branches.resize(branches_id.size());

	for (int count = 0; count < branches.size(); count++) {
		branches.at(count).branch_ind = count;
		branches.at(count).objects.resize(branches_id.at(count).size());
		for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++) {
			branches.at(count).objects.at(count1).ID = branches_id.at(count).at(count1);
			for (int aux = 0; aux < objects.size(); aux++) {
				if (objects.at(aux).ID == branches.at(count).objects.at(count1).ID) {
					objects.at(aux).branch = count;
					branches.at(count).objects.at(count1).Class = objects.at(aux).Class;
					branches.at(count).objects.at(count1).volume = objects.at(aux).volume;
					branches.at(count).objects.at(count1).name = objects.at(aux).name;
					branches.at(count).objects.at(count1).adjacent = objects.at(aux).adjacent;

				}
			}
		}
	}


	n_branch = branches.size();

}

void HydroNet_model::HydroNet_Executions() {

	bool solver_flag = false;
	double head_aux, pos_aux, temperature_aux, obj_inlet_pos_aux;
	

	//// HEAD LOSS IN BRANCHES
	// Calculates head loss in every branch

	if (flag_first_exec || flag_valve_change || flag_thermostat_changes)

		HydroNet_HeadLoss();




	//// FLOWS IN BRANCHES and HEAD IN VOLUMETRIC PUMPS

	// Sets boundary conditions: flow = 0
	// - All branches that are not inside a closed cycle
	// - All branches that contain closed valves or thermostats

	bound_flows.assign( n_branch,NAN); // vector filled with NaN to store boundary conditions
	for (int count = 0; count < bound_flows.size(); count++)
		if (!branch_cycle.at(count))
			bound_flows.at(count) = 0;
	
	for (int count = 1, obj_aux, branch_aux; count < valve_var.size(); count++) {
		if (valve_var.at(count).opening == 0) {
			obj_aux = valve_var.at(count).valve_var_id;
			for (int count1 = 0; count1 < branches_ind.size(); count1++)
				for (int count2 = 0; count2 < branches_ind.at(count1).size(); count2++)
					if (branches_ind.at(count1).at(count2) == obj_aux) {
						branch_aux = count1;
						break;
						break;
					}
			
			bound_flows.at(branch_aux) = 0;
		}
	}
	for (int count = 1, obj_aux, branch_aux; count < thermostat.size(); count++) {
		if (thermostat.at(count).opening == 0) {
			obj_aux = thermostat.at(count).themostat_id;
			for (int count1 = 0; count1 < branches_ind.size(); count1++)
				for (int count2 = 0; count2 < branches_ind.at(count1).size(); count2++)
					if (branches_ind.at(count1).at(count2) == obj_aux) {
						branch_aux = count1;
						break;
						break;
					}
			
			bound_flows.at(branch_aux) = 0;
		}
	}


	// Loop to solve flows
	if (flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes) {
		solver_flag = true;

		while (solver_flag) {

			solver_flag = false;

			HydroNet_FlowSolver();

			// Checks whether the obtained head of volumetric pumps are higher than the maximum
			// head of each volumetric pump. In such case, assign flow = 0 as boundary
			// condition and call the solver again
			for (int count = 0, branch_aux; count < pump_volum.size(); count++) {// counter of volumetric pumps
				branch_aux = objects.at(pump_volum.at(count).Pump_volum).branch; // pump branch
				pump_volum.at(count).pump_head = head_vol_pump.at(branch_aux); // store in table
				if (abs(pump_volum.at(count).pump_head) > pump_volum.at(count).head_max) {
					// calculated head is higher than maximum head
					bound_flows.at(branch_aux) = 0; // assign flow = 0 as boundary condition
					solver_flag = true; // re-calculate flows
				}
			}
		}



		// When needed, branch directions are flipped so that flow and branch directions are the same
		for (int count_branch = 0; count_branch< flows.size(); count_branch++)
			if (flows.at(count_branch) < 0) {
				branches_ind.at(count_branch) = flip(branches_ind.at(count_branch));
				branches_id.at(count_branch) = flip(branches_id.at(count_branch));
				branch_temp_pos.at(count_branch) = flip(branch_temp_pos.at(count_branch));
				for (int count1 = 0; count1<branch_temp_pos.at(count_branch).size(); count1++)
					branch_temp_pos.at(count_branch).at(count1) = 100 - branch_temp_pos.at(count_branch).at(count1);
				branch_temperature.at(count_branch) = flip(branch_temperature.at(count_branch));
				branch_htx.at(count_branch) = flip(branch_htx.at(count_branch));
				/*branch_htx_Tout.at(count_branch)
					= flip(branch_htx_Tout.at(count_branch));*/

				for (int count = 0, count_obj; count < branches_ind.at(count_branch).size(); count++) {
					std::vector<double> obj_inlet_pos_aux;
					count_obj = branches_ind.at(count_branch).at(count);
					obj_inlet_pos_aux = obj_inlet_pos.at(count_obj); // saves this vector before overwriting it
					for (int count1 = 0; count1<obj_inlet_pos.at(count_obj).size(); count1++)
						obj_inlet_pos.at(count_obj).at(count1) = 100 - obj_outlet_pos.at(count_obj).at(count1); // oulets are now inlets
					for (int count1 = 0; count1<obj_outlet_pos.at(count_obj).size(); count1++)
						obj_outlet_pos.at(count_obj).at(count1) = 100 - obj_inlet_pos_aux.at(count1); // inlets are now outlets
				}
			}
		for (int count = 0; count < flows.size(); count++)
			flows.at(count) = abs(flows.at(count)); // All flows are positive now

	}



	//// PUMP POWER

	if (flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes) {

		// Volumetric pumps
		for (int count = 0, branch_aux, object_aux; count < pump_volum.size(); count++) { // counter of volumetric pumps
			object_aux = pump_volum.at(count).Pump_volum; // index of pump in the object's list
			branch_aux = objects.at(object_aux).branch; // branch that contains the pump (only one)
			pos_aux = obj_inlet_pos.at(branch_aux).at(object_aux) + obj_outlet_pos.at(branch_aux).at(object_aux) / 2;
			temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos.at(branch_aux), branch_temperature.at(branch_aux));
			// temperature at the pump's middle

			pump_volum.at(count).pump_power = density(fluid_type, temperature_aux) * gravity_acc * head_vol_pump.at(branch_aux) * flows.at(branch_aux) / efficiency(head_vol_pump.at(branch_aux), flows.at(branch_aux));
		}

		// Turbopumps
		for (int count = 0, object_aux, branch_aux; count < pump_turbo.size(); count++) { // counter of turbopumps
			object_aux = pump_turbo.at(count).Pump_turbo; // index of pump in the object's list
			branch_aux = objects.at(object_aux).branch; // branch that contains the pump (only one)

			pos_aux = obj_inlet_pos.at(branch_aux).at(object_aux) + obj_outlet_pos.at(branch_aux).at(object_aux) / 2;
			temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos.at(branch_aux), branch_temperature.at(branch_aux));
			// temperature at the pump's middle

			// Head curve of pump: head(speed, flow)
			head_aux = (pump_turbo.at(count).coef_N0 + pump_turbo.at(count).coef_N1 * pump_turbo.at(count).pump_speed + pump_turbo.at(count).coef_N2 *Math_wamH::pow2(pump_turbo.at(count).pump_speed)) + flows.at(branch_aux) * (pump_turbo.at(count).Q_coef_N0 + pump_turbo.at(count).Q_coef_N1 * pump_turbo.at(count).pump_speed + pump_turbo.at(count).Q_coef_N2 *Math_wamH::pow2(pump_turbo.at(count).pump_speed)) + Math_wamH::pow2(flows.at(branch_aux)) * (pump_turbo.at(count).Q2_coef_N0 + pump_turbo.at(count).Q2_coef_N1 * pump_turbo.at(count).pump_speed + pump_turbo.at(count).Q2_coef_N2 * Math_wamH::pow2(pump_turbo.at(count).pump_speed));

			pump_turbo.at(count).pump_power = density(fluid_type, temperature_aux) * gravity_acc * head_aux * flows.at(branch_aux) / efficiency(head_vol_pump.at(branch_aux), flows.at(branch_aux));
		}
	}



	//// TEMPERATURES

	HydroNet_Temperature();


}

double HydroNet_model::HydroNet_GetObjTemperature(double obj_pos, std::vector<double> branch_temp_pos, std::vector<double> branch_temperature) {

	double distance = -1000.0;
	double temperature = -500.0;
	int count = 1;


	//Returns temperature at the specified position

	// Index of the temperature position closest to the specified positionand the difference of position between them

	if (branch_temp_pos.size() == 2)
		temperature = branch_temperature.at(0);
	else {

		if (obj_pos == branch_temp_pos.at(1))
			temperature = branch_temperature.at(1);
		else if (obj_pos == branch_temp_pos.at(branch_temp_pos.size() - 1))
			temperature = branch_temperature.at(branch_temperature.size() - 1);
		else {
			while (distance < 0.0) {
				distance = branch_temp_pos.at(count + 1) - obj_pos;		// count+1 because the fert positin is checked before
				count++;
			}
			temperature = branch_temperature.at(count - 1);				// -1 because vector estarts at 0
		}
	}
	return temperature;


}

HydroNet_model::strInsertVolum HydroNet_model::HydroNet_InsertVolume(std::vector<double> positions, std::vector<double> temperatures, double start_pos, double end_pos, int temp_action, std::string fluid_type, double value, double volume, int max_divisions, int size_position, int size_temperature) {

	strInsertVolum InsertVolum;

	int start_pos_ind;
	int end_pos_ind;
	int last_element;
	int count = 0;
	int count_pos;
	double temp_to_insert;
	double mass_start;
	double mass_end;
	double temp_mass_end;
	double mass_sum = 0;
	double temp_mass_sum = 0;
	double temp_mass_start;
	bool flag_start_exists = false;
	bool flag_end_exists = false;
	bool flag_inserted = false;
	std::vector<double> _temperatures;
	std::vector<double> _positions;
	strSize _size;
	strSize size;
	size.position = size_position;
	size.temperature = size_temperature;
	InsertVolum.temperature = temperatures;
	InsertVolum.position = positions;
	InsertVolum.temp_size = size.temperature;
	InsertVolum.pos_size = size.position;


	_temperatures.resize(2 * max_divisions);
	_positions.resize(2 * max_divisions);

	_size.position = 0;
	_size.temperature = 0;



	// Obtains range of positions where a volume is to be placed inside a branch and replaces current volumes by the new volume.In addition, a temperature for the volume is calculated.
	// Fluid_type is only mandatory for temp_action = 0 and 1 (must calculate outlet temperature).
	// Value is only mandatory for temp_action = 1 and 2 (not just averaging).
	// Volume is only mandatory for temp_action = 1 (impose heat).

	if (start_pos == end_pos)
		return InsertVolum;					// to avoid dividing by zero


	// Finds start position inside positions vector

	while ((count < size.position) && (start_pos > positions.at(count))) {

		count++;
	}
	if (start_pos == positions.at(count)) {      // Start position allready exixt in possitions
		flag_start_exists = true;
		start_pos_ind = count;
	}
	else
		start_pos_ind = count - 1;				// Position of the highest value in the vector position lower than start_position

												// Finds end position inside positions vector

	count = 0;

	while ((count < size.position) && (end_pos > positions.at(count))) {

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
			mass_start = density(fluid_type, temperatures.at(start_pos_ind )) * (end_pos - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind) * mass_start;
		}
		else {
			mass_start = density(fluid_type, temperatures.at(start_pos_ind )) * (positions.at(start_pos_ind) - start_pos);
			// Not a real mass.Volume is not included because it appears both in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind) * mass_start;
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
			mass_end = density(fluid_type, temperatures.at(end_pos_ind )) * (end_pos - positions.at(end_pos_ind ));
			temp_mass_end = temperatures.at(end_pos_ind) * mass_end;
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
			mass_start = density(fluid_type, temperatures.at(start_pos_ind )) * (end_pos - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind ) * mass_start;
		}

		else {
			mass_start = density(fluid_type, temperatures.at(start_pos_ind )) * (positions.at(start_pos_ind) - start_pos);
			// Not a real mass.Volume is not included because it appears both
			// in the numerator and the denominator multiplying all terms.
			temp_mass_start = temperatures.at(start_pos_ind) * mass_start;
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
			mass_end = density(fluid_type, temperatures.at(end_pos_ind )) * (end_pos - positions.at(end_pos_ind ));
			temp_mass_end = temperatures.at(end_pos_ind ) * mass_end;
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
	for (int count = 0, aux = 0; flag_start_exists == false && flag_end_exists == false && count <  size.temperature + 2
		|| flag_start_exists == true && flag_end_exists == false && count < size.temperature + 1
		|| flag_start_exists == false && flag_end_exists == true && count <  size.temperature + 1
		|| flag_start_exists == true && flag_end_exists == true && count <  size.temperature; count++) {

		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count--;
			aux--;
		}
		if (positions.at(count) < start_pos) {										// Saves values of temperatures until the new temperature
			_temperatures.at(_size.temperature) = temperatures.at(aux);
			_size.temperature++;
			aux++;
		}

		else if (flag_inserted == true) {											// Saves values of temperatures after the new temperature

			_temperatures.at(_size.temperature) = temperatures.at(aux);
			_size.temperature++;
			aux++;
		}

		else {																//Saves value of new temperature betwin two previusly known
			_temperatures.at(_size.temperature) = temp_to_insert;
			_size.temperature++;
			flag_inserted = true;
			if (flag_start_exists == true && flag_end_exists == true)
				aux++;

			if (flag_start_exists == false && flag_end_exists == false) {
				_temperatures.at(_size.temperature) = temperatures.at(aux - 1);
				_size.temperature++;
			}
		}
		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count++;
			aux++;
		}

	}
	

	//Inserts specified positions and removes intermediate positions in the branch

	flag_inserted = false;

	for (int count = 0, aux = 0; flag_start_exists == false && flag_end_exists == false && count <  size.position + 2
		|| flag_start_exists == true && flag_end_exists == false && count <  size.position
		|| flag_start_exists == false && flag_end_exists == true && count <  size.position
		|| flag_start_exists == true && flag_end_exists == true && count <  size.position; count++) {

		if (flag_start_exists == false && flag_end_exists == false && flag_inserted == true) {
			count--;

		}

		if (flag_inserted != true && positions.at(count) < start_pos) {										// Saves values of positions until the new temperature
			_positions.at(_size.position) = positions.at(aux);
			_size.position++;
			aux++;
		}

		else if (flag_inserted == true) {											// Saves values of positions after the new temperature


			_positions.at(_size.position) = positions.at(aux);
			_size.position++;
			aux++;
		}

		else {																//Saves value of new temperature betwin two previusly known


			_positions.at(_size.position) = start_pos;
			_size.position++;
			_positions.at(_size.position) = end_pos;
			_size.position++;
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
	

	InsertVolum.temperature = _temperatures;
	InsertVolum.position = _positions;
	InsertVolum.temp_size = _size.temperature;
	InsertVolum.pos_size = _size.position;

	return InsertVolum;
}

void HydroNet_model::HydroNet_FlowSolver() {

	Eigen::RowVectorXi solved(n_branch);
	std::vector<bool> trues;
	Eigen::RowVectorXi volum_pump(n_branch);
	Eigen::RowVectorXi solv_nodes(nodes_id.size());
	bool while_flag;
	Eigen::RowVectorXi solv_prev1(n_branch);
	Eigen::RowVectorXi solv_prev2(n_branch);
	int count_uns;
	int ind_aux;
	double flow_sum;
	double total_head;
	int last_id;
	bool change_sign;
	bool sign_aux;
	std::vector<double> _sol_aux;
	double sol_aux;
	int n_unknown;
	std::vector<std::vector<double>> coef1_mat, coef2_mat, const_mat;
	std::vector<double> const_vec;
	std::vector<int> unk_ind_aux;
	std::vector<bool> volpump_ind_aux, flows_ind_aux;
	int last_row=0;
	int n_idp_mesh;
	std::vector<int> idp_mesh;
	std::vector<bool> remain_branch;
	int best_mesh;
	int most_uns;
	Eigen::VectorXd x0;
	std::vector<std::vector<int>> idp_mesh_branches;
	bool flag;
	std::vector<double> Y;
	std::string p;
	strSize size;

	trues.assign(n_branch, true);


	
	//Solver for hydraulic networks.
	//Receives information about meshes, branches and nodes in the the network.
	//Receives known flows imposed as boundary conditions - a vector with the flow values in the corresponding branches and NaN in the rest.
	//Receives head losses in every branch.Data must be divided in three vectors: one for head losses that do not depend on the flow, another one for head losses proportional to the flow and the last one for losses proportional to the flow squared.
	//This function attempts to solve the network with algebraic methods.If some flows cannot be solved, it uses a numerical method.
	//Receives information about volumetric pumps.If there are volumetric pumps in the circuit, the fucntion returns its pump head.



	//SOLVING PROCEDURE

	//First, solves branches where flow is known, is zero or can be calculated independently because there is a volumetric pump(flow depends on pump speed).Check if all branches are solved.

	//Second, enters a loop :
	//NODES STEP : Solves nodes where there is only one unknown flow.Check if all branches are solved.Loop until no more nodes are solved.
	//If it is not the first iteration and no branches have been solved in the nodes step, break main loop.
	//MESH STEP : solves meshes that depend only on one flow(simple meshes; only one branch remains to be solved).Check if all branches are solved.Loop until no more meshes are solved.
	//If no branches have been solved in the mesh step, break main loop.
	//Loop

	//Third, sets up and solves the system of non - linear equations.Equations  depend on flow(nodes, etc.), on flow ^ 2 (pipes, etc.) and on head of volumetric pumps.All branches are solved.

	// Creates result vectors

	flows.assign(n_branch, 0);						// stores flow values
	solved.setConstant(false);					// stores solving status
	unk_ind_aux.resize(n_branch);
	remain_branch.resize(n_branch);
	size.remain_branch = 0;
	size.unk_ind_aux = 0;


	// 1. INDEPENDENT FLOWS

	//Assigns flows entered as a boundary condition

	for (int count = 0; count < n_branch; count++)
		if (0 == isnan(bound_flows.at(count))) {
			flows.at(count) = bound_flows.at(count);
			solved(count) = true;
		}


	if (solved.isOnes())																// all branches are solved
		return;																			// end function


																						// Solves closed branches containing VOLUMETRIC PUMPS

	volum_pump.setConstant(false);													// marks where are there volumetric pumps for later use

	for (int branch_aux = 0; branch_aux < branches.size(); branch_aux++)																// Enter every branch
		for (int object_aux = 0, aux = 0; object_aux < branches.at(branch_aux).objects.size(); object_aux++) {		 //Looks objects in branch
																													 //	branch_aux = objects.branch{ pump_volum.obj_index(branch_auxbranch_aux) };			
																													 // pump branch
			p = "pump_volum";
			if ((branches.at(branch_aux).objects.at(object_aux).Class == p) && (solved(branch_aux) == false)) {		//Enters when the object is a pump_volum										// the branch is not solved yet
				volum_pump(branch_aux) = true;																										// there is a volumetric pump
				aux = 0;
				while ((aux + 1 < pump_volum.size()) && (branches.at(branch_aux).objects.at(object_aux).ID != pump_volum.at(aux).Pump_volum))				// Looks up for the pump at the branch in the pump_volum vector
					aux++;
				//calculates flow

				flows.at(branch_aux) = (pump_volum.at(aux).coef0 + pump_volum.at(aux).coef1 * pump_volum.at(aux).pump_speed + pump_volum.at(aux).coef2 * pump_volum.at(aux).pump_speed * pump_volum.at(aux).pump_speed + pump_volum.at(aux).coef3 * pump_volum.at(aux).pump_speed * pump_volum.at(aux).pump_speed * pump_volum.at(aux).pump_speed)*1E-3;
				solved(branch_aux) = true; // set as solved
			}
		}

	// There is no check of whether all branches are solved because volumetric pump head must still be found
	head_vol_pump.resize(n_branch);															// vector with the head value of every volumetric pump in its corresponding branch

																							// 2. LOOP FOR NODES WITH ONE UNKNOWN FLOW

	n_nodes = nodes_id.size();
	solv_nodes.setConstant(false);														// vector to mark nodes whose all connected branches are solved

																						// Loop for solving NODES and MESHES that have only one unknown flow.
																						// At the moment the MESH sub - loop is not finished since its implementation is very complex and it seems that solving a system of equations will be faster

	while_flag = false;																	// flag to know if it is the first iteration of the while loop
	solv_prev1.resize(solved.size());												// flows solved in the previous step(main loop)

	while ((solved != solv_prev1)) {												// in the previous step, some flow was solved(NODES + MESHES)

		solv_prev1 = solved;														// flows solved after the meshes step or at start
		solv_prev2.resize(solved.size());											// flows solved in the previous step(sub - loop)

		while (solved != solv_prev2) {												// in the previous step, some flow was solved(NODES)

			solv_prev2 = solved;															// flows solved after the nodes sub loop

			for (int count = 0; count < n_nodes; count++) {
				count_uns = 0;																	 // counter of unsolved branches
				for (int count1 = 0; count1 < node_branches.at(count).size(); count1++) {			// branches connected to node
					if (solved(node_branches.at(count).at(count1)) == false) {								// branch unsolved
						ind_aux = count1;														// save index of unsolved branch
						count_uns = count_uns + 1;													// counter of unsolved branches
						if (count_uns == 2) {													 // too many unknowns
							break;													// go to next node
						}
					}
				}

				if (count_uns == 0)													// all connected branches are solved
					solv_nodes(count) = true;												// set as solved

				else if (count_uns == 1) {									// one unknown; node can be solved directly

																			// Continuity equation : sum of flows in a node = 0
					flow_sum = 0;

					for (int count1 = 0, aux; count1 < node_branches.at(count).size(); count1++) {
						aux = node_branches.at(count).at(count1);
						// branches connected to node
						// Node can be at the end or at the begining of the branch
						// Criterium: positive flow--> enters node
						if (nodes_id.at(count) == branches_id.at(aux).at(branches_id.at(aux).size() - 1))			// at the end, keep sign
																													// flow enters the node
							flow_sum = flow_sum + flows.at(aux);
						else if (nodes_id.at(count) == branches_id.at(aux).at(0))								// at the beginning, change sign
																												// flow leaves the node
							flow_sum = flow_sum - flows.at(aux);												// change sign
					}


					ind_aux = node_branches.at(count).at(ind_aux);									// change ind_aux to refer directly to the branch or flow index

																									// Node can be at the end or at the begining of the branch
					if (nodes_id.at(count) == branches_id.at(ind_aux).at(0))									// at the beginning
						flows.at(ind_aux) = flow_sum;
					else
						flows.at(ind_aux) = -flow_sum;															// change sign

					solved(ind_aux) = true;																// set branch as solved

					if (solved.isOnes() && volum_pump.isOnes())														// all branches are solved(flows and volumetric pump head)
						return;

					solv_nodes(count) = 1/*true*/;																 // set node as solved
				}
			}
		} // while for nodes

		if (while_flag == true && solved == solv_prev1)
			// it is not the first iteration and no flows have been solved during the nodes step
			// no more flows will be solved in the meshes step because there have been no changes
			break;																							// go to next solving method


		solv_prev1 = solved;																		// flows solved after the nodes step

				//////////////////////////// middle of loop

		solv_prev2.resize(solved.size());																	// flows solved in the previous step(sub - loop)

		while (solved != solv_prev2) {																		// in the previous step, some flow was solved(MESH)

			solv_prev2 = solved; // flows solved after the mesh sub - loop

			for (int count = 1; count < n_mesh; count++) {														 // meshes
				count_uns = 0;																				// counter of unsolved branches
				for (int count1 = 0; count1 < mesh_branches.at(count).size(); count1++) {										 // branches that form the mesh
					if (solved(count1) != true || volum_pump(count1)) {
						// branch unsolved or contains volumetric pump(pump head must be calculated)
						ind_aux = count1;																				// save unsolved branch
						count_uns = count_uns + 1;																			// counter of unsolved branches
						if (count_uns == 2)																		// too many unknowns
							break;																					// go to next mesh

					}
				}
				if (count_uns == 1) {																			 // one unknown; mesh can be solved directly

						 // II Kirchhoff's Law: sum of head increments in a mesh = 0
					total_head = 0;
					last_id = 0;
					for (int count1 = 0; count1 < mesh_branches.at(count).size(); count++) {								// branches that form the mesh
																															// Unsolved one must be included in the loop to find out mesh flow direction
																															// Branches are already ordered along the mesh
																															// Criterium: positive mesh flow--> direction of the first branch in the mesh
																															// If there is a pump or a heat exchanger of type 'T out' in the branch, flow and branch direction are the same
																															// If flow direction is unknown, assume same direction as branch
																															// Flow may be solved or not

																															// Check whether the branch has the same direction as the mesh flow
						change_sign = false;
						if (last_id == 0 || branches_id.at(count1).at(1) == last_id)
							// it is the first branch in the mesh OR the branch is in the same direction as the mesh flow
							last_id = branches_id.at(count1).at(branches_id.at(count1).size());								 // refresh other end of branch
						else {// it is not the first branch in the mesh AND the branch is in the opposite direction
							change_sign = true;																			 // invert sign
							last_id = branches_id.at(count1).at(1);																// refresh other end of branch
						}

						if (count1 != ind_aux || volum_pump(count1) == true) {									// flow in the branch is known
							if (flows.at(count1) < 0)																 // flow direction is opposite to branch direction
								change_sign = ~change_sign;																			// invert sign

																																	// Add head increment
							total_head = total_head + (head_loss.at(count1) + hydr_resist1.at(count1) * fabs(flows.at(count1)) + hydr_resist2.at(count1) * flows.at(count1) * flows.at(count1))* (double)(-2 * change_sign + 1);
							// add head if change_sign = false; subtract if change_sign = true
						}
						else																							// branch unknown is a flow : save sign - false : positive, true : negative,
							sign_aux = change_sign;



						// Solve unknown
						if (volum_pump(ind_aux) == true) {																	// the unknown branch has a volumetric pump
																															// Solve unknown head increment
							head_vol_pump.at(ind_aux) = abs(total_head);												 // calculate pump head
							head_loss.at(ind_aux) = head_loss.at(ind_aux) - abs(total_head);
							// consider pump head as a negative head loss in the branch
							volum_pump(ind_aux) = false;																			// remove branch from branches with volumetric pumps
						}
						else {																						// the unknown is a flow; solve it
							if (total_head < 0)																	// branch flow is in the opposite direction of the mesh flow
								sign_aux = ~sign_aux;															// invert branch flow sign

																												// Solve branch
							_sol_aux = roots(hydr_resist2.at(ind_aux), hydr_resist1.at(ind_aux), head_loss.at(ind_aux) - abs(total_head));
							// absolute value of head increment in the branch and
							// absolute value of total head must be equal
							// If there are two solutions, choose the positive one
							if (_sol_aux.size() == 2)
								if (_sol_aux.at(1) > 0)
									sol_aux = _sol_aux.at(1);
								else
									sol_aux = _sol_aux.at(2);



							// Save result
							flows.at(ind_aux) = sol_aux * sign_aux*(-1);													 // add sign now
							solved(ind_aux) = true;																			// set as solved
						}
					}

					if (solved.isOnes() && volum_pump.isOnes())															 // all branches are solved(flows and volumetric pump head)
						return;
				}
			}
		} // while for meshes

		if (solved == solv_prev1)
			// no flows have been solved during the meshes step
			// no more flows will be solved in the nodes step because there have been no changes
			break;																								// go to next solving method


		while_flag = true;																				 // first iteration finished
	}

	//// 3. SYSTEM OF EQUATIONS
	// A * X* | X | +B * X + C * X / | X | +D = 0 --> A = coef2_mat, B = coef1_mat, C = const_mat, D = const_vec
	// C is for constant head loss whose sign is unknown

	n_unknown = n_branch - solved.sum() + volum_pump.sum(); // number of unknown variables :
															// flows - solved flows + volumetric pump head
															// in branches with volumetric pumps, pump head is unknown

															// Coefficient matrices
	coef2_mat.resize(n_unknown); // coefficients for unknowns ^ 2 in each equation(row)
	coef1_mat.resize(n_unknown); // coefficients that multiply unknowns in each equation(row)
	const_mat.resize(n_unknown); // matrix for independent head loss in unsolved branches

	for (int count = 0; count < n_unknown; count++) {
			for (int count1 = 0; count1 < n_unknown; count1++) {
			coef2_mat.at(count).assign(n_unknown, 0);
			coef1_mat.at(count).assign(n_unknown, 0);
			const_mat.at(count).assign(n_unknown, 0);
		}
	}
	const_vec.assign(n_unknown, 0); // vector for constant terms of the equations :
									// independent head loss and head calculated with known flows

									// Auxiliary vector that contains the indices of the unsolved branches
	for (int count = 0; count < solved.size(); count++) {
		if (volum_pump(count) == true || solved(count) == false) {
			unk_ind_aux.at(size.unk_ind_aux) = (count);
			size.unk_ind_aux++;
		}

	}


	// Auxiliary vector that contains 1 at the indices where the vector
	// 'unknowns' or 'X' refer to branches that have volumetric pumps
	volpump_ind_aux.assign(size.unk_ind_aux, false);
	flows_ind_aux.assign(size.unk_ind_aux, true);
	for (int count = 0; count < size.unk_ind_aux; count++) {
		if ((volum_pump(unk_ind_aux.at(count)) == true)) {
			volpump_ind_aux.at(count) = true;
			flows_ind_aux.at(count) = false;
		}

	}

	

	// NODE EQUATIONS
	// All unsolved nodes are considered for the system but one of them is
	// discarded since it will be dependent on the rest.
	// Number of node equations = number of unsolved nodes - 1
	// Exception: there are tanks in the network.

	// Set first unsolved node as solved.Unsolved nodes will be included in the
	// system of equations.
	if (n_tanks == 0)
		for (int count = 0; count < n_nodes; count++) {
			if (solv_nodes(count) != true) // node not solved
				solv_nodes(count) = true; // set node as solved
			break; // done; exit

		}

	// Formulates node equations
	
	for (int count = 0; count < (solv_nodes.size() - solv_nodes.sum()); count++) { // unsolved nodes loop


																				   // Continuity equation : sum of flows in a node = 0
		for (int count1 = 0, aux; count1 < node_branches.at(count).size(); count1++) {// branches connected to node(branch index)

			aux = node_branches.at(count).at(count1);

			// Flow may be solved or not
			if (solved(aux) == true) // the flow in this branch has been solved
				const_vec.at(count) = flows.at(aux);
			else if (unk_ind_aux.at(count1) == aux)// flow is unknown
				coef1_mat.at(count).at(count1) = 1;


			// Node can be at the end or at the begining of the branch
			// Criterium: positive branch direction--> enters node(node is at the end of the branch)
			if (nodes_id.at(count) == branches_id.at(count1).at(1)) // node is at the beginning of the branch,
																	// branch exits node, change sign
				const_vec.at(count) = -const_vec.at(count);

		}
		last_row = count + 1; // save row to continue entering the equations
	}
	
	// MESH EQUATIONS

	// Number of independent meshes : branches to solve - nodes to solve + 1 = branches to solve - independent nodes to solve
	
	n_idp_mesh = (solved.size() - solved.sum()) - (solv_nodes.size() - solv_nodes.sum()) + volum_pump.sum() - n_tanks;
	
	// head of volumetric pumps must be solved as well
	// additionaly for every tank, there is a fake branch and two known pressure nodes : -1 mesh per tank

	// Selects meshes for the system of equations.Criterion : mesh must contain
	// the maximum possible number of unsolved branches.Next selected meshes
	// must contain the maximum number of unsolved branches excluding the
	// unsolved branches contained in meshes previously chosen.This way mesh
	// independence and solution of all branches are guaranteed.

	// Vector that stores the indices of the meshes that are included in the system of equations
	idp_mesh.assign(n_idp_mesh, 0);

	// Vector that marks the branches that are not covered by the current
	// selection of mesh equations.Remaining branches to be added to the system.
	
	for (int count = 0; count < solved.size(); count++) {
		remain_branch.at(size.remain_branch) = (!solved(count) || volum_pump(count));
		size.remain_branch++;
	}
	

	if (n_mesh > n_idp_mesh) { // there are more meshes than necessary->selection
		for (int count = 0; count < n_idp_mesh; count++) { // selected mesh loop
			best_mesh = 0; // selected mesh for the system of equations
			most_uns = 0; // number of unsolved branches included in the selected mesh
			for (int count1 = 0; count1 < n_mesh; count1++) { // mesh loop
				count_uns = 0; // counter of unsolved branches included in the selected mesh
				for (int count2 = 0; count2 < mesh_branches.at(count1).size(); count2++) { // branches contained in the mesh(vector index)
					if (remain_branch.at(mesh_branches.at(count1).at(count2)) == true) // branch must be solved
						count_uns = count_uns + 1; // add one point

				}
				if (count_uns > most_uns) { // solves more branches than previous meshes
					most_uns = count_uns; // refresh most unsolved branches
					best_mesh = count1; // refresh best mesh
				}
			}
			if (best_mesh == 0) {
				flag = false;
				for (int count1 = 0; count < n_mesh; count++) {
					for (int aux = 0; aux < idp_mesh.size(); aux++)
						if (count1 == idp_mesh.at(aux))
							flag = true;
					if (count1 == n_mesh - 1 && flag == false)
						idp_mesh.at(count) = count1;
					break;
				}

			}
			idp_mesh.at(count) = best_mesh; // save best mesh for the system
											// Remove all branches contained in the selected mesh from the remaining branches
			for (int count2 = 0; count2 < mesh_branches.at(best_mesh).size(); count2++) { // branches contained into mesh(vector index)
				remain_branch.at(mesh_branches.at(best_mesh).at(count2)) = false; // branch does not have to be solved
			}
		}
	}

	else // all meshes are necessary
		 // save all meshes as independent
		idp_mesh.assign(n_mesh, 0);

	idp_mesh_branches.resize(idp_mesh.size());

	for (int count = 0; count< idp_mesh.size(); count++)
		idp_mesh_branches.at(count) = mesh_branches.at(idp_mesh.at(count));


	// Initial values for the non - linear solver
	x0.setConstant(n_unknown, 0.0001);
	

	// Formulates mesh equations

	for (int count = 0; count < n_idp_mesh; count++) { // independent mesh loop

													   // II Kirchhoff's Law: sum of head increments in a mesh = 0
		last_id = 0;
		for (int count1 = 0, aux, pos; count1 < idp_mesh_branches.at(count).size(); count1++) {// branches in the mesh(branches index)

			aux = idp_mesh_branches.at(count).at(count1);

			for (int position = 0; position < size.unk_ind_aux; position++)
				if (unk_ind_aux.at(position) == aux)
					pos = position;

			// Branches are already ordered along the mesh
			// Criterium: positive mesh flow-- > direction of the first branch in the mesh
			// If there is a pump or a heat exchanger of type 'T out' in the branch, flow and branch direction are the same
			// If flow direction is unknown, assume same direction as branch
			// Flow may be solved or not

			// Check whether the branch has the same direction as the mesh flow
			change_sign = false;
			if (last_id == 0)
				// it is the first branch in the mesh OR the branch is in the same direction as the mesh flow
				last_id = branches_id.at(aux).at(branches_id.at(aux).size() - 1); // refresh other end of branch
			else { // it is not the first branch in the mesh AND the branch is in the opposite direction
				change_sign = true; // invert sign
				last_id = branches_id.at(aux).at(0); // refresh other end of branch
			}


			if (solved(aux) == true) {// the flow in this branch has been solved
				if (flows.at(aux) < 0.0)// flow direction is opposite to branch direction
					change_sign = ~change_sign; // invert sign

				const_vec.at(last_row + count) = const_vec.at(last_row + count) + (head_loss.at(aux) + hydr_resist1.at(aux) * abs(flows.at(aux)) + hydr_resist2.at(aux) * flows.at(aux)* flows.at(aux)) * (double)(-2 * change_sign + 1);
				// add head if change_sign = false; subtract if change_sign = true
				if (volum_pump(aux) == true)   //&& (unk_ind_aux.at(count) == aux)) // if there is a volumetric pump, include pump head as an unknown variable of the system
					coef1_mat.at(last_row + count).at(pos) = -(-2 * change_sign + 1) * 10000; // include pump head
																							  //unk_ind_aux.at(count) = aux;
																							  // change sign since head losses are positive and pump head must be positive too
																							  // multiplies by 10000 so the solution's magnitude order is closer to flow values

			}

			else { 
				   // assume flow mesh and branch flow have the same direction


				const_mat.at(last_row + count).at(pos) = head_loss.at(aux) * (-2 * change_sign + 1);
				coef1_mat.at(last_row + count).at(pos) = hydr_resist1.at(aux) * (-2 * change_sign + 1);
				coef2_mat.at(last_row + count).at(pos) = hydr_resist2.at(aux) * (-2 * change_sign + 1);

			}
		}
	}
	// The system is complete

	// Non-linear solver
	const_vec;
	const_mat;
	coef1_mat;
	coef2_mat;


	Y = NonLinearSolver(x0, coef1_mat, coef2_mat, const_mat, const_vec);


																		 // Save 
	for (int count = 0; count < solved.size(); count++) {
		for (int count1 = 0; count1 < Y.size(); count1++) {
			if (solved(count) == false && count == unk_ind_aux.at(count1))
				flows.at(count) = Y.at(count1);
		}
	}


	// Save head in volumetric pumps 
	for (int count = 0; count < volum_pump.size(); count++) {
		for (int count1 = 0; count1 < Y.size(); count1++) {
			if (volum_pump(count) == true && volpump_ind_aux.at(count1) == true)
				head_vol_pump.at(count) = Y.at(count1) * 10000;
		}
	}


	// previously pump head was divided by 10000 so its magnitude order was
	// close to that of the flows



}

void HydroNet_model::HydroNet_Temperature() {

	// Finds positions where temperature changes and temperatures in those ranges.
	// Moves positions and temperatures according to flows and, simultaneously,
	// enters heat from heat exchangers. Returns refreshed ranges and temperatures
	// in every branch as well as temperatures for the inlet of heat exchangers.
	// Merges ranges (volumes) and its temperatures if there are too many of them in a branch.

	int  extra_div, position_count = 0, node_count, end_node_ind, start_node_ind, node_aux, temp_action, pos_aux;
	std::vector<std::vector<double>>new_pos, new_temp;
	std::vector<std::vector<int>>  node_inlet_branches, node_outlet_branches;
	std::vector<double>  overflow_temperature, overflow_temp_volpos_aux, overflow_temperature_aux;
	double delta_pos, start_pos, end_pos, value = -9595959, sum_aux, temp_aux, flow_sum, sum_mass, sum_temp_mass, position, obj_pos;
	bool flows_flag = false;
	std::vector<bool>flag_inlet, to_remove;
	std::vector<double> node_volume_in, node_mass, node_temperature, volume_share, new_pos_aux, aux_vec;
	std::vector<int>branch_ind_order_aux;
	strInsertVolum InsertVolum;
	strSize size;
	Eigen::RowVectorXd overflow(n_branch);






	// Maximum divisions inside each branch; if there are more, some are merged
	//max_divisions = 20;



	/// VOLUME THAT STAYS INSIDE THE SAME BRANCH AND OVERFLOWING VOLUME

	new_pos.resize(n_branch); // temporary cell to store positions in the current instant
	new_temp.resize(n_branch); // temporary cell to store temperatures in the current instant
	overflow.setZero(); // vector to store volumes that move out from its original branch
	for (int count = 0; count < branch_volume.size(); count++) {
		if (branch_volume.at(count) > 0)
			overflow(count) = flows.at(count) * dt;
	}

	overflow_temperature.assign(n_branch, 0); // vector to store temperatures of overflowed volumes

	size.node_inlet_branches.resize(n_nodes);
	size.node_outlet_branches.resize(n_nodes);
	size.new_temp.resize(n_branch);
	size.new_pos.resize(n_branch);
	for (int count = 0; count < n_branch; count++) {
		size.new_temp.at(count) = 2 * max_divisions;
		size.new_pos.at(count) = 2 * max_divisions;
	}
	for (int count = 0; count < n_nodes; count++) {
		size.node_inlet_branches.at(count) = 0;
		size.node_outlet_branches.at(count) = 0;
	}
	size.overflow_temperature_aux = 0;
	size.overflow_temp_volpos_aux = 0;
	size.branch_ind_order_aux = 0;
	size.aux_vec = 0;
	size.to_remove = 0;
	size.position = 0;
	size.temperature = 0;

	overflow_temperature_aux.resize(2 * max_divisions);
	overflow_temp_volpos_aux.resize(2 * max_divisions);
	branch_ind_order_aux.resize(2 * max_divisions);
	aux_vec.resize(2 * max_divisions);
	node_inlet_branches.resize(n_nodes);
	node_outlet_branches.resize(n_nodes);
	for (int count = 0; count < n_branch; count++) {
		new_pos.at(count).resize(2 * max_divisions);
		new_temp.at(count).resize(2 * max_divisions);
	}
	for (int count = 0; count < n_nodes; count++) {
		node_inlet_branches.at(count).resize(n_branch);
		node_outlet_branches.at(count).resize(n_branch);
	}
	// Branches without volume
	for (int count_branch = 0; count_branch < branch_volume.size(); count_branch++) {
		if (branch_volume.at(count_branch) == 0) {
			new_pos.at(count_branch).at(0) = 0;
			new_pos.at(count_branch).at(1) = 100;
			size.new_pos.at(count_branch) = 2;
			for (int count = 0; count < branch_temperature.at(count_branch).size(); count++)
				new_temp.at(count_branch).at(count) = branch_temperature.at(count_branch).at(count);
			size.new_temp.at(count_branch) = branch_temperature.at(count_branch).size();

		}
	}


	// Branch loop excluding branches without volume
	for (int count_branch = 0; count_branch < branch_volume.size(); count_branch++) {
		if (branch_volume.at(count_branch) > 0) {
			delta_pos = flows.at(count_branch) * dt / branch_volume.at(count_branch) * 100;
			// % of branch volume that flow moved

			// NO FLOW MOVEMENT: refreshes temperatures inside the branch
			if (delta_pos == 0) {// flow is zero, no movement

				for (int count = 0; count < branch_temp_pos.at(count_branch).size(); count++)
					new_pos.at(count_branch).at(count) = branch_temp_pos.at(count_branch).at(count);
				size.new_pos.at(count_branch) = branch_temp_pos.at(count_branch).size();
				for (int count = 0; count < branch_temperature.at(count_branch).size(); count++)
					new_temp.at(count_branch).at(count) = branch_temperature.at(count_branch).at(count);
				size.new_temp.at(count_branch) = branch_temperature.at(count_branch).size();
				// If there are heat exchangers or pumps, inserts volume for it and refreshes temperature
				//// Heat exchanger of type T_out
				//for (int count_obj_htx = 0, ind_aux; count_obj_htx < branch_htx_Tout.at(count_branch).size(); count_obj_htx++) {

				//	ind_aux = heat_exch_Tout.at(count_obj_htx).heat_exch_Tout;
				//	// index of the heat exchanger in the table objects

				//	start_pos = obj_inlet_pos.at(ind_aux).at(0); // object's inlet position
				//	end_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position
				//	temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
				//	value = heat_exch_Tout.at(count_obj_htx).T_out; // final temperature

				//	InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
				//	for (int count = 0; count< InsertVolum.position.size(); count++)
				//		new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
				//	size.new_pos.at(count_branch) = InsertVolum.pos_size;
				//	for (int count = 0; count< InsertVolum.temperature.size(); count++)
				//		new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
				//	size.new_temp.at(count_branch) = InsertVolum.temp_size;
				//}

				// Heat exchanger 
				for (int count_obj_htx = 0, ind_aux; count_obj_htx < branch_htx.at(count_branch).size(); count_obj_htx++) {

					ind_aux = heat_exch.at(branch_htx.at(count_branch).at(count_obj_htx)).heat_exch;
					// index of the heat exchanger in the table objects

					start_pos = obj_inlet_pos.at(ind_aux).at(0); // object's inlet position
					end_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position
					if (heat_exch.at(branch_htx.at(count_branch).at(count_obj_htx)).type == "fix") {
						temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
						value = heat_exch.at(branch_htx.at(count_branch).at(count_obj_htx)).heat; // heat to add
					}
					else {
						temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
						value = heat_exch.at(branch_htx.at(count_branch).at(count_obj_htx)).T_out; // New volume's temperature
					}

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.position.size(); count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temperature.size(); count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}

				// Voluimetric Pump
				for (int count_obj_pump = 0, ind_aux; count_obj_pump < branch_pump_volum.at(count_branch).size(); count_obj_pump++) {//count_obj_htx = branch_htx_fix{ count_branch }

					ind_aux = pump_volum.at(branch_pump_volum.at(count_branch).at(count_obj_pump)).Pump_volum;
					// index of the heat exchanger in the table objects

					start_pos = obj_inlet_pos.at(ind_aux).at(0); // object's inlet position
					end_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position
					temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = pump_volum.at(branch_pump_volum.at(count_branch).at(count_obj_pump)).heat; // heat to add

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.position.size(); count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temperature.size(); count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}

				// Turbo pump
				for (int count_obj_pump = 0, ind_aux; count_obj_pump < branch_pump_turbo.at(count_branch).size(); count_obj_pump++) {

					ind_aux = pump_turbo.at(branch_pump_turbo.at(count_branch).at(count_obj_pump)).Pump_turbo;
					// index of the heat exchanger in the table objects

					start_pos = obj_inlet_pos.at(ind_aux).at(0); // object's inlet position
					end_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position
					temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = pump_turbo.at(branch_pump_turbo.at(count_branch).at(count_obj_pump)).heat; // heat to add

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.position.size(); count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temperature.size(); count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}
			}




			////////////////////////////////////
			// THERE IS FLOW MOVEMENT INSIDE THE BRANCH
			else {

				position_count = 1; // position counter

				new_pos.at(count_branch).at(0);

				size.new_pos.at(count_branch) = 1;
				size.new_temp.at(count_branch) = 0;

				// REFRESHES TEMPERATURES INSIDE THE BRANCH WITHOUT HEAT EXCHANGE
				while (((position_count) < (branch_temp_pos.at(count_branch).size())) && ((branch_temp_pos.at(count_branch).at(position_count) + delta_pos) <= 100)) {
					// Positions are moved along the branch according to flowing flow

					position = branch_temp_pos.at(count_branch).at(position_count);
					size.position = 1;

					// Flow must be positive (same direction as branch)
					new_pos.at(count_branch).at(size.new_pos.at(count_branch)) = (position + delta_pos); // add new position at the end of the vector
					size.new_pos.at(count_branch)++;
					new_temp.at(count_branch).at(size.new_temp.at(count_branch)) = (branch_temperature.at(count_branch).at(position_count - 1));
					size.new_temp.at(count_branch)++;

					position_count = position_count + 1;
				}

				new_pos.at(count_branch).at(size.new_pos.at(count_branch)) = 100; // add vector end
				size.new_pos.at(count_branch)++;
				new_temp.at(count_branch).at(size.new_temp.at(count_branch)) = (branch_temperature.at(count_branch).at(branch_temperature.at(count_branch).size() - 1));
				size.new_temp.at(count_branch)++;
				// add new temperature at the end of the vector


				////////////////////////////////////
				// VOLUME OVERFLOW: makes volume buffer without taking into account heat exchangers
				// First step to find the temperature of the buffer volume
				if (position_count <= branch_temp_pos.at(count_branch).size()) {
					// old position exists inside the vector, thus after the previous in-branch
					// while loop the new position has to be in another branch: OVERFLOW

					// Makes a new pair of vectors like branch_temp_pos and branch_temperature
					// for the volume that goes out of the branch
					overflow_temp_volpos_aux.at(0) = 100; // initialization: branch end
					size.overflow_temp_volpos_aux = 1;
					// position specified as % of the branch volume
					size.overflow_temperature_aux = 0;

					if (branch_temp_pos.at(count_branch).size() == 2) {
						overflow_temperature_aux.at(size.overflow_temperature_aux) = (branch_temperature.at(count_branch).at(0));
						size.overflow_temperature_aux++;
						overflow_temp_volpos_aux.at(0) = 0;
						overflow_temp_volpos_aux.at(1) = 100;
						size.overflow_temp_volpos_aux = 2;
					}

					for (int count = position_count; count < branch_temp_pos.at(count_branch).size(); count++) {

						overflow_temp_volpos_aux.at(size.overflow_temp_volpos_aux) = (branch_temp_pos.at(count_branch).at(position_count) + delta_pos);
						size.overflow_temp_volpos_aux++;
						overflow_temperature_aux.at(size.overflow_temperature_aux) = (branch_temperature.at(count_branch).at(position_count - 1));
						size.overflow_temperature_aux++;
					}
				}


				////////////////////////////////////
				// HEAT EXCHANGERS & PUMPS
				// If there are heat exchangers or pumps, checks whether the volume out of them remains inside the branch or it
				// goes out, inserts volume for it, refreshes temperature inside the branch and stores overflowed volume.
				// There can only be heat exchangers if branch volume is > 0.
				// Second step to find the temperature of the buffer volume

				/*// Heat exchanger of type T_out
				for (int count = 0, count_obj_htx, ind_aux; count < branch_htx_Tout.at(count_branch).size(); count++) {// loop starting from the closest exchanger to the branch's end
					count_obj_htx = branch_htx_Tout.at(count_branch).at(count);

					for (int aux = 0; aux<heat_exch_Tout.size(); aux++)
						if (count_obj_htx == heat_exch_Tout.at(aux).heat_exch_Tout) {
							ind_aux = aux;
							break;
						}
					// index of the heat exchanger in the table objects

					temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = heat_exch_Tout.at(count_obj_htx).T_out; // new temperature

					start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
					end_pos = start_pos + delta_pos;

					if (end_pos > 100) {
						start_pos = 100;

						InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
						for (int count = 0; count< InsertVolum.position.size(); count++)
							overflow_temp_volpos_aux.at(count) = InsertVolum.position.at(count);
						size.overflow_temp_volpos_aux = InsertVolum.pos_size;
						for (int count = 0; count< InsertVolum.temperature.size(); count++)
							overflow_temperature_aux.at(count) = InsertVolum.temperature.at(count);
						size.overflow_temperature_aux = InsertVolum.temp_size;
						// For the volume inside the branch
						start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
						end_pos = 100;
					}

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count< InsertVolum.position.size(); count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count< InsertVolum.temperature.size(); count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}*/

				// Heat exchanger of type fixed heat
				for (int count = 0, count_obj_htx, ind_aux; count < branch_htx.at(count_branch).size(); count++) { // loop starting from the closest exchanger to the branch's end
					count_obj_htx = branch_htx.at(count_branch).at(count);
					for (int aux = 0; aux < heat_exch.size(); aux++)
						if (count_obj_htx == heat_exch.at(aux).heat_exch) {
							ind_aux = aux;
							break;
						}

					// index of the heat exchanger in the table objects
					if (heat_exch.at(ind_aux).type == "fix") {
						temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
						value = heat_exch.at(ind_aux).heat;

						start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
						end_pos = start_pos + delta_pos;

						if (end_pos > 100) {

							value = heat_exch.at(ind_aux).heat* (end_pos - 100) / (end_pos - start_pos);
							// heat to add averaged for the volume that goes out of the branch

							start_pos = 100;

							InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
							for (int count = 0; count < InsertVolum.position.size(); count++)
								overflow_temp_volpos_aux.at(count) = InsertVolum.position.at(count);
							size.overflow_temp_volpos_aux = InsertVolum.pos_size;
							for (int count = 0; count < InsertVolum.temperature.size(); count++)
								overflow_temperature_aux.at(count) = InsertVolum.temperature.at(count);
							size.overflow_temperature_aux = InsertVolum.temp_size;

							// For the volume inside the branch
							start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
							end_pos = 100;
							value = heat_exch.at(ind_aux).heat* (100 - start_pos) / (end_pos - start_pos);
							// heat to add averaged for the volume that stays inside the branch
						}
					}
					else {
						temp_action = 2; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
						value = heat_exch.at(count_obj_htx).T_out; // new temperature

						start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
						end_pos = start_pos + delta_pos;

						if (end_pos > 100) {
							start_pos = 100;

							InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
							for (int count = 0; count < InsertVolum.position.size(); count++)
								overflow_temp_volpos_aux.at(count) = InsertVolum.position.at(count);
							size.overflow_temp_volpos_aux = InsertVolum.pos_size;
							for (int count = 0; count < InsertVolum.temperature.size(); count++)
								overflow_temperature_aux.at(count) = InsertVolum.temperature.at(count);
							size.overflow_temperature_aux = InsertVolum.temp_size;
							// For the volume inside the branch
							start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
							end_pos = 100;
						}
					}

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.pos_size; count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temp_size; count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}

				// Volumetric pumps
				for (int count = 0, count_obj_pmp, ind_aux; count < branch_pump_volum.at(count_branch).size(); count++) { // loop starting from the closest exchanger to the branch's end
					count_obj_pmp = branch_pump_volum.at(count_branch).at(count);
					for (int aux = 0; aux < pump_volum.size(); aux++)
						if (count_obj_pmp == pump_volum.at(aux).Pump_volum) {
							ind_aux = aux;
							break;
						}
					// index of the heat exchanger in the table objects

					temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = pump_volum.at(count_obj_pmp).heat;

					start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
					end_pos = start_pos + delta_pos;

					if (end_pos > 100) {

						value = pump_volum.at(count_obj_pmp).heat* (end_pos - 100) / (end_pos - start_pos);
						// heat to add averaged for the volume that goes out of the branch

						start_pos = 100;

						InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
						for (int count = 0; count < InsertVolum.position.size(); count++)
							overflow_temp_volpos_aux.at(count) = InsertVolum.position.at(count);
						size.overflow_temp_volpos_aux = InsertVolum.pos_size;
						for (int count = 0; count < InsertVolum.temperature.size(); count++)
							overflow_temperature_aux.at(count) = InsertVolum.temperature.at(count);
						size.overflow_temperature_aux = InsertVolum.temp_size;

						// For the volume inside the branch
						start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
						end_pos = 100;
						value = pump_volum.at(count_obj_pmp).heat* (100 - start_pos) / (end_pos - start_pos);
						// heat to add averaged for the volume that stays inside the branch
					}

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.pos_size; count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temp_size; count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}
				// Turbo pump
				for (int count = 0, count_obj_pmp, ind_aux; count < branch_pump_turbo.at(count_branch).size(); count++) { // loop starting from the closest exchanger to the branch's end
					count_obj_pmp = branch_pump_turbo.at(count_branch).at(count);
					for (int aux = 0; aux < pump_turbo.size(); aux++)
						if (count_obj_pmp == pump_turbo.at(aux).Pump_turbo) {
							ind_aux = aux;
							break;
						}
					// index of the heat exchanger in the table objects

					temp_action = 1; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
					value = pump_turbo.at(count_obj_pmp).heat;

					start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
					end_pos = start_pos + delta_pos;

					if (end_pos > 100) {

						value = pump_turbo.at(count_obj_pmp).heat* (end_pos - 100) / (end_pos - start_pos);
						// heat to add averaged for the volume that goes out of the branch

						start_pos = 100;

						InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
						for (int count = 0; count < InsertVolum.position.size(); count++)
							overflow_temp_volpos_aux.at(count) = InsertVolum.position.at(count);
						size.overflow_temp_volpos_aux = InsertVolum.pos_size;
						for (int count = 0; count < InsertVolum.temperature.size(); count++)
							overflow_temperature_aux.at(count) = InsertVolum.temperature.at(count);
						size.overflow_temperature_aux = InsertVolum.temp_size;

						// For the volume inside the branch
						start_pos = obj_outlet_pos.at(ind_aux).at(0); // object's outlet position is the start position
						end_pos = 100;
						value = pump_volum.at(count_obj_pmp).heat* (100 - start_pos) / (end_pos - start_pos);
						// heat to add averaged for the volume that stays inside the branch
					}

					InsertVolum = HydroNet_InsertVolume(new_pos.at(count_branch), new_temp.at(count_branch), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.new_pos.at(count_branch), size.new_temp.at(count_branch));
					for (int count = 0; count < InsertVolum.pos_size; count++)
						new_pos.at(count_branch).at(count) = InsertVolum.position.at(count);
					size.new_pos.at(count_branch) = InsertVolum.pos_size;
					for (int count = 0; count < InsertVolum.temp_size; count++)
						new_temp.at(count_branch).at(count) = InsertVolum.temperature.at(count);
					size.new_temp.at(count_branch) = InsertVolum.temp_size;
				}

				////////////////////////////////////
				// TEMPERATURE OF THE BUFFER VOLUME (OVERFLOW)
				// Final averaging of temperatures in the buffer volume

				start_pos = 100;
				end_pos = overflow_temp_volpos_aux.at(size.overflow_temp_volpos_aux - 1);
				temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature

				InsertVolum = HydroNet_InsertVolume(overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(count_branch), max_divisions, size.overflow_temp_volpos_aux, size.overflow_temperature_aux);
				overflow_temperature.at(count_branch) = InsertVolum.temperature.at(InsertVolum.temp_size - 1);
			}
			size.overflow_temperature_aux = 0;
			size.overflow_temp_volpos_aux = 0;
		}
	}





	//// SEPARES BRANCHES AT THE INLET AND OUTLET OF EVERY NODE
	// Separates branches that are inlets to nodes and those that are outlets to nodes from node_branches
	// All branches should start and end at a node!!
	for (int count = 0; count < flows.size(); count++)
		if (flows.at(count) > 0)
			flows_flag = true;

	if (flows_flag == true) { // if there is no movement in the branch, there is no need for this

							  // stores the sum of all volumes leaving every node to check if outlet branches are overflowed
							  // no used because, later, volume_share is more precise

		for (int node_count = 0, count; node_count < n_nodes; node_count++) { // node loop
			flag_inlet.assign(node_branches.at(node_count).size(), false); // marks inlet branches to this node
			count = -1; // counter of branches in every node
			for (int count1 = 0, count_branch; count1 < node_branches.at(node_count).size(); count1++) {

				count_branch = node_branches.at(node_count).at(count1);
				count++;
				if (branches_id.at(count_branch).at(branches_id.at(count_branch).size() - 1) == nodes_id.at(node_count))
					// flow goes into the node from the branch: it is an inlet branch to the node
					flag_inlet.at(count) = true; // mark branch as inlet to node
												 

			}

			// Separates inlet and outlet branches to each node
			for (int aux = 0; aux < flag_inlet.size(); aux++)
				if (flag_inlet.at(aux) == true) {
					node_inlet_branches.at(node_count).at(size.node_inlet_branches.at(node_count)) = node_branches.at(node_count).at(aux);
					size.node_inlet_branches.at(node_count)++;
				}

			for (int aux = 0; aux < flag_inlet.size(); aux++)
				if (flag_inlet.at(aux) == false) {
					node_outlet_branches.at(node_count).at(size.node_outlet_branches.at(node_count)) = node_branches.at(node_count).at(aux);
					size.node_outlet_branches.at(node_count)++;
				}
			
		}
	}


	//// MANAGES OVERFLOWED VOLUMES
	// All branches should start and end at a node!!
	// The number of branches without volume must be minimized!!

	while (overflow.isZero()) { // while there are overflows

											// Makes an enthalpy balance of all volumes mixing in every node
		node_volume_in.assign(n_nodes, 0); // stores the sum of all volumes mixing in every node
		node_mass.assign(n_nodes, 0); // temporary vector to store mass in every node
		node_temperature.assign(n_nodes, 0); // stores temperatures in every node after mixing

											 // Volume and mass in nodes
		for (int count = 0, count_branch; count < overflow.size(); count++) { // 1 : n_branch // branch loop
			if (overflow(count) > 0) {
				count_branch = count;
				// do not consider branches that do not have overflows (first iteration: do not considerer nodes without volume)
				// branch_count = 1 : n_branch would be valid too (overflow = 0 does not add volume or mass)
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size() - 1))
						end_node_ind = aux; // index of the node at the end of the branch

				node_volume_in.at(end_node_ind) = node_volume_in.at(end_node_ind) + overflow(count_branch);
				node_mass.at(end_node_ind) = node_mass.at(end_node_ind) + density(fluid_type, overflow_temperature.at(count_branch)) * overflow(count_branch);
			}
		}

		// Temperature in nodes, provisional
		// Processes nodes that have overflow buffer volumes > 0

		for (int count = 0, count_node; count < node_mass.size(); count++) { // to avoid dividing by zero
			if (node_mass.at(count) > 0) {
				count_node = count;
				for (int count1 = 0, count_branch; count1 < size.node_inlet_branches.at(count_node); count1++) {
					count_branch = node_inlet_branches.at(count_node).at(count1);
					for (int aux = 0; aux < nodes_id.size(); aux++)
						if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size() - 1))
							end_node_ind = aux; // index of the node at the end of the branch
					node_temperature.at(end_node_ind) = node_temperature.at(end_node_ind) + overflow_temperature.at(count_branch) * density(fluid_type, overflow_temperature.at(count_branch)) *overflow(count_branch) / node_mass.at(end_node_ind);
					// T branch to node * branch mass to node / total mass to node
				}
			}
		}

		// Control goes to buffer volumes in the nodes node_volume_in
		for (int count = 0; count < overflow.size(); count++) {
			overflow(count) = 0;
			overflow_temperature.at(count) = 0;
		}


		// Passes buffer volumes that are at the inlet node of branches with
		// volume = 0 to the end node

		sum_aux = 0; // sum of volumes at the start of branches without volume
		for (int count = 0, count_branch; count < flows.size(); count++) {
			if ((branch_volume.at(count) == 0) && (flows.at(count) > 0)) {
				count_branch = count;
				// loop of branches without volume but with flow movement

				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = aux;

				sum_aux = sum_aux + node_volume_in.at(start_node_ind);
			}
		}

		while (sum_aux > 0) {
			// there are buffer volumes at the start of branches without volume; they must be past through

			// Priority: the node at the start of the branch has no inlet branches
			// without volume (they also would get passed a volume from behind)

			for (int count = 0; count < flows.size(); count++)
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0)) {
					branch_ind_order_aux.at(size.branch_ind_order_aux) = (count);// branches without volume but with flow movement
					size.branch_ind_order_aux++;
				}

			for (int count = 0, count_branch; count < flows.size(); count++) {
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0)) {

					count_branch = count;
					for (int aux = 0; aux < nodes_id.size(); aux++)
						if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
							start_node_ind = aux;
					for (int aux = 0; aux < size.node_inlet_branches.at(start_node_ind); aux++)
						if (branch_volume.at(node_inlet_branches.at(start_node_ind).at(aux)) == 0) {
							// any of the inlet branches to the start node has no volume
							for (int count1 = 0; count < size.branch_ind_order_aux; count1++) {
								if (branch_ind_order_aux.at(count1) == count_branch) {
									for (int aux1 = count1; aux1 < size.branch_ind_order_aux && aux1<2 * max_divisions - 1; aux1++) {
										branch_ind_order_aux.at(aux1) = branch_ind_order_aux.at(aux1 + 1);
									}
									size.branch_ind_order_aux--;
									branch_ind_order_aux.at(size.branch_ind_order_aux) = (count_branch); // put the branch at the bottom of the list
									size.branch_ind_order_aux++;
								}
							}
						}
				}
			}

			// Moves volumes according to priority
			for (int count = 0, count_branch; count < size.branch_ind_order_aux; count++) {
				count_branch = branch_ind_order_aux.at(count);
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
						start_node_ind = aux;// index of the node at the start of the branch
				for (int aux = 0; aux < nodes_id.size(); aux++)
					if (nodes_id.at(aux) == branches_id.at(count_branch).at(branches_id.at(count_branch).size() - 1))
						end_node_ind = aux; // index of the node at the end of the branch

				node_volume_in.at(end_node_ind) = node_volume_in.at(end_node_ind) + node_volume_in.at(start_node_ind);
				temp_aux = ((node_temperature.at(end_node_ind) * node_mass.at(end_node_ind)) + (node_temperature.at(start_node_ind) * node_mass.at(start_node_ind))) / (node_mass.at(end_node_ind) + node_mass.at(start_node_ind));
				// balance between volume already present and incoming volume
				node_temperature.at(end_node_ind) = temp_aux;
				node_mass.at(end_node_ind) = node_mass.at(end_node_ind) + node_mass.at(start_node_ind);
				node_volume_in.at(start_node_ind) = 0;
				node_mass.at(start_node_ind) = 0;
				node_temperature.at(start_node_ind) = 0;
				new_temp.at(count_branch).at(size.new_temp.at(count_branch)) = (temp_aux); // assigns temperature to the branch
				size.new_temp.at(count_branch)++;
				if (size.new_temp.at(count_branch) <= 2) {
					for (int aux = 0; aux < size.new_temp.at(count_branch) && aux < 2 * max_divisions - 1; aux++)
						new_temp.at(count_branch).at(aux) == new_temp.at(count_branch).at(aux + 1);
					size.new_temp.at(count_branch)--;
				}
				
			}

			sum_aux = 0; // sum of volumes at the start of branches without volume
			for (int count = 0, count_branch; count < flows.size(); count++) {
				if ((branch_volume.at(count) == 0) && (flows.at(count) > 0)) {
					count_branch = count;
					// loop of branches without volume but with flow movement
					for (int aux = 0; aux < nodes_id.size(); aux++)
						if ((nodes_id.at(aux) == branches_id.at(count_branch).at(0)))
							start_node_ind = aux; // index of the node at the start of the branch
					sum_aux = sum_aux + node_volume_in.at(start_node_ind);
				}
			}
			size.branch_ind_order_aux = 0;
		}



		// Distributes volume among branches
		for (int count = 0; count < node_volume_in.size(); count++)
			if (node_volume_in.at(count) > 0) { // loop of nodes that have a buffer volume
				node_count = count;
				volume_share.resize(size.node_outlet_branches.at(node_count));
				flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count), size.node_outlet_branches.at(node_count)); // sum of flows leaving the node
				for (int aux = 0; aux < size.node_outlet_branches.at(node_count); aux++)
					volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
				// volume that goes to every branch that is an outlet to the node (proportional to flow)
				for (int count_branch = 0, branch_aux; count_branch < size.node_outlet_branches.at(node_count); count_branch++) { // loop of branches that are outlets to the node
					branch_aux = node_outlet_branches.at(node_count).at(count_branch); // branch index
																					 
					if (volume_share.at(count_branch) > 0) {// processes only branches where volume actually goes in
						if (volume_share.at(count_branch) <= branch_volume.at(branch_aux)) {
							// volume in the branch is bigger than incoming volume
							// also avoids dividing by zero

							// Adds position and temperature at the beginning
							for (int aux = 0; aux < size.new_pos.at(branch_aux) && aux < 2 * max_divisions - 1; aux++)
								new_pos.at(branch_aux).at(size.new_pos.at(branch_aux) - aux) = new_pos.at(branch_aux).at(size.new_pos.at(branch_aux) - aux - 1);
							new_pos.at(branch_aux).at(0) = volume_share.at(count_branch) / branch_volume.at(branch_aux) * 100;
							for (int aux = 0; aux < size.new_temp.at(branch_aux) && aux < 2 * max_divisions - 1; aux++)
								new_temp.at(branch_aux).at(size.new_temp.at(branch_aux) - aux) = new_temp.at(branch_aux).at(size.new_temp.at(branch_aux) - aux - 1);
							new_temp.at(branch_aux).at(0) = node_temperature.at(node_count);
							size.new_pos.at(branch_aux)++;
							size.new_temp.at(branch_aux)++;


							// Removes positions lower than or equal to the new one
							to_remove.assign(2 * max_divisions, false);
							size.to_remove = size.new_pos.at(branch_aux);
							for (int pos_aux = 1; pos_aux < size.new_pos.at(branch_aux); pos_aux++) {
								if (new_pos.at(branch_aux).at(pos_aux) <= new_pos.at(branch_aux).at(0)) {
									to_remove.at(pos_aux) = true;
								}
							}
							for (int aux = 0; aux < size.to_remove; aux++)
								if (to_remove.at(aux) == true) {
									for (int aux1 = aux; aux1 < size.new_pos.at(branch_aux) && aux1 < 2 * max_divisions - 2; aux1++)
										new_pos.at(branch_aux).at(aux1) = new_pos.at(branch_aux).at(aux1 + 1);
									size.new_pos.at(branch_aux)--;
								}

							// Does not remove temperature that is after the last position to remove
							for (int count1 = 0; count1 < size.to_remove; count1++)  {
								if (to_remove.at(count1) == true) {
									to_remove.at(count1) = false;
									break;
								}
							}
							//to_remove.erase(to_remove.begin());
							for (int aux = 0; aux < size.to_remove && aux < 2 * max_divisions - 2; aux++)
								to_remove.at(aux) = to_remove.at(aux + 1);
							size.to_remove--;
							for (int aux = 0; aux < size.to_remove; aux++)
								if (to_remove.at(aux) == true) {
									for (int aux1 = aux; aux1 < size.new_temp.at(branch_aux) && aux1 < 2 * max_divisions - 2; aux1++) {
										new_temp.at(branch_aux).at(aux1) = new_temp.at(branch_aux).at(aux1 + 1);
									}
									size.new_temp.at(branch_aux)--;
								}

							// Adds position 0 at the beginning
							for (int aux = 0; aux < size.new_pos.at(branch_aux) && aux < 2 * max_divisions - 1; aux++)
								new_pos.at(branch_aux).at(size.new_pos.at(branch_aux) - aux) = new_pos.at(branch_aux).at(size.new_pos.at(branch_aux) - aux - 1);

							size.new_pos.at(branch_aux)++;

							new_pos.at(branch_aux).at(0) = 0;
						}

						else { // volume in the branch is smaller than incoming volume

							   // Fills branch
							new_pos.at(branch_aux) = { 0, 100 };
							new_temp.at(branch_aux).at(0) = node_temperature.at(node_count);
							size.new_pos.at(branch_aux) = 2;
							size.new_temp.at(branch_aux) = 1;

							// Stores the excess volume
							overflow(branch_aux) = volume_share.at(branch_aux) - branch_volume.at(branch_aux);
							overflow_temperature.at(branch_aux) = node_temperature.at(node_count);
						}
					}
				}
			}

		for (int count = 0; count < node_volume_in.size(); count++)
			node_volume_in.at(count) = 0; // just in case it is used after the while loop
										  // control passes to overflow in branches
	}




	//// MERGES VOLUMES INSIDE A BRANCH
	// If there are more volumes than the value specified in max_divisions

	for (int count_branch = 0; count_branch < n_branch; count_branch++) {
		extra_div = size.new_pos.at(count_branch) - 1 - max_divisions; // number of divisions over max_divisions
		while (extra_div > 0) {
			new_pos_aux = new_pos.at(count_branch);

			// Finds smaller combination of adjacent volumes
			for (int aux = 0; aux < new_pos.size() - 3; aux++)
				aux_vec.at(size.aux_vec) = (new_pos_aux.at(3 + aux) - new_pos_aux.at(aux));
			size.aux_vec++;
			pos_aux = MinComponent(aux_vec);

			// Obtains temperature by means of an enthalpy balance
			temp_aux = (new_temp.at(count_branch).at(pos_aux) * density(fluid_type, new_temp.at(count_branch).at(pos_aux)) * (new_pos.at(count_branch).at(pos_aux + 1) - new_pos.at(count_branch).at(pos_aux)) + new_temp.at(count_branch).at(pos_aux + 1) * density(fluid_type, new_temp.at(count_branch).at(pos_aux + 1)) * (new_pos.at(count_branch).at(pos_aux + 2) - new_pos.at(count_branch).at(pos_aux + 1))) / (density(fluid_type, new_temp.at(count_branch).at(pos_aux)) * (new_pos_aux.at(pos_aux + 1) - new_pos_aux.at(pos_aux)) + density(fluid_type, new_temp.at(count_branch).at(pos_aux + 1)) * (new_pos_aux.at(pos_aux + 2) - new_pos_aux.at(pos_aux + 1)));
			// Merges volumes
			for (int aux1 = pos_aux + 1; aux1 < size.new_pos.at(count_branch) && aux1 < 2 * max_divisions - 2; aux1++) {// remove position in the middle
				new_pos.at(count_branch).at(aux1) = new_pos.at(count_branch).at(aux1 + 1);
			}
			size.new_pos.at(count_branch)--;
			for (int aux1 = pos_aux + 1; aux1 < size.new_temp.at(count_branch) && aux1 < 2 * max_divisions - 2; aux1++) {// remove second temperature
				new_temp.at(count_branch).at(aux1) = new_temp.at(count_branch).at(aux1 + 1);
			}
			size.new_temp.at(count_branch)--;
			new_temp.at(count_branch).at(pos_aux) = temp_aux; // enter temperature from balance

			extra_div = extra_div - 1; // next division to merge
		}
		size.aux_vec = 0;

	}



	//// CALCULATES INLET TEMPERATURE TO HEAT EXCHANGERS & PUMPS
	// Asumption: the flow that will pass through the heat exchanger and pumps in the next
	// instant is equal to the branch flow in the current instant.
	// Exception: if flow = 0 in the current instant, inserts a volume for the
	// heat exchanger or pump and obtains its temperature.
	// Backwards direction from the outlet of the heat exchanger or pump.

	/*// Heat exchangers of type Tout
	for (int count_htx = 0, obj_aux, branch_aux; count_htx < heat_exch_Tout.size(); count_htx++) {

		obj_aux = heat_exch_Tout.at(count_htx).heat_exch_Tout; // index of the heat exchanger in objects table
															   //branch_aux = objects.branch.at(obj_aux); // branch that contains the heat exchanger

		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}
		if (flows.at(branch_aux) == 0) { // no movement in the branch

										 // Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux).at(0);
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the heat exchanger
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_Tout.at(count_htx).T_in = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux).at(0) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch_Tout.at(count_htx).T_in = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);
		}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			  // First, obtains temperature of the volume between the start of the
			  // branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);
			heat_exch_Tout.at(count_htx).T_in = temp_aux;

			// Initializes enthalpy balance
			sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

												 // Stores volume that will come from previous branches
			overflow(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (overflow.isZero()) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow(count) > 0) {
						count_branch = count;

						node_aux = branches_id.at(count_branch).at(1); // node at the start of the branch
						flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count), size.node_outlet_branches.at(node_count)); // sum of flows going into the branch

						for (int aux = 0; aux < size.node_outlet_branches.at(node_count); aux++)
							volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
						// volume that comes from every inlet branch (proportional to flow)
						for (int count_branch1 = 0, branch_aux1; count_branch1 < size.node_inlet_branches.at(node_aux); count_branch1++) { // loop of inlet branches
							branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
							if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
								// volume in the branch is bigger than outcoming volume
								// also avoids dividing by zero

								// Inserts a volume for the affected area and obtains its temperature
								start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

								// Reads the temperature in the middle of the affected volume
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);

																															   // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch1);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);
							}

							else {// volume in the branch is smaller than outcoming volume

								  // Inserts a volume for the entire branch and obtains its temperature
								start_pos = 0;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

								// Reads the temperature in the middle of the branch
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);

																															  // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch1);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);

								// Stores the excess volume
								overflow(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
							}
						}
					}
				}
			}
			for (int count = 0; count < overflow.size(); count++)
				overflow(count) = 0;// leaves everything as found

									// Stores resulting inlet temperature
			heat_exch_Tout.at(count_htx).T_in = sum_temp_mass / sum_mass;
		}
	}*/


	// Heat exhangers
	for (int count_htx = 0, branch_aux, obj_aux; count_htx < heat_exch.size(); count_htx++) {


		obj_aux = heat_exch.at(count_htx).heat_exch; // index of the heat exchanger in objects table
		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}// branch that contains the heat exchanger

		if (flows.at(branch_aux) == 0) { // no movement in the branch

										 // Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux).at(0);
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the heat exchanger
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch.at(count_htx).T_in = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);//branch_temp_pos_aux, branch_temp_aux);
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux).at(0) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			heat_exch.at(count_htx).T_in = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);
		}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			  // First, obtains temperature of the volume between the start of the
			  // branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Reads and stores the temperature in the middle of the affected volume
			obj_pos = (start_pos + end_pos) / 2;
			temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);
			heat_exch.at(count_htx).T_in = temp_aux;

			// Initializes enthalpy balance
			sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

												 // Stores volume that will come from previous branches
			overflow(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (overflow.isZero()) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow(count) > 0) {
						count_branch = count;
						for (int aux = 0; aux<nodes_id.size(); aux++)
							if (branches_id.at(count_branch).at(1) == nodes_id.at(aux))
								node_aux = aux; // node at the start of the branch
						flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count), size.node_outlet_branches.at(node_count)); // sum of flows going into the branch

						for (int aux = 0; aux < size.node_outlet_branches.at(node_count); aux++)
							volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
						// volume that comes from every inlet branch (proportional to flow)
						for (int count_branch1 = 0, branch_aux1; count_branch1 < size.node_inlet_branches.at(node_aux); count_branch1++) { // loop of inlet branches
							branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
							if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
								// volume in the branch is bigger than outcoming volume
								// also avoids dividing by zero

								// Inserts a volume for the affected area and obtains its temperature
								start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the affected volume
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);

																															   // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch1);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);
							}

							else {// volume in the branch is smaller than outcoming volume

								  // Inserts a volume for the entire branch and obtains its temperature
								start_pos = 0;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the branch
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);

																															  // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);

								// Stores the excess volume
								overflow(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
							}
						}
					}
				}
			}

			for (int count = 0; count < overflow.size(); count++)
				overflow(count) = 0;// leaves everything as found

									// Stores resulting inlet temperature
			heat_exch.at(count_htx).T_in = sum_temp_mass / sum_mass;

		}
	}

	// Volumetric pump
	for (int count_pmp = 0, branch_aux, obj_aux; count_pmp < pump_volum.size(); count_pmp++) {


		obj_aux = pump_volum.at(count_pmp).Pump_volum; // index of the heat exchanger in objects table
		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}// branch that contains the heat exchanger

		if (flows.at(branch_aux) == 0) { // no movement in the branch

										 // Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux).at(0);
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux).at(0) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));
		}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			  // First, obtains temperature of the volume between the start of the
			  // branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			// Initializes enthalpy balance
			sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

												 // Stores volume that will come from previous branches
			overflow(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (overflow.isZero()) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow(count) > 0) {
						count_branch = count;
						for (int aux = 0; aux<nodes_id.size(); aux++)
							if (branches_id.at(count_branch).at(1) == nodes_id.at(aux))
								node_aux = aux; // node at the start of the branch
						flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count), size.node_outlet_branches.at(node_count)); // sum of flows going into the branch

						for (int aux = 0; aux < size.node_outlet_branches.at(node_count); aux++)
							volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
						// volume that comes from every inlet branch (proportional to flow)
						for (int count_branch1 = 0, branch_aux1; count_branch1 < size.node_inlet_branches.at(node_aux); count_branch1++) { // loop of inlet branches
							branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
							if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
								// volume in the branch is bigger than outcoming volume
								// also avoids dividing by zero

								// Inserts a volume for the affected area and obtains its temperature
								start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the affected volume
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);

																															   // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch1);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);
							}

							else {// volume in the branch is smaller than outcoming volume

								  // Inserts a volume for the entire branch and obtains its temperature
								start_pos = 0;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the branch
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);

																															  // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);

								// Stores the excess volume
								overflow(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
							}
						}
					}
				}
			}

			for (int count = 0; count < overflow.size(); count++)
				overflow(count) = 0;// leaves everything as found

								
		}
	}


	// Turbo pump
	for (int count_pmp = 0, branch_aux, obj_aux; count_pmp < pump_turbo.size(); count_pmp++) {


		obj_aux = pump_turbo.at(count_pmp).Pump_turbo; // index of the heat exchanger in objects table
		for (int count = 0; count < branches.size(); count++)
			for (int count1 = 0; count1 < branches.at(count).objects.size(); count1++)
				if (obj_aux == branches.at(count).objects.at(count1).ID) {
					branch_aux = count;
					break;
				}// branch that contains the heat exchanger

		if (flows.at(branch_aux) == 0) { // no movement in the branch

										 // Inserts a volume for the heat exchanger and obtains its temperature
			start_pos = obj_inlet_pos.at(obj_aux).at(0);
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			
		}

		else if ((flows.at(branch_aux) * dt) <= (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux))) {
			// there is movement in the branch and the affected volume is inside the branch

			// Inserts a volume for the affected area and obtains its temperature
			start_pos = obj_outlet_pos.at(obj_aux).at(0) - flows.at(branch_aux) * dt / branch_volume.at(branch_aux) * 100;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			}

		else {// there is movement in the branch and the affected volume reaches the previous branches

			  // First, obtains temperature of the volume between the start of the
			  // branch and the heat exchanger outlet, inserting a volume
			start_pos = 0;
			end_pos = obj_outlet_pos.at(obj_aux).at(0);
			temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
			InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux), new_temp.at(branch_aux), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux), max_divisions, size.new_pos.at(branch_aux), size.new_temp.at(branch_aux));

			
			// Initializes enthalpy balance
			sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));
			sum_temp_mass = temp_aux * sum_mass; // sum of temperature times mass

												 // Stores volume that will come from previous branches
			overflow(branch_aux) = (flows.at(branch_aux) * dt) - (obj_outlet_pos.at(obj_aux).at(0) / 100 * branch_volume.at(branch_aux));

			// Loop to obtain temperatures of the incoming flows
			while (overflow.isZero()) {// while there are overflows
				for (int count = 0, count_branch; count < overflow.size(); count++) { // loop of branches with overflow
					if (overflow(count) > 0) {
						count_branch = count;
						for (int aux = 0; aux<nodes_id.size(); aux++)
							if (branches_id.at(count_branch).at(1) == nodes_id.at(aux))
								node_aux = aux; // node at the start of the branch
						flow_sum = vec_sum_elements(flows, node_outlet_branches.at(node_count), size.node_outlet_branches.at(node_count)); // sum of flows going into the branch

						for (int aux = 0; aux < size.node_outlet_branches.at(node_count); aux++)
							volume_share.at(aux) = flows.at(node_outlet_branches.at(node_count).at(aux)) / flow_sum*node_volume_in.at(node_count);
						// volume that comes from every inlet branch (proportional to flow)
						for (int count_branch1 = 0, branch_aux1; count_branch1 < size.node_inlet_branches.at(node_aux); count_branch1++) { // loop of inlet branches
							branch_aux1 = node_inlet_branches.at(node_aux).at(count_branch1); // branch index
							if (volume_share.at(count_branch1) < branch_volume.at(branch_aux1)) {
								// volume in the branch is bigger than outcoming volume
								// also avoids dividing by zero

								// Inserts a volume for the affected area and obtains its temperature
								start_pos = (1 - volume_share.at(count_branch1) / branch_volume.at(branch_aux1)) * 100;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the affected volume
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature); //branch_temp_pos_aux, branch_temp_aux);

																															   // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch1);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);
							}

							else {// volume in the branch is smaller than outcoming volume

								  // Inserts a volume for the entire branch and obtains its temperature
								start_pos = 0;
								end_pos = 100;
								temp_action = 0; // 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
								InsertVolum = HydroNet_InsertVolume(new_pos.at(branch_aux1), new_temp.at(branch_aux1), start_pos, end_pos, temp_action, fluid_type, value, branch_volume.at(branch_aux1), max_divisions, size.new_pos.at(branch_aux1), size.new_temp.at(branch_aux1));

								// Reads the temperature in the middle of the branch
								obj_pos = (start_pos + end_pos) / 2;
								temp_aux = HydroNet_GetObjTemperature(obj_pos, InsertVolum.position, InsertVolum.temperature);// branch_temp_pos_aux, branch_temp_aux);

																															  // Enters temperature and mass in the enthalpy balance
								sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share.at(count_branch);
								sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share.at(count_branch1);

								// Stores the excess volume
								overflow(branch_aux1) = volume_share.at(count_branch1) - branch_volume.at(branch_aux1);
							}
						}
					}
				}
			}

			for (int count = 0; count < overflow.size(); count++)
				overflow(count) = 0;// leaves everything as found

									// Stores resulting inlet temperature
			

		}
	}
}




int main() {
	HydroNet_model objeto;

	objeto.HydroNet_Main();

	return 0;
}
