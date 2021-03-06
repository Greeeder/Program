/*// HydroNet_FlowSolver.cpp : Defines the entry point for the console application.
//*/

#include "stdafx.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Dense>
#include <string>


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

struct strObjects
{
	int ID;							// Object's idintificator
	std::string name;				// Object's name
	std::string Class;				// Object's class
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

std::vector<double> roots(double a, double b, double c) {			//Roots a*x^2+b*x+c=0
	double r1, r2;
	std::vector<double> r(2);
	r1 = abs((-b + sqrt(b*b - 4 * a*c)) / (2 * a));
	r2 = abs((-b - sqrt(b*b - 4 * a*c)) / (2 * a));

	r.at(0) = std::max(r1, r2);
	r.at(1) = std::min(r1, r2);

	return r;
}

class Solver {
	// Functor
private:
	std::vector<double> A, B, C, D;
	int n;



public:
	Solver(Eigen::VectorXd x0, Eigen::VectorXd data) {

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



std::vector<double> NonLinearSolver(Eigen::VectorXd _sol,std::vector<std::vector<double>> coef1_mat, std::vector<std::vector<double>> coef2_mat, std::vector<std::vector<double>> const_mat,std::vector<double> const_vec)
{
	/*
	A is an independent constant vector
	B is a vector dependent on the solution sing
	C ia a vector containing the lineal coefficients
	D ia a vector containing the quadratic coefficients
	Y=D*X*|X|+C*X+B*X/|X|+A n-dimensional
	*/

	int n= _sol.size();
	std::vector <double> sol;
	int info;
	Eigen::VectorXd data(10 * n);
	//int count = 0;
	std::vector<double> _data;
	_data.resize(10 * n);


	//The matrix an vector's coeficients are stored in a single Eigen vector
	for (int count=0; count < n; count++) {
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
			_sol[j] = _sol[j]/abs(_sol[j]*1E-9);

	//Computation
	Solver functor(_sol, data);
	/*Eigen::LevenbergMarquardt<Solver,double> lm(functor);
	lm.parameters.xtol = 1e-16;
	lm.parameters.ftol = 1e-16;
	info = lm.minimize(x);*/
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

struct strSize
{
	int unk_ind_aux;
	int remain_branch;
};


void HydroNet_FlowSolver(std::vector<strBranches> branches, std::vector<strPump_volum> pump_volum, std::vector<std::vector<int>>branches_id, std::vector<double> bound_flows, std::vector<int> nodes_id, std::vector<std::vector<int>> mesh_branches, std::vector<std::vector<int>> node_branches, int n_mesh, int n_branch, std::vector<double> head_loss, std::vector<double> hydr_resist1, std::vector<double> hydr_resist2, int n_tanks) {
	
	std::vector<double> flows;
	//std::vector<bool> solved;
	Eigen::RowVectorXi solved(n_branch);
	std::vector<bool> trues;
	//std::vector<bool> volum_pump;
	Eigen::RowVectorXi volum_pump(n_branch);
	std::vector<double> head_vol_pump;
	int n_nodes;
	//std::vector<bool> solv_nodes;
	Eigen::RowVectorXi solv_nodes(nodes_id.size());
	bool while_flag;
	//std::vector<bool> solv_prev1;
	//std::vector<bool> solv_prev2;
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
	std::vector<bool> volpump_ind_aux,flows_ind_aux;
	int last_row;
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


	/*
	 Solver for hydraulic networks.
		 Receives information about meshes, branches and nodes in the the network.
		 Receives known flows imposed as boundary conditions - a vector with the flow values in the corresponding branches and NaN in the rest.
		 Receives head losses in every branch.Data must be divided in three vectors: one for head losses that do not depend on the flow, another one for head losses proportional to the flow and the last one for losses proportional to the flow squared.
		 This function attempts to solve the network with algebraic methods.If some flows cannot be solved, it uses a numerical method.
		 Receives information about volumetric pumps.If there are volumetric pumps in the circuit, the fucntion returns its pump head.



		 SOLVING PROCEDURE

		 First, solves branches where flow is known, is zero or can be calculated independently because there is a volumetric pump(flow depends on pump speed).Check if all branches are solved.

		 Second, enters a loop :
	 NODES STEP : Solves nodes where there is only one unknown flow.Check if all branches are solved.Loop until no more nodes are solved.
		 If it is not the first iteration and no branches have been solved in the nodes step, break main loop.
		 MESH STEP : solves meshes that depend only on one flow(simple meshes; only one branch remains to be solved).Check if all branches are solved.Loop until no more meshes are solved.
		 If no branches have been solved in the mesh step, break main loop.
		 Loop

		 Third, sets up and solves the system of non - linear equations.Equations  depend on flow(nodes, etc.), on flow ^ 2 (pipes, etc.) and on head of volumetric pumps.All branches are solved.

		*/

		// Creates result vectors

	flows.assign(n_branch, 0);						// stores flow values
	solved.setConstant(false);					// stores solving status
	unk_ind_aux.resize(n_branch);
	remain_branch.resize(n_branch);
	size.remain_branch = 0;
	size.unk_ind_aux = 0;


// 1. INDEPENDENT FLOWS

	//Assigns flows entered as a boundary condition

	for (int count=0; count < n_branch; count++)
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
			if ((branches.at(branch_aux).objects.at(object_aux).Class == p) && (solved(branch_aux) == false) ){		//Enters when the object is a pump_volum										// the branch is not solved yet
				volum_pump(branch_aux) = true;																										// there is a volumetric pump
				aux = 0;
				while ((aux+1 < pump_volum.size())&&(branches.at(branch_aux).objects.at(object_aux).ID != pump_volum.at(aux).Pump_volum))				// Looks up for the pump at the branch in the pump_volum vector
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

					for (int count1=0,aux; count1 < node_branches.at(count).size(); count1++) {
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
							total_head = total_head + (head_loss.at(count1) + hydr_resist1.at(count1) * fabs(flows.at(count1)) + hydr_resist2.at(count1) * flows.at(count1) * flows.at(count1) )* (double)(-2 * change_sign + 1);
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
							flows.at(ind_aux) = sol_aux * sign_aux*(-1);																		 // add sign now
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
	
	n_unknown = n_branch - solved.sum()+ volum_pump.sum(); // number of unknown variables :
// flows - solved flows + volumetric pump head
	// in branches with volumetric pumps, pump head is unknown

	// Coefficient matrices
	coef2_mat.resize(n_unknown); // coefficients for unknowns ^ 2 in each equation(row)
	coef1_mat.resize(n_unknown); // coefficients that multiply unknowns in each equation(row)
	const_mat.resize(n_unknown); // matrix for independent head loss in unsolved branches

	for (int count = 0; count < n_unknown; count++) {
		/*coef2_mat.resize(n_unknown); 
		coef1_mat.resize(n_unknown); 
		const_mat.resize(n_unknown);*/
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
	
	// Auxiliary vector that contains 1 at the indices where the vector
		// 'unknowns' or 'X' refer to flows
	//volpump_ind_aux.flip;

	//flows_ind_aux = volpump_ind_aux;
	
	//volpump_ind_aux.flip;
	
			


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
	//solv_nodes.flip();
	for (int count = 0; count < (solv_nodes.size() - solv_nodes.sum()); count++) { // unsolved nodes loop
	

		// Continuity equation : sum of flows in a node = 0
		for (int count1 = 0, aux; count1 < node_branches.at(count).size() ; count1++) {// branches connected to node(branch index)

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
		last_row = count+1; // save row to continue entering the equations
	}
	//solv_nodes.flip();

	// MESH EQUATIONS
	
		// Number of independent meshes : branches to solve - nodes to solve + 1 = branches to solve - independent nodes to solve
	//solv_nodes.flip();
	//solved.flip();
	n_idp_mesh = (solved.size() - solved.sum()) - (solv_nodes.size() - solv_nodes.sum()) + volum_pump.sum() - n_tanks;
	//solved.flip();
	//solv_nodes.flip();

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
	//solved.flip();
	for (int count = 0; count < solved.size(); count++) {
		remain_branch.at(size.remain_branch) = (!solved(count) || volum_pump(count));
		size.remain_branch++;
	}
	//solved.flip();


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
		idp_mesh.assign(n_mesh,1);

	idp_mesh_branches.resize(idp_mesh.size());

	for (int count=0; count< idp_mesh.size() ; count++)
			idp_mesh_branches.at(count)=mesh_branches.at(idp_mesh.at(count));


	// Initial values for the non - linear solver
	x0.setConstant(n_unknown, 0.0001);
	/*x0.resize(flows.size()) ; // size equal to number of unknowns
	for (int count = 0; count < solved.size(); count++) {
		if (solved.at(count) == false)
			x0.at(count) = 0.001;
		else
			x0.at(count) == flows.at(count);

	}*/
	// flow 1 liter per second
	// x0 = zeros(n_unknown, 1); // size equal to number of unknowns
	// x0(flows_ind_aux) = 0.001; // flow 1 liter per second
	// x0(volpump_ind_aux) = 10; // head  10 m.c.f
	// in the future, this vector can be feed with values obtained in a
	// previous call or from tables; remember to divide pump head by 10000 or so

	// Formulates mesh equations

	for (int count = 0; count < n_idp_mesh; count++) { // independent mesh loop

		// II Kirchhoff's Law: sum of head increments in a mesh = 0
		last_id = 0;
		for (int count1 = 0,aux, pos; count1 < idp_mesh_branches.at(count).size(); count1++) {// branches in the mesh(branches index)
			
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
				last_id = branches_id.at(aux).at(branches_id.at(aux).size()-1); // refresh other end of branch
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

			else { //if (unk_ind_aux.at(count) = aux){// flow is unknown
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


	//f = @(x) coef2_mat*(x.*abs(x)) + coef1_mat*x + const_mat*(x. / abs(x)) + const_vec;
	// oldopts = optimoptions('fsolve');
	// opts = optimoptions(oldopts,'MaxFunEvals',1000000), 'Display','off');
	Y = NonLinearSolver(x0 , coef1_mat,  coef2_mat,  const_mat,  const_vec); //,opts);


	// Save 
	for (int count = 0; count < solved.size(); count++) {
		for (int count1 = 0; count1 < Y.size(); count1++) {
			if (solved(count) == false && count== unk_ind_aux.at(count1))
				flows.at( count) = Y.at(count1);
		}
	}
	
	///flows(not(solved)) = Y(flows_ind_aux);

	// Save head in volumetric pumps 
	for (int count = 0; count < volum_pump.size(); count++) {
		for (int count1 = 0; count1 < Y.size(); count1++) {
			if (volum_pump(count) == true && volpump_ind_aux.at(count1) == true)
				head_vol_pump.at(count) = Y.at(count1) * 10000;
		}
	}

	///head_vol_pump(volum_pump) = Y(volpump_ind_aux) * 10000;
	// previously pump head was divided by 10000 so its magnitude order was
	// close to that of the flows



	}

	int main()
{
	std::vector<strBranches> branches;
	std::vector<strPump_volum> pump_volum;
	std::vector<std::vector<int>>branches_id;
	std::vector<double> bound_flows;
	std::vector<int> nodes_id;
	std::vector<std::vector<int>> mesh_branches;
	std::vector<std::vector<int>> node_branches;
	int n_mesh;
	int n_branch;
	std::vector<double> head_loss;
	std::vector<double> hydr_resist1;
	std::vector<double> hydr_resist2;
	int n_tanks;

	branches_id.resize(6);

	branches_id.at(0) = { 12,13,1,2,3 };
	branches_id.at(1) = { 12,11,10,8,7};
	branches_id.at(2) = { 7,6 };
	branches_id.at(3) = { 6,4,3 };
	branches_id.at(4) = { 6,14,5,3 };
	branches_id.at(5) = { 12,5,9,7 };

	bound_flows.assign(6,NAN);
	bound_flows.at(1) = 0;

	nodes_id = { 3,6,7,12 };

	mesh_branches.resize(6);

	mesh_branches.at(0) = { 0,1,2,3 };
	mesh_branches.at(1) = { 0,1,2,4 };
	mesh_branches.at(2) = { 1,5 };
	mesh_branches.at(3) = { 0,5,2,3 };
	mesh_branches.at(4) = { 0,5,2,4 };
	mesh_branches.at(5) = { 3,4 };

	node_branches.resize(4);

	node_branches.at(0) = { 0,3,4 };
	node_branches.at(1) = { 2,3,4};
	node_branches.at(2) = { 1,2,5 };
	node_branches.at(3) = { 0,1,5 };

	n_mesh = 6;

	n_branch = 6;

	head_loss = { 0,0,0,0,0,1 };

	hydr_resist1 = { 0,0,0,0,0,0 };

	hydr_resist2 = { 2.120926015593668E06,1.016417857504270E06,0,50000,100000,5.1642E5};

	n_tanks = 0;

	pump_volum.resize(6);

	pump_volum.at(0).Pump_volum = 1;
	pump_volum.at(0).coef0 = 1;
	pump_volum.at(0).coef1 = 0;
	pump_volum.at(0).coef2 = 0;
	pump_volum.at(0).coef3 = 0;
	pump_volum.at(0).head_max = 1000000;
	pump_volum.at(0).inlet_object = 1;
	pump_volum.at(0).outlet_oject = 2;
	pump_volum.at(0).pump_speed = 105;

	branches.resize(6);
	
	
	branches.at(0).branch_ind = 1;
	branches.at(1).branch_ind = 2;
	branches.at(2).branch_ind = 3;
	branches.at(4).branch_ind = 4;
	branches.at(4).branch_ind = 5;
	branches.at(5).branch_ind = 6;

	
	branches.at(0).objects.resize(3);
	branches.at(1).objects.resize(3);
	branches.at(2).objects.resize(1);
	branches.at(3).objects.resize(1);
	branches.at(4).objects.resize(2);
	branches.at(5).objects.resize(2);

	branches.at(0).objects.at(0).Class = "pump_volum";
	branches.at(0).objects.at(0).ID = 1;
	branches.at(0).objects.at(0).name = "Main_Pump";
	branches.at(0).objects.at(0).volume = 1E-3;
	branches.at(0).objects.at(0).inlet_object = 13;
	branches.at(0).objects.at(0).outlet_oject = 2;

	branches.at(0).objects.at(1).Class = "pipe";
	branches.at(0).objects.at(1).ID = 2;
	branches.at(0).objects.at(1).name = "Enigine_Inlet";
	branches.at(0).objects.at(1).volume = 8.8357E-5;
	branches.at(0).objects.at(1).inlet_object = 1;
	branches.at(0).objects.at(1).outlet_oject = 3;

	branches.at(0).objects.at(2).Class = "pipe";
	branches.at(0).objects.at(2).ID = 13;
	branches.at(0).objects.at(2).name = "Return_Pipe";
	branches.at(0).objects.at(2).volume = 6.2832E-4;
	branches.at(0).objects.at(2).inlet_object = 12;
	branches.at(0).objects.at(2).outlet_oject = 1;

	branches.at(1).objects.at(0).Class = "thermostat";
	branches.at(1).objects.at(0).ID = 8;
	branches.at(1).objects.at(0).name = "Thermostat_rad";
	branches.at(1).objects.at(0).volume = 0;
	branches.at(1).objects.at(0).inlet_object = 7;
	branches.at(1).objects.at(0).outlet_oject = 10;

	branches.at(1).objects.at(1).Class = "heat_exch_fix";
	branches.at(1).objects.at(1).ID = 10;
	branches.at(1).objects.at(1).name = "Radiator";
	branches.at(1).objects.at(1).volume = 5E-4;
	branches.at(1).objects.at(1).inlet_object = 8;
	branches.at(1).objects.at(1).outlet_oject = 11;

	branches.at(1).objects.at(2).Class = "pipe";
	branches.at(1).objects.at(2).ID = 11;
	branches.at(1).objects.at(2).name = "Radiator_outlet";
	branches.at(1).objects.at(2).volume = 3.1416E-4;
	branches.at(1).objects.at(2).inlet_object = 10;
	branches.at(1).objects.at(2).outlet_oject = 12;

	branches.at(2).objects.at(0).Class = "node";
	
	
	

	branches.at(3).objects.at(0).Class = "heat_exch_fix";
	branches.at(3).objects.at(0).ID = 4;
	branches.at(3).objects.at(0).name = "Block_cyl1";
	branches.at(3).objects.at(0).volume = 5E-4;
	branches.at(3).objects.at(0).inlet_object = 3;
	branches.at(3).objects.at(0).outlet_oject = 6;

	branches.at(4).objects.at(0).Class = "heat_exch_fix";
	branches.at(4).objects.at(0).ID = 5;
	branches.at(4).objects.at(0).name = "Head_cyl1";
	branches.at(4).objects.at(0).volume = 5E-4;
	branches.at(4).objects.at(0).inlet_object = 3;
	branches.at(4).objects.at(0).outlet_oject = 14;

	branches.at(4).objects.at(1).Class = "heat_exch_fix";
	branches.at(4).objects.at(1).ID = 14;
	branches.at(4).objects.at(1).name = "Head_cyl2";
	branches.at(4).objects.at(1).volume = 5E-4;
	branches.at(4).objects.at(1).inlet_object = 5;
	branches.at(4).objects.at(1).outlet_oject = 6;

	branches.at(5).objects.at(0).Class = "thermostat";
	branches.at(5).objects.at(0).ID = 9;
	branches.at(5).objects.at(0).name = "Thermostat_byp";
	branches.at(5).objects.at(0).volume = 0;
	branches.at(5).objects.at(0).inlet_object = 7;
	branches.at(5).objects.at(0).outlet_oject = 15;

	branches.at(5).objects.at(1).Class = "pipe";
	branches.at(5).objects.at(1).ID = 5;
	branches.at(5).objects.at(1).name = "Bypass_Pipe";
	branches.at(5).objects.at(1).volume = 3.1416E-4;
	branches.at(5).objects.at(1).inlet_object = 9;
	branches.at(5).objects.at(1).outlet_oject = 12;

	HydroNet_FlowSolver(branches, pump_volum, branches_id, bound_flows, nodes_id, mesh_branches, node_branches, n_mesh, n_branch, head_loss, hydr_resist1, hydr_resist2, n_tanks);

    return 0;
}

