#include "Globales.h"
#include "TPipe.hpp"
#include "AllBoundaryConditions.hpp"
#include "LaxWendroff.hpp"
#include "AllExtraSourceTerms.hpp"
#include "AllPipeMethods.hpp"
#include <chrono>

namespace __Ambient {
	double p_Pa = 10000;	   //!< The ambient pressure [Pa]
	double p_bar = 1;		   //!< The ambient pressure [bar]
	double T_K = 300;		   //!< The ambient temperature [K]
	double T_degC = 27;		   //!< The ambient temperature [degC]
	double HR = 50;			   //!< The humidity [%]
}
;

void LW() {
	RowVector x(2);
	RowVector A(2);
	double dx;
	Pipe_ptr pipe = make_shared<TPipe>();
	int n_nodes = 100;
	RowVector p(n_nodes);
	RowVector T(n_nodes);
	RowVector u(n_nodes);
	IOFormat format(StreamPrecision, DontAlignCols, "\t", "\n", "", "", "", "");

	x(0) = 0.;
	x(1) = 1.;
	A.setConstant(1.);
	dx = x(1) / (n_nodes - 1);

	p.setConstant(5E5);
	T.setConstant(1200.);
	u.setConstant(0.);
	p.tail(n_nodes / 2).setConstant(1E5);
	T.tail(n_nodes / 2).setConstant(300.);

	pipe->setGeometry(x, dx, A);
	use_laxwendroff(pipe);
	pipe->setPTU(p, T, u);
	close_pipe_end(pipe, nmLeft);
	close_pipe_end(pipe, nmRight);
	set_null_source(pipe);
// 	auto start = chrono::steady_clock::now();
	cout << pipe->getCurrentTime() << "\t"
		 << pipe->getPressure().format(format) << endl;
	while(pipe->getCurrentTime() < 5E-3) {
		pipe->setTimeStep(pipe->getMaxTimeStep());
		pipe->Solve();
		cout << pipe->getCurrentTime() << "\t"
			 << pipe->getPressure().format(format) << endl;
	}
// 	auto end = chrono::steady_clock::now();
// 	auto diff = end - start;
// 	cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}


void G() {
	RowVector x(2);
	RowVector A(2);
	double dx;
	double dt;
	Pipe_ptr pipe = std::make_shared<TPipe> ();
	int n_nodes = 100;
	RowVector p(n_nodes);
	RowVector T(n_nodes);
	RowVector u(n_nodes);
	IOFormat format(StreamPrecision, DontAlignCols, "\t", "\n",
					"", "", "", "");

	x(0) = 0.;
	x(1) = 1.;
	A.setConstant(1.);
	dx = (x(1) - x(0)) / n_nodes;

	p.setConstant(5E5);
	T.setConstant(1200.);
	u.setConstant(0.);
	p.tail(n_nodes / 2).setConstant(1E5);
	T.tail(n_nodes / 2).setConstant(300.);

	pipe->setGeometry(x, dx, A);
	use_godunov(pipe, &KT);
	pipe->setPTU(p, T, u);
	close_pipe_end(pipe, nmLeft);
	close_pipe_end(pipe, nmRight);
	set_null_source(pipe);

	cout << pipe->getCurrentTime() << "\t"
		 << pipe->getPressure().format(format) << endl;
// 	auto start = chrono::steady_clock::now();
	while(pipe->getCurrentTime() < 5E-3) {
		pipe->setTimeStep(pipe->getMaxTimeStep());
		pipe->Solve();
		cout << pipe->getCurrentTime() << "\t"
			 << pipe->getPressure().format(format) << endl;
	}
// 	auto end = chrono::steady_clock::now();
// 	auto diff = end - start;
// 	cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

void G2P() {
	RowVector x(2);
	RowVector A(2);
	double dx;
	double dt;
	Pipe_ptr left_pipe = std::make_shared<TPipe>();
	Pipe_ptr right_pipe = std::make_shared<TPipe>();
	BoundaryCondition_ptr left_wall;
	BoundaryCondition_ptr right_wall;
	int n_nodes = 50;
	RowVector p(n_nodes);
	RowVector T(n_nodes);
	RowVector u(n_nodes);
	IOFormat format(StreamPrecision, DontAlignCols, "\t", "\n", "", "", "", "");

	x(0) = 0.;
	x(1) = 0.5;
	A.setConstant(1.);
	dx = (x(1) - x(0)) / (n_nodes);

	p.setConstant(5E5);
	T.setConstant(1200.);
	u.setConstant(0.);
	left_pipe->setGeometry(x, dx, A);
	use_godunov(left_pipe, &KT);
	left_pipe->setPTU(p, T, u);

	x(0) = 0.5;
	x(1) = 1;
	A.setConstant(1.);
	dx = (x(1) - x(0)) / (n_nodes);

	p.setConstant(1E5);
	T.setConstant(300.);
	u.setConstant(0.);
	right_pipe->setGeometry(x, dx, A);
	use_godunov(right_pipe, &KT);
	right_pipe->setPTU(p, T, u);

	close_pipe_end(left_pipe, nmLeft);
	close_pipe_end(right_pipe, nmLeft);
	attach_pipes(left_pipe, nmRight, right_pipe, nmRight);
	set_null_source(left_pipe);
	set_null_source(right_pipe);

	cout << left_pipe->getCurrentTime() << "\t" << left_pipe->getPressure().format(format) << "\t" <<
		 right_pipe->getPressure().format(format) << endl;
	while(left_pipe->getCurrentTime() < 5E-4) {
		dt = min(left_pipe->getMaxTimeStep(), right_pipe->getMaxTimeStep());
		left_pipe->setTimeStep(dt);
		right_pipe->setTimeStep(dt);
		left_pipe->Solve();
		right_pipe->Solve();
		cout << left_pipe->getCurrentTime() << "\t" << left_pipe->getPressure().format(format) << "\t" <<
			 right_pipe->getPressure().format(format) << endl;
	}
}


void MUSCL(Limiter_pt lim) {
	RowVector x(2);
	RowVector A(2);
	double dx;
	double dt;
	Pipe_ptr pipe = std::make_shared<TPipe> ();
	TPipe * foo = pipe.get();
	int n_nodes = 100;
	RowVector p(n_nodes);
	RowVector T(n_nodes);
	RowVector u(n_nodes);
	IOFormat format(StreamPrecision, DontAlignCols, "\t", "\n",
					"", "", "", "");

	x(0) = 0.;
	x(1) = 1.;
	A.setConstant(1.);
	dx = (x(1) - x(0)) / n_nodes;

	p.setConstant(5E5);
	T.setConstant(1200.);
	u.setConstant(0.);
	p.tail(n_nodes / 2).setConstant(1E5);
	T.tail(n_nodes / 2).setConstant(300.);

	pipe->setGeometry(x, dx, A);
	use_muscl(pipe, &KT, lim);
	pipe->setPTU(p, T, u);
	close_pipe_end(pipe, nmLeft);
	close_pipe_end(pipe, nmRight);
	set_null_source(pipe);

// 	cout << pipe->getCurrentTime() << "\t"
// 		<< pipe->getPressure().format(format) << endl;
	auto start = chrono::steady_clock::now();
	while(pipe->getCurrentTime() < 5E-3) {
		pipe->setTimeStep(pipe->getMaxTimeStep());
		pipe->Solve();
// 		cout << pipe->getCurrentTime() << "\t"
// 			<< pipe->getPressure().format(format) << endl;
	}
	auto end = chrono::steady_clock::now();
	auto diff = end - start;
	cout << chrono::duration <double, milli> (diff).count() << " ms" << endl;
}

int main() {
// 	LW();
// 	G();
// 	G2P();
	MUSCL(&Minmod);
// 	MUSCL(&VanLeer);
// 	MUSCL(&Superbee);

	return EXIT_SUCCESS;
}
