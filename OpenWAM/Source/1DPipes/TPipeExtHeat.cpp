/*
 * TPipeExtHeat.cpp
 *
 *  Created on: 3 de feb. de 2016
 *      Author: farnau
 */

#include "TPipeExtHeat.h"

TPipeExtHeat::TPipeExtHeat() {
	// TODO Auto-generated constructor stub

}

TPipeExtHeat::~TPipeExtHeat() {
	// TODO Auto-generated destructor stub
}

void TPipeExtHeat::setText(double T) {

	Text = ArrayXd::Ones(ncells) * T;
}

void TPipeExtHeat::setText(ArrayXd T) {

	Text = T;

}
