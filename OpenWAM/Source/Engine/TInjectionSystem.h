/**
 * @file TInjectionSystem.h
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
 * @date 11 de mar. de 2016
 *
 * @section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenWAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 * TODO Insert here the description.
 */
#ifndef SOURCE_ENGINE_TINJECTIONSYSTEM_H_
#define SOURCE_ENGINE_TINJECTIONSYSTEM_H_

#include "TInjector.h"

class TInjectionSystem {
private:

	double InjectionPressure;
	double SOI;
	double FuelMass;

	std::vector<TInjector_ptr> Injector;
public:
	TInjectionSystem();
	virtual ~TInjectionSystem();

	void addInjector(TInjector_ptr inj);

	void setInjectionPressure(double IP){
		InjectionPressure = IP;
	}

	double getFuelInjected(int i){
		return Injector[i]->getFuelInjected();
	}

	double* getInjectionPressure_ptr() {
		return &InjectionPressure;
	}

	double* getSOI_ptr() {
		return &SOI;
	}

	double* getFuelMass_ptr() {
		return &FuelMass;
	}
};

typedef unique_ptr<TInjectionSystem> InjectionSystem_ptr;

#endif /* SOURCE_ENGINE_TINJECTIONSYSTEM_H_ */
