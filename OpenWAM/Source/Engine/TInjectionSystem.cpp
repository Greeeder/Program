/**
 * @file TInjectionSystem.cpp
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
#include "TInjectionSystem.h"

TInjectionSystem::TInjectionSystem() {
	// TODO Auto-generated constructor stub

}

TInjectionSystem::~TInjectionSystem() {
	// TODO Auto-generated destructor stub
}

void TInjectionSystem::addInjector(TInjector_ptr inj){
	Injector.push_back(inj);
}
