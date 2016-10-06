/* --------------------------------------------------------------------------------*\
 ==========================|
 \\   /\ /\   // O pen     | OpenWAM: The Open Source 1D Gas-Dynamic Code
  \\ |  X  | //  W ave     |
   \\ \/_\/ //   A ction   | CMT-Motores Termicos / Universidad Politecnica Valencia
    \\/   \//    M odel    |
 ----------------------------------------------------------------------------------
 License

 This file is part of OpenWAM.

 OpenWAM is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenWAM is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.


 \*--------------------------------------------------------------------------------*/

/**
 * @file TIntegrable.cpp
 * @author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
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
 * This file defines general integrable objects.
 */

#include "TIntegrable.hpp"
#include "iostream"
#include "cstdlib"

TIntegrable::TIntegrable() {
	FInsOutput = false;
	FCurrentTime = 0;
	FMaxTimeStep = 1000;
	FTimeStep = 1;
}

TIntegrable::~TIntegrable() {}

double TIntegrable::getCurrentTime() const {
	return FCurrentTime;
}

bool TIntegrable::getInsOutputStatus() const {
	return FInsOutput;
}


std::string TIntegrable::getName() const {
	return FName;
}

double TIntegrable::getTimeStep() const {
	return FTimeStep;
}

void TIntegrable::ReadInsResults(const pugi::xml_node& node) {}

void TIntegrable::setName(const std::string & name) {
	FName = name;
}

void TIntegrable::setTimeStep(double dt) {
	FTimeStep = dt;
}

void TIntegrable::WriteInsHeader(std::stringstream& output) const {
	output << "\t" << getName();
}

void TIntegrable::WriteInsHeader(std::vector<std::string>& output) const
{
	output.push_back(getName());
}

void TIntegrable::WriteInsResults(std::stringstream& output) const {
	output << "\t" << 0;
}

void TIntegrable::WriteInsResults(std::vector<float>& output) const
{
	output.push_back(0.);
}

TIntegrableGroup::TIntegrableGroup() {}

TIntegrableGroup::TIntegrableGroup(const std::vector<Integrable_ptr>& members):
	TIntegrableGroup() {
	FMembers = members;
}

void TIntegrableGroup::IntegrateWithoutUpdating() {
	for (auto member : FMembers) {
		member->IntegrateWithoutUpdating();
	}
}

double TIntegrableGroup::getMaxTimeStep() {
	FMaxTimeStep = 1E10;
	for (auto member : FMembers) {
		FMaxTimeStep = std::min(FMaxTimeStep, member->getMaxTimeStep());
	}
	return FMaxTimeStep;
}

std::vector< Integrable_ptr > TIntegrableGroup::getMembers() const {
	return FMembers;
}

void TIntegrableGroup::setTimeStep(double dt) {
	TIntegrable::setTimeStep(dt);
	for (auto member : FMembers) {
		member->setTimeStep(dt);
	}
}

void TIntegrableGroup::Solve() {
	IntegrateWithoutUpdating();
	UpdateStateVector();
	FCurrentTime += FTimeStep;
}

void TIntegrableGroup::Solve(double t) {
	while ((t - FCurrentTime) > 0.) {
		FMaxTimeStep = getMaxTimeStep();
		if (FMaxTimeStep > (t - FCurrentTime)) {
			setTimeStep(t - FCurrentTime);
		} else {
			setTimeStep(FMaxTimeStep);
		}
		Solve();
	}
}

void TIntegrableGroup::UpdateStateVector() {
	for (auto member : FMembers) {
		member->UpdateStateVector();
	}
}

TParallelIntegrableGroup::TParallelIntegrableGroup() {}

TParallelIntegrableGroup::TParallelIntegrableGroup(
	const std::vector<Integrable_ptr>& members) {
	FMembers = members;
}

void TParallelIntegrableGroup::IntegrateWithoutUpdating() {
	std::vector<std::thread> workers;
	auto fun = [&] (Integrable_ptr member) {
		member->IntegrateWithoutUpdating();};
	for (auto member: FMembers) {
		workers.push_back(std::thread(fun, member));
	}
	for (auto i = 0; i < workers.size(); i++) {
		workers[i].join();
	}
}
