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
 * @file AllBoundaryConditions.hpp
 * @author Francisco Jose Arnau <farnau@mot.upv.es>
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
 * This file defines a BC creator function.
 */

#include "AllBoundaryConditions.hpp"

void create_BCs(const xml_node & openwam, const vector<Pipe_ptr> & pipes, const vector<BasicPlenum_ptr>& plenum, 
	const TComponentArray_ptr& fluid, const ControlUnit_ptr& ControlUnit, const EngineBlock_ptr& engine) {
	try {
		std::string ConnectionType;
		auto block = GetNodeChild(openwam, "BlockOfConnections");
		for (auto node = GetNodeChild(block, "Bob_NewConnection");
				node; node = node.next_sibling("Bob_NewConnection")) {
			ConnectionType = node.attribute("BoundaryType").value();
			auto ID = IDtoInt(node.attribute("Boundary_ID").as_string());
			if(ConnectionType == "ConstantConditions") {
				auto node_con = GetNodeChild(node, "Con_ConstantConditions");
				auto node_prop = GetNodeChild(node_con, "GasProperties");

				double p, T, u;
				RowVector Y;
				Y.setZero(fluid->size());
				ReadGasProperties(node_prop, p, T, u, Y, fluid);
				auto FluidBC = make_shared<TFluid>(fluid);
				FluidBC->SetComposition(Y);

				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_constant_BC(pipe, nmLeft, p, T, u, FluidBC);
					} else if (pipe->getRightNode() == ID) {
						attach_to_constant_BC(pipe, nmRight, p, T, u, FluidBC);
					}
				}
			} else if (ConnectionType == "ClosedEnd") {
				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						close_pipe_end(pipe, nmLeft);
					} else if (pipe->getRightNode() == ID) {
						close_pipe_end(pipe, nmRight);
					}
				}
			} else if (ConnectionType == "IncidentPressure") {
				auto node_con = GetNodeChild(node, "Con_IncidentPressureBC");
				auto node_prop = GetNodeChild(node_con, "GasProperties");
				auto node_units = GetNodeChild(node_prop, "Units");
				std::string data_file = node_con.attribute("DataFile").as_string();
				RowVector t, Y, p, A_A;
				double dummy;
				if (data_file != "") {
					std::string t_unit = node_units.attribute("Time").as_string();
					std::string p_unit = node_units.attribute("Pressure").as_string();
					load_incident_pressure_data(data_file, '\t', "t", "p", "A_A",
						t_unit, p_unit, &t, &p, &A_A);
				}
				Y.setZero(fluid->size());
				ReadGasProperties(node_prop, dummy, dummy, dummy, Y, fluid);
				auto FluidBC = make_shared<TFluid>(fluid);
				FluidBC->SetComposition(Y);

				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_incident_pressure_BC(pipe, nmLeft, t, p, A_A, FluidBC);
					} else if (pipe->getRightNode() == ID) {
						attach_to_incident_pressure_BC(pipe, nmRight, t, p, A_A, FluidBC);
					}
				}
			} else if (ConnectionType == "Junction") {
				std::vector<Pipe_ptr> junction_pipes;
				std::vector<nmPipeEnd> junction_ends;
				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						junction_pipes.push_back(pipe);
						junction_ends.push_back(nmLeft);
					} else if (pipe->getRightNode() == ID) {
						junction_pipes.push_back(pipe);
						junction_ends.push_back(nmRight);
					}
				}
				if (junction_pipes.size() == 2) {
					attach_pipes(junction_pipes[0], junction_ends[0],
						junction_pipes[1], junction_ends[1]);
				} else {
					std::stringstream message;
					message << "Issues at junction id " <<
						ID << ": number of pipes attached: "
						<< junction_pipes.size();
					throw Exception(message.str());
				}
			} else if (ConnectionType == "ExternalPlenumConnection") {
				auto node_con = GetNodeChild(node, "Con_ExternalPlenumConnection");
				auto node_prop = GetNodeChild(node_con, "GasProperties");
				auto con_ID = GetAttributeAsInt(node_con, "Connection_ID");

				double p, T, u;
				RowVector Y;
				Y.setZero(fluid->size());
				ReadGasProperties(node_prop, p, T, u, Y, fluid);
				auto FluidBC = make_shared<TFluid>(fluid);
				FluidBC->SetComposition(Y);

				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_external_connection(pipe, nmLeft, p, T, FluidBC, u, con_ID);
					} else if (pipe->getRightNode() == ID) {
						attach_to_external_connection(pipe, nmRight, p, T, FluidBC, u, con_ID);
					}
				}
			} else if (ConnectionType == "PressureBC") {
				auto node_con = GetNodeChild(node, "Con_PressureBC");
				auto node_prop = GetNodeChild(node_con, "GasProperties");

				double p, T, u;
				RowVector Y;
				Y.setZero(fluid->size());
				ReadGasProperties(node_prop, p, T, u, Y, fluid);
				auto FluidBC = make_shared<TFluid>(fluid);
				FluidBC->SetComposition(Y);

				for (auto pipe: pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_pressure_BC(pipe, nmLeft, p, T, FluidBC);
					} else if (pipe->getRightNode() == ID) {
						attach_to_pressure_BC(pipe, nmRight, p, T, FluidBC);
					}
				}
			}
			else if (ConnectionType == "PipeToPlenum"){
				auto node_con = GetNodeChild(node, "Con_PipeToPlenum");
				auto plenum_id = IDtoInt(node_con.attribute("PlenumID").as_string());
				auto plm = FindObjectByID(plenum, node_con.attribute("PlenumID").as_string());

				auto valve_id = node_con.attribute("ValveID").as_string();

				auto node_valve = openwam.child("BlockOfValves").find_child_by_attribute("Valve_ID", valve_id);
				string ValveType = node_valve.attribute("Valve_type").as_string();
				TTipoValvula_ptr valve;

				if (ValveType == "FixDC") {
					valve = make_unique<TCDFijo>();
				}
				else if (ValveType == "Stator") {
					valve = make_unique<TEstatorTurbina>();
				}
				else if (ValveType == "Rotor") {
					valve = make_unique<TRotorTurbina>();
				}
				else if (ValveType == "Throttle") {
					valve = make_unique<TMariposa>();
					auto node_act = node_con.child("Actuator");
					if (node_act) {
						int id = IDtoInt(node_act.attribute("ID").as_string());
						string param = node_act.attribute("Parameter").as_string();
						if (param == "Lift") {
							ControlUnit->AppendOutput_ptr(id, dynamic_cast<TMariposa*>(valve.get())->GetLift_ptr());
						}
					}
				}
				else{
					cout << "ERROR: Valve type " << ValveType << " not allowed to be connected to a plenum" << endl;
					cout << "       Plenum ID: " << plenum_id << endl;
				}
				valve->LeeDatosInicialesXML(node_valve, 0, false, NULL);

				for (auto pipe : pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_0Dto1Dconnection(ID, pipe, nmLeft, plm, move(valve));
					}
					else if (pipe->getRightNode() == ID) {
						attach_to_0Dto1Dconnection(ID, pipe, nmRight, plm, move(valve));
					}
				}

			}
			else if (ConnectionType == "IntakeValve") {
				auto node_con = GetNodeChild(node, "Con_IntakeValve");
				auto cyl = FindObjectByID(engine->getCylinders(), node_con.attribute("CylinderID").as_string());

				auto valve_id = node_con.attribute("ValveID").as_string();

				auto node_valve = openwam.child("BlockOfValves").find_child_by_attribute("Valve_ID", valve_id);
				string ValveType = node_valve.attribute("Valve_type").as_string();
				TTipoValvula_ptr valve;

				if (ValveType == "CamValve") {
					valve = make_unique<TValvula4T>();
				}
				else {
					cout << "ERROR: Valve type " << ValveType << " not allowed to be connected to a cylinder" << endl;
					cout << "       Cylinder ID: " << cyl->getID() << endl;
				}
				dynamic_cast<TValvula4T*>(valve.get())->setCylinder(cyl);
				valve->LeeDatosInicialesXML(node_valve, 0, false, NULL);

				for (auto pipe : pipes) {
					if (pipe->getLeftNode() == ID) {
						attach_to_0Dto1Dconnection(ID, pipe, nmLeft, cyl, move(valve));
					}
					else if (pipe->getRightNode() == ID) {
						attach_to_0Dto1Dconnection(ID, pipe, nmRight, cyl, move(valve));
					}
				}

				cyl->appendIntakeValve();
			}
		}
	} catch(exception & N) {
		std::cout << " ERROR : create_BCs" << std::endl;
		std::cout << " Error message : " << N.what() << std::endl;
		throw Exception(N.what());
	}
}
