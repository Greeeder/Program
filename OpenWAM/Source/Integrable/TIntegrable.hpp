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

/*!
 * \file TIntegrable.hpp
 * \author Luis Miguel Garcia-Cuevas Gonzalez <luiga12@mot.upv.es>
 *
 * \section LICENSE
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
 * \section DESCRIPTION
 * This file declares general integrable objects.
 */

#ifndef TIntegrable_hpp
#define TIntegrable_hpp

#include <memory>
#include <vector>
#include <sstream>

#include <algorithm>
#include <mutex>
#include <thread>

#include "pugixml.hpp"

/*!
 * \brief A time-integrable object.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 */
class TIntegrable {
  protected:
	double FCurrentTime; ///< Current time for the pipe. [s]
	double FMaxTimeStep; ///< Maximum allowable time-step. [s]
	double FTimeStep; ///< Time step. [s]
	std::string FName; //!< Object name.
	mutable std::mutex mtx; //!< Mutex, for critical sections in multithreading.
	bool FInsOutput; //!< Instantaneous output switch.

  public:
	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 */
	TIntegrable();

	/*!
	 * \brief Destructor.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/22
	 * 
	 * Destructs the object.
	 *
	 */
	virtual ~TIntegrable();

	/*!
	 * \brief Gets the current time.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Current time. [s]
	 */
	virtual double getCurrentTime() const;

	/*!
	 * \brief Returns whether the object writes instantaneous results.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return True if the object writes instantaneous results.
	 */
	bool getInsOutputStatus() const;

	/*!
	 * \brief Returns the maximum allowable time-step.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Returns the maximum allowable time-step due to stability criteria.
	 *
	 * \return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep() = 0;

	/*!
	 * \brief Returns the object name.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/17
	 * 
	 * \return Object name;
	 */
	std::string getName() const;

	/*!
	 * \brief Gets the current time-step.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \return Current time-step. [s]
	 */
	double getTimeStep() const;

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void IntegrateWithoutUpdating() = 0;

	/*!
	 * \brief Reads the instantaneous results setup from an XML node.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param node XML node with the object data.
	 */
	virtual void ReadInsResults(const pugi::xml_node& node);

	/*!
	 * \brief Sets the object name.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/17
	 * 
	 * \param Object name;
	 */
	void setName(const std::string & name);

	/*!
	 * \brief Sets the time-step for the integration of this pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param dt Time-step. [s]
	 */
	virtual void setTimeStep(double dt);

	/*!
	 * \brief Integrate the flow, updating the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct, updating the state vector
	 * at the end. Uses an already-set time-step.
	 */
	virtual void Solve() = 0;

	/*!
	 * \brief Integrate the flow, updating the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct, updating the state
	 * vector at the end. Computes until the time passed is reached.
	 *
	 * \param t Time at the end of the integration. [s]
	 */
	virtual void Solve(double t) = 0;

	/*!
	 * \brief Updates the state vector with the already-computed results.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * After solving the system and computing the state vector update, this
	 * function is called to update it.
	 */
	virtual void UpdateStateVector() = 0;

	/*!
	 * \brief Writes the instantaneous results header.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsHeader(std::stringstream & output) const;

	/*!
	* \brief Writes the instantaneous results header.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*/
	virtual void WriteInsHeader(std::vector<std::string>& output) const;

	/*!
	 * \brief Writes the instantaneous results.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void WriteInsResults(std::stringstream & output) const;

	/*!
	* \brief Writes the instantaneous results.
	*
	* \author L.M. García-Cuevas <luiga12@mot.upv.es>
	*/
	virtual void WriteInsResults(std::vector<float>& output) const;
};

/*!
 * \brief Shared pointer to a TIntegrable object.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 */
typedef std::shared_ptr<TIntegrable> Integrable_ptr;

class TIntegrableGroup: public TIntegrable {
  protected:
	std::vector<Integrable_ptr> FMembers; ///< Vector of members of the group.

  public:

	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 */
	TIntegrableGroup();

	/*!
	 * \brief Constructs the group with a vector of members.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 * 
	 * \param members Vector of members.
	 */
	TIntegrableGroup(const std::vector<Integrable_ptr> & members);
	
	/*!
	 * \brief Returns the maximum allowable time-step.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Returns the maximum allowable time-step due to stability criteria.
	 *
	 * \return Maximum allowable time-step. [s]
	 */
	virtual double getMaxTimeStep();

	/*!
	 * \brief Returns the group members.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/17
	 * 
	 * \return The group members.
	 */
	std::vector<Integrable_ptr> getMembers() const;

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 */
	virtual void IntegrateWithoutUpdating();

	/*!
	 * \brief Sets the time-step for the integration of this pipe.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * \param dt Time-step. [s]
	 */
	virtual void setTimeStep(double dt);

	/*!
	 * \brief Integrate the flow, updating the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct, updating the state vector
	 * at the end. Uses an already-set time-step.
	 */
	virtual void Solve();

	/*!
	 * \brief Integrate the flow, updating the state vector.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * Integrates the flow evolution inside the duct, updating the state
	 * vector at the end. Computes until the time passed is reached.
	 *
	 * \param t Time at the end of the integration. [s]
	 */
	virtual void Solve(double t);

	/*!
	 * \brief Updates the state vector with the already-computed results.
	 *
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * 
	 * After solving the system and computing the state vector update, this
	 * function is called to update it.
	 */
	virtual void UpdateStateVector();
};

/*!
 * \brief Shared pointer to a TIntegrableGroup object.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 */
typedef std::shared_ptr<TIntegrableGroup> IntegrableGroup_ptr;

/*!
 * \brief A group of TIntegrables that can be computed in parallel.
 * 
 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
 * \date 2016/03/23
 */
class TParallelIntegrableGroup: public TIntegrableGroup {
  public:
	/*!
	 * \brief Default constructor.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 */
	TParallelIntegrableGroup();

	/*!
	 * \brief Constructs the group with a vector of members.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 * 
	 * \param members Vector of members.
	 */
	TParallelIntegrableGroup(const std::vector<Integrable_ptr> & members);

	/*!
	 * \brief Integrates one time-step without updating the state vector.
	 * 
	 * \author L.M. García-Cuevas <luiga12@mot.upv.es>
	 * \date 2016/03/23
	 * 
	 * Integrates without updating the state vector, computing everything
	 * in parallel.
	 */
	virtual void IntegrateWithoutUpdating();
};

#endif
