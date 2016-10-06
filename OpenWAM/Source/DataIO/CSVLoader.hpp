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
 * \file CSVLoader.hpp
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
 * This file declares functions for data loading from csv-like files.
 */

#ifndef CSVLoader_hpp
#define CSVLoader_hpp

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include "Math_wam.h"

/*!
 * \brief A map of data.
 */
typedef std::map<std::string, RowVector> DataMap;

/*!
 * \brief Gets the number of points in a datafile.
 * 
 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
 * \date 30/05/2016
 *
 * \param file_name File name.
 * \return The number of points in the data file.
 */
int get_datafile_number_of_points(const std::string& file_name);

/*!
 * \brief Loads data from a comma-separated-values file.
 * 
 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
 * \date 30/05/2016
 * 
 * Loads data from a file with fields separated by commas.
 * 
 * \param file_name File name.
 * \return A DataMap with the file data.
 */
DataMap load_csv(const std::string& file_name);

/*!
 * \brief Loads data from a file.
 * 
 * \author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
 * \date 30/05/2016
 * 
 * Loads data from a file with fields separated by a custom separator.
 * 
 * \param file_name File name.
 * \param separator Field separator.
 * \return A DataMap with the file data.
 */
DataMap load_table(const std::string& file_name, char separator);

#endif
