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
 * \file CSVLoader.cpp
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
 * This file defines functions for data loading from csv-like files.
 */

#include "CSVLoader.hpp"

int get_datafile_number_of_points(const std::string& file_name)
{
	std::ifstream file(file_name);
	std::string line;
	int number_of_points = -2;
	if (file.is_open())
	{
		while(file.good())
		{
			getline(file, line);
			number_of_points++;
		}
		file.close();
	}
	return number_of_points;
}

DataMap load_csv(const std::string& file_name) {
	return load_table(file_name, ',');
}

DataMap load_table(const std::string& file_name, char separator) {
	DataMap data;
	int number_of_points = get_datafile_number_of_points(file_name);
	std::ifstream file(file_name);
	std::string line, field;
	std::stringstream line_ss, x_ss;
	vector<std::string> header;
	double x;
	int f, i = 0;
	vector<std::string>::size_type j = 0;
	if (file.is_open())
	{
		getline(file, line);
		line_ss.str(line);
		while (line_ss.good())
		{
			getline(line_ss, field, separator);
			f = field.find('"');
			while (f > -1)
			{
				field.erase(field.begin() + f);
				f = field.find('"');
			}
			f = field.find('\r');
			while (f > -1)
			{
				field.erase(field.begin() + f);
				f = field.find('\r');
			}
			data[field] = RowVector(number_of_points);
			header.push_back(field);
		}
		while (file.good() & (i < number_of_points))
		{
			getline(file, line);
			line_ss.clear();
			line_ss.str(line);
			j = 0;
			while (line_ss.good() & (j < header.size()))
			{
				getline(line_ss, field, separator);
				x_ss.clear();
				x_ss.str(field);
				x_ss >> x;
				data.at(header[j])[i] = x;
				j++;
			}
			i++;
		}
	}
	return data;
}
