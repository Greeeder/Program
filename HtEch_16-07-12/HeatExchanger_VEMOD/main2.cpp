#include <iostream>
#include<fstream>

using namespace std;

int main() {
	
	int block;
	int count;
	int lines;
	std::ifstream HeatExchanger_file;
	string line_scan;
	
	HeatExchanger_file.open("HeatExchanger_efficiency_table.txt", std::ios::out);
	
	line = 0;
	
	while (!HeatExchanger_file.eof()) {
		getline(HeatExchanger_file,line_scan);
		lines++;	// lines counter
	}
	
	std::string efficiency_data[line];
	std::string NTu_data[line];
	std::string Cr_data[lines];
	
	block = 0;
	for (count = 0; count < lines; count++) {
		getline(HeatExchanger_file, NTu_data[block],'\t');
		
		getline(HeatExchanger_file, Cr_data[block],'\t');
		
		getline(HeatExchanger_file, efficiency_data[block],'\n');
	
		block++;
	}
	HeatExchanger_file.close();
	
	
	cout <<"\n";
	cout <<line;
	cout <<"\n";
	cout << "fin";

	return 0;
}
