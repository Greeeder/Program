
#include <vector>
#include "TFluid.h"

class HeatExchanger_model
{
public:
	
	// METHODS

		/*!
		* \brief Reads heat exchanger data from the XML file.
		*
		* \author J. Salvador
		* \date 10/06/2016
		*
		*/
		void ReadHeatExchangerXML();


		/*!
		* \brief Obtains exchanged heat and outlet temperatures of the fluids in the heat exchanger using the efficiency method.
		*
		* \author J. Salvador
		* \date 09/06/2016
		*/
		void calculate_heat_exchange();
	


private:

	// VARIABLES
		
		// General data
		int HeatExchanger_ID;				//!< ID of heat exchanger
		std::string name;					//!< Heat exchanger name
		TFluid::TFluid_ptr fluid_type_1;	//!< Fluid: air, water, oil, etc.
		TFluid::TFluid_ptr fluid_type_2;	//!< Fluid: air, water, oil, etc.
		double UA;							//!< Overall heat transfer coefficient * heat transfer area [W/K]

		// Inputs
		double dt;							//!< Time step [s]
		double Temperature_inlet_1;			//!< Temperature of the first fluid at the inlet of the heat exchanger [K]
		double Temperature_inlet_2;			//!< Temperature of the second fluid at the inlet of the heat exchanger [K]
		double mass_flow_1;					//!< Mass flow of the first fluid [kg/s]
		double mass_flow_2;					//!< Mass flow of the second fluid [kg/s]

		// Results
		double efficiency;					//!< Heat exchanger efficiency
		double exchanged_heat_1;			//!< Heat being exchanged from the second fluid to the first fluid [W]
		double exchanged_heat_2;			//!< Heat being exchanged from the first fluid to the second fluid [W]
		double heat_1;						//!< Heat being exchanged from the second fluid to the first fluid [J]
		double heat_2;						//!< Heat being exchanged from the first fluid to the second fluid [J]
		double Temperature_outlet_1;		//!< Temperature of the first fluid at the outlet of the heat exchanger [K]
		double Temperature_outlet_2;		//!< Temperature of the second fluid at the outlet of the heat exchanger [K]


		// Characteristics for efficiency
		int exchanger_type;					//!< 1:Simple co-current, 2:Simple counter-current, 3:Shell and tube, 4:Crossed flows (simple pass, fluids not mixed), 5:Crossed flows (simple pass, fluid 1 mixed, fluid 2 not mixed), 6:Crossed flows (simple pass, fluid 1 not mixed, fluid 2 mixed), 7:Crossed flows (simple pass, both fluids mixed)
		int n_shell_passes;					//!< Passes through shell

		// Efficiency table
		bool flag_efficiency_table;				//!< True if efficiency is given by a table as a function of NTU and Cr; False if efficiency is calculated from the heat exchanger characteristics
		std::vector <double> NTU_col;			//!< Column NTU for the efficiency table
		std::vector <double> Cr_col;			//!< Column Cr for the efficiency table
		std::vector <double> efficiency_col;	//!< Efficiency values for every combination of NTU and Cr in the table
		int n_rows_table;						//!< Rows in the table
	


	// METHODS

		/*!
		* \brief Looks up a value in a bidimensional table
		*
		* \author J. Salvador
		* \date 08/06/2016
		*
		* \param n_rows			Number of rows of x, y and z columns
		* \param x   			Value of parameter x
		* \param y   			Value of parameter y
		* \param *x_column		Values of parameter x in the table
		* \param *y_column		Values of parameter y in the table
		* \param *z_column		Values of queried variable z in the table
		*
		* \return			Queried value of the variable z 
		*/
		double lookup_table_2dim (int n_rows, double x, double y, std::vector<double> *x_column, std::vector<double> *y_column, std::vector<double> *z_column);


		/*!
		* \brief Heat capacity. Dummy fuction to be ignored when the model is integrated in VEMOD
		*
		* \author J. Salvador
		* \date 08/06/2016
		*
		* \param fluid_type		Fluid: air, water, oil, etc.
		* \param temperature	Fluid temperature [K]
		*
		* \return			Fixed value of heat capacity 4181.3 J/kg/K (water)
		*/
		// double cp(std::string fluid_type, double temperature);

};