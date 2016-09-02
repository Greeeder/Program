
#include <vector>


class HeatExchanger_model
{
public:

	struct stResults {	//!< Liner nodes data (geometry, etc.)
		int HeatExchanger_ID;	//!< ID of heat exchanger
		double Temperature_outlet_1;	//!< Temperature of the first fluid at the outlet of the heat exchanger
		double Temperature_outlet_2;	//!< Temperature of the second fluid at the outlet of the heat exchanger
		double exchanged_heat_1;	// Heat being exchanged from the second fluid to the first fluid
		double exchanged_heat_2;	// Heat being exchanged from the first fluid to the second fluid
	};

	std::vector<stResults> results;


	/*!
	* \brief Reads heat exchanger data and fills a stHeatExchanger structure for every heat exchanger
	*
	* \author J. Salvador
	* \date 10/06/2016
	*
	*/
	void create_HeatExchangers();

	/*!
	* \brief Refreshes the Heat Exchanger data
	
	
	*/






	void HeatExchanger_refresh(int HeatExchanger_Id, std::vector<stHeatExchanger> HeatExchangers);


private:

	struct stHeatExchanger {	//!< Liner nodes data (geometry, type, etc.)
		int ID_htx;	//!< ID of heat exchanger
		int ID_htx_HydroNet;	// ID inside the model of hydraulic circuits	[HydroNet]
		std::string name;	//!< Heat exchanger name
		std::string fluid_type_1;	//!< Air, water, oil, etc.
		std::string fluid_type_2;	//!< Air, water, oil, etc.
		double hydraulic_resistance;	//!< (Head loss / Flow^2) for a flow passing through the heat exchanger [m.c.f / (m3/s)]	[HydroNet]
		double head_loss;	//!< Fixed head loss for a flow passing through the heat exchanger [m.c.f]	[HydroNet]
		double UA;	//!< Overall heat transfer coefficient * heat transfer area
		int exchanger_type;	//!< 1:Simple co-current, 2:Simple counter-current, 3:Shell and tube, 4:Crossed flows (simple pass, fluids not mixed)
		int n_shell_passes;	//!< Passes through shell
		bool flag_efficiency_table;	//!< True if efficiency is given by a table as a function of NTU and Cr;
		//!< False if efficiency is calculated from the heat exchanger characteristics
		std::vector <double> NTU;	//!< column NTU for the efficiency table
		std::vector <double> Cr;	//!< column Cr for the efficiency table
		std::vector <double> efficiency;	//!< Efficiency values for every combination of NTU and Cr in the table
		int n_rows_table; //!< Rows in the table
	};

	
	std::vector<stHeatExchanger>HeatExchangers; // **** stHeatExchanger Heatexchager;


	/*!
	* \brief Looks up a value in a bidimensional table
	*
	* \author J. Salvador
	* \date 08/06/2016
	*
	* \param n_rows		Number of rows of x, y and z columns
	* \param x   		Value of parameter x
	* \param y   		Value of parameter y
	* \param x_column   Values of parameter x in the table
	* \param y_column  	Values of parameter y in the table
	* \param z_column  	Values of queried variable z in the table
	*
	* \return			Queried value of the variable z 
	*/
	double lookup_table_2dim(int n_rows, double x, double y, double x_column[], double y_column[], double z_column[]);

	/*!
	* \brief Obtains exchanged heat and outlet temperatures of the fluids in the heat exchanger using the efficiency method and writes them in Results
	*
	* \author J. Salvador
	* \date 09/06/2016
	*
	* \param HeatExchanger_ID		ID of heat exchanger
	* \param Temperature_inlet_1   	Temperature of the first fluid at the inlet of the heat exchanger
	* \param Temperature_inlet_2   	Temperature of the second fluid at the inlet of the heat exchanger
	* \param mass_flow_1			Mass flow of the first fluid through the heat exchanger
	* \param mass_flow_2			Mass flow of the second fluid through the heat exchanger
	*
	*/
	void calculate_heat_exchange(int HeatExchanger_ID, double Temperature_inlet_1, double Temperature_inlet_2, double mass_flow_1, double mass_flow_2);

};
