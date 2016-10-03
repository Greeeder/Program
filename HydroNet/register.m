function   [temperature_after_pump,temperature_after_engine,temperature_after_radiator,time] ...
=register( branch_temp_pos, branch_temperature,dt,time,flag_first_exec,temperature_after_pump,...
temperature_after_engine,temperature_after_radiator)
if flag_first_exec 
temperature_after_pump(1) = HydroNet_GetObjTemperature(0, branch_temp_pos{5}, branch_temperature{5}) ;

temperature_after_engine(1) = HydroNet_GetObjTemperature( 0, branch_temp_pos{3}, branch_temperature{3} );

temperature_after_radiator(1) = HydroNet_GetObjTemperature( 0, branch_temp_pos{1}, branch_temperature{1} );

time(1) = 0;
else
temperature_after_pump(end +1) = HydroNet_GetObjTemperature( 0, branch_temp_pos{5}, branch_temperature{5}) ;

temperature_after_engine(end +1) = HydroNet_GetObjTemperature(0, branch_temp_pos{3}, branch_temperature{3} );

temperature_after_radiator(end +1) = HydroNet_GetObjTemperature( 0, branch_temp_pos{1}, branch_temperature{1} );

time(end+1)= time(end)+dt;
end

end