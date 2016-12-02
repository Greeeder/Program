function   [temperature_thermostat,temperature_after_engine,temperature_after_radiator] ...
=register( objects, obj_inlet_pos, obj_outlet_pos,branch_temp_pos, branch_temperature, ...
temperature_thermostat, temperature_after_engine, temperature_after_radiator, count_exec)


%[pos, branch] = HydroNet_GetPosition(6, objects, obj_inlet_pos, obj_outlet_pos, 'end');
[pos, branch] = HydroNet_GetPosition(1, objects, obj_inlet_pos, obj_outlet_pos, 'end');
temperature_thermostat(count_exec) = HydroNet_GetPosTemperature(pos, branch_temp_pos{branch}, branch_temperature{branch}) ;

% %[pos, branch] = HydroNet_GetPosition(14, objects, obj_inlet_pos, obj_outlet_pos, 'end');
% [pos, branch] = HydroNet_GetPosition(4, objects, obj_inlet_pos, obj_outlet_pos, 'startoo');
% temperature_after_engine(count_exec) = HydroNet_GetObjTemperature( pos, branch_temp_pos{branch}, branch_temperature{branch} );
% 
% %[pos, branch] = HydroNet_GetPosition(1, objects, obj_inlet_pos, obj_outlet_pos, 'end');
% [pos, branch] = HydroNet_GetPosition(8, objects, obj_inlet_pos, obj_outlet_pos, 'startoo');
% temperature_after_radiator(count_exec) = HydroNet_GetObjTemperature( pos, branch_temp_pos{branch}, branch_temperature{branch} );


end