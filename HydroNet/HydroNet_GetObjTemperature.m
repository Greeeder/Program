function [ temperature ] = HydroNet_GetObjTemperature( object_index, branch_temp_pos, branch_temperature, obj_inlet_pos, obj_outlet_pos, fluid_type)
    % Returns temperature in the middle of the specified object
    % Inserts a volume for the object and obtains its mean temperature
    
    start_pos = obj_inlet_pos{object_index};
    end_pos = obj_outlet_pos{object_index};
    
    temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
    
    [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(branch_temp_pos, branch_temperature, ...
        start_pos, end_pos, temp_action, fluid_type,NaN,NaN, NaN);
    
    temperature = HydroNet_GetPosTemperature((start_pos + end_pos)/2, branch_temp_pos_aux, branch_temp_aux );

end
