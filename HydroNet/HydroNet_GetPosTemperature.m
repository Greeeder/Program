function [ temperature ] = HydroNet_GetPosTemperature( obj_pos, branch_temp_pos, branch_temperature )
    % Returns temperature at the specified position
    
    if obj_pos == branch_temp_pos(1) % first position: first temperature
		temperature = branch_temperature(1);
    elseif or(obj_pos == branch_temp_pos(end),obj_pos>100)  % last position: last temperature
		temperature = branch_temperature(end);
    else % find closest position
        distance = -999;
        count = 1; 
        while distance < 0 % the object position is before the current vector position
            count = count + 1;
			distance = branch_temp_pos(count) - obj_pos;
        end
		temperature = branch_temperature(count - 1);
    end
    
end
