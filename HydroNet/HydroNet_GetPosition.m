function [ position, branch ] = HydroNet_GetPosition( ind_obj, objects, obj_inlet_pos, obj_outlet_pos, location)
%location: 'start', 'middle' (default), 'end'
    
    branch = objects.branch{ind_obj};
   
    if strcmpi(location, 'start')
       position = obj_inlet_pos{ind_obj};
    elseif strcmpi(location, 'end')
        position = obj_outlet_pos{ind_obj};
    else % middle
        position = (obj_inlet_pos{ind_obj} + obj_outlet_pos{ind_obj}) / 2;
    end

end

