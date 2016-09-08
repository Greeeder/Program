function [ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, ...
    heat_exch_fix, heat_exch_Tout, tank ] = HydroNet_ReadObj( circuit_file )
% Reads circuit from configuration file and prepares the data in arrays.
%   Input is the name of the file that contains the information of the
%   circuit. Outputs are an object's list with id, name, class and joints
%   (table  'objects') and one table for all objects of each class (except
%   nodes). In addition an analysis is performed to check if all objects
%   are connected and identifiers are right.


%circuit_file = 'Coolant_basic_nodes.txt';
circuit_text = fileread(circuit_file); % read file

circuit_obj = strsplit(circuit_text, '#'); % looks for objects and splits text
fluid_type = sscanf(circuit_obj{1}, '%s %*'); % read fisrt string
circuit_obj = circuit_obj(2:end); % remove first string (fluid type)
circuit_obj = cellstr(circuit_obj); % converts into cell array of strings



%% CREATES TABLE OF OBJECTS

n_obj = numel(circuit_obj); % number of objects
if n_obj == 0 % no objects
    handle = msgbox('No elements defined in circuit');
    uiwait(handle); % keep msgbox in the foreground
end

obj_id = zeros(n_obj, 1); % vector that contants id of objects
obj_name = cell(n_obj, 1); % cell because it contains strings
obj_class = cell(n_obj, 1); % cell because it contains strings
obj_volume = zeros(n_obj, 1); % vector that contants volume of objects
obj_adjac = cell(n_obj, 1); % cell because it contains vectors (adjacent objects)

for count = 1 : n_obj % read id, name, class inlet and outlet of objects
    
   obj_id(count) = sscanf(circuit_obj{count}, '%d %*'); % id: first integer at first line
   obj_name{count} = sscanf(circuit_obj{count}, '%*[^\n] %s %*'); % name: first string at 2nd line
   obj_class{count} = sscanf(circuit_obj{count}, '%*[^\n] %*[^\n] %s %*'); % class
       
	if strcmp(obj_class{count}, 'node') % is node
        n_aux = sscanf(circuit_obj{count}, '%*[^\n] %*[^\n] %*[^\n] %d %*'); % number of joints
        adjac = zeros(n_aux, 1); % vector of connections
        for count1 = 1 : n_aux % reads every connection and stores them in an array
            str_aux = '%*[^\n] %*[^\n] %*[^\n] %d %*'; % to read line 4
            for count2 = 1 : count1  % line breaks to add
                str_aux = strjoin({'%*[^\n]', str_aux}); % add one line break for every connection
            end
            adjac(count1) = sscanf(circuit_obj{count}, str_aux); % stores connection
        end
        obj_adjac{count} = adjac; % contains vector with all connections
       
	else % one inlet and one outlet (arbitratry)
        adjac = zeros(2, 1);
        adjac(1) = sscanf(circuit_obj{count}, '%*[^\n] %*[^\n] %*[^\n] %d %*'); % inlet
        adjac(2) = sscanf(circuit_obj{count}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %d %*'); % outlet
        obj_adjac{count} = adjac; % contains vector with inlet and outlet
	end
end

objects = table(obj_id, obj_name, obj_class, obj_volume, obj_adjac, 'VariableNames', {'id' 'name' 'class' 'volume' 'adjacent'});



%% CHECKS THAT OBJECT'S ID ARE NOT 0 AND ARE NOT DUPLICATED

if any(obj_id == 0) % there is some 0
    handle = msgbox('Identifier cannot be 0');
    display(objects);
    for count = find((obj_id == 0))
        fprintf('Id of %s is 0\n', objects.name{count});
    end
    uiwait(handle);
end

if ne(numel(unique(obj_id)), numel(obj_id))
    % vector of unique values must have the same size as original
    handle = msgbox('Some identifier is duplicated');
    display(objects);
    uiwait(handle);
end



%% CHECKS THAT CLASSES ARE RECOGNISED

for count = 1 : n_obj
    switch obj_class{count}
        case 'node'
        case 'pipe'
        case 'valve_fix'
        case 'valve_var'
        case 'thermostat'
        case 'pump_volum'
        case 'pump_turbo'
        case 'heat_exch_fix'
        case 'heat_exch_UA'
        case 'heat_exch_Tout'
        case 'tank'
        otherwisewise
            handle = msgbox('Unrecognised class');
            fprintf('Class of #%d %s is %s\n', objects.id(count), objects.name{count}, objects.class{count});
            uiwait(handle);
    end
end



%% CHECKS: CONNECTIVITY AMONG OBJECTS, NO SELF-CONNECTIONS AND NO REPETITIONS

for count = 1 : n_obj % objects loop
    
    id_aux = objects.id(count); % id of current object
    
    if any(objects.adjacent{count} == id_aux)
    	% the object is connected to itself
        handle = msgbox('Object self-connected');
        display(objects);
        fprintf('Object #%d is connected to itself', id_aux);
        uiwait(handle);
    end
    
    if ne(numel(objects.adjacent{count}), numel(unique(objects.adjacent{count})))
    	% there are repeated connections
        handle = msgbox('Connection between objects repeated');
        display(objects);
        fprintf('Object #%d has a repeated connection', id_aux);
        uiwait(handle);
    end
    
    for count1 = 1 : numel(objects.adjacent{count}) % connections loop
        id_aux1 = objects.adjacent{count}(count1); % id of adjacent object
        vec_aux = objects.adjacent(objects.id == id_aux1); % connections of the adjacent object
        if not(any(vec_aux{1} == id_aux))
            % the current object is not present among the connections of the adjacent object
            handle = msgbox('Wrong connection between objects');
            display(objects);
            fprintf('Object #%d is linked to #%d but #%d is not linked to #%d', ...
                    id_aux, id_aux1, id_aux1, id_aux);
            uiwait(handle);
        end
    end
end



%% CREATES TABLES FOR CLASSES

% PIPE
% Looks for objects of class 'pipe'
aux_pos = strcmp(obj_class, 'pipe'); % is pipe? (0,1)
aux_ind = find(aux_pos); % indexes (row of object)
pipe = zeros(sum(aux_pos), 4); % to store parameters
for count = 1 : sum(aux_pos) % only the pipes
    pipe(count, 1) = objects.id(aux_ind(count)); % identifier
    pipe(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % length
    pipe(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % diameter
    pipe(count, 4) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % friction coef.
    objects.volume(aux_ind(count)) = pipe(count, 2) * pi / 4 * pipe(count, 3)^2; % calculate volume = length * pi / 4 * diameter^2
end


% VALVE_FIX
aux_pos = strcmp(obj_class, 'valve_fix');
aux_ind = find(aux_pos);
valve_fix = zeros(sum(aux_pos), 2);
for count = 1 : sum(aux_pos)
    valve_fix(count, 1) = objects.id(aux_ind(count)); % identifier
    valve_fix(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % head loss
end


% VALVE_VAR
aux_pos = strcmp(obj_class, 'valve_var');
aux_ind = find(aux_pos);
valve_var = zeros(sum(aux_pos), 6);
for count = 1 : sum(aux_pos)
    valve_var(count, 1) = objects.id(aux_ind(count)); % identifier
    valve_var(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % opening
    valve_var(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #0 of curve
    valve_var(count, 4) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #1 of curve
    valve_var(count, 5) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #2 of curve
    valve_var(count, 6) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #3 of curve 
end


% THERMOSTAT
aux_pos = strcmp(obj_class, 'thermostat');
aux_ind = find(aux_pos);
thermostat = zeros(sum(aux_pos), 10);
for count = 1 : sum(aux_pos)
    thermostat(count, 1) = objects.id(aux_ind(count)); % identifier
    % Rest of parameters (coefficients) - use a loop because there are many coefficients (8)
    str_aux = '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'; % to read line 5 (if it were a float)
    for count1 = 2 : size(thermostat, 2) % line breaks to add
        str_aux = strjoin({'%*[^\n]', str_aux}); % add one line break for every new coefficient
        thermostat(count, count1) = sscanf(circuit_obj{aux_ind(count)}, str_aux); % stores coefficients
    end
end


% PUMP_VOLUM
aux_pos = strcmp(obj_class, 'pump_volum');
aux_ind = find(aux_pos);
pump_volum = zeros(sum(aux_pos), 10);
for count = 1 : sum(aux_pos)
    pump_volum(count, 1) = objects.id(aux_ind(count)); % identifier
    pump_volum(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %d %*'); % outlet object
    objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
    pump_volum(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % pump speed
    pump_volum(count, 4) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % max. head
    pump_volum(count, 5) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #0 of curve
    pump_volum(count, 6) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #1 of curve
    pump_volum(count, 7) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #2 of curve
    pump_volum(count, 8) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % coef. #3 of curve
    pump_volum(count, 9) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % pump head % initialization of instantaneous variable
    pump_volum(count, 10) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % pump power % initialization of instantaneous variable
end


% PUMP_TURBO
aux_pos = strcmp(obj_class, 'pump_turbo');
aux_ind = find(aux_pos);
pump_turbo = zeros(sum(aux_pos), 14);
for count = 1 : sum(aux_pos)
    pump_turbo(count, 1) = objects.id(aux_ind(count)); % identifier
    pump_turbo(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %d %*'); % outlet object
    objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
    pump_turbo(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % pump speed
    pump_turbo(count, 4) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % max head
    % Rest of parameters (coefficients) - use a loop because there are many coefficients (9)
    str_aux = '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'; % to read line 7
    for count1 = 2 : size(pump_turbo, 2) % line breaks to add
        str_aux = strjoin({'%*[^\n]', str_aux}); % add one line break for every new coefficient
        pump_turbo(count, count1) = sscanf(circuit_obj{aux_ind(count)}, str_aux); % stores coefficients
    end
    str_aux = strjoin({'%*[^\n]', str_aux}); % add one more line break
    pump_turbo(count, size(pump_turbo, 2)) = sscanf(circuit_obj{aux_ind(count)}, str_aux); % pump power % initialization of instantaneous variable
end


% HEAT_EXCH_FIX
aux_pos = strcmp(obj_class, 'heat_exch_fix');
aux_ind = find(aux_pos);
heat_exch_fix = zeros(sum(aux_pos), 4);
for count = 1 : sum(aux_pos)
    heat_exch_fix(count, 1) = objects.id(aux_ind(count)); % identifier
    objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
    heat_exch_fix(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % heat
    heat_exch_fix(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % hydr. resistance
    heat_exch_fix(count, 4) = -999; % inlet temperature % initialization of instantaneous variable
end


% HEAT_EXCH_UA % will be of type heat_exch_Tout
% aux_pos = strcmp(obj_class, 'heat_exch_UA');
% aux_ind = find(aux_pos);
% heat_exch_UA = zeros(sum(aux_pos), 4);
% for count = 1 : sum(aux_pos)
%     heat_exch_UA(count, 1) = objects.id(aux_ind(count)); % identifier
%     objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
%     heat_exch_UA(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % UA
%     heat_exch_UA(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % mCp_min
%     heat_exch_UA(count, 4) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % hydr. resistance
% end


% HEAT_EXCH_Tout
aux_pos = strcmp(obj_class, 'heat_exch_Tout');
aux_ind = find(aux_pos);
heat_exch_Tout = zeros(sum(aux_pos), 5);
for count = 1 : sum(aux_pos)
    heat_exch_Tout(count, 1) = objects.id(aux_ind(count)); % identifier
    heat_exch_Tout(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %d %*'); % outlet object
    objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
    heat_exch_Tout(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % hydr. resistance
    heat_exch_Tout(count, 4) = -999; % outlet temperature % initialization of instantaneous variable
    heat_exch_Tout(count, 5) = -999; % inlet temperature % initialization of instantaneous variable
end


% TANK
% aux_pos = strcmp(obj_class, 'tank');
aux_ind = find(aux_pos);
tank = zeros(sum(aux_pos), 3);
for count = 1 : sum(aux_pos)
    tank(count, 1) = objects.id(aux_ind(count)); % identifier
    tank(count, 2) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %d %*'); % outlet object (to be feed by the tank)
    objects.volume(aux_ind(count)) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*'); % volume
    tank(count, 3) = sscanf(circuit_obj{aux_ind(count)}, '%*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %*[^\n] %f %*');
        % head drop between inlet and outlet   
end



% Convert to tables or create empty tables

pipe = array2table(pipe, 'VariableNames',{'id' 'length' 'diameter' 'friction_coef'});
    pipe.Properties.VariableUnits{'length'} = 'm';
    pipe.Properties.VariableUnits{'diameter'} = 'm';
    pipe.Properties.VariableDescriptions{'diameter'} = 'Inner diameter';
    pipe.Properties.VariableUnits{'friction_coef'} = ' ';
    
    
valve_fix = array2table(valve_fix, 'VariableNames',{'id' 'head_loss'});
    valve_fix.Properties.VariableUnits{'head_loss'} = 'm.c.f.';
    
    
valve_var = array2table(valve_var, 'VariableNames',{'id' 'opening' 'coef0' 'coef1' 'coef2' 'coef3'});
    valve_var.Properties.VariableUnits{'opening'} = '[0 - 1]';
    valve_var.Properties.VariableDescriptions{'opening'} = '0 closed; 1 wide open. Instantaneous input';
    valve_var.Properties.VariableUnits{'coef0'} = 'm.c.f.';
    valve_var.Properties.VariableUnits{'coef1'} = 'm.c.f.';
    valve_var.Properties.VariableUnits{'coef2'} = 'm.c.f.';
    valve_var.Properties.VariableUnits{'coef3'} = 'm.c.f.';
    valve_var.Properties.VariableDescriptions{'coef0'} = 'Head loss: constant';
    valve_var.Properties.VariableDescriptions{'coef1'} = 'Head loss: coefficient that multiplies opening';
    valve_var.Properties.VariableDescriptions{'coef2'} = 'Head loss: coefficient that multiplies opening^2';
    valve_var.Properties.VariableDescriptions{'coef3'} = 'Head loss: coefficient that multiplies opening^3';
    
    
thermostat = array2table(thermostat, 'VariableNames',{'id' 'T_coef0' 'T_coef1' 'T_coef2' 'T_coef3' 'opening' 'h_coef0' 'h_coef1' 'h_coef2' 'h_coef3'});
    thermostat.Properties.VariableUnits{'T_coef0'} = ' - ';
    thermostat.Properties.VariableUnits{'T_coef1'} = 'K^-1';
    thermostat.Properties.VariableUnits{'T_coef2'} = 'K^-2';
    thermostat.Properties.VariableUnits{'T_coef3'} = 'K^-3';
    thermostat.Properties.VariableUnits{'opening'} = '[0 - 1]';
    thermostat.Properties.VariableDescriptions{'T_coef0'} = 'Opening: constant';
    thermostat.Properties.VariableDescriptions{'T_coef1'} = 'Opening: coefficient that multiplies temperature';
    thermostat.Properties.VariableDescriptions{'T_coef2'} = 'Opening: coefficient that multiplies temperature^2';
    thermostat.Properties.VariableDescriptions{'T_coef3'} = 'Opening: coefficient that multiplies temperature^3';
    thermostat.Properties.VariableDescriptions{'opening'} = '0 closed; 1 wide open';
    thermostat.Properties.VariableUnits{'h_coef0'} = 'm.c.f.';
    thermostat.Properties.VariableUnits{'h_coef1'} = 'm.c.f.';
    thermostat.Properties.VariableUnits{'h_coef2'} = 'm.c.f.';
    thermostat.Properties.VariableUnits{'h_coef3'} = 'm.c.f.';
    thermostat.Properties.VariableDescriptions{'h_coef0'} = 'Head loss: constant';
    thermostat.Properties.VariableDescriptions{'h_coef1'} = 'Head loss: coefficient that multiplies opening';
    thermostat.Properties.VariableDescriptions{'h_coef2'} = 'Head loss: coefficient that multiplies opening^2';
    thermostat.Properties.VariableDescriptions{'h_coef3'} = 'Head loss: coefficient that multiplies opening^3';
    
    

pump_volum = array2table(pump_volum, 'VariableNames',{'id' 'outlet' 'pump_speed' 'head_max' 'coef0' 'coef1' 'coef2' 'coef3' ' pump_head' 'pump_power'});
    pump_volum.Properties.VariableUnits{'pump_speed'} = 'rad/s';
    pump_volum.Properties.VariableUnits{'head_max'} = 'm.c.f.';
    pump_volum.Properties.VariableUnits{'coef0'} = 'm^3/s';
    pump_volum.Properties.VariableUnits{'coef1'} = 'm^3/rad';
    pump_volum.Properties.VariableUnits{'coef2'} = 'm^3/rad^2*s';
    pump_volum.Properties.VariableUnits{'coef3'} = 'm^3/rad^3*s^2';
    pump_volum.Properties.VariableUnits{'pump_head'} = 'm.c.f.';
    pump_volum.Properties.VariableUnits{'pump_power'} = 'W';
    pump_volum.Properties.VariableDescriptions{'outlet'} = 'Outlet object';
    pump_volum.Properties.VariableDescriptions{'pump_speed'} = 'Instantaneous input';
    pump_volum.Properties.VariableDescriptions{'head_max'} = 'Maximum head: if head loss is higher, flow is zero';
    pump_volum.Properties.VariableDescriptions{'coef0'} = 'Flow: constant';
    pump_volum.Properties.VariableDescriptions{'coef1'} = 'Flow: coefficient that multiplies speed';
    pump_volum.Properties.VariableDescriptions{'coef2'} = 'Flow: coefficient that multiplies speed^2';
    pump_volum.Properties.VariableDescriptions{'coef3'} = 'Flow: coefficient that multiplies speed^3';
    pump_volum.Properties.VariableDescriptions{'pump_head'} = 'Instantaneous result';
    pump_volum.Properties.VariableDescriptions{'pump_power'} = 'Instantaneous result';
    
    
pump_turbo = array2table(pump_turbo, 'VariableNames',{'id' 'outlet' 'pump_speed' 'head_max' 'coef_N0' 'coef_N1' 'coef_N2' 'Q_coef_N0' 'Q_coef_N1' 'Q_coef_N2' 'Q2_coef_N0' 'Q2_coef_N1' 'Q2_coef_N2' 'pump_power'});
    pump_turbo.Properties.VariableUnits{'pump_speed'} = 'rad/s';
    pump_turbo.Properties.VariableUnits{'head_max'} = 'm.c.f.';
    pump_turbo.Properties.VariableUnits{'coef_N0'} = 'm.c.f.';
    pump_turbo.Properties.VariableUnits{'coef_N1'} = 'm.c.f./rad*s';
    pump_turbo.Properties.VariableUnits{'coef_N2'} = 'm.c.f./rad^2*s^2';
    pump_turbo.Properties.VariableUnits{'Q_coef_N0'} = 'm.c.f./m^3*s';
    pump_turbo.Properties.VariableUnits{'Q_coef_N1'} = 'm.c.f./m^3/rad*s^2';
    pump_turbo.Properties.VariableUnits{'Q_coef_N2'} = 'm.c.f./m^3/rad^2*s^3';
    pump_turbo.Properties.VariableUnits{'Q2_coef_N0'} = 'm.c.f./m^6*s^2';
    pump_turbo.Properties.VariableUnits{'Q2_coef_N1'} = 'm.c.f./m^6/rad*s^3';
    pump_turbo.Properties.VariableUnits{'Q2_coef_N2'} = 'm.c.f./m^6/rad^2*s^4';
    pump_turbo.Properties.VariableDescriptions{'outlet'} = 'Outlet object';
    pump_turbo.Properties.VariableDescriptions{'pump_speed'} = 'Instantaneous input';
    pump_turbo.Properties.VariableDescriptions{'head_max'} = 'Maximum head: if head loss is higher, flow is zero';
    pump_turbo.Properties.VariableDescriptions{'coef_N0'} = 'Head loss: constant';
    pump_turbo.Properties.VariableDescriptions{'coef_N1'} = 'Head loss: coefficient that multiplies speed';
    pump_turbo.Properties.VariableDescriptions{'coef_N2'} = 'Head loss: coefficient that multiplies speed^2';
    pump_turbo.Properties.VariableDescriptions{'Q_coef_N0'} = 'Head loss: coefficient that multiplies flow';
    pump_turbo.Properties.VariableDescriptions{'Q_coef_N1'} = 'Head loss: coefficient that multiplies flow and speed';
    pump_turbo.Properties.VariableDescriptions{'Q_coef_N2'} = 'Head loss: coefficient that multiplies flow and speed^2';
    pump_turbo.Properties.VariableDescriptions{'Q2_coef_N0'} = 'Head loss: coefficient that multiplies flow^2';
    pump_turbo.Properties.VariableDescriptions{'Q2_coef_N1'} = 'Head loss: coefficient that multiplies flow^2 and speed';
    pump_turbo.Properties.VariableDescriptions{'Q2_coef_N2'} = 'Head loss: coefficient that multiplies flow^2 and speed^2';
    
    
heat_exch_fix = array2table(heat_exch_fix, 'VariableNames',{'id' 'heat' 'hydr_resist' 'inlet_temp'});
    heat_exch_fix.Properties.VariableUnits{'heat'} = 'W';
    heat_exch_fix.Properties.VariableUnits{'hydr_resist'} = 'm.c.f./m^6*s^2';
    heat_exch_fix.Properties.VariableUnits{'inlet_temp'} = 'K';
    heat_exch_fix.Properties.VariableDescriptions{'heat'} = 'Heat power: (+) energy into circuit; (-) energy out of circuit';
    heat_exch_fix.Properties.VariableDescriptions{'hydr_resist'} = 'Hydraulic resistance = head loss / flow^2';
    heat_exch_fix.Properties.VariableDescriptions{'inlet_temp'} = 'Temperature at the inlet of the heat exchanger. Instantaneous result';
    
    
% will be of type heat_exch_Tout
% heat_exch_UA = array2table(heat_exch_UA, 'VariableNames',{'id' 'UA' 'mCp_min' 'hydr_resist'});
%     heat_exch_UA.Properties.VariableUnits{'UA'} = 'W/K';
%     heat_exch_UA.Properties.VariableUnits{'mCp_min'} = 'W/K';
%     heat_exch_UA.Properties.VariableUnits{'hydr_resist'} = 'm.c.f./m^6*s^2';
%     heat_exch_UA.Properties.VariableDescriptions{'UA'} = 'Constant parameter UA of exchanger: heat / temperature difference';
%     heat_exch_UA.Properties.VariableDescriptions{'mCp_min'} = 'Minimum value from both fluids of the product of mass flow times heat capacity';
%     heat_exch_UA.Properties.VariableDescriptions{'hydr_resist'} = 'Hydraulic resistance = head loss / flow^2';
    
    
heat_exch_Tout = array2table(heat_exch_Tout, 'VariableNames',{'id' 'outlet' 'hydr_resist' 'T_out' 'inlet_temp'});
    heat_exch_Tout.Properties.VariableUnits{'hydr_resist'} = 'm.c.f./m^6*s^2';
    heat_exch_Tout.Properties.VariableUnits{'T_out'} = 'K';
    heat_exch_Tout.Properties.VariableUnits{'inlet_temp'} = 'K';
    heat_exch_Tout.Properties.VariableDescriptions{'outlet'} = 'Outlet object';
    heat_exch_Tout.Properties.VariableDescriptions{'hydr_resist'} = 'Hydraulic resistance = head loss / flow^2';
    heat_exch_Tout.Properties.VariableDescriptions{'T_out'} = 'Fluid temperature at the outlet. Instantaneous input'; % instantaneous variable
    heat_exch_Tout.Properties.VariableDescriptions{'inlet_temp'} = 'Temperature at the inlet of the heat exchanger. Instantaneous result'; % instantaneous variable
    
    
tank = array2table(tank, 'VariableNames',{'id' 'outlet' 'head_drop'});



% display(pipe);
% display(valve_fix);
% display(valve_var);
% display(thermostat);
% display(pump_volum);
% display(pump_turbo);
% display(heat_exch_fix);
% display(heat_exch_UA);
% display(heat_exch_Tout);
% display(tank);

