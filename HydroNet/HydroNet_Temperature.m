
function [ new_pos, new_temp ] = HydroNet_Temperature( fluid_type, objects, obj_inlet_pos, obj_outlet_pos, ...
    heat_exch_Tout, heat_exch_fix, nodes_ind, branches_ind, node_branches, n_nodes, n_branch, branch_volume, ...
    branch_temp_pos, branch_temperature, branch_htx_Tout, branch_htx_fix, flows, dt )
% Finds positions where temperature changes and temperatures in those ranges.
% Moves positions and temperatures according to flows and, simultaneously,
% enters heat from heat exchangers. Returns refreshed ranges and temperatures
% in every branch as well as temperatures for the inlet of heat exchangers.
% Merges ranges (volumes) and its temperatures if there are too many of them in a branch.


% Maximum divisions inside each branch; if there are more, some are merged
max_divisions = 20;



%% VOLUME THAT STAYS INSIDE THE SAME BRANCH AND OVERFLOWING VOLUME

new_pos = cell(n_branch, 1); % temporary cell to store positions in the current instant
new_temp = cell(n_branch, 1); % temporary cell to store temperatures in the current instant
overflow = zeros(n_branch, 1); % vector to store volumes that move out from its original branch
overflow(branch_volume > 0) = flows(branch_volume > 0) * dt;
overflow_temperature = zeros(n_branch, 1); %ones(n_branch, 1) * (-999999); % vector to store temperatures of overflowed volumes


% Branches without volume
for count_branch = transpose(find(branch_volume == 0))
    new_pos{count_branch} = [0, 100]; 
    new_temp{count_branch} = branch_temperature{count_branch};
end


% Branch loop excluding branches without volume
for count_branch = transpose(find(branch_volume > 0))
    
    delta_pos = flows(count_branch) * dt / branch_volume(count_branch) * 100;
        % \% of branch volume that flow moved
    
    % NO FLOW MOVEMENT: refreshes temperatures inside the branch
    if delta_pos == 0 % flow is zero, no movement
        
        new_pos{count_branch} = branch_temp_pos{count_branch};
        new_temp{count_branch} = branch_temperature{count_branch};
        
        % If there are heat exchangers, inserts volume for it and refreshes temperature
        % Heat exchanger of type T_out
        for count_obj_htx = branch_htx_Tout{count_branch}
            
            ind_aux = heat_exch_Tout.obj_index(count_obj_htx);
                % index of the heat exchanger in the table objects
            
            start_pos = obj_inlet_pos{ind_aux}; % object's inlet position
            end_pos = obj_outlet_pos{ind_aux}; % object's outlet position
            temp_action = 2; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
            value = heat_exch_Tout.T_out(count_obj_htx); % final temperature
            
            [new_pos{count_branch}, new_temp{count_branch}] = HydroNet_InsertVolume...
                (new_pos{count_branch}, new_temp{count_branch}, start_pos, end_pos, temp_action, fluid_type, value);
        end
        % Heat exchanger of type fixed heat
        for count_obj_htx = branch_htx_fix{count_branch}
            
            ind_aux = heat_exch_fix.obj_index(count_obj_htx);
                % index of the heat exchanger in the table objects
            
            start_pos = obj_inlet_pos{ind_aux}; % object's inlet position
            end_pos = obj_outlet_pos{ind_aux}; % object's outlet position
            temp_action = 1; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature
            value = heat_exch_fix.heat(count_obj_htx); % heat to add
            
            [new_pos{count_branch}, new_temp{count_branch}] = HydroNet_InsertVolume(new_pos{count_branch}, ...
                new_temp{count_branch}, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch));
        end

        
    %%%%%%%%%%%%%%%%%%
    % THERE IS FLOW MOVEMENT INSIDE THE BRANCH
    else 
        
        position_count = 2; % position counter
        
        % REFRESHES TEMPERATURES INSIDE THE BRANCH WITHOUT HEAT EXCHANGE
        while and((branch_temp_pos{count_branch}(position_count) + delta_pos) <= 100, ... % new position is inside the same branch
                (position_count) < numel(branch_temp_pos{count_branch})) % old position exists inside the vector
            % Positions are moved along the branch according to flowing flow
            
            position = branch_temp_pos{count_branch}(position_count);
            
            % Flow must be positive (same direction as branch)
            new_pos{count_branch} = [new_pos{count_branch}, (position + delta_pos)]; % add new position at the end of the vector
            new_temp{count_branch} = [new_temp{count_branch}, branch_temperature{count_branch}(position_count - 1)];
            
            position_count = position_count + 1;
        end
        % New position is in another branch OR there are no more positions in the branch
        % Closes branch
        new_pos{count_branch} = [0, new_pos{count_branch}, 100]; % add vector start and end
        new_temp{count_branch} = [new_temp{count_branch}, branch_temperature{count_branch}(position_count - 1)];
            % add new temperature at the end of the vector
        
        
        %%%%%%%%%%%%%%%%%%
        % VOLUME OVERFLOW: makes volume buffer without taking into account heat exchangers
        % First step to find the temperature of the buffer volume
        if (position_count) <= numel(branch_temp_pos{count_branch})
            % old position exists inside the vector, thus after the previous in-branch
            % while loop the new position has to be in another branch: OVERFLOW
            
            % Makes a new pair of vectors like branch_temp_pos and branch_temperature
            % for the volume that goes out of the branch
            overflow_temp_volpos_aux = 100; % initialization: branch end
                % position specified as \% of the branch volume
            overflow_temperature_aux = [];
            
            for position_count = position_count : numel(branch_temp_pos{count_branch})
                
                overflow_temp_volpos_aux = [overflow_temp_volpos_aux, branch_temp_pos{count_branch}(position_count) + delta_pos];
                    % position specified as \% of the branch volume (> 100%)
                overflow_temperature_aux = [overflow_temperature_aux, branch_temperature{count_branch}(position_count - 1)];
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%
        % HEAT EXCHANGERS
        % If there are heat exchangers, checks whether the volume out of them remains inside the branch or it
        % goes out, inserts volume for it, refreshes temperature inside the branch and stores overflowed volume.
        % There can only be heat exchangers if branch volume is > 0.
        % Second step to find the temperature of the buffer volume
        
        % Heat exchanger of type T_out
        for count_obj_htx = branch_htx_Tout{count_branch} % loop starting from the closest exchanger to the branch's end
            
            ind_aux = heat_exch_Tout.obj_index(count_obj_htx);
                % index of the heat exchanger in the table objects
            
            temp_action = 2; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature    
            value = heat_exch_Tout.T_out(count_obj_htx); % new temperature
            
            start_pos = obj_outlet_pos{ind_aux}; % object's outlet position is the start position
            end_pos = start_pos + delta_pos;  
            
            if end_pos > 100
                start_pos = 100;
                 
                [overflow_temp_volpos_aux, overflow_temperature_aux] = HydroNet_InsertVolume...
                    (overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value);
                
                % For the volume inside the branch
                start_pos = obj_outlet_pos{ind_aux}; % object's outlet position is the start position
                end_pos = 100;
            end
            
            [new_pos{count_branch}, new_temp{count_branch}] = HydroNet_InsertVolume...
                (new_pos{count_branch}, new_temp{count_branch}, start_pos, end_pos, temp_action, fluid_type, value);
        end
        
         % Heat exchanger of type fixed heat
        for count_obj_htx = branch_htx_fix{count_branch} % loop starting from the closest exchanger to the branch's end
            
            ind_aux = heat_exch_fix.obj_index(count_obj_htx);
                % index of the heat exchanger in the table objects
            
            temp_action = 1; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature            
            value = heat_exch_fix.heat(count_obj_htx);
            
            start_pos = obj_outlet_pos{ind_aux}; % object's outlet position is the start position
            end_pos = start_pos + delta_pos;
            
            if end_pos > 100
                
                value = heat_exch_fix.heat(count_obj_htx) * (end_pos - 100)/(end_pos - start_pos);
                    % heat to add averaged for the volume that goes out of the branch
                
                start_pos = 100;
                
                [overflow_temp_volpos_aux, overflow_temperature_aux] = HydroNet_InsertVolume(overflow_temp_volpos_aux, ...
                    overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch));
                
                % For the volume inside the branch
                start_pos = obj_outlet_pos{ind_aux}; % object's outlet position is the start position
                end_pos = 100;     
                value = heat_exch_fix.heat(count_obj_htx) * (100 - start_pos)/(end_pos - start_pos);
                    % heat to add averaged for the volume that stays inside the branch
            end
            
            [new_pos{count_branch}, new_temp{count_branch}] = HydroNet_InsertVolume(new_pos{count_branch}, ...
                new_temp{count_branch}, start_pos, end_pos, temp_action, fluid_type, value, branch_volume(count_branch));
        end
        
        
        %%%%%%%%%%%%%%%%%%
        % TEMPERATURE OF THE BUFFER VOLUME (OVERFLOW)
        % Final averaging of temperatures in the buffer volume
        
        start_pos = 100;
        end_pos = overflow_temp_volpos_aux(end);
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        
        [~, overflow_temperature(count_branch)] = HydroNet_InsertVolume...
            (overflow_temp_volpos_aux, overflow_temperature_aux, start_pos, end_pos, temp_action, fluid_type);
    end
end



%% SEPARES BRANCHES AT THE INLET AND OUTLET OF EVERY NODE
% Separates branches that are inlets to nodes and those that are outlets to nodes from node_branches
% All branches should start and end at a node!!

if any(flows > 0) % if there is no movement in the branch, there is no need for this
    
    node_inlet_branches = node_branches;
    node_outlet_branches = node_branches;
    %node_volume_out = zeros(n_nodes, 1);
        % stores the sum of all volumes leaving every node to check if outlet branches are overflowed
        % no used because, later, volume_share is more precise
    for node_count = 1 : n_nodes % node loop
        flag_inlet = false(numel(node_branches{node_count}), 1); % marks inlet branches to this node
        count = 0; % counter of branches in every node
        for count_branch = node_branches{node_count} % alternative: branch_count = 1 : numel(node_branches{node_count})
            count = count + 1;
            if branches_ind{count_branch}(end) == nodes_ind(node_count)
                % flow goes into the node from the branch: it is an inlet branch to the node
                flag_inlet(count) = true; % mark branch as inlet to node
    %         else % flow goes into the branch from the node: it is an outlet branch from the node
    %             node_volume_out(node_count) = node_volume_out(node_count) + branch_volume(branch_count);
            end
        end
        % Separates inlet and outlet branches to each node
        node_inlet_branches{node_count} = node_inlet_branches{node_count}(flag_inlet);
        node_outlet_branches{node_count} = node_outlet_branches{node_count}(not(flag_inlet));
    end
end


%% MANAGES OVERFLOWED VOLUMES
% All branches should start and end at a node!!
% The number of branches without volume must be minimized!!

while any(overflow ~= 0) % while there are overflows

    % Makes an enthalpy balance of all volumes mixing in every node
    node_volume_in = zeros(n_nodes, 1); % stores the sum of all volumes mixing in every node
    node_mass = zeros(n_nodes, 1); % temporary vector to store mass in every node
    node_temperature = zeros(n_nodes, 1); % stores temperatures in every node after mixing
    
    % Volume and mass in nodes
    for count_branch = transpose(find(overflow > 0)) % 1 : n_branch % branch loop
        % do not consider branches that do not have overflows (first iteration: do not considerer nodes without volume)
        % branch_count = 1 : n_branch would be valid too (overflow = 0 does not add volume or mass)
        end_node_ind = find(nodes_ind == branches_ind{count_branch}(end)); % index of the node at the end of the branch
        node_volume_in(end_node_ind) = node_volume_in(end_node_ind)  + overflow(count_branch);
        node_mass(end_node_ind) = node_mass(end_node_ind) + ...
            density(fluid_type, overflow_temperature(count_branch)) * overflow(count_branch);
    end
    
    % Temperature in nodes, provisional
    % Processes nodes that have overflow buffer volumes > 0   
    for count_node = transpose(find(node_mass > 0)) % to avoid dividing by zero
        for count_branch = node_inlet_branches{count_node}
            end_node_ind = find(nodes_ind == branches_ind{count_branch}(end)); % index of the node at the end of the branch
            node_temperature(end_node_ind) = node_temperature(end_node_ind) + ...
                overflow_temperature(count_branch) * density(fluid_type, overflow_temperature(count_branch)) * ...
                overflow(count_branch) / node_mass(end_node_ind);
                % T branch to node * branch mass to node / total mass to node
        end
    end
    
    % Control goes to buffer volumes in the nodes node_volume_in
    overflow(:) = 0;
    overflow_temperature(:) = 0;
    
    
    % Passes buffer volumes that are at the inlet node of branches with
    % volume = 0 to the end node
    
    sum_aux = 0; % sum of volumes at the start of branches without volume
    for count_branch = transpose(find((branch_volume == 0) .* (flows > 0)));
        % loop of branches without volume but with flow movement
        start_node_ind = find(nodes_ind == branches_ind{count_branch}(1));
        sum_aux = sum_aux + node_volume_in(start_node_ind);
    end 
    
    while sum_aux > 0
        % there are buffer volumes at the start of branches without volume; they must be past through
        
        % Priority: the node at the start of the branch has no inlet branches
        % without volume (they also would get passed a volume from behind)
        branch_ind_order_aux = find((branch_volume == 0) .* (flows > 0)); % branches without volume but with flow movement
        for count_branch = transpose(find((branch_volume == 0) .* (flows > 0)))
            start_node_ind = find(nodes_ind == branches_ind{count_branch}(1)); % index of the node at the start of the branch
            if any(branch_volume(node_inlet_branches{start_node_ind}) == 0)
                % any of the inlet branches to the start node has no volume
                branch_ind_order_aux(branch_ind_order_aux == count_branch) = []; % remove the branch
                branch_ind_order_aux = [branch_ind_order_aux; count_branch]; % put the branch at the bottom of the list 
            end
        end
        
        % Moves volumes according to priority
        for count_branch = transpose(branch_ind_order_aux)
            start_node_ind = find(nodes_ind == branches_ind{count_branch}(1)); % index of the node at the start of the branch
            end_node_ind = find(nodes_ind == branches_ind{count_branch}(end)); % index of the node at the end of the branch

            node_volume_in(end_node_ind) = node_volume_in(end_node_ind) + node_volume_in(start_node_ind);
            temp_aux = (node_temperature(end_node_ind) * node_mass(end_node_ind) + ...
                 node_temperature(start_node_ind) * node_mass(start_node_ind)) / (node_mass(end_node_ind) + node_mass(start_node_ind)); 
                % balance between volume already present and incoming volume
            node_temperature(end_node_ind) = temp_aux;
            node_mass(end_node_ind) = node_mass(end_node_ind) + node_mass(start_node_ind);
            node_volume_in(start_node_ind) = 0;
            node_mass(start_node_ind) = 0;
            node_temperature(start_node_ind) = 0;
            new_temp{count_branch} = temp_aux; % assigns temperature to the branch 
        end
        
        sum_aux = 0; % sum of volumes at the start of branches without volume
        for count_branch = transpose(find((branch_volume == 0) .* (flows > 0)))
            % loop of branches without volume but with flow movement
            start_node_ind = find(nodes_ind == branches_ind{count_branch}(1)); % index of the node at the start of the branch
            sum_aux = sum_aux + node_volume_in(start_node_ind);
        end
    end
    
    
    % Distributes volume among branches
    for node_count = transpose(find(node_volume_in > 0)) % loop of nodes that have a buffer volume
        flow_sum = sum(flows(node_outlet_branches{node_count})); % sum of flows leaving the node
        volume_share = flows(node_outlet_branches{node_count}) / flow_sum .* node_volume_in(node_count);
            % volume that goes to every branch that is an outlet to the node (proportional to flow)
        for count_branch = 1 : numel(node_outlet_branches{node_count}) % loop of branches that are outlets to the node
            branch_aux = node_outlet_branches{node_count}(count_branch); % branch index
            if volume_share(count_branch) > 0 % processes only branches where volume actually goes in
                if volume_share(count_branch) <= branch_volume(branch_aux)
                    % volume in the branch is bigger than incoming volume
                    % also avoids dividing by zero

                    % Adds position and temperature at the beginning
                    new_pos{branch_aux} = [(volume_share(count_branch) / branch_volume(branch_aux) * 100), new_pos{branch_aux}];
                    new_temp{branch_aux} = [node_temperature(node_count), new_temp{branch_aux}];

                    % Removes positions lower than or equal to the new one
                    to_remove = false(numel(new_pos{branch_aux}), 1);
                    for pos_aux = 2 : numel(new_pos{branch_aux})
                        if new_pos{branch_aux}(pos_aux) <= new_pos{branch_aux}(1)
                            to_remove(pos_aux) = true;
                        end
                    end
                    new_pos{branch_aux} = new_pos{branch_aux}(not(to_remove));

                    % Does not remove temperature that is after the last position to remove
                    for count = numel(to_remove) : -1 : 1
                        if to_remove(count)
                            to_remove(count) = false;
                            break;
                        end
                    end
                    to_remove = to_remove(1:(end-1));
                    new_temp{branch_aux} = new_temp{branch_aux}(not(to_remove));

                    % Adds position 0 at the beginning
                    new_pos{branch_aux} = [0, new_pos{branch_aux}];

                else % volume in the branch is smaller than incoming volume

                    % Fills branch
                    new_pos{branch_aux} = [0, 100];
                    new_temp{branch_aux} = node_temperature(node_count);

                    % Stores the excess volume
                    overflow(branch_aux) = volume_share(count_branch) - branch_volume(branch_aux);
                    overflow_temperature(branch_aux) = node_temperature(node_count);
                end
            end
        end
    end
    node_volume_in(:) = 0; % just in case it is used after the while loop
        % control passes to overflow in branches
end




%% MERGES VOLUMES INSIDE A BRANCH
% If there are more volumes than the value specified in max_divisions

for count_branch = 1 : n_branch
    extra_div = (numel(new_pos{count_branch}) - 1) - max_divisions; % number of divisions over max_divisions
    while extra_div > 0
        new_pos_aux = new_pos{count_branch};
        
        % Finds smaller combination of adjacent volumes
        pos_aux = min(new_pos_aux(3:end) - new_pos_aux(1:end-2));
        
        % Obtains temperature by means of an enthalpy balance
        temp_aux = (new_temp(pos_aux) * density(fluid_type, new_temp(pos_aux)) * (new_pos(pos_aux + 1) - new_pos(pos_aux)) + ...
            new_temp(pos_aux + 1) * density(fluid_type, new_temp(pos_aux + 1)) * (new_pos(pos_aux + 2) - new_pos(pos_aux + 1))) / ...
            (density(fluid_type, new_temp(pos_aux)) * (new_pos_aux(pos_aux + 1) - new_pos_aux(pos_aux)) + ...
            density(fluid_type, new_temp(pos_aux + 1)) * (new_pos_aux(pos_aux + 2) - new_pos_aux(pos_aux + 1))); 
        
        % Merges volumes
        new_pos{count_branch}(new_pos_aux + 1) = []; % remove position in the middle
        new_temp{count_branch}(new_pos_aux + 1) = []; % remove second temperature
        new_temp{count_branch}(new_pos_aux) = temp_aux; % enter temperature from balance
        
        extra_div = extra_div - 1; % next division to merge
    end
end



%% CALCULATES INLET TEMPERATURE TO HEAT EXCHANGERS
% Asumption: the flow that will pass through the heat exchanger in the next
% instant is equal to the branch flow in the current instant.
% Exception: if flow = 0 in the current instant, inserts a volume for the
% heat exchanger and obtains its temperature.
% Backwards direction from the outlet of the heat exchanger.

% Heat exchangers of type Tout
for count_htx = 1 : size(heat_exch_Tout, 1)
    
    obj_aux = heat_exch_Tout.obj_index(count_htx); % index of the heat exchanger in objects table
    branch_aux = objects.branch{obj_aux}; % branch that contains the heat exchanger
    
    if flows(branch_aux) == 0 % no movement in the branch
        
        % Inserts a volume for the heat exchanger and obtains its temperature
        start_pos = obj_inlet_pos{obj_aux};
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);
        
        % Reads and stores the temperature in the middle of the heat exchanger
        obj_pos = (start_pos + end_pos) / 2;
        heat_exch_Tout.inlet_temp(count_htx) = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        
    elseif (flows(branch_aux) * dt) <= (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux))
        % there is movement in the branch and the affected volume is inside the branch
        
        % Inserts a volume for the affected area and obtains its temperature
        start_pos = obj_outlet_pos{obj_aux} - flows(branch_aux) * dt / branch_volume(branch_aux) * 100;
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);        
        
        % Reads and stores the temperature in the middle of the affected volume
        obj_pos = (start_pos + end_pos) / 2;
        heat_exch_Tout.inlet_temp(count_htx) = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        
    else % there is movement in the branch and the affected volume reaches the previous branches
        
        % First, obtains temperature of the volume between the start of the
        % branch and the heat exchanger outlet, inserting a volume
        start_pos = 0;
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);        
        
        % Reads and stores the temperature in the middle of the affected volume
        obj_pos = (start_pos + end_pos) / 2;
        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        heat_exch_Tout.inlet_temp(count_htx) = temp_aux;
        
        % Initializes enthalpy balance
        sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux));
        sum_temp_mass = temp_aux * sum_mass; % sum of temperature times mass
        
        % Stores volume that will come from previous branches
        overflow(branch_aux) = (flows(branch_aux) * dt) - (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux));
        
        % Loop to obtain temperatures of the incoming flows
        while any(overflow ~= 0) % while there are overflows
            for count_branch = transpose(find(overflow > 0)) % loop of branches with overflow
                node_aux = branches_ind{count_branch}(1); % node at the start of the branch
                flow_sum = sum(flows(node_inlet_branches{node_aux})); % sum of flows going into the branch
                volume_share = flows(node_inlet_branches{node_aux}) / flow_sum .* overflow(count_branch);
                    % volume that comes from every inlet branch (proportional to flow)
                for count_branch1 = 1 : numel(node_inlet_branches{node_aux}) % loop of inlet branches
                    branch_aux1 = node_inlet_branches{node_aux}(count_branch1); % branch index
                    if volume_share(count_branch1) < branch_volume(branch_aux1)
                        % volume in the branch is bigger than outcoming volume
                        % also avoids dividing by zero
                        
                        % Inserts a volume for the affected area and obtains its temperature
                        start_pos = (1 - volume_share(count_branch1) / branch_volume(branch_aux1)) * 100;
                        end_pos = 100;
                        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
                        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux1}, new_temp{branch_aux1}, ...
                            start_pos, end_pos, temp_action, fluid_type);  
                        
                        % Reads the temperature in the middle of the affected volume
                        obj_pos = (start_pos + end_pos) / 2;
                        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
                        
                        % Enters temperature and mass in the enthalpy balance
                        sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share(count_branch1);
                        sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share(count_branch1);

                    else % volume in the branch is smaller than outcoming volume

                        % Inserts a volume for the entire branch and obtains its temperature
                        start_pos = 0;
                        end_pos = 100;
                        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
                        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux1}, new_temp{branch_aux1}, ...
                            start_pos, end_pos, temp_action, fluid_type);  
                        
                        % Reads the temperature in the middle of the branch
                        obj_pos = (start_pos + end_pos) / 2;
                        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
                        
                        % Enters temperature and mass in the enthalpy balance
                        sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share(count_branch1);
                        sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share(count_branch1);

                        % Stores the excess volume
                        overflow(branch_aux1) = volume_share(count_branch1) - branch_volume(branch_aux1);
                    end
                end
            end
        end
        overflow(:) = 0; % leaves everything as found
        
        % Stores resulting inlet temperature
        heat_exch_Tout.inlet_temp(count_htx) = sum_temp_mass / sum_mass;
    end
end


% Heat exhangers of type fixed heat
for count_htx = 1 : size(heat_exch_fix, 1)
    
    obj_aux = heat_exch_fix.obj_index(count_htx); % index of the heat exchanger in objects table
    branch_aux = objects.branch{obj_aux}; % branch that contains the heat exchanger
    
    if flows(branch_aux) == 0 % no movement in the branch
        
        % Inserts a volume for the heat exchanger and obtains its temperature
        start_pos = obj_inlet_pos{obj_aux};
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);
        
        % Reads and stores the temperature in the middle of the heat exchanger
        obj_pos = (start_pos + end_pos) / 2;
        heat_exch_fix.inlet_temp(count_htx) = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        
    elseif (flows(branch_aux) * dt) <= (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux))
        % there is movement in the branch and the affected volume is inside the branch
        
        % Inserts a volume for the affected area and obtains its temperature
        start_pos = obj_outlet_pos{obj_aux} - flows(branch_aux) * dt / branch_volume(branch_aux) * 100;
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);        
        
        % Reads and stores the temperature in the middle of the affected volume
        obj_pos = (start_pos + end_pos) / 2;
        heat_exch_fix.inlet_temp(count_htx) = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        
    else % there is movement in the branch and the affected volume reaches the previous branches
        
        % First, obtains temperature of the volume between the start of the
        % branch and the heat exchanger outlet, inserting a volume
        start_pos = 0;
        end_pos = obj_outlet_pos{obj_aux};
        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux}, new_temp{branch_aux}, ...
            start_pos, end_pos, temp_action, fluid_type);        
        
        % Reads and stores the temperature in the middle of the affected volume
        obj_pos = (start_pos + end_pos) / 2;
        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
        heat_exch_fix.inlet_temp(count_htx) = temp_aux;
        
        % Initializes enthalpy balance
        sum_mass = density(fluid_type, temp_aux) * (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux));
        sum_temp_mass = temp_aux * sum_mass; % sum of temperature times mass
        
        % Stores volume that will come from previous branches
        overflow(branch_aux) = (flows(branch_aux) * dt) - (obj_outlet_pos{obj_aux} / 100 * branch_volume(branch_aux));
        
        % Loop to obtain temperatures of the incoming flows
        while any(overflow ~= 0) % while there are overflows
            for count_branch = transpose(find(overflow > 0)) % loop of branches with overflow
                node_aux = branches_ind{count_branch}(1); % node at the start of the branch
                flow_sum = sum(flows(node_inlet_branches{node_aux})); % sum of flows going into the branch
                volume_share = flows(node_inlet_branches{node_aux}) / flow_sum .* overflow(count_branch);
                    % volume that comes from every inlet branch (proportional to flow)
                for count_branch1 = 1 : numel(node_inlet_branches{node_aux}) % loop of inlet branches
                    branch_aux1 = node_inlet_branches{node_aux}(count_branch1); % branch index
                    if volume_share(count_branch1) < branch_volume(branch_aux1)
                        % volume in the branch is bigger than outcoming volume
                        % also avoids dividing by zero
                        
                        % Inserts a volume for the affected area and obtains its temperature
                        start_pos = (1 - volume_share(count_branch1) / branch_volume(branch_aux1)) * 100;
                        end_pos = 100;
                        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
                        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux1}, new_temp{branch_aux1}, ...
                            start_pos, end_pos, temp_action, fluid_type);  
                        
                        % Reads the temperature in the middle of the affected volume
                        obj_pos = (start_pos + end_pos) / 2;
                        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
                        
                        % Enters temperature and mass in the enthalpy balance
                        sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share(count_branch1);
                        sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share(count_branch1);

                    else % volume in the branch is smaller than outcoming volume

                        % Inserts a volume for the entire branch and obtains its temperature
                        start_pos = 0;
                        end_pos = 100;
                        temp_action = 0; % 0: average temperature; 1: calculate temperature with heat; 2: impose temperature 
                        [branch_temp_pos_aux, branch_temp_aux] = HydroNet_InsertVolume(new_pos{branch_aux1}, new_temp{branch_aux1}, ...
                            start_pos, end_pos, temp_action, fluid_type);  
                        
                        % Reads the temperature in the middle of the branch
                        obj_pos = (start_pos + end_pos) / 2;
                        temp_aux = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos_aux, branch_temp_aux);
                        
                        % Enters temperature and mass in the enthalpy balance
                        sum_mass = sum_mass + density(fluid_type, temp_aux) * volume_share(count_branch1);
                        sum_temp_mass = sum_temp_mass + temp_aux * density(fluid_type, temp_aux) * volume_share(count_branch1);

                        % Stores the excess volume
                        overflow(branch_aux1) = volume_share(count_branch1) - branch_volume(branch_aux1);
                    end
                end
            end
        end
        overflow(:) = 0; % leaves everything as found
        
        % Stores resulting inlet temperature
        heat_exch_fix.inlet_temp(count_htx) = sum_temp_mass / sum_mass;
    end
end


end
