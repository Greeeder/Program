function [ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, ...
    heat_exch_fix, tank, nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, ...
    node_branches, n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, ...
    branch_htx_fix] = HydroNet_Create( circuit_file )
    % Creates a model for the specified circuit
     

%% ANALYSIS OF CIRCUIT OBJECTS
% Reads circuit configuration and specifications of objects and stores data in tables

[ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, ...
    heat_exch_fix, heat_exch_Tout, tank ] = HydroNet_ReadObj(circuit_file);

n_obj = size(objects, 1);



%% NODES
% List of all nodes

nodes_ind = find(strcmp(objects.class, 'node'));
nodes_id = objects.id(nodes_ind);
n_nodes = size(nodes_ind, 1);

% Using a loop
% nodes_ind = zeros(n_obj, 1); % stores position in object list
% nodes_id = zeros(n_obj, 1); % stores id
%     % maximum size: all objects are nodes
%     
% node_count = 0;
% for count = 1 : n_obj % objects loop
%     if strcmp(objects.class{count}, 'node')
%         node_count = node_count + 1;
%         nodes_ind(node_count) = count;
%         nodes_id(node_count) = objects.id(count);
%     end
% end
% nodes_ind = nodes_ind(1 : node_count); % adjust size
% nodes_id = nodes_id(1 : node_count);



%% TANKS
% List of all tanks

tanks_ind = zeros(n_obj, 1); % stores tank ids
    % maximum size: all objects are tanks
    
tank_count = 0;
for count = 1 : n_obj % objects loop
    if strcmp(objects.class{count}, 'tank')
        tank_count = tank_count + 1;
        tanks_ind(tank_count) = count;
    end
end
tanks_ind = tanks_ind(1 : tank_count); % adjust size

n_tanks = length(tanks_ind);



%% ADD OBJECT'S POSITION INTO CERTAIN TABLES
% Index of the object in the object's table

valve_var.obj_index = find(strcmp(objects.class, 'valve_var'));
thermostat.obj_index = find(strcmp(objects.class, 'thermostat'));
pump_volum.obj_index = find(strcmp(objects.class, 'pump_volum'));
pump_turbo.obj_index = find(strcmp(objects.class, 'pump_turbo'));
heat_exch_fix.obj_index = find(strcmp(objects.class, 'heat_exch_fix'));
heat_exch_Tout.obj_index = find(strcmp(objects.class, 'heat_exch_Tout'));
tank.obj_index = find(strcmp(objects.class, 'tank'));



%% CIRCUIT COMPONENTS
% Obtains nodes, branches and meshes in the circuit.
% Determines branches connected to nodes and branches that are part of a mesh.
% All valves and thermostats are assumed open.

[ objects, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, n_mesh, n_branch ] = HydroNet_Components...
    ( objects, n_obj, nodes_ind, pump_volum, pump_turbo, heat_exch_Tout, tank );



%% BRANCH VOLUME

branch_volume = zeros(n_branch, 1);
for count = 1 : n_branch % branch loop
    sum_aux = 0;
    for count1 = branches_ind{count} % loop of objects in the branch
        sum_aux = sum_aux + objects.volume(count1);
    end
    branch_volume(count) = sum_aux;
end



%% INLET AND OUTLET OF OBJECTS
% Stores volumetric position of the inlet and outlet of objects, as a \% of the total branch volume 

% Inlet positions
obj_inlet_pos_cell = cell(n_branch, 1); % branch / objects format
obj_inlet_pos = cell(n_obj, 1); % "vector" of objects format
for count_branch = 1 : n_branch % branch loop
    inlet_aux = zeros(1, numel(branches_ind{count_branch})); % auxiliary vector to store inlets in branch
    for count_obj = branches_ind{count_branch} % initialization loop
        obj_inlet_pos{count_obj} = [obj_inlet_pos{count_obj} 0]; % all objects of the branch have the inlet at 0%
    end
    if branch_volume(count_branch) ~= 0 % to avoid dividing by zero if branch volume is 0
        for count_obj = 2 : numel(branches_ind{count_branch}) % loop of objects in the branch
            inlet_aux(count_obj) = inlet_aux(count_obj - 1) + ...
                objects.volume(branches_ind{count_branch}(count_obj - 1)) / branch_volume(count_branch) * 100;
                % \% of the previous object inlet + volume previous object / branch volume * 100
            obj_inlet_pos{branches_ind{count_branch}(count_obj)} = ...
                [obj_inlet_pos{branches_ind{count_branch}(count_obj)}(1:(end-1)) inlet_aux(count_obj)];
        end
    end
    obj_inlet_pos_cell{count_branch} = inlet_aux; % stores inlets of the current branch
end


% Outlet positions
obj_outlet_pos = cell(n_obj, 1);
for count = 1 : n_obj
	if (objects.volume(count) == 0)
    	obj_outlet_pos{count} = obj_inlet_pos{count};
	else % volume > 0
    	obj_outlet_pos{count} = obj_inlet_pos{count} + objects.volume(count) / branch_volume(objects.branch{count}) * 100;
    end
end



%% INDEXES OF HEAT EXCHANGERS IN EVERY BRANCH

branch_htx_fix = cell(n_branch, 1);
branch_htx_Tout = cell(n_branch, 1);

for count_branch = 1 : n_branch
	branch_htx_fix{count_branch} = transpose(find(cell2mat(objects.branch(heat_exch_fix.obj_index)) == count_branch));
        % index in the object's table of heat exchangers of type fix heat that are present in the branch
    % Note: a heat exchanger cannot be in two branches but a branch can contain several heat exchangers    
    if isempty(branch_htx_fix{count_branch})
        branch_htx_fix{count_branch} = []; % to avoid problems with [0x1 double]
    end
    
    branch_htx_Tout{count_branch} = transpose(find(cell2mat(objects.branch(heat_exch_Tout.obj_index)) == count_branch));
        % index in the object's table of heat exchangers of type T_out that are present in the branch
    if isempty(branch_htx_Tout{count_branch})
        branch_htx_Tout{count_branch} = []; % to avoid problems with [0x1 double]
    end
    
    % Sort according to closeness to the end of the branch 
    [~, ind_aux] = sort(obj_inlet_pos_cell{count_branch}(branch_htx_fix{count_branch}),'descend');
        % sorted indexes of heat exchangers in obj_inlet_pos
    branch_htx_fix{count_branch} = branch_htx_fix{count_branch}(ind_aux);
        % reorder heat exchangers
    % Repeat for branch_htx_Tout
    [~, ind_aux] = sort(obj_inlet_pos_cell{count_branch}(branch_htx_Tout{count_branch}),'descend');
    branch_htx_Tout{count_branch} = branch_htx_Tout{count_branch}(ind_aux);
end


end

