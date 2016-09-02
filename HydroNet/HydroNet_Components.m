function [ objects, branches_ind, branches_id, branch_cycle, ...
    mesh_branches, node_branches, n_mesh, branch_count ] = HydroNet_Components...
    ( objects, n_obj, nodes_ind, pump_volum, pump_turbo, heat_exch_Tout, tank ) 
% Obtains all information about a network and its components - meshes,
% nodes and branches.

% Stores all nodes in the network (index and id)
% Finds all paths using a Depth-first search algorithm.
% Checks whether paths are closed cycles (mesh).
% Stores all objects in each mesh (index and id)
% Finds all branches (paths without bifurcations crossed by a single flow)
% Stores all objects in each branch (index and id)
% Checks whether branches form part of closed cycles (meshes).
% Finds branches that form part of each mesh
% Finds branches connected to each node

% Inputs are the table of objects, tables for classes valve and thermostat
% (for checking opening) and temperature at the inlet of thermostats.

% Outputs are:
%   Nodes in the network (object indices)
%   Nodes in the network (object id)
%   Objects in every mesh (object indices)
%   Objects in every mesh (object id)
%   Objects in every branch (object indices)
%   Objects in every branch (object id)
%   Logical vector that marks whether branches form part of closed cycles
%   Branches that form part of each mesh
% 	Branches connected to each node
%   Number of meshes
%   Number of branches


% Initialization for testing
% circuit_file =  'Coolant_basic_nodes.txt';
% %circuit_file = 'Test_nodes.txt';
% [ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, ...
%     heat_exch_fix, heat_exch_Tout, tank ] = HydroNet_ReadObj(circuit_file);
% temperature = 0;



%% EDGES
% Looks for all connections among objects (indirected).

edge_count = 0; % counter of edges
edges = false(n_obj); % matrix that contains false if there is no connection between
    % objects and true if they are connected (rows and columns are the indices of objects)

for count = 1 : n_obj - 1 % objects loop (no need to check last object)
    % counter is equal to index
    for count1 = 1 : numel(objects.adjacent{count}) % connections loop
        id_aux = objects.adjacent{count}(count1); % id of adjacent object
        ind_aux = find(objects.id == id_aux);
        if ind_aux > count % the adjacent object has not been visited
            % the edge has not been stored
            edge_count = edge_count + 1; % next edge
            edges(count, ind_aux) = true; % stores edge
        end
    end
end % edges is a triangular matrix

edges = or(edges, transpose(edges)); % full symmetrical matrix (false OR true = true)


% Alternative: Store list in a cell
% edges = cell((obj_num * (obj_num - 1))/2, 1); % stores all edges
    % maximum size: all objects are nodes and are connected with each other   
% obj_id = zeros(obj_num, 1);
% obj_count = 0; % counter of visited objects
% for count = 1 : obj_num - 1 % objects loop (no need to check last object)
%     id_aux = objects.id(count); % id of current object
%     obj_count = obj_count + 1; % next visited object
%     obj_id(obj_count) = id_aux; % store object as visited
%     for count1 = 1 : numel(objects.adjacent{count}) % connections loop
%         id_aux1 = objects.adjacent{count}(count1); % id of adjacent object
%         if not(any((id_aux1 == obj_id)) % the adjacent object is not in the list of visited objects
%             % the edge has not been stored
%             edge_count = edge_count + 1; % next edge
%             edges{edge_count} = [id_aux id_aux1]; % stores edge
%         end
%     end
% end
% edges = edges(1 : edge_count); % adjust size



%% PATHS

% Selects starting object: first node of the list
if isempty(nodes_ind) % no nodes found
    start_obj = 1; % index of first object in the list
    paths_ind = cell(1); % stores all paths (indices) % only one path
else
    start_obj = nodes_ind(1); % index of first node
    paths_ind = cell(edge_count, 1); % stores all paths (indices)
        % maximum size: same as number of edges
end


% Finds all paths using a Depth-first search algorithm.

path_count = 0; % counter of paths
mat_aux = edges; % auxiliary matrix to mark visited objects
path = zeros(1, n_obj + 1); % stores indices of objects in a path
    % +1 because in a closed cycle the starting object is also the end object
path(1) = start_obj; % initialize path
obj_count = 1; % counter of objects along the path
row_aux = start_obj; % initialize row (starting object)
column = 0;

while not(and(row_aux == start_obj, column == n_obj))
    % while there are unexplored paths from the starting object (did not reach the end of the starting row)
    
    column = column + 1; % next column
    
    if mat_aux(row_aux, column) % there is a connection
        
        mat_aux(row_aux, column) = false; % mark edge as explored
        prev_row = row_aux;
        row_aux = column; % move to connection's row
        flag_path = false; % path not finished
        obj_count = obj_count + 1; % next object along the path
        path(obj_count) = row_aux; % store index of object in the path
        mat_aux(row_aux, prev_row) = false; % going back the same way is not allowed

        if not(any(mat_aux(row_aux, :))) % there are no more ways to go
            flag_path = true; % end path
        elseif any(path(1 : obj_count-2) == row_aux) % the new object was already in the path - the cycle is closed
            flag_path = true; % end path
%         elseif strcmp(objects.class(row_aux), 'valve_var') % is a valve and may be closed
%             if  valve_var.opening(valve_var.id == objects.id(row_aux)) == 0 % it is closed
%                 flag_path = true; % end path
%             end
%         elseif strcmp(objects.class(row_aux), 'thermostat') % is a thermostat and may be closed
%             pos_aux = find(thermostat.id == objects.id(row_aux)); % position in thermostat table
%             if thermostat.T_coef0(pos_aux) + thermostat.T_coef1(pos_aux) * temperature ...
%                     + thermostat.T_coef2(pos_aux) * temperature^2 == 0 % it is closed
%                 flag_path = true; % end path
%             end
        end

        if flag_path % path has ended
            path = path(1 : obj_count); % adjust size
            path_count = path_count + 1; % next path
            paths_ind{path_count} = path; % save path
            path = path(1 : end - 1); % remove last object from path - step back
            obj_count = obj_count - 1; % previous object along the path
            if row_aux == start_obj % the starting row never gets fully re-initialized
                mat_aux(row_aux, prev_row) = true; % just allow the same way back for future paths
            else
                mat_aux(row_aux, :) = edges(row_aux, :); % re-initialize connections of object
            end
            row_aux = path(end); % return to upper object
        else % path goes further
            column = 0; % stay in the new row and start from the beginning
        end
    end
    
    while and(column >= n_obj, ne(row_aux, start_obj)) % loop to step back in a path
        % end of the row and it is not the starting object % skip fully explored objects
        mat_aux(row_aux, :) = edges(row_aux, :); % re-initialize connections of object
        column = path(end); % continue the upper object where it was left
        path = path(1 : end - 1); % remove last object from path - step back
        obj_count = obj_count - 1; % previous object along the path
        row_aux = path(end); % return to upper object
        if numel(path) > 1
            mat_aux(row_aux, path(end-1)) = false; % the previous object along the path is not allowed
        end
    end
end
n_paths = path_count; % number of paths
paths_ind = paths_ind(1 : n_paths); % adjust size




% Find and isolate closed paths (cycles)
path_cycle = false(n_paths, 1);
    % vector that indicates whether the path forms a mesh (closed cycle)
for count_path = 1 : n_paths % path loop
    % Criterium to know if the paths are cycles: the last object and one of the other objects in the path are the same
    % Two cases: pure cycle (first object and last object are the same)
    %            branched cycle (first objects are not a part of the cycle)
    if paths_ind{count_path}(1) == paths_ind{count_path}(end) % pure cycle
        path_cycle(count_path) = true;
    else %any(paths_ind{count}(2 : end - 2) == paths_ind{count}(end)) % branched cycle
        for count_obj = 2 : numel(paths_ind{count_path}) - 3 % minus first, last and next-to-last
            if paths_ind{count_path}(count_obj) == paths_ind{count_path}(end) % cycle starts here
                path_cycle(count_path) = true;
                % Create new path for the branch that is out of the cycle
                paths_ind{size(paths_ind, 1) + 1} = paths_ind{count_path}(1 : count_obj);
                % Convert branched cycle in a pure cycle
                paths_ind{count_path} = paths_ind{count_path}(count_obj : end);
                break % path processed; go to next
            end
        end
    end
end
n_paths = size(paths_ind, 1); % refresh number of paths
path_cycle = [path_cycle; false(n_paths - numel(path_cycle), 1)];


% Removes repeated paths and cycles (all cycles are duplicated once - two directions)
to_remove = false(n_paths, 1); % stores paths to be removed
for count = 1 : n_paths % path loop
    if not(to_remove(count)) % the path is not marked to be deleted
        for count1 = count + 1 : n_paths % compare with next paths
            if numel(paths_ind{count1}) == numel(paths_ind{count})
                % to compare for equality, vectors must have the same length
%                 if or(paths_ind{count1} == paths_ind{count}, paths_ind{count1} == flip(paths_ind{count}))
%                     % paths are equal OR path is equal to the inverse path
%                     to_remove(count1) = true; % mark to be removed
%                 end
                % all objects of one branch are present in the other one
                % (even if start/end node is not the same)
                vec_aux = false(numel(paths_ind{count}), 1);
                for count_obj = 1 :  numel(paths_ind{count})
                    if any(paths_ind{count}(count_obj) == paths_ind{count1})
                        vec_aux(count_obj) = true;
                    end
                end
                if all(vec_aux) % all objects are present
                    to_remove(count1) = true; % mark to be removed
                end
            end
        end
    end
end
paths_ind = paths_ind(not(to_remove)); % remove repeated paths and cycles
path_cycle = path_cycle(not(to_remove));
n_paths = size(paths_ind, 1); % refresh number of paths


% % Convert paths to change indices by id's
% paths_id = cell(path_count, 1); % store all paths (id)
% for count = 1 : path_count % path loop
%     for count1 = 1 : numel(paths_ind{count}) % objects in a path
%         paths_id{count}(count1) = objects.id(paths_ind{count}(count1)); % assign id
%     end
% end



%% MESHES (cycles)
mesh_ind = paths_ind(path_cycle); % paths that are cycles (indices)
%mesh_id = paths_id(path_cycle); % (id)

% Convert meshes to change indices by id's
% mesh_id = cell(size(mesh_ind, 1), 1); % store all meshes (id)
% for count = 1 : size(mesh_id, 1) % mesh loop
%     for count1 = 1 : numel(mesh_ind{count}) % objects in a mesh
%         mesh_id{count}(count1) = objects.id(mesh_ind{count}(count1)); % assign id
%     end
% end



%% BRANCHES
% Paths without bifurcations crossed by a single flow.
% Branches that are not included in a mesh will be obtained too because they may be
% traversed by a flow if there is a boundary condition at one end (flow or head).


branches_ind = cell(edge_count * n_paths, 1); % stores unique branches of all paths (indices)
    % maximum size: number of edges * number of paths
branch_cycle = false(edge_count * n_paths, 1); % stores whether the branch is inside a mesh
branch_mesh = cell(edge_count * n_paths, 1); % stores the mesh where the branch is contained (indices)
    % is a cell because it will store all meshes that contain that branch
    % if the branch is not inside a mesh, the path is not stored
    % it is used to make a cell later that stores branches in each mesh
    
    
% Find branches in paths by splitting them at nodes and tanks
branch_count = 0; % branch counter
for count = 1 : n_paths % path loop
    start_obj = 1; % starting index for the first branch
    for count1 = 2 : (numel(paths_ind{count}) - 1) % counter of objects in a path excluding first and last
        if or(strcmp(objects.class{paths_ind{count}(count1)}, 'node'), ... % is a node
                strcmp(objects.class{paths_ind{count}(count1)}, 'tank'))
            branch_count = branch_count + 1; % next branch
            branches_ind{branch_count} = paths_ind{count}(start_obj:count1); % save branch
            start_obj = count1; % last object of the present branch is the first object of the next one
            branch_cycle(branch_count) = path_cycle(count); % save whether the branch is inside a mesh
            if path_cycle(count) % path is a mesh
                branch_mesh{branch_count} = find(find(path_cycle) == count); % save origin mesh (index)
            end
        end
    end
    if numel(paths_ind{count}) == 2 % branches with only 2 objects
        count1 = 1; % needed because the loop above was avoided 
    end
    count1 = count1 + 1; % next object % count1 = end of path
    branch_count = branch_count + 1; % next branch
    branches_ind{branch_count} = paths_ind{count}(start_obj:count1); % save branch
    branch_cycle(branch_count) = path_cycle(count); % save whether the branch is inside a mesh
    if path_cycle(count) % path is a mesh
        branch_mesh{branch_count} = find(find(path_cycle) == count); % save origin mesh (index)
    end
end
branches_ind = branches_ind(1 : branch_count); % adjust size
branch_cycle = branch_cycle(1 : branch_count);
branch_mesh = branch_mesh(1 : branch_count);


% Removes repeated branches
to_remove = false(branch_count, 1); % stores branches to be removed
for count = 1 : branch_count % branch loop
    if not(to_remove(count)) % the branch is not marked to be deleted
        for count1 = count + 1 : branch_count % compare with next branches
            if numel(branches_ind{count1}) == numel(branches_ind{count})
                % to compare for equality, vectors must have the same length
                if or(branches_ind{count1} == branches_ind{count}, branches_ind{count1} == flip(branches_ind{count}))
                    % branches are equal OR branch is equal to the inverse branch
                    to_remove(count1) = true; % mark to be removed
                    if and(not(branch_cycle(count)), branch_cycle(count1))
                        % first branch is not in a cycle but second (and same) branch is
                        branch_cycle(count) = true; % branch is inside a cycle
                    end
                    branch_mesh{count} = [branch_mesh{count} branch_mesh{count1}]; % add mesh
                end
            end
        end
    end
end
branches_ind = branches_ind(not(to_remove)); % remove repeated branches
branch_cycle = branch_cycle(not(to_remove));
branch_mesh = branch_mesh(not(to_remove));
branch_count = size(branches_ind, 1); % refresh number of branches


% Finds and merges split branches.
% If there are no nodes, the starting object may be inside a mesh. In that case,
% it is in the middle of a branch but the algorithms have splitted it in two parts.
if isempty(nodes_ind) % no nodes found (first condition)
    for count = 1 : branch_count % branch loop
        if or(branches_ind{count}(1) == 1, branches_ind{count}(end) == 1)
            % starting object is in the branch (its index is 1)
            if branch_cycle(count) % is in a mesh (second condition)
                for count1 = count : branch_count % looks for the other side
                    if or(branches_ind{count1}(1) == 1, branches_ind{count1}(end) == 1)
                        % Join branches
                        if and(branches_ind{count}(end) == 1, branches_ind{count}(1) == 1)
                            % starting object is the last of first branch and the first of the second branch
                            branches_ind{count} = [branches_ind{count}, branches_ind{count1}(2:end)]; % join directly % do not repeat starting object
                        elseif and(branches_ind{count}(1) == 1, branches_ind{count}(end) == 1)
                            % starting object is the last of first branch and the last of the second branch
                            branches_ind{count} = [branches_ind{count}, flip(branches_ind{count1}(1 : end-1))]; % flip second branch
                        elseif and(branches_ind{count}(1) == 1, branches_ind{count}(end) == 1)
                            % starting object is the first of first branch and the first of the second branch
                            branches_ind{count} = [flip(branches_ind{count}(2:end)), branches_ind{count1}]; % flip first branch
                        elseif and(branches_ind{count}(1) == 1, branches_ind{count}(end) == 1)
                            % starting object is the first of first branch and the last of the second branch
                            branches_ind{count} = [branches_ind{count1}, branches_ind{count}(2:end)]; % join second branch and then first
                        end
                        branches_ind = branches_ind(1 : count1 - 1, count1 + 1 : end); % remove second branch
                        branch_cycle = branch_cycle(1 : count1 - 1, count1 + 1 : end);
                        branch_mesh = branch_mesh(1 : count1 - 1, count1 + 1 : end);
                        branch_count = branch_count - 1; % refresh number of branches
                        break % all done; exit algorithm
                    end
                end 
            end
            break % exit algorithm
        end
    end
end


% Stores branches connected to each node
node_branches = cell(numel(nodes_ind), 1);
for count = 1 : numel(nodes_ind) % counter of nodes
    for count1 = 1 : branch_count % branch loop
        if or(nodes_ind(count) == branches_ind{count1}(1), nodes_ind(count) == branches_ind{count1}(end))
            % branch starts or ends at the node
            node_branches{count} = [node_branches{count} count1]; % add index of branch
        end
    end
end


% Stores branches that form each mesh
n_mesh = size(mesh_ind, 1);
mesh_branches = cell(n_mesh, 1);
for count = 1 : n_mesh % mesh loop
    for count1 = 1 : branch_count % branch loop
        if any(count == branch_mesh{count1}) % the branch belongs to the current mesh
            mesh_branches{count} = [mesh_branches{count}; count1]; % add index of accordant branch
        end
    end
end
% Re-order branches
for count_mesh = 1 : n_mesh % mesh loop
    if numel(mesh_branches{count_mesh}) > 2 % if there are only one or two branches, the current order is accepted
        vec_aux = zeros(numel(mesh_branches{count_mesh}), 1); % stores the ordered branches in the mesh
        vec_aux(1) = mesh_branches{count_mesh}(1); % the first branch it finds is chosen as the first of the mesh
        % Looks for the following branches
        for count_branch1 = 1 : numel(vec_aux) - 1 % after last branch comes first branch again (cycle)
            ind_aux = vec_aux(count_branch1); % current branch index
            for count_branch2 = 2 : numel(mesh_branches{count_mesh}) % branches to compare with % first branch cannot be
                ind_aux1 = mesh_branches{count_mesh}(count_branch2); % index of branch to compare
                if ind_aux ~= ind_aux1 % does not ckeck itself
                    if branches_ind{ind_aux1}(1) == branches_ind{ind_aux}(end)
                        % second branch starts with the last object of first branch
                        vec_aux(count_branch1 + 1) = ind_aux1; % is the next element
                        break % next branch
                    elseif branches_ind{ind_aux1}(end) == branches_ind{ind_aux}(end)
                        % second branch ends with the last object of first branch
                        branches_ind{ind_aux1} = flip(branches_ind{ind_aux1}); % inverts branch
                        vec_aux(count_branch1 + 1) = ind_aux1; % is the next element
                        break % next branch
                    end
                end
            end
        end
    mesh_branches{count_mesh} = vec_aux; % ordered branches in the mesh
    end
end


% Adds every object's branch to the objects table
objects.branch = cell(n_obj, 1);
for count = 1 : branch_count % branch loop
    for count1 = branches_ind{count} % objects in a branch loop
    	objects.branch{count1} = [objects.branch{count1} count];
    end
end


% Converts branches to change indices by id's
branches_id = cell(branch_count, 1); % store all branches (id)
for count = 1 : branch_count % branch loop
    for count1 = 1 : numel(branches_ind{count}) % objects in a branch
        branches_id{count}(count1) = objects.id(branches_ind{count}(count1)); % assign id
    end
end


% Inverts branches that are in the wrong direction
% Flow direction is determined in some objects (and their branches):
% pumps and heat exchangers of type 'T out'

ind_outlet_aux = [pump_volum.obj_index pump_volum.outlet; pump_turbo.obj_index pump_turbo.outlet; ...
    heat_exch_Tout.obj_index heat_exch_Tout.outlet; tank.obj_index tank.outlet];
    % column 1: object index of all volumetric pumps, turbopumps, tanks and heat exchangers of type 'T out'
    % column 2: id ot the outlets of those objects

for count = 1 : size(ind_outlet_aux, 1) % counter of objects that have a direction
    branch_aux = objects.branch{ind_outlet_aux(count, 1)}; % branch that contains the object
    ind_aux = find(branches_ind{branch_aux} == ind_outlet_aux(count, 1)); % position of the object in the branch
    ind_aux1 = find(branches_id{branch_aux} == ind_outlet_aux(count, 2)); % position of the outlet object in the branch
    % Compare direction of branch and object
    if ind_aux1 < ind_aux % the outlet is located earlier in the branch % opposite direction
        % this works for tanks in the return branch because (empty < x) == false
        branches_id{branch_aux} = flip(branches_id{branch_aux}); % invert branch
        branches_ind{branch_aux} = flip(branches_ind{branch_aux});
    end
end


end
