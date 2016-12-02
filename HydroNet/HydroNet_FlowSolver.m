function [ flows, head_vol_pump ] = HydroNet_FlowSolver( objects, pump_volum, branches_id, flows_prev, bound_flows, ...
    nodes_id, mesh_branches, node_branches, n_mesh, n_branch, head_loss, hydr_resist1, hydr_resist2, n_tanks )
% Solver for hydraulic networks.
% Receives information about meshes, branches and nodes in the network.
% Receives known flows imposed as boundary conditions - a vector with the flow
% values in the corresponding branches and NaN in the rest.
% Receives head losses in every branch. Data must be divided in three
% vectors: one for head losses that do not depend on the flow, another one
% for head losses proportional to the flow and the last one for losses
% proportional to the flow squared.
% This function attempts to solve the network with algebraic methods. If
% some flows cannot be solved, it uses a numerical method.
% Receives information about volumetric pumps. If there are volumetric pumps
% in the circuit, the fucntion returns its pump head.


% % % % % %
% SOLVING PROCEDURE
%
% First, solves branches where flow is known, is zero or can be calculated
%   independently because there is a volumetric pump (flow depends on pump
%   speed). Check if all branches are solved.
%
% Second, enters a loop:
    % NODES STEP: Solves nodes where there is only one unknown flow. Check
        % if all branches are solved. Loop until no more nodes are solved.
    % If it is not the first iteration and no branches have been solved in
    % the nodes step, break main loop.
    % MESH STEP: solves meshes that depend only on one flow (simple meshes;
        % only one branch remains to be solved). Check if all branches are
        % solved. Loop until no more meshes are solved.
    % If no branches have been solved in the mesh step, break main loop.
    % Loop
%
% Third, sets up and solves the system of non-linear equations. Equations
%   depend on flow (nodes, etc.), on flow^2 (pipes, etc.) and on head of
%   volumetric pumps. All branches are solved.
% % % % % %


% Creates result vectors
flows = NaN(n_branch, 1); % stores flow values
solved = false(n_branch, 1); % stores solving status



%% 1. INDEPENDENT FLOWS

% Assigns flows entered as a boundary condition
flows(not(isnan(bound_flows))) = bound_flows(not(isnan(bound_flows)));
solved(not(isnan(bound_flows))) = true;
if solved % all branches are solved 
    return % end function
end


% Solves closed branches containing VOLUMETRIC PUMPS
volum_pump = false(n_branch, 1); % marks where are there volumetric pumps for later use
for count = 1 : size(pump_volum, 1) % counter of volumetric pumps
    branch_aux = objects.branch{pump_volum.obj_index(count)}; % pump branch
     if not(solved(branch_aux)) % the branch is not solved yet
        volum_pump(branch_aux) = true; % there is a volumetric pump
        % calculates flow
        flows(branch_aux) = pump_volum.coef0(count) ...
            + pump_volum.coef1(count) * pump_volum.pump_speed(count) ...
            + pump_volum.coef2(count) * pump_volum.pump_speed(count)^2 ...
            + pump_volum.coef3(count) * pump_volum.pump_speed(count)^3;
        if find(branches_id{branch_aux} == pump_volum.id(count)) > find(branches_id{branch_aux} == pump_volum.outlet(count))
            % branch is in the opposite direction of the pump
            flows(branch_aux) = - flows(branch_aux);
        end
        solved(branch_aux) = true; % set as solved
     end
end
% There is no check of whether all branches are solved because volumetric
% pump head must still be found
head_vol_pump = nan(n_branch, 1);
    % vector with the head value of every volumetric pump in its corresponding branch




%% 2. LOOP FOR NODES WITH ONE UNKNOWN FLOW

n_nodes = numel(nodes_id);
solv_nodes = false(n_nodes, 1); % vector to mark nodes whose all connected branches are solved


% Loop for solving NODES and MESHES that have only one unknown flow.

while_flag = false; % flag to know if it is the first iteration of the while loop
solv_prev1 = nan(length(solved), 1);  % flows solved in the previous step (main loop)

% (NODES + MESHES) LOOP
while not(all(solved == solv_prev1)) % in the previous step, some flow was solved (mesh sub-loop)

     % NODE SUB-LOOP
    
    solv_prev1 = solved; % flows solved after the meshes step or at start
    solv_prev2 = nan(length(solved), 1); % flows solved in the previous step (sub-loop)
    
    while not(all(solved == solv_prev2)) % in the previous step, some flow was solved

        solv_prev2 = solved; % flows solved after the nodes sub loop
        for count = 1 : n_nodes
            count_uns = 0; % counter of unsolved branches
            for count1 = 1 : numel(node_branches{count}) % branches connected to node
                if solved(node_branches{count}(count1)) == false % branch unsolved
                    ind_aux = count1; % save index of unsolved branch
                    count_uns = count_uns + 1; % counter of unsolved branches
                    if count_uns == 2 % too many unknowns
                        break % go to next node
                    end
                end
            end
            
            if count_uns == 0 % all connected branches are solved
                solv_nodes(count) = true; % set as solved
                
            elseif count_uns == 1 % one unknown; node can be solved directly
                
                % Continuity equation: sum of flows in a node = 0
                flow_sum = 0;
                for count1 = node_branches{count}
                    % branches connected to node
                    % Node can be at the end or at the begining of the branch
                    % Criterium: positive flow --> enters node
                    if and(nodes_id(count) == branches_id{count1}(end), not(isnan(flows(count1)))) % at the end, keep sign
                        % flow enters the node
                        flow_sum = flow_sum + flows(count1);
                    elseif and(nodes_id(count) == branches_id{count1}(1), not(isnan(flows(count1)))) % at the beginning, change sign
                        % flow leaves the node
                        flow_sum = flow_sum - flows(count1); % change sign
                    end
                end
                
                ind_aux = node_branches{count}(ind_aux); % change ind_aux to refer directly to the branch or flow index
                
                % Node can be at the end or at the begining of the branch
                if nodes_id(count) == branches_id{ind_aux}(1) % at the beginning
                    flows(ind_aux) = flow_sum;
                else
                    flows(ind_aux) = - flow_sum; % change sign
                end
                
                solved(ind_aux) = true; % set branch as solved
                if and(solved, not(volum_pump)) % all branches are solved (flows and volumetric pump head)
                    return
                end
                solv_nodes(count) = true; % set node as solved
            end
        end
    end % while for nodes

    if and(while_flag == true, solved == solv_prev1)
       % it is not the first iteration and no flows have been solved during the nodes step
       % no more flows will be solved in the meshes step because there have been no changes
       break % go to next solving method
    end

    solv_prev1 = solved; % flows solved after the nodes step
    
    %%%%%%%%%%%%%% middle of loop
    
    % MESH SUB-LOOP
    
    solv_prev2 = nan(length(solved), 1);  % flows solved in the previous step ( node sub-loop)
    
    while not(solved == solv_prev2) % in the previous step, some flow was solved (during node sub-loop)

        solv_prev2 = solved; % flows solved after the mesh sub-loop
        
        for count = 1 : n_mesh % meshes
            count_uns = 0; % counter of unsolved branches
            for count1 = transpose(mesh_branches{count}) % branches that form the mesh
                if not(solved(count1))|| volum_pump(count1)||flows(count1)==0
                    % branch unsolved or contains volumetric pump (pump head must be calculated)
                    ind_aux = count1; % save unsolved branch
                    count_uns = count_uns + 1; % counter of unsolved branches
                    if count_uns == 2 % too many unknowns
                        break % go to next mesh
                    end
                end
            end
            if count_uns == 1 % one unknown; mesh can be solved directly
                
                % II Kirchhoff's Law: sum of head increments in a mesh = 0
                total_head = 0;
                last_id = 0;
                for count1 = transpose(mesh_branches{count}) % branches that form the mesh
                    % Unsolved one must be included in the loop to find out mesh flow direction
                    % Branches are already ordered along the mesh
                    % Criterium: positive mesh flow --> direction of the first branch in the mesh
                    % If there is a pump or a heat exchanger of type 'T out' in the branch, flow and branch direction are the same
                    % If flow direction is unknown, assume same direction as branch
                    % Flow may be solved or not

                    % Check whether the branch has the same direction as the mesh flow
                    change_sign = false;
                    if or(last_id == 0, branches_id{count1}(1) == last_id)
                        % it is the first branch in the mesh OR the branch is in the same direction as the mesh flow
                        last_id = branches_id{count1}(end); % refresh other end of branch
                    else % it is not the first branch in the mesh AND the branch is in the opposite direction
                        change_sign = true; % invert sign
                        last_id = branches_id{count1}(1); % refresh other end of branch
                    end

                    if or(count1 ~= ind_aux, volum_pump(count1)) % flow in the branch is known
                        if flows(count1) < 0 % flow direction is opposite to branch direction
                            change_sign = not(change_sign); % invert sign
                        end
                        % Add head increment
                        total_head = total_head + (head_loss(count1) + hydr_resist1(count1) ...
                            * abs(flows(count1)) + hydr_resist2(count1) * flows(count1)^2) * (-2 * change_sign + 1);
                            % add head if change_sign = false; subtract if change_sign = true
                    else % branch unknown is a flow: save sign - false: positive, true: negative, 
                        sign_aux = change_sign;
                    end
                end
                
                % Solve unknown
                if volum_pump(ind_aux) % the unknown branch has a volumetric pump
                    % Solve unknown head increment
                    head_vol_pump(ind_aux) = abs(total_head); % calculate pump head
                    head_loss(ind_aux) =  head_loss(ind_aux) - abs(total_head);
                        % consider pump head as a negative head loss in the branch
                    volum_pump(ind_aux) = false; % remove branch from branches with volumetric pumps
                else % the unknown is a flow; solve it
                    if total_head < 0 % branch flow is in the opposite direction of the mesh flow
                         sign_aux = not(sign_aux); % invert branch flow sign
                    end
                    
                    % Solve branch
                    sol_aux = roots([hydr_resist2(ind_aux) hydr_resist1(ind_aux) (head_loss(ind_aux)-abs(total_head))]);
                        % absolute value of head increment in the branch and
                        % absolute value of total head must be equal
                    % If there are two solutions, choose the positive one
                    if numel(sol_aux) == 2
                        if sol_aux(1) > 0
                            sol_aux = sol_aux(1);
                        else
                            sol_aux = sol_aux(2);
                        end
                    end
                    
                    % Save result
                    flows(ind_aux) = sol_aux * (1+2*sign_aux*(-1)); % add sign now
                    solved(ind_aux) = true; % set as solved
                end
            end
            
            if and(solved, not(volum_pump)) % all branches are solved (flows and volumetric pump head)
                return
            end
        end
    end % while for meshes
    
    if solved == solv_prev1
        % no flows have been solved during the meshes step
        % no more flows will be solved in the nodes step because there have been no changes
        break % go to next solving method
    end
    
    while_flag = true; % first iteration finished
end

% display(solved);


%% 3. SYSTEM OF EQUATIONS
% A * X*|X| + B * X + C * X/|X| + D = 0 --> A = coef2_mat, B = coef1_mat, C = const_mat, D = const_vec
% C is for constant head loss whose sign is unknown

n_unknown = n_branch - sum(solved) + sum(volum_pump); % number of unknown variables:
    % flows - solved flows + volumetric pump head
    % in branches with volumetric pumps, pump head is unknown

% Coefficient matrices
coef2_mat = zeros(n_unknown); % coefficients for unknowns^2 in each equation (row)
coef1_mat = zeros(n_unknown); % coefficients that multiply unknowns in each equation (row)
const_mat = zeros(n_unknown); % matrix for independent head loss in unsolved branches
const_vec = zeros(n_unknown, 1); % vector for constant terms of the equations:
    % independent head loss and head calculated with known flows

% Auxiliary vector that contains the indices of the unsolved branches
unk_ind_aux = [find(volum_pump); find(not(solved))];

% Auxiliary vector that contains 1 at the indices where the vector
% 'unknowns' or 'X' refer to branches that have volumetric pumps
if size(pump_volum,1)==0
    volpump_ind_aux = unk_ind_aux*0;
else
    volpump_ind_aux = unk_ind_aux==find(volum_pump);
end
% Auxiliary vector that contains 1 at the indices where the vector
% 'unknowns' or 'X' refer to flows
flows_ind_aux = not(volpump_ind_aux);


% NODE EQUATIONS
% All unsolved nodes are considered for the system but one of them is
% discarded since it will be dependent on the rest.
% Number of node equations = number of unsolved nodes - 1
% Exception: there are tanks in the network.

%Set last unsolved node as solved. Unsolved nodes will be included in the
%system of equations.
if n_tanks == 0
    for count = 1 : n_nodes
        if not(solv_nodes(n_nodes - count + 1)) % node not solved
            solv_nodes(n_nodes - count + 1) = true; % set node as solved
            break % done; exit
        end
    end
end

% Formulates node equations
last_row=0;
for count = find(transpose(not(solv_nodes))) % unsolved nodes loop
    last_row=last_row+1;
    % Continuity equation: sum of flows in a node = 0
    for count1 = node_branches{count} % branches connected to node (branch index)
        change_sign = false;
        if branches_id{count1}(1) ~= nodes_id(count)
            change_sign=true;
        end
        
        % Flow may be solved or not
        if solved(count1) % the flow in this branch has been solved
            const_mat(last_row, (count1))  = flows(count1)*(1)*(-2*change_sign+1);
        else % flow is unknown
            coef1_mat(last_row, (unk_ind_aux==count1)) = 1*(-2*change_sign+1);
        end
       
        
        % Node can be at the end or at the begining of the branch
        % Criterium: positive branch direction --> enters node (node is at the end of the branch)
%         if nodes_id(count)  ~= branches_id{count1}(1) % node is at the beginning of the branch,
%             % branch exits node, change sign
%         	const_mat(last_row, (unk_ind_aux==count1)) = - const_mat(last_row, (unk_ind_aux==count1));
%         end
    end
end
%last_row = count; % save row to continue entering the equations


% MESH EQUATIONS

% Number of independent meshes: branches to solve - nodes to solve + 1 = branches to solve - independent nodes to solve
n_idp_mesh = sum(not(solved)) - sum(not(solv_nodes)) + sum(volum_pump) - n_tanks;
    % head of volumetric pumps must be solved as well
    % additionaly for every tank, there is a fake branch and two known pressure nodes: -1 mesh per tank

% Selects meshes for the system of equations. Criterion: mesh must contain
% the maximum possible number of unsolved branches. Next selected meshes
% must contain the maximum number of unsolved branches excluding the
% unsolved branches contained in meshes previously chosen. This way mesh
% independence and solution of all branches are guaranteed.

% Vector that stores the indices of the meshes that are included in the system of equations
idp_mesh = zeros(n_idp_mesh, 1);

% Vector that marks the branches that are not covered by the current
% selection of mesh equations. Remaining branches to be added to the system.
remain_branch = not(solved) + volum_pump;
branch_merit = solved * (- n_idp_mesh); % solved branches get a penalization
if n_mesh > n_idp_mesh % there are more meshes than necessary -> selection
    for count = 1 : n_idp_mesh % selected mesh loop
        best_mesh = 0; % selected mesh for the system of equations
        most_uns = 0; % number of unsolved branches included in the selected mesh
        for count1 = 1 : n_mesh % mesh loop
            count_uns = 0; % counter of unsolved branches included in the selected mesh
            for count2 = transpose(mesh_branches{count1}) % branches contained in the mesh (vector index)
                if remain_branch(count2) % branch must be solved
                    count_uns = count_uns + 1; % add one point
                end
            end
            if count_uns > most_uns % solves more branches than previous meshes
                most_uns = count_uns; % refresh most unsolved branches
                best_mesh = count1; % refresh best mesh
            end
        end
        if best_mesh == 0 % any mesh will work
            highest_merit = -inf;
            for count1 = 1 : n_mesh % mesh loop
                count_merit = 1.2 * numel(mesh_branches{count1}); % the merit of having more branches is bonused by 20%
                if not(any(count1 == idp_mesh)) % looks for a mesh not chosen yet
                    for count2 = transpose(mesh_branches{count1}) % branches contained into mesh (vector index)
                        count_merit = count_merit + branch_merit(count2);
                    end
                    if count_merit > highest_merit % the merit of this mesh is higher than the previous ones
                        highest_merit = count_merit; % refresh highest merit
                        best_mesh = count1; % refresh best mesh
                    end
                end
            end
            idp_mesh(count) = best_mesh;
            for count2 = transpose(mesh_branches{count1}) % branches contained into mesh (vector index)
                branch_merit(count2) = branch_merit(count2) - 1; % penalize branches in the chosen mesh
            end
        else
            idp_mesh(count) = best_mesh; % save best mesh for the system
            % Remove all branches contained in the selected mesh from the remaining branches
            for count2 = transpose(mesh_branches{best_mesh}) % branches contained into mesh (vector index)
                remain_branch(count2) = false; % branch does not have to be solved
                branch_merit(count2) = branch_merit(count2) - 1; % penalize branches in the chosen mesh
            end
        end
    end
    
else % all meshes are necessary
    % save all meshes as independent
    idp_mesh = transpose(1:n_mesh);
end
idp_mesh_branches = mesh_branches(idp_mesh);


% Initial values for the non-linear solver
% if not(isnan(flows_prev))
    x0 = zeros(n_unknown, 1) ;
    for count_unkown = 1:n_unknown % same number of unknowns as branches
        if flows_ind_aux(count_unkown) % flow in the branch must be solved
            x0(count_unkown) = flows_prev(unk_ind_aux(count_unkown));
            if x0(count_unkown) == 0
                x0(count_unkown) = 0.0000001;
            end
        else % head of the volumetric pump must be solved
            pump_branch_pos = find(branches_id{unk_ind_aux(count_unkown)} == pump_volum.id);
            x0(count_unkown) = pump_volum.pump_head(branches_id{unk_ind_aux(count_unkown)}(pump_branch_pos))/10000;
                % pump head divided by 10000 so its magnitude order is close to that of the flows
        end
    end
% else
% x0 = ones(n_unknown, 1) *.0005; % size equal to number of unknowns
% end
    % flow 1 liter per second
    
% x0 = zeros(n_unknown, 1); % size equal to number of unknowns
% x0(flows_ind_aux) = 0.001; % flow 1 liter per second
% x0(volpump_ind_aux) = 10; % head  10 m.c.f
    % in the future, this vector can be feed with values obtained in a
    % previous call or from tables; remember to divide pump head by 10000 or so

    
% Formulates mesh equations
% flag_changed=false;
for count = 1 : n_idp_mesh % independent mesh loop

    % II Kirchhoff's Law: sum of head increments in a mesh = 0
    last_id = 0;
    for count1 = transpose(idp_mesh_branches{count}) % branches in the mesh (branches index)
        % Branches are already ordered along the mesh
        % Criterium: positive mesh flow --> direction of the first branch in the mesh
        % If there is a pump or a heat exchanger of type 'T out' in the branch, flow and branch direction are the same
        % If flow direction is unknown, assume same direction as branch
        % Flow may be solved or not

        % Check whether the branch has the same direction as the mesh flow
        change_sign = false;
        
        if last_id == 0 % is the first branch in the mesh
            last_id = branches_id{count1}(end); % initialization: assume that the last aboject according to mesh direction is the object at the end of the branch
        %if count1 < numel(idp_mesh_branches{count}) % is is not the last branch in the mesh
            next_branch_aux = idp_mesh_branches{count}(2); % determines next branch in the mesh

            if and(branches_id{count1}(end) ~= branches_id{next_branch_aux}(1), branches_id{count1}(end) ~= branches_id{next_branch_aux}(end))
                % last element of the present branch is not in the next branch (at the ends)
                change_sign = true; % invert sign
                last_id = branches_id{count1}(1);   % object at the end of the branch according to the mesh direction            
            end
        
        else % is not the first branch in the mesh
            if last_id == branches_id{count1}(end)
                change_sign = true; % invert sign
                last_id = branches_id{count1}(1); % object at the end of the branch according to the mesh direction
            else
                last_id = branches_id{count1}(end); % object at the end of the branch according to the mesh direction
            end
            
            
        end
        
%         if or(last_id == 0, branches_id{count1}(1) == last_id)
%             % it is the first branch in the mesh OR the branch is in the same direction as the mesh flow
%             last_id = branches_id{count1}(end); % refresh other end of branch
%         else % it is not the first branch in the mesh AND the branch is in the opposite direction
%             change_sign = true; % invert sign
%             last_id = branches_id{count1}(end); % refresh other end of branch
%         end

        if solved(count1) % the flow in this branch has been solved
            if flows(count1) <= 0 % flow direction is opposite to branch direction
                change_sign = not(change_sign); % invert sign
            end
%             x0(count1) = x0(count1) * (-2 * change_sign + 1);
%             const_mat(last_row + count, (unk_ind_aux==count1)) = const_vec(last_row + count) + (head_loss(count1) + hydr_resist1(count1) ...
%                 * abs(flows(count1)) + hydr_resist2(count1) * flows(count1)^2) * (-2 * change_sign + 1);
            const_vec(last_row + count) = const_vec(last_row + count) + (head_loss(count1) + hydr_resist1(count1) ...
                 * abs(flows(count1)) + hydr_resist2(count1) * flows(count1)^2) * (-2 * change_sign + 1);
                % add head if change_sign = false; subtract if change_sign = true
            if volum_pump(count1) % if there is a volumetric pump, include pump head as an unknown variable of the system
                coef1_mat(last_row + count, (unk_ind_aux==count1)) = - (-2 * change_sign + 1) * 10000; % include pump head
                    % change sign since head losses are positive and pump head must be positive too
                    % multiplies by 10000 so the solution's magnitude order is closer to flow values
            end

        else % flow is unknown
            % assume flow mesh and branch flow have the same direction
          
            const_mat(last_row + count, (unk_ind_aux==count1)) = head_loss(count1) *  (-2 * change_sign + 1);
            coef1_mat(last_row + count, (unk_ind_aux==count1)) = hydr_resist1(count1) * (-2 * change_sign + 1);
            coef2_mat(last_row + count, (unk_ind_aux==count1)) = hydr_resist2(count1) * (-2 * change_sign + 1);
        end
    end
end
% The system is complete


% Non-linear solver
f = @(x) coef2_mat*(x.*abs(x)) + coef1_mat*x + const_mat*(x./abs(x)) + const_vec;
opts = optimoptions('fsolve','MaxFunEvals',80000,'MaxIter',40000);%, 'Display','off');
Y = fsolve(f,x0,opts);


% Save flows
flows(not(solved)) = Y(flows_ind_aux);

% Save head in volumetric pumps 
if size(pump_volum,1) > 0 % there are volumetric pumps
	head_vol_pump(volum_pump) = Y(volpump_ind_aux) * 10000;
        % previously pump head was divided by 10000 so its magnitude order was
        % close to that of the flows
end


end

