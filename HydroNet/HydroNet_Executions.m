function ...
    [ flows, branches_ind, branches_id, branch_temp_pos, branch_temperature, branch_htx, branch_pump_volum, branch_pump_turbo, ...
	heat_exch, obj_inlet_pos, obj_outlet_pos,pump_volum, pump_turbo, head_loss, hydr_resist1, hydr_resist2 , thermostat] ...
    = HydroNet_Executions ( fluid_type, objects, ...
	pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch, tank, ...
	nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, ...
	n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx, branch_pump_volum, branch_pump_turbo, ...
	branch_temp_pos, branch_temperature, gravity_acc, dt, flag_pump_speed_change, flag_valve_change, ...
	flag_thermostat_changes, flag_first_exec, head_loss, hydr_resist1, hydr_resist2, flows, weight_fraction )
    % Main execution of the HydroNet model


%% HEAD LOSS IN BRANCHES
% Calculates head loss in every branch

if flag_first_exec || flag_valve_change || flag_thermostat_changes || flag_pump_speed_change

    [ head_loss, hydr_resist1, hydr_resist2 ] = HydroNet_HeadLoss( objects, pipe, valve_fix, ...
        valve_var, thermostat, pump_turbo, heat_exch, tank, n_branch);                  
end



%% FLOWS IN BRANCHES and HEAD IN VOLUMETRIC PUMPS

% Sets boundary conditions: flow = 0
% - All branches that are not inside a closed cycle
% - All branches that contain closed valves or thermostats

bound_flows = nan(n_branch, 1); % vector filled with NaN to store boundary conditions
bound_flows(not(branch_cycle)) = 0; % outside a closed cycle
for count = 1 : size(valve_var, 1)
    if valve_var.opening(count) == 0
        obj_aux = valve_var.obj_index;
        branch_aux = objects.branch{obj_aux};
        bound_flows(branch_aux) = 0;
    end
end
% In thermostats, this additionally avoids the calculation of very small flows
% (which may cause problems for the solver and for heatings/coolings)
for count = 1 : size(thermostat, 1)
%     if thermostat.opening(count) < 0.01
    if thermostat.opening(count) < 0.01
        obj_aux = thermostat.obj_index;
        branch_aux = objects.branch{obj_aux(count)};
        bound_flows(branch_aux) = 0;
    end
end


% Loop to solve flows
if flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes

    solver_flag = true;

    while solver_flag

        solver_flag = false;

        [ flows, head_vol_pump ] = HydroNet_FlowSolver( objects, pump_volum, branches_id, flows, bound_flows, ...
        nodes_id, mesh_branches, node_branches, n_mesh, n_branch, head_loss, hydr_resist1, hydr_resist2, n_tanks );

        % Checks whether the obtained head of volumetric pumps are higher than the maximum
        % head of each volumetric pump. In such case, assign flow = 0 as boundary
        % condition and call the solver again
        for count = 1 : size(pump_volum, 1) % counter of volumetric pumps
            branch_aux = objects.branch{pump_volum.obj_index(count)}; % pump branch
            pump_volum.pump_head(count) = head_vol_pump(branch_aux); % store in table
            if abs(pump_volum.pump_head(count)) > pump_volum.head_max(count)
                % calculated head is higher than maximum head
                bound_flows(branch_aux) = 0; % assign flow = 0 as boundary condition
                solver_flag = true; % re-calculate flows
            end
            
        end
    end


    % When needed, branch directions are flipped so that flow and branch directions are the same
    for count_branch = transpose(find(flows < 0))
        branches_ind{count_branch} = flip(branches_ind{count_branch});
        branches_id{count_branch} = flip(branches_id{count_branch});
        branch_temp_pos{count_branch} = 100 - flip(branch_temp_pos{count_branch});
        branch_temperature{count_branch} = flip(branch_temperature{count_branch});
        branch_htx{count_branch} = flip(branch_htx{count_branch});
        branch_pump_volum{count_branch} = flip(branch_pump_volum{count_branch});
        branch_pump_turbo{count_branch} = flip(branch_pump_turbo{count_branch});

        for count_obj = branches_ind{count_branch}
            if any(nodes_ind==count_obj) % object is a node
                branch_pos_in_node = find(node_branches{nodes_ind==count_obj} == count_branch);
                    % index in a node of node_branches of the current branch
                obj_inlet_pos_aux = obj_inlet_pos{count_obj}; % saves this vector before overwriting it
                obj_inlet_pos{count_obj}(branch_pos_in_node) = 100 - obj_outlet_pos{count_obj}(branch_pos_in_node); % oulets are now inlets
                obj_outlet_pos{count_obj}(branch_pos_in_node) = 100 - obj_inlet_pos_aux(branch_pos_in_node); % inlets are now outlets 
            else % object is not a node
                obj_inlet_pos_aux = obj_inlet_pos{count_obj}; % saves this vector before overwriting it
                obj_inlet_pos{count_obj} = 100 - obj_outlet_pos{count_obj}; % oulets are now inlets
                obj_outlet_pos{count_obj} = 100 - obj_inlet_pos_aux; % inlets are now outlets     
            end
        end
    end
    flows = abs(flows); % All flows are positive now

end

%   display(flows);

%% PUMP POWER

if flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes
    
    % Volumetric pumps
    for count = 1 : size(pump_volum, 1) % counter of volumetric pumps
        object_aux = pump_volum.obj_index(count); % index of pump in the object's list
        branch_aux = objects.branch{object_aux}; % branch that contains the pump (only one)
        pos_aux = obj_inlet_pos{object_aux} + obj_outlet_pos{object_aux} / 2;
        temperature_aux = HydroNet_GetPosTemperature(pos_aux, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
            % temperature at the pump's middle
        
        [density,~]=fluidsproperties(fluid_type,temperature_aux,weight_fraction);
       

        pump_volum.pump_power(count) = density * gravity_acc * head_vol_pump(branch_aux) * ...
            flows(branch_aux) / efficiency(head_vol_pump(branch_aux), flows(branch_aux));
    end

    % Turbopumps
    for count = 1 : size(pump_turbo, 1) % counter of turbopumps
        object_aux = pump_turbo.obj_index(count); % index of pump in the object's list
        branch_aux = objects.branch{object_aux}; % branch that contains the pump (only one)
        pos_aux = obj_inlet_pos{object_aux} + obj_outlet_pos{object_aux} / 2;
        temperature_aux = HydroNet_GetPosTemperature(pos_aux, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
            % temperature at the pump's middle

        % Head curve of pump: head(speed, flow)
        head_aux = (pump_turbo.coef_N0(count) + pump_turbo.coef_N1(count) * pump_turbo.pump_speed + ...
            pump_turbo.coef_N2(count) * pump_turbo.pump_speed^2) + ...
            flows(branch_aux) * (pump_turbo.Q_coef_N0(count) + pump_turbo.Q_coef_N1(count) * ...
            pump_turbo.pump_speed + pump_turbo.Q_coef_N2(count) * pump_turbo.pump_speed^2) + ...
            flows(branch_aux)^2 * (pump_turbo.Q2_coef_N0(count) + pump_turbo.Q2_coef_N1(count) * ...
            pump_turbo.pump_speed + pump_turbo.Q2_coef_N2(count) * pump_turbo.pump_speed^2);
        
        [density,~]=fluidsproperties(fluid_type,temperature_aux,weight_fraction);

        pump_turbo.pump_head(count) = head_aux;
        pump_turbo.pump_power(count) = density * gravity_acc * head_aux * ...
            flows(branch_aux) / efficiency(head_vol_pump(branch_aux), flows(branch_aux));
    end
end




%% TEMPERATURES
% for count_htxch = 1:size(heat_exch,1)
if and(heat_exch.heat(2)<=0,heat_exch.inlet_temp(2) < 365)
    heat_exch.heat(2) = 0;
else
    heat_exch.heat(2) = -110000;
end
% end


[ branch_temp_pos, branch_temperature, heat_exch, pump_volum, pump_turbo, thermostat ] = HydroNet_Temperature(fluid_type, objects, obj_inlet_pos, ...
    obj_outlet_pos, heat_exch, pump_volum, pump_turbo, nodes_ind, branches_ind, node_branches, n_nodes, n_branch, branch_volume, ...
    branch_temp_pos, branch_temperature, branch_htx, branch_pump_volum, branch_pump_turbo, flows, dt, weight_fraction, thermostat);


end
