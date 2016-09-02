function ...
    [ flows, branches_ind, branches_id, branch_temp_pos, branch_temperature, branch_htx_fix, branch_htx_Tout, ...
	heat_exch_Tout, heat_exch_fix, obj_inlet_pos, obj_outlet_pos ] = HydroNet_Executions ( fluid_type, objects, ...
	pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, heat_exch_fix, tank, ...
	nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, ...
	n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, branch_htx_fix, ...
	branch_temp_pos, branch_temperature, gravity_acc, dt, flag_pump_speed_change, flag_valve_change, ...
	flag_thermostat_changes, flag_first_exec, head_loss, hydr_resist1, hydr_resist2, flows )
    % Main execution of the HydroNet model


%% HEAD LOSS IN BRANCHES
% Calculates head loss in every branch

if flag_first_exec || flag_valve_change || flag_thermostat_changes

    [ head_loss, hydr_resist1, hydr_resist2 ] = HydroNet_HeadLoss( objects, pipe, valve_fix, ...
        valve_var, thermostat, pump_turbo, heat_exch_fix, heat_exch_Tout, tank, branches_ind, ...
        branch_cycle, n_branch, obj_inlet_pos, obj_outlet_pos, branch_temp_pos, branch_temperature);
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
for count = 1 : size(thermostat, 1)
    if thermostat.opening(count) == 0
        obj_aux = thermostat.obj_index;
        branch_aux = objects.branch{obj_aux};
        bound_flows(branch_aux) = 0;
    end
end


% Loop to solve flows
if flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes

    solver_flag = true;

    while solver_flag

        solver_flag = false;

        [ flows, head_vol_pump ] = HydroNet_FlowSolver( objects, pump_volum, branches_id, bound_flows, ...
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
        branch_htx_fix{count_branch} = flip(branch_htx_fix{count_branch});
        branch_htx_Tout{count_branch} = flip(branch_htx_Tout{count_branch});

        for count_obj = branches_ind{count_branch}
            obj_inlet_pos_aux = obj_inlet_pos{count_obj}; % saves this vector before overwriting it
            obj_inlet_pos{count_obj} = 100 - obj_outlet_pos{count_obj}; % oulets are now inlets
            obj_outlet_pos{count_obj} = 100 - obj_inlet_pos_aux; % inlets are now outlets        
        end
    end
    flows = abs(flows); % All flows are positive now

end



%% PUMP POWER

if flag_first_exec || flag_pump_speed_change || flag_valve_change || flag_thermostat_changes
    
    % Volumetric pumps
    for count = 1 : size(pump_volum, 1) % counter of volumetric pumps
        object_aux = pump_volum.obj_index(count); % index of pump in the object's list
        branch_aux = objects.branch{object_aux}; % branch that contains the pump (only one)
        pos_aux = obj_inlet_pos{object_aux} + obj_outlet_pos{object_aux} / 2;
        temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
            % temperature at the pump's middle

        pump_volum.pump_power(count) = density(fluid_type, temperature_aux) * gravity_acc * head_vol_pump(branch_aux) * ...
            flows(branch_aux) * efficiency(head_vol_pump(branch_aux), flows(branch_aux));
    end

    % Turbopumps
    for count = 1 : size(pump_turbo, 1) % counter of turbopumps
        object_aux = pump_turbo.obj_index(count); % index of pump in the object's list
        branch_aux = objects.branch{object_aux}; % branch that contains the pump (only one)
        pos_aux = obj_inlet_pos{object_aux} + obj_outlet_pos{object_aux} / 2;
        temperature_aux = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
            % temperature at the pump's middle

        % Head curve of pump: head(speed, flow)
        head_aux = (pump_turbo.coef_N0(count) + pump_turbo.coef_N1(count) * pump_turbo.pump_speed + ...
            pump_turbo.coef_N2(count) * pump_turbo.pump_speed^2) + ...
            flows(branch_aux) * (pump_turbo.Q_coef_N0(count) + pump_turbo.Q_coef_N1(count) * ...
            pump_turbo.pump_speed + pump_turbo.Q_coef_N2(count) * pump_turbo.pump_speed^2) + ...
            flows(branch_aux)^2 * (pump_turbo.Q2_coef_N0(count) + pump_turbo.Q2_coef_N1(count) * ...
            pump_turbo.pump_speed + pump_turbo.Q2_coef_N2(count) * pump_turbo.pump_speed^2);

        pump_turbo.pump_power(count) = density(fluid_type, temperature_aux) * gravity_acc * head_aux * ...
            flows(branch_aux) * efficiency(head_vol_pump(branch_aux), flows(branch_aux));
    end
end



%% TEMPERATURES

[ branch_temp_pos, branch_temperature ] = HydroNet_Temperature(fluid_type, objects, obj_inlet_pos, ...
    obj_outlet_pos, heat_exch_Tout, heat_exch_fix, nodes_ind, branches_ind, node_branches, n_nodes, n_branch, branch_volume, ...
    branch_temp_pos, branch_temperature, branch_htx_Tout, branch_htx_fix, flows, dt);


end
