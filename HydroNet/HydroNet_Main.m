function [ ] = HydroNet_Main( )
% Controls data flow


%% MODEL CREATION
% All valves and thermostats are assumed open.

circuit_file = 'Coolant_basic_nodes.txt';
%circuit_file = 'Circuit_test_1.txt';

[ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, ...
    heat_exch_fix, tank, nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, ...
    node_branches, n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, ...
    branch_htx_fix] = HydroNet_Create( circuit_file );



%% INITIALIZATION

gravity_acc = 9.81; % gravity acceleration [m/s^2] for pump power

% Cell to store volumetric position of the points where temperature changes,
% as a % of the total volume of the branch % temperature between those points is constant
% branch_temp_pos = cell(n_branch, 1);
branch_temp_pos = {[0,10,100]; [0, 100]; [0, 100]; [0,20,100]; [0,30,100]; [0, 100]; [0,10,100]};

% Cell to store temperatures of every volume section in the branch
% branch_temperature = cell(n_branch, 1);
branch_temperature = {[300,250]; 350; 400; [450,300]; [500,400]; 550; [600,550]};

% Output variables related to branches
head_loss = zeros(n_branch, 1);
hydr_resist1 = zeros(n_branch, 1);
hydr_resist2 = zeros(n_branch, 1);
flows = zeros(n_branch, 1);



%% EXECUTIONS
% Some variables must be recalculated only if other change

n_exec = 5;

dt_vec = [0.1, 0.1, 0.1, 0.1, 0.2];

volum_pump_speed = {105; 210; 314; 419; 419};
turbo_pump_speed = {105; 210; 314; 419; 419};
valve_opening = {20; 50; 70; 80; 30};

thst_sensitivity = ones(n_exec, 1) * 5; % thermostat sensitivity [%]


% Execution loop
flag_first_exec = true;
for count_exec = 1 : n_exec
    
    % TIME SPAN
    dt = dt_vec(count_exec);
    
    
    % PUMP SPEED
    % If the speed of any pump changes, flows must be recalculated
    
    flag_pump_speed_change = false;
    if or(any(not((volum_pump_speed{count_exec} == pump_volum.pump_speed(:)))), ...
        any(not((turbo_pump_speed{count_exec} == pump_turbo.pump_speed(:)))))
        flag_pump_speed_change = true;
    end
    
    pump_volum.pump_speed(:) = volum_pump_speed{count_exec};
    pump_turbo.pump_speed(:) = turbo_pump_speed{count_exec};
    
    
    % VALVES
    % If the valve opening varies, vectors for head losses and thus flows must be recalculated.
    % If the valve is closed, flow in its branch is zero.
    
    flag_valve_change = false;
    if any(not((valve_opening{count_exec} == valve_var.opening(:))))
        flag_valve_change = true;
    end
    
    valve_var.opening(:) = valve_opening{count_exec};
    
    
    % THERMOSTATS
    % If the thermostat opening varies because there is a change in the inlet
    % temperature, vectors for head losses and thus flows must be recalculated.
    % Head losses and flows will be recalculated only if the opening variation
    % of any thermostat is larger than the sensitivity of the thermostat[%]
    % If the thermostat is closed, flow in its branch is zero.
    
    thermostat_opening_new = zeros(size(thermostat, 1), 1);
    for count_thst = 1 : size(thermostat, 1);
        obj_aux = thermostat.obj_index(count_thst);
        branch_aux = objects.branch{obj_aux};
        obj_pos = (obj_inlet_pos{obj_aux} + obj_outlet_pos{obj_aux}) / 2;
        thermostat_temp = HydroNet_GetObjTemperature(obj_pos, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
        thermostat_opening_new(count_thst) = thermostat.T_coef0(count_thst) + thermostat.T_coef1(count_thst) * thermostat_temp ...
            + thermostat.T_coef2(count_thst) * thermostat_temp^2;
    end
    
    flag_thermostat_changes = false;
    if any(abs(thermostat_opening_new - thermostat.opening(:)) ./ thermostat.opening(:) * 100 > thst_sensitivity(count_exec))
        flag_thermostat_changes = true;
    end
    
    thermostat.opening(:) = thermostat_opening_new;
    
    
    % MAIN EXECUTION FUNCTION CALL
    [ flows, branches_ind, branches_id, branch_temp_pos, branch_temperature, branch_htx_fix, branch_htx_Tout, ...
        heat_exch_Tout, heat_exch_fix, obj_inlet_pos, obj_outlet_pos ] = HydroNet_Executions ( fluid_type, objects, ...
        pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch_Tout, heat_exch_fix, tank, ...
        nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, ...
        n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx_Tout, branch_htx_fix, ...
        branch_temp_pos, branch_temperature, gravity_acc, dt, flag_pump_speed_change, flag_valve_change, ...
        flag_thermostat_changes, flag_first_exec, head_loss, hydr_resist1, hydr_resist2, flows ); 
    
    flag_first_exec = false;
    
    display(flows);
end


end
