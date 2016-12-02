function [ ] = HydroNet_Main()
% Controls data flow
flag_save=0;

%% MODEL CREATION
% All valves and thermostats are assumed open.

% circuit_file='Cyl&pump.txt';
%circuit_file='2Cyl&pump.txt';
%circuit_file='Cyl&Rad&pump.txt';
%circuit_file='Cyl&RadValve&pump.txt';
%circuit_file='Cyl&RadThermost&pump.txt';
%circuit_file ='Coolant_basic_nodes - copia.txt';
%circuit_file = 'Coolant_basic_nodes.txt';
%circuit_file = 'Circuit_test_1.txt';
%circuit_file = 'Cyl&pump&rad-valves.txt';
%circuit_file = 'test.txt';
%circuit_file='cyl-2cyl-3cyl-rad-2rad-pump.txt';
%circuit_file='TestEngine.txt';
%circuit_file='no_nodes _cyl_pump.txt';
circuit_file='Engine_Renault.txt';
[ fluid_type, objects, pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, ...
    heat_exch, tank, nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, ...
    node_branches, n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx, ...
    branch_pump_volum, branch_pump_turbo] = HydroNet_Create( circuit_file );



%% INITIALIZATION

T_initial = 320;

gravity_acc = 9.81; % gravity acceleration [m/s^2] for pump power

% Cell to store volumetric position of the points where temperature changes,
% as a % of the total volume of the branch % temperature between those points is constant
 branch_temp_pos =cell(n_branch, 1);
 for count=1:n_branch
     branch_temp_pos{count}=[0,100];
 end
%branch_temp_pos = {[0,10,100]; [0, 100]; [0, 100]; [0,20,100]; [0,30,100]; [0, 100]; [0,10,100]};
%branch_temp_pos = {[0,100];[0,100]};
%branch_temp_pos = {[0,100];[0,100];[0,100]};
%branch_temp_pos = {[0,100]; [0, 100]; [0,100]; [0,100]; [0,100]; [0, 100]; [0, 100]};
%branch_temp_pos = {[0,10,100]};
%branch_temp_pos = {[0,100]; [0, 100]; [0,100]; [0,100]; [0,100]; [0, 100]; [0, 100]};

% Cell to store temperatures of every volume section in the branch
 branch_temperature = cell(n_branch, 1);
 for count=1:n_branch
     branch_temperature{count} = T_initial;
 end
%branch_temperature = {[300,300]; 300; 300; [300,300]; [300,300]; 300; [300,300]};
%branch_temperature = {300;300};
%branch_temperature = {360;360;300};
%branch_temperature = {300;300;300};
%branch_temperature = {20;20;20;20;20;20;20};
%branch_temperature = {[20,20]};
%branch_temperature = {300;300;300;300;300;300;300};


% Output variables related to branches
head_loss = zeros(n_branch, 1);
hydr_resist1 = zeros(n_branch, 1);
hydr_resist2 = zeros(n_branch, 1);
flows = ones(n_branch, 1) * 0.0005;
pump_volum.pump_head(:) = 5;
% pump_volum.inlet_temp(:) = T_initial;
% pump_turbo.inlet_temp(:) = T_initial;
% heat_exch.inlet_temp(:) = T_initial;
thermostat.Ref_temp(:) = T_initial;
weight_fraction=0;
thermostat_opening_old=thermostat.opening(:);


% Some variables must be recalculated only if other change

n_exec = 1600;
dt = 0.05;
dt_vec = ones(n_exec, 1) * dt;


turbo_pump_speed = cell(n_exec, 1);
for count = 1 : 5000%n_exec/2-
    turbo_pump_speed{count} = 100;
end
for count = 5000 + 1 :n_exec 
    turbo_pump_speed{count} = 100;
end

volum_pump_speed = cell(n_exec, 1);
for count = 1 : n_exec/2
    volum_pump_speed{count} = 100;
end
for count = n_exec/2 + 1 :n_exec 
    volum_pump_speed{count} = 200;
end
  
%valve_opening = {0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;10;20;30;40;50;70;100;150;210;280;320;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360;360};
%valve_opening = {1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0};
% valve_opening = {1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0};

    
Tmin = thermostat.T_closed(:);
Tmax = thermostat.T_open(:);
y_end = thermostat.Shape_factor(:);
x0 = (Tmax+Tmin)/2;
k = - log(y_end) / (Tmin - x0);
thrm_speed = 0.01; % Thermostas maximum opening change per second


% Temperature register
temperature_thermostat = zeros(1, n_exec);
temperature_after_engine = zeros(1, n_exec);
temperature_after_radiator = zeros(1, n_exec);
time = dt : dt : n_exec*dt;

% Values for plotting
temp_inlet_cyl = zeros(1, n_exec);
temp_inlet_rad = zeros(1, n_exec);
flowss = zeros(n_branch,n_exec);
pump_power = zeros(1, n_exec);
speed = zeros(1, n_exec);
pump_head = zeros(1, n_exec);
thermos = zeros(2,n_exec);
[n_obj,~]=size(objects);
branchesTempInlet = ones(n_obj,n_exec)*-999;
branchesTempOutlet = ones(n_obj,n_exec)*-999;
branchesTemps = cell(1,n_branch);
resol=0.5; % Resolution in branch volume percentatage 
n_pos=100/resol+1;
for count = 1:n_branch
    branchesTemps{count} = ones(n_pos,n_exec)*-999;
end



%% EXECUTION LOOP

flag_first_exec = true;
for count_exec = 1 : n_exec
    
    % TIME SPAN
    dt = dt_vec(count_exec);
    
    
    % PUMP SPEED
     %If the speed of any pump changes, flows must be recalculated
    
    flag_pump_speed_change = false;
    if  or(any(not((turbo_pump_speed{count_exec} == pump_turbo.pump_speed(:)))), ...
        any(not((volum_pump_speed{count_exec} == pump_volum.pump_speed(:)))))
        flag_pump_speed_change = true;
    end
    
    pump_volum.pump_speed(:) = volum_pump_speed{count_exec};
    pump_turbo.pump_speed(:) = turbo_pump_speed{count_exec};
    
    
    % VALVES
    % If the valve opening varies, vectors for head losses and thus flows must be recalculated.
    % If the valve is closed, flow in its branch is zero.
    
     flag_valve_change = false;
%    if any(not((valve_opening{count_exec} == valve_var.opening(:))))
%       flag_valve_change = true;
%    end
%     
%     valve_var.opening(:) = valve_opening{count_exec};
    if and(heat_exch.inlet_temp(1)>365,valve_var.opening(1)~=1)
        valve_var.opening(1)=1;
        flag_valve_change = true;
    end
    
    
    % THERMOSTATS
    % If the thermostat opening varies because there is a change in the inlet
    % temperature, vectors for head losses and thus flows must be recalculated.
    % Head losses and flows will be recalculated only if the opening variation
    % of any thermostat is larger than the sensitivity of the thermostat[%]
    % If the thermostat is closed, flow in its branch is zero.
    thermostat_opening_new = zeros(size(thermostat, 1), 1);
    for count_thst = 1 : size(thermostat, 1);
%         branch_aux = thermostat.Node_Branch_ref(count_thst); % inlet branch to the thermostat
%         if branches_id{branch_aux}(1) == thermostat.Node(count_thst)
%             pos_aux = 0;
%         else
%             pos_aux = 100;
%         end
        %thermostat_temp = HydroNet_GetObjTemperature(pos_aux, branch_temp_pos{branch_aux}, branch_temperature{branch_aux});
        thermostat_temp = thermostat.Ref_temp(count_thst);
        thermostat_opening_new(count_thst) = 1./(1 + exp(-k(count_thst) * (thermostat_temp - x0(count_thst))));
            %thermostat.to_open(count_thst)=0;
        
        % Checks if the thermostat is openinig faster than it can
        if abs(thermostat_opening_new(count_thst)-thermostat.opening(count_thst)) > dt*thrm_speed
            
            thermostat_opening_new(count_thst) = thermostat_opening_old(count_thst) ...
                + dt * thrm_speed * (thermostat_opening_new(count_thst) - thermostat.opening(count_thst)) ...
                / abs(thermostat_opening_new(count_thst) - thermostat.opening(count_thst));
                % old opening + delta time * max opening per second * sign
        end
    end
  
    flag_thermostat_changes = false;
    for count_thst = 1 : size(thermostat, 1);
        if and(thermostat_opening_new(count_thst) ~= thermostat.opening(count_thst), ...
                thermostat_opening_new(count_thst) > 0.01)
            flag_thermostat_changes = true;
            thermostat.opening(count_thst) = thermostat_opening_new(count_thst); % assign opening
            thermostat_opening_old(count_thst) = thermostat_opening_new(count_thst);
%         end
        else
            thermostat_opening_old(count_thst) = thermostat_opening_new(count_thst);
            thermostat.opening(count_thst) = 0;
        end
    end
    
%     if size(thermostat,1)~=0
%     thermos(1, count_exec)=thermostat.opening(1);
% %     thermos(2, count_exec)=thermostat.opening(2);
%     end
    
%     display(thermostat_temp);
%     display(thermostat.opening(1));
%     display(thermostat_opening_new);
    % MAIN EXECUTION FUNCTION CALL
    [ flows, branches_ind, branches_id, branch_temp_pos, branch_temperature, branch_htx, branch_pump_volum, branch_pump_turbo, ...
        heat_exch, obj_inlet_pos, obj_outlet_pos,pump_volum, pump_turbo, head_loss, hydr_resist1, hydr_resist2 , thermostat] = HydroNet_Executions ( fluid_type, objects, ...
        pipe, valve_fix, valve_var, thermostat, pump_volum, pump_turbo, heat_exch, tank, ...
        nodes_ind, nodes_id, n_nodes, branches_ind, branches_id, branch_cycle, mesh_branches, node_branches, ...
        n_mesh, n_branch, n_tanks, branch_volume, obj_inlet_pos, obj_outlet_pos, branch_htx, branch_pump_volum, branch_pump_turbo, ...
        branch_temp_pos, branch_temperature, gravity_acc, dt, flag_pump_speed_change, flag_valve_change, ...
        flag_thermostat_changes, flag_first_exec, head_loss, hydr_resist1, hydr_resist2, flows,weight_fraction ); 
    
    [temperature_thermostat,temperature_after_engine,temperature_after_radiator] ...
    = register( objects, obj_inlet_pos, obj_outlet_pos,branch_temp_pos, branch_temperature, ...
    temperature_thermostat, temperature_after_engine, temperature_after_radiator, count_exec);



    temp_inlet_cyl(count_exec) = heat_exch.inlet_temp(1);
%     if(size(heat_exch,1)>2)
     temp_inlet_rad(count_exec) = heat_exch.inlet_temp(2);
%     end
%     if and(size(pump_volum,1)==0, size(pump_turbo,1)~=0)
%     	pump_power(count_exec) = pump_turbo.pump_power(1);
%         speed(count_exec) = pump_turbo.pump_speed(1);
%         pump_head(count_exec) = pump_turbo.pump_head(1);
%     else not(and(size(pump_volum,1)==0, size(pump_turbo,1)~=0))
%         pump_power(count_exec) = pump_volum.pump_power(1);
%         speed(count_exec) = pump_volum.pump_speed(1);
%         pump_head(count_exec) = pump_volum.pump_head(1);
%     end
    
    for count_branc =1:n_branch
        flowss(count_branc,count_exec)=flows(count_branc); 
    end 
%     
%     [branchesTempInlet,branchesTempOutlet]=brances_tem_reg(branchesTempInlet,branchesTempOutlet,branch_temp_pos,branch_temperature, ...
%     obj_inlet_pos,obj_outlet_pos,objects,n_branch,branches_id,count_exec);
% 
% 
[branchesTemps] = brances_temps_reg(branchesTemps,branch_temp_pos,branch_temperature, ...
        n_branch,n_pos,count_exec);

    flag_first_exec = false;

%     display(flows);
    display(count_exec);

    
end




%% PLOTTING

display(flows);

if flag_save
    fig = figure('visible','off');
else
    fig = figure;
end
hold on

plot(time,temperature_thermostat,'b','LineWidth',0.7)
% plot(time,temperature_after_engine,'r','LineWidth',0.7)
% plot(time,temperature_after_radiator,'k','LineWidth',0.7)
plot(time,temp_inlet_cyl,'r','LineWidth',0.7)
plot(time,temp_inlet_rad,'k','LineWidth',0.7)
% ttl=strcat('Temperature ');
ttl=strcat('Temperatures ',' - ','Thermostat max speed - ', num2str(thrm_speed));
title_txt = title (ttl);
% legend('Cylinder inlet', )
% L = legend('Cylinder inlet');
L = legend('Thermostat','Cylinder inlet','Radiator inlet');
% L = legend('Thermostat','Engine inlet','Radiator inlet');
L.FontSize = 18;
L.Location ='southeast';
x_label = xlabel('Time [s]'); 
y_label = ylabel('Temperature [K]'); 
% axis([0 1000 280 inf]);
% axis([0 1000 345 375]);
if flag_save
    title_txt.FontSize = 26;
    x_label.FontSize = 22;
    y_label.FontSize = 22;
    set(gca,'fontsize',18); % numbers
    set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Temperatures2.png');
    name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\barrido Thermostat Speedn\','Thermostat_max_speed_', num2str(thrm_speed),'_',circuit_file(1 : end-4),'_Temperatures.png');
    print(name_save,'-dpng','-r300');
    fprintf('File saved: %s.png\n', name_save);
end


% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% hold on
% 
% plot(time,transpose(flowss))
% 
% title_txt = title (strcat('Thermostat Speed - ',num2str(thrm_speed),' - ' ,'Flows'));
% legend('Pump','Cylinder','Radiator')
% legend_=[]
% for count_branch=1:n_branch
%         branche_name = strcat('Branch',' ',num2str(count_branch));
%         legend_=[legend_; branche_name];
%         
% end
% legend(legend_)
% legend('Pump','Cylinder 1','Cylinder 2','Cylinder 3','Cylinder 4','Cylinder 5','Cylinder 6','Radiator 1','Radiator 2','Radiator 3')
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Flow [m^{3}/s]'); 
% bar(flows)
% % axis([0 1000 300 380])
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Flows.png');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end


% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% hold on
% title_txt = title ('Power');
% plot(time,pump_power,'LineWidth',1)
% legend('Power')
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Power [W]');
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Power');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end
% 
% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% hold on
% title_txt = title ('Speed');
% plot(time,speed,'k')
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Pump speed [rad/s]');
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Speed');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end
% 
% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% hold on
% title_txt = title ('Pump Head');
% plot(time,pump_head,'LineWidth',1)
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Pump head [m.c.f]'); 
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Pump head');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end
% 
% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% % title_txt = title (strcat('Thermostats opening'));
% title_txt = title (strcat('Thermostats opening',' - ','Thermostat max speed - ', num2str(thrm_speed)));
% hold on 
% plot(time,transpose(thermos),'LineWidth',1)
% legend('Rad','Cyl')
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Opening'); 
% axis([0 1000 0 1.01]);
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
% %     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),' Openings.png');
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\barrido Thermostat Speedn\','Thermostat_max_speed_', num2str(thrm_speed),'_',circuit_file(1 : end-4),'_Openings.png');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end

% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% title_txt = title (strcat('Thermostat Speed - ',num2str(thrm_speed),' - ' ,'Objects inlet temperature'));
% hold on 
% for object=1:n_obj
%     if not(strcmp(objects.class(object),'node'))
%         plot(time,branchesTempInlet(object,:),'LineWidth',1)
%     end
% end
% hold on
% 
% legend_=[]
% for object=1:n_obj
%     if not(strcmp(objects.class(object),'node'))
%         legend_=[legend_ objects.name(object)];
%         hold on;
%     end
% end
% legend(legend_) 
% % legend(legend_)
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Temperature [K]'); 
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\', ...
%         'Thermostat Speed',num2str(thrm_speed),circuit_file(1 : end-4),' inlets_temps.png');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end  
% 
% if flag_save
%     fig = figure('visible','off');
% else
%     fig = figure;
% end
% title_txt = title (strcat('Branches temperature'));
% hold on 
% [n_obj,~]=size(objects)
% for object=1:n_obj
%     if not(strcmp(objects.class(object),'node'))
%         plot(time,branchesTempOutlet(object,:),'LineWidth',1)
%         hold on
%     end
% end
% legend_=[]
% for object=1:n_obj
%     if not(strcmp(objects.class(object),'node'))
%         legend_=[legend_; objects.name(object)];
%         hold on;
%     end
% end
% legend(legend_)
% hold on;
% x_label = xlabel('Time [s]'); 
% y_label = ylabel('Temperature [K]'); 
% axis([0 1000 300 inf]);
% if flag_save
%     title_txt.FontSize = 26;
%     x_label.FontSize = 22;
%     y_label.FontSize = 22;
%     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\', ...
%        circuit_file(1 : end-4),' branches_temps.png');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end
% 
pos_axis=0:resol:100;
% % pos_axis = 100 - pos_axis;

% for branch_count = 1:n_branch
%     if flag_save
%     fig = figure('visible','off');
%     else
%         fig = figure;
%     end
%     title_txt = title (strcat('Branch ',num2str(branch_count),' - ' ,'Branch temperature'));
%     hold on 
%     imagesc(time,pos_axis,branchesTemps{branch_count});
%     colormap(jet);
%     x_label = xlabel('Time [s]'); 
%     y_label = ylabel('Branch Volume [%]'); 
%     axis([0 1000 0 100])
%     colorbar;
%     if flag_save
%         title_txt.FontSize = 26;
%         x_label.FontSize = 22;
%         y_label.FontSize = 22;
%         set(gca,'fontsize',18); % numbers
%         set(gcf, 'PaperPosition', [0 0 42 21]);
%         name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\', ...
%             circuit_file(1 : end-4),'Branch-',num2str(branch_count),' Branch temperature.png');
%         print(name_save,'-dpng','-r300');
%         fprintf('File saved: %s.png\n', name_save);
%     end
% end
% 
% 
% if flag_save
% fig = figure('visible','off');
% else
%     fig = figure;
% end
% % % %     title_txt = title (strcat('Branch ',num2str(branch_count),' Thermostat Speed - ',num2str(thrm_speed),' - ' ,'Branch temperature at the end'));
% hold on 
% for branch_count=1:n_branch
% % for branch_count=[1,2,3]
%     s(branch_count) = subplot(n_branch,1,branch_count);
%  %    s(branch_count) = subplot(3,1,find(branch_count==[1,2,3]));
%      plot(pos_axis,branchesTemps{branch_count}(:,7000),'LineWidth',1)
% end
% for branch_count=1:n_branch
% % for branch_count=[1,2,3]
%    title_txt = title (s(branch_count),strcat('Branch ',num2str(branch_count)));
%    x_label = xlabel(s(branch_count),'Branch Volume [%]'); 
%     x_label = xlabel(s(branch_count),'Branch Volume [%]'); 
%     y_label = ylabel(s(branch_count),'Temperature [K]'); 
% %     axis(s(branch_count),[0 100 372 373])
% %     set(s(branch_count),'ytick',linspace(372,373,5));
%     if flag_save
%         title_txt.FontSize = 22;
%         if branch_count~=3
%             axis(s(branch_count) ,[0 100 366 369])
%             set(s(branch_count),'ytick',linspace(366,369,4));
%         else
%             axis(s(branch_count) ,[0 100 295 305])
%             set(s(branch_count),'ytick',linspace(295,305,5));
%         end
%         x_label.FontSize = 18;
%         y_label.FontSize = 18;
%         set(s(branch_count),'fontsize',18); % numbers
%     end
% end
% 
% 
% if flag_save
% %     title_txt.FontSize = 26;
% %     
% %         x_label.FontSize = 16;
% %         y_label.FontSize = 16;
% %     
% %     set(gca,'fontsize',18); % numbers
%     set(gcf, 'PaperPosition', [0 0 42 21]);
%     name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\', ...
%         circuit_file(1 : end-4),' Branches temperature at 7500 seconds.png');
%     print(name_save,'-dpng','-r300');
%     fprintf('File saved: %s.png\n', name_save);
% end

if flag_save
    fig = figure('visible','off');
else
    fig = figure;
end
hold on
s(1)=subplot(2,1,1);
[y_axis_1,line_1_1,line_1_2]=plotyy(time,transpose(flowss),time,pump_head);
s(2)=subplot(2,1,2);
[y_axis_2,line_2_1,line_2_2]=plotyy(time,speed,time,pump_power);
title_txt_1 = title (s(1),'Flows - Pump Head');
title_txt_2 = title (s(2),'Speed - Pump Power');
x_label_1 = xlabel(s(1),'Time [s]'); 
y_label_1_1 = ylabel(y_axis_1(1),'Flow [m^{3}/s]'); 
y_label_1_2 = ylabel(y_axis_1(2),'Pump head [m.c.f]'); 
x_label_2 = xlabel(s(2),'Time [s]'); 
y_label_2_1 = ylabel(y_axis_2(1),'Speed [rad/s]'); 
y_label_2_2 = ylabel(y_axis_2(2),'Pump power [w]'); 
set(line_1_1,'LineWidth',1);
set(line_1_2,'LineWidth',1);
set(line_2_1,'LineWidth',1);
set(line_2_2,'LineWidth',1);
legend_=cell(1,n_branch+1);
 for branch=1:n_branch
            legend_{branch}= strcat('Branch', num2str(branch));
            hold on;
 end
legend_{n_branch+1}=  'Pump Head';
L1 = legend(s(1),legend_);
L1.FontSize = 14;
L1.Location ='best';
L2 = legend(s(2),'Pump speed', 'Pump power');
L2.FontSize = 14;
L2.Location ='best';
line_1_2.LineStyle='--';
line_2_2.LineStyle='--';
% axis(y_axis_1(1) , [0 1000 0 0.0015]);
% axis(y_axis_1(2) , [0 1000 0 4]);
% axis(y_axis_2(1) , [0 1000 0 150]);
% axis(y_axis_2(2) , [0 1000 0 35]);
% set(y_axis_1(1),'ytick',linspace(0,0.0015,4));
% set(y_axis_1(2),'ytick',linspace(0,4,5));
% set(y_axis_2(1),'ytick',[0,50,100,150]);
% set(y_axis_2(2),'ytick',linspace(0,35,8));
if flag_save
    title_txt_1.FontSize = 26;
    title_txt_2.FontSize = 26;
    x_label_1.FontSize = 22;
    x_label_2.FontSize = 22;
    y_label_1_1.FontSize = 22;
    y_label_1_2.FontSize = 22;
    y_label_2_1.FontSize = 22;
    y_label_2_2.FontSize = 22;
    set([y_axis_1(1) y_axis_1(2) y_axis_2(1) y_axis_2(2)],'fontsize',18); % numbers
    set(gcf, 'PaperPosition', [0 0 42 21]);
    name_save = strcat('C:\Users\Jaime\Documents\& PFC-CMT\Program\HydroNet\Graph\',circuit_file(1 : end-4),'_Flows-Pump head_Speed-Power');
    print(name_save,'-dpng','-r300');
    fprintf('File saved: %s.png\n', name_save);
end
% display (flows);
end
