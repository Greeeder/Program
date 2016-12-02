function [ head_loss, hydr_resist1, hydr_resist2 ] = HydroNet_HeadLoss( objects, pipe, valve_fix, ...
    valve_var, thermostat, pump_turbo, heat_exch, tank, n_branch)
% Calculates head loss in every branch of an hydraulic circuit


head_loss = zeros(n_branch, 1); % head loss (m.c.f.) - independent of flow
hydr_resist1 = zeros(n_branch, 1); % hydraulic resistance (m.c.f. / m^3 * s)
    % proportional to flow - for turbopumps
hydr_resist2 = zeros(n_branch, 1); % hydraulic resistance (m.c.f. / m^6 * s^2)
    % proportional to flow^2

    
% Check out every object that generates head loss, find out its branch and add
% its head loss to the branch's total head loss
    
for count_obj = 1 : size(pipe, 1)
    % pipe: hydraulic resistance depends on characteristics
    branch_aux = objects.branch{pipe.obj_index(count_obj)};
	hydr_resist2(branch_aux) = hydr_resist2(branch_aux) + 8 * pipe.friction_coef(count_obj) ...
        * pipe.length(count_obj) / (pi^2) / 9.81 / (pipe.diameter(count_obj)^5);
end

for count_obj = 1 : size(valve_fix, 1)
    % valve with fixed head loss
    branch_aux = objects.branch{valve_fix.obj_index(count_obj)};
    head_loss(branch_aux) = head_loss(branch_aux) + valve_fix.head_loss(count_obj);
end

for count_obj = 1 : size(valve_var, 1)
    % valve with head loss dependent on opening
    branch_aux = objects.branch{valve_var.obj_index(count_obj)};
    hydr_resist2(branch_aux) = hydr_resist2(branch_aux) + valve_var.coef0(count_obj) ...
        + valve_var.coef1(count_obj) * valve_var.opening(count_obj) ...
        + valve_var.coef2(count_obj) * valve_var.opening(count_obj)^2 ...
        + valve_var.coef3(count_obj) * valve_var.opening(count_obj)^3;
end

for count_obj = 1 : size(thermostat, 1)
    % valve with head loss dependant on opening dependant on temperature
    branch_aux = objects.branch{thermostat.obj_index(count_obj)};
    opening = thermostat.opening(count_obj);
     % To avoid dividing by 0
    if opening == 0
        opening=1E-16;
    end
    hydr_resist2(branch_aux) = hydr_resist2(branch_aux) + thermostat.h_coef0(count_obj) ...
        + thermostat.h_coef1(count_obj) / opening ...  % OJOO! change curve
        + thermostat.h_coef2(count_obj) / opening^2 ...
        + thermostat.h_coef3(count_obj) / opening^3;
end

for count_obj = 1 : size(pump_turbo, 1)
    % turbopump: positive head -> negative loss or resistance
    branch_aux = objects.branch{pump_turbo.obj_index(count_obj)};
    speed = pump_turbo.pump_speed(count_obj);
    head_loss(branch_aux) = head_loss(branch_aux) - (pump_turbo.coef_N0(count_obj) + ...
        pump_turbo.coef_N1(count_obj) * speed + pump_turbo.coef_N2(count_obj) * speed^2);
    hydr_resist1(branch_aux) = hydr_resist1(branch_aux) - (pump_turbo.Q_coef_N0(count_obj) + pump_turbo.Q_coef_N1(count_obj) * ...
        speed + pump_turbo.Q_coef_N2(count_obj) * speed^2);
    hydr_resist2(branch_aux) = hydr_resist2(branch_aux) - (pump_turbo.Q2_coef_N0(count_obj) + pump_turbo.Q2_coef_N1(count_obj) * ...
        speed + pump_turbo.Q2_coef_N2(count_obj) * speed^2);
end

for count_obj = 1 : size(heat_exch, 1)
    % heat exchanger: hydraulic resistance given
    branch_aux = objects.branch{heat_exch.obj_index(count_obj)};
    hydr_resist2(branch_aux) = hydr_resist2(branch_aux) + heat_exch.hydr_resist(count_obj);
end

% for count_obj = 1 : size(tank, 1) %%%% TO DO!
% end

end

