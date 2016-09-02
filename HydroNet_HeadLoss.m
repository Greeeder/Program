function [ head_loss, hydr_resist1, hydr_resist2 ] = HydroNet_HeadLoss( objects, pipe, valve_fix, ...
    valve_var, thermostat, pump_turbo, heat_exch_fix, heat_exch_Tout, tank, branches_ind, ...
    branch_cycle, n_branch, temperature )

% Calculates head loss in every branch of an hydraulic circuit


head_loss = zeros(n_branch, 1); % head loss (m.c.f.) - independent of flow
hydr_resist1 = zeros(n_branch, 1); % hydraulic resistance (m.c.f. / m^3 * s)
    % proportional to flow - for turbopumps
hydr_resist2 = zeros(n_branch, 1); % hydraulic resistance (m.c.f. / m^6 * s^2)
    % proportional to flow^2

    
% Branch loop
for count = transpose(find(branch_cycle)) % goes only to branches that are part of a cycle
%for count = n_branch % analyzes all branches
    for count1 = branches_ind{count} % objects in branch (indices)
        id_aux = objects.id(count1); % id of current object
        class_aux = objects.class{count1}; % class of current object
        
        % I can't use a switch because of the heat exchangers
        if strcmp(class_aux, 'pipe') % pipe: hydraulic resistance depends on characteristics
            hydr_resist2(count) = hydr_resist2(count) + 8 * pipe.friction_coef(pipe.id == id_aux) ...
                * pipe.length(pipe.id == id_aux) / (pi^2) / 9.81 / (pipe.diameter(pipe.id == id_aux)^5);
        
        elseif strcmp(class_aux, 'valve_fix') % valve with fixed head loss
            head_loss(count) = head_loss(count) + valve_fix.head_loss(valve_fix.id == id_aux);
        
        elseif strcmp(class_aux, 'valve_var') % valve with head loss dependant on opening
            ind_aux = find(valve_var.id == id_aux); % index of current object in class table
            head_loss(count) = head_loss(count) + valve_var.coef0(ind_aux) ...
                + valve_var.coef1(ind_aux) * valve_var.opening(ind_aux) ...
                + valve_var.coef2(ind_aux) * valve_var.opening(ind_aux)^2 ...
                + valve_var.coef3(ind_aux) * valve_var.opening(ind_aux)^3;
        
        elseif strcmp(class_aux, 'thermostat')
            % valve with head loss dependant on opening dependant on temperature
            ind_aux = find(thermostat.id == id_aux);
            opening = thermostat.T_coef0(ind_aux) ...
                + thermostat.T_coef1(ind_aux) * temperature ...
                + thermostat.T_coef2(ind_aux) * temperature^2 ...
                + thermostat.T_coef3(ind_aux) * temperature^3;
            head_loss(count) = head_loss(count) + thermostat.h_coef0(ind_aux) ...
                + thermostat.h_coef1(ind_aux) * opening ...
                + thermostat.h_coef2(ind_aux) * opening^2 ...
                + thermostat.h_coef3(ind_aux) * opening^3;
        
        elseif strcmp(class_aux, 'pump_turbo')
            % turbopump: positive head -> negative loss or resistance
            ind_aux = find(pump_turbo.id == id_aux);
            head_loss(count) = head_loss(count) - (pump_turbo.coef_N0(ind_aux) + pump_turbo.coef_N1(ind_aux) * ...
                pump_turbo.pump_speed + pump_turbo.coef_N2(ind_aux) * pump_turbo.pump_speed^2);
            hydr_resist1(count) = hydr_resist1(count) - (pump_turbo.Q_coef_N0(ind_aux) + pump_turbo.Q_coef_N1(ind_aux) * ...
                pump_turbo.pump_speed + pump_turbo.Q_coef_N2(ind_aux) * pump_turbo.pump_speed^2);
            hydr_resist2(count) = hydr_resist2(count) - (pump_turbo.Q2_coef_N0(ind_aux) + pump_turbo.Q2_coef_N1(ind_aux) * ...
                pump_turbo.pump_speed + pump_turbo.Q2_coef_N2(ind_aux) * pump_turbo.pump_speed^2);
        
        elseif strncmp(class_aux, 'heat_exch', 9) % class starts with 'heat_exch'
            % heat exchanger: hydraulic resistance given
            ind_aux = find(eval(strcat(class_aux, '.id')) == id_aux); % index of current object in class table
            hydr_resist2(count) = hydr_resist2(count) + eval(strcat(class_aux, '.hydr_resist(', num2str(ind_aux), ')'));

%         elseif strcmp(class_aux, 'tank')
%             % tank: %%%%%%%%%%%%%%%%%% TO DO!

        end
    end
end


end

