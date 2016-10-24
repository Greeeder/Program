function [ T_out_1, T_out_2, heat ] = Heat_Exchanger( fluid_1, mass_flow_1, T_in_1, fluid_2, mass_flow_2, T_in_2, heat_exch_props,weight_fraction_1, weight_fraction_2 )
%( fluid_1, flow_1, T_in_1, fluid_2, flow_2, T_in_2, heat_exch_props )
% Inputs: type of fluid (air, water, oil...), mass flow and inlet
%   temperatures of both flows; table with heat exchanger properties.
% Outputs: outlet temperatures of 1 and 2 and heat from 1 to 2
% Calculates efficiency from the UA value of the heat exchanger and from
%   the flows' characteristics
[~,Cp_1]=fluidproperties(fluid_1,T_in_1,weight_fraction_1);
[~,Cp_2]=fluidproperties(fluid_2,T_in_2,weight_fraction_2);

mCp_1 = mass_flow_1 * Cp_1;
mCp_2 = mass_flow_2 * Cp_2;
% mCp_1 = density(fluid_1, T_in_1) * flow_1 * cp(fluid_1, T_in_1);
% mCp_2 = density(fluid_2, T_in_2) * flow_2 * cp(fluid_2, T_in_2);
mCp_min = min(mCp_1, mCp_2);
mCp_max = max(mCp_1, mCp_2);

NTU = heat_exch_props.UA / mCp_min;

C = mCp_min / mCp_max;


% Efficiency is found in a matrix stored at heat_exch_props.efficiency(NTU, C)
% NTU values are stored in heat_exch_props.eff_NTU
% C values are stored in heat_exch_props.eff_C

% If NTU or C are outside the table, assign extreme values
if NTU > heat_exch_props.eff_NTU(end)
    NTU =  heat_exch_props.eff_NTU(end);
elseif NTU < heat_exch_props.eff_NTU(1)
    NTU = heat_exch_props.eff_NTU(1);
end

if C > heat_exch_props.eff_C(end)
    C =  heat_exch_props.eff_C(end);
elseif C < heat_exch_props.eff_C(1)
    C = heat_exch_props.eff_C(1);
end

efficiency = interp2(heat_exch_props.eff_NTU, heat_exch_props.eff_C, heat_exch_props.efficiency, NTU, C);

%%% Alternative way
% [NTU_diff, NTU_ind] = min(NTU - heat_exch_props.eff_NTU); % get index of closer NTU and difference
% 
% if and(NTU_diff < 0, NTU_ind == 1) % searched value is outside of the table (lower)
%     % NTU minimum
%     NTU_range = [heat_exch_props.eff_NTU(NTU_ind) heat_exch_props.eff_NTU(NTU_ind + 1)];
% elseif and(NTU_diff >= 0, NTU_ind == numel(heat_exch_props.eff_NTU)) % searched value is outside of the table (higher)
%     % NTU maximum
%     NTU_range = [heat_exch_props.eff_NTU(NTU_ind - 1) heat_exch_props.eff_NTU(NTU_ind1)];
% else  % searched value is inside table
%     if NTU_diff >= 0
%         % get index of closer NTU and next one (higher)
%         NTU_range = [heat_exch_props.eff_NTU(NTU_ind) heat_exch_props.eff_NTU(NTU_ind + 1)];
%     else
%         % get index of closer NTU and previous one (lower)
%         NTU_range = [heat_exch_props.eff_NTU(NTU_ind - 1) heat_exch_props.eff_NTU(NTU_ind)];
%     end
% end
% 
% [C_diff, C_ind] = min(C - heat_exch_props.eff_C); % get index of closer C and difference
% 
% if and(C_diff < 0, C_ind == 1) % searched value is outside of the table (lower)
%     % C minimum
%     C_range = [heat_exch_props.eff_C(C_ind) heat_exch_props.eff_C(C_ind + 1)];
% elseif and(C_diff >= 0, C_ind == numel(heat_exch_props.eff_C)) % searched value is outside of the table (higher)
%     % C maximum
%     C_range = [heat_exch_props.eff_C(C_ind - 1) heat_exch_props.eff_C(C_ind1)];
% else  % searched value is inside table
%     if C_diff >= 0
%         % get index of closer C and next one (higher)
%         C_range = [heat_exch_props.eff_C(C_ind) heat_exch_props.eff_C(C_ind + 1)];
%     else
%         % get index of closer C and previous one (lower)
%         C_range = [heat_exch_props.eff_C(C_ind - 1) heat_exch_props.eff_C(C_ind)];
%     end
% end
% 
% efficiency_range = heat_exch_props.efficiency(NTU_range, C_range);
% 
% % Linear interpolation
% efficiency = interp2(NTU_range, C_range, efficiency_range, NTU, C);


%%% Alternative: Efficiency given by a correlation
% efficiency = heat_exch_props.coef_NTU0_Cr0 + heat_exch_props.coef_NTU1_Cr0 * NTU + heat_exch_props.coef_NTU0_Cr1 * Cr + ...
% 	heat_exch_props.coef_NTU1_Cr1 * NTU * Cr + heat_exch_props.coef_NTU2_Cr1 * NTU^2 * Cr + heat_exch_props.coef_NTU1_Cr2 * NTU * Cr^2 + ...
% 	heat_exch_props.coef_NTU2_Cr2 * NTU^2 * Cr^2;
%
% if efficiency > 1
% 	efficiency = 1;
% elseif efficiency < 1
% 	efficiency = 0;
% end



heat = efficiency * mCp_min * (T_in_1 - T_in_2); % (+) heat from 1 to 2; (-) heat from 2 to 1


T_out_1 = T_in_1 - heat / mCp_1;
T_out_2 = T_in_2 + heat / mCp_2;

end

