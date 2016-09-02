function [ positions, temperatures ] = HydroNet_InsertVolume( positions, temperatures, ...
    start_pos, end_pos, temp_action, fluid_type, value, volume )
% Obtains range of positions where a volume is to be placed inside a branch and replaces
% current volumes by the new volume. In addition, a temperature for the volume is calculated.
% Fluid_type is only mandatory for temp_action = 0 and 1 (must calculate outlet temperature).
% Value is only mandatory for temp_action = 1 and 2 (not just averaging).
% Volume is only mandatory for temp_action = 1 (impose heat).

if start_pos == end_pos
    return % to avoid dividing by zero
end

% Finds start position inside positions vector
ind_aux = find(start_pos == positions, 1); % checks if the position already exists
flag_start_exists = false;
if isempty(ind_aux) % position does not exist in the vector yet
    count_pos = 1; % position counter
    while positions(count_pos) <= start_pos
            % position is lower or equal and counter is inside the vector
        count_pos = count_pos + 1;
        if count_pos > numel(positions) % counter is outside the vector
            break % to avoid exceeding the dimension of "positions"
        end
    end
    start_pos_ind = count_pos; % index of the start position in positions
else % position already exists in the vector
    start_pos_ind = ind_aux; % index of the start position in positions
    flag_start_exists = true;
    count_pos = start_pos_ind;
end


% Finds end position inside positions vector
ind_aux = find(end_pos == positions, 1); % checks if the position already exists
flag_end_exists = false;
if isempty(ind_aux) % position does not exist in the vector yet
    while positions(count_pos) <= end_pos
            % position is lower or equal and counter is inside the vector
        count_pos = count_pos + 1;
        if count_pos > numel(positions) % counter is outside the vector
            break % to avoid exceeding the dimension of "positions"
        end
    end
    end_pos_ind = count_pos; % index of the end position in positions
else % position already exists in the vector
    end_pos_ind = ind_aux; % index of the end position in positions
    flag_end_exists = true;
end



% Obtains volume temperature

if temp_action == 2 % impose temperature given by 'value'
    
    temp_to_insert = value;
    
else % 0 or 1: calculates mean temperature with an enthalpy balance
    % Branch volume will not be used now because it is not mandatory since
    % it is a constant that multiplies and divides everything
    
    % Temperatures at both ends of the volume, multiplied by density
    % and by the percentage range that they occupy inside the volume
    % Those volumes are taken incomplete because they are cut by the
    % new bounds start_pos and end_pos
    
    if flag_start_exists % start_pos already exists in the vector
        mass_start = 0;
        temp_mass_start = 0;
    elseif start_pos_ind == end_pos_ind % both are inside the same volume initially
        mass_start = density(fluid_type, temperatures(start_pos_ind - 1)) * (end_pos - start_pos);
            % Not a real mass. Volume is not included because it appears both
            % in the numerator and the denominator multiplying all terms. 
        temp_mass_start = temperatures(start_pos_ind - 1) * mass_start;
    else
        mass_start = density(fluid_type, temperatures(start_pos_ind - 1)) * (positions(start_pos_ind) - start_pos);
            % Not a real mass. Volume is not included because it appears both
            % in the numerator and the denominator multiplying all terms. 
        temp_mass_start = temperatures(start_pos_ind - 1) * mass_start;
    end
    
    if flag_end_exists % end_pos already exists in the vector
        mass_end = 0;
        temp_mass_end = 0;
    elseif start_pos_ind == end_pos_ind % both are inside the same volume initially
         mass_end = 0;
         temp_mass_end = 0;
    else
        mass_end = density(fluid_type, temperatures(end_pos_ind - 1)) * (end_pos - positions(end_pos_ind - 1));
        temp_mass_end = temperatures(end_pos_ind - 1) * mass_end;
    end
    
    mass_sum = mass_start + mass_end; % denominator
    temp_mass_sum = temp_mass_start + temp_mass_end; % numerator
    
    
    % Volumes completely inside new volume
    
    if flag_end_exists
        last_element = (end_pos_ind - 1);
    else
        last_element = (end_pos_ind - 2);
    end

    for count_pos = start_pos_ind : last_element
        mass_sum = mass_sum + density(fluid_type, temperatures(count_pos)) * (positions(count_pos + 1) - positions(count_pos));
        temp_mass_sum = temp_mass_sum + temperatures(count_pos) * density(fluid_type, temperatures(count_pos)) * ...
            (positions(count_pos + 1) - positions(count_pos));
    end
    
    
    % Enter averaged temperature
    temp_to_insert = temp_mass_sum / mass_sum; 
    
    if temp_action == 1 % calculate temperature imposing exchanged heat
        
        mass_sum = mass_sum * volume; % volume needed to calculate real mass
        temp_to_insert = temp_to_insert + value / mass_sum / cp(fluid_type, temp_to_insert);    
    end
end



% Inserts calculated temperature and removes intermediate temperatures in the branch

if start_pos == positions(1) % start position is the first position of the vector
    lower_temp_group = [];
else % start position is not the first one of the vector
    lower_temp_group = temperatures(1 : start_pos_ind - 1);
end

if end_pos == positions(end) % end position is the last position of the vector
    upper_temp_group = [];
elseif flag_end_exists % end position is an existing position of the vector
    upper_temp_group = temperatures(end_pos_ind : end);
else % end position is a new position of the vector
    upper_temp_group = temperatures(end_pos_ind - 1 : end);
end

temperatures = [lower_temp_group, temp_to_insert, upper_temp_group];


% Inserts specified positions and removes intermediate positions in the branch

if start_pos == positions(1) % start position is the first position of the vector
    lower_pos_group = [];
else % start position is not the first one of the vector
    lower_pos_group = positions(1 : start_pos_ind - 1);
end

if end_pos == positions(end) % end position is the last position of the vector
    upper_temp_group = [];
elseif flag_end_exists % end position is an existing position of the vector
    upper_temp_group = positions(end_pos_ind + 1 : end);
else % end position is not equal to the last one of the vector (but the index can be the last)
    % and is not present in the original vector
    upper_temp_group = positions(end_pos_ind : end);
end

positions = [lower_pos_group, start_pos, end_pos, upper_temp_group];

end
