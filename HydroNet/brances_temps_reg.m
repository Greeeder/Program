function [branchesTemps]=brances_temps_reg(branchesTemps,branch_temp_pos,branch_temperature, ...
    n_branch,n_pos,count_exec)
    
    for branch_count = 1:n_branch
        for pos_count = 1:n_pos
            if or(branch_count~=3,count_exec>=1883)
            position = pos_count*100/n_pos;
            else
            position =100 - pos_count*100/n_pos;
            end
            branchesTemps{branch_count}(pos_count,count_exec) = HydroNet_GetPosTemperature(position, branch_temp_pos{branch_count}, branch_temperature{branch_count}) ;
            
        end
    end
end