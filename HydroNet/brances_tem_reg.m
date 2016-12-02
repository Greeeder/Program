function [branchesTempInlet,branchesTempOutlet]=brances_tem_reg(branchesTempInlet,branchesTempOutlet,branch_temp_pos,branch_temperature, ...
    obj_inlet_pos,obj_outlet_pos,objects,n_branch,branches_id,count_exec)

for count_branc=1:n_branch
    for obj= branches_id{count_branc}(2:1:end-1)
        if obj ~=objects.id(obj)
            obj=find(objects.id==obj);
        end
        [pos, branch] = HydroNet_GetPosition(obj, objects, obj_inlet_pos, obj_outlet_pos, 'start');
        if pos<0
            pos=0;
        end
        branchesTempInlet(obj,count_exec) = HydroNet_GetPosTemperature(pos, branch_temp_pos{branch}, branch_temperature{branch}) ;
        [pos, branch] = HydroNet_GetPosition(obj, objects, obj_inlet_pos, obj_outlet_pos, 'end');
        if pos>100
            pos=100;
        end
        branchesTempOutlet(obj,count_exec) = HydroNet_GetObjTemperature(pos, branch_temp_pos{branch}, branch_temperature{branch}) ;
    end
end