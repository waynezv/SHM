function Trans_Dist_Mat2 = MSF_Data_Pert(Dist_Mat, tree_num, pert_strength)
for i = 1:tree_num
    if i == 1
        Trans_Dist_Mat2 = MST_Data(Dist_Mat);
    else
        Trans_Dist_Mat2 = squeeze( max( cat(3, Trans_Dist_Mat2, MST_Data_Per(Dist_Mat, pert_strength)) , [], 3) );
    end
end