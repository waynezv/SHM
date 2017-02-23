function label2 = eight_connectivity(label1)
[height, width] = size(label1);
label2 = zeros(height, width);
sup_num1 = max(max(label1));
count = 0;
for i = 1:sup_num1
    flag = zeros(height, width);
    idx1 = label1 == i;
    flag(idx1) = 1;
    [label_tmp num] = bwlabel(flag, 8);
    for j = 1:num
        idx2 = find(label_tmp==j);
        if flag(idx2(1)) == 1
            count = count+1;
            label2(idx2) = count;
        end
    end
end