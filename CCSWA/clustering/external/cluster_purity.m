function score = cluster_purity(label_gt, label_in)
N = length(label_in);
label_in2 = zeros(N,1);
label_in_unique = unique(label_in);
num_label_in = length(label_in_unique);
for i = 1:num_label_in
    idx = label_in == label_in_unique(i);
    label_in2(idx) = mode(label_gt(idx));
end
score = sum(label_gt==label_in2)/N;