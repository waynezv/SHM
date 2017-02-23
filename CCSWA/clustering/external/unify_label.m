function label_out = unify_label(label_in)
N = length(label_in);
label_out = zeros(N, 1);
label_in_unique = unique(label_in);
N_label = length(label_in_unique);
for i = 1:N_label
    idx = label_in == label_in_unique(i);
    label_out(idx) = i;
end