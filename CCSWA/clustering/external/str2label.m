function label = str2label(str_label)
str_label_unique = unique(str_label);
N = length(str_label);
label = zeros(N,1);
NL = length(str_label_unique);
for i = 1:N
    i
    label(i,1) = find(strcmp(str_label{i},str_label_unique));
end