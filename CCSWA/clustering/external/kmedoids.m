function [label, energy, index] = kmedoids(D,k)
n = size(D,1);
[~, label] = min(D(randsample(n,k),:),[],1);
last = 0;
while any(label ~= last)
    [~, index] = min(D*sparse(1:n,label,1,n,k,n),[],1);
    last = label;
    [val, label] = min(D(index,:),[],1);
end
energy = sum(val);
