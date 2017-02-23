function IDX = spectral_clust2(W, k1, k2)
tic; [V,~] = evecs(W,k1); ttt = toc;
disp(['evecs took ' num2str(ttt) ' seconds']);
%% normalize rows of V
count = 0;
for m=1:size(V,1)
    if norm(V(m,:))==0
        count=count+1;
        count
    end
    V(m,:)=(V(m,:)+eps)/(norm(V(m,:)+eps));
end
%% call kmeans
tic
[IDX,~] = kmeans(V,k2,'distance','cosine','emptyaction','singleton');
ttt = toc;
disp(['kmeans took ' num2str(ttt) ' seconds']);