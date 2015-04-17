function distance=dist(x,y)
z=repmat(y',size(x,1),1);
distance=sqrt(sum((x-z).^2,2));