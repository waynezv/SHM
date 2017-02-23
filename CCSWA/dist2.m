function distance=dist2(x,y)
z=y';
distance=sqrt(sum((x-z).^2,2));