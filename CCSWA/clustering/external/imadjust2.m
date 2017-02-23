function img_out = imadjust2(img_in, map)
img_out = img_in;
for i=1:length(map)
    idx = img_in==i;
    img_out(idx)=255*map(i);
end