function[img_bound] = disp_bound(label, img)
[height width] = size(label);
img_bound = img;
for rr = 1 : height
    for cc = 1 : width
        for r = -1 : 1
            for c = -1 : 1
                if rr + r >= 1 && rr + r <= height && cc + c >= 1 && cc + c <= width && (r~=0 || c~=0)
                    if label(rr, cc) ~= label(rr + r, cc + c)
                        img_bound(rr, cc, :) = 0;
                    end
                end
            end
        end
    end
end