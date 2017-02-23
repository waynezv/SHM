%David Fouhey
%CV Fall 2012 - Provided Code
%Code to get a reasonable filter bank

function [filterBank] = createFilterBank()
    filterBank = {};
    %3 scales
    for scale=1:1:3
        scaleMultiply = sqrt(2)^scale;
        %Some arbitrary scales
        gaussianSigmas = [1, 2, 4];
        logSigmas = [1, 2, 4, 8];
        dGaussianSigmas = [2, 4];
        %Gaussians
        for s=gaussianSigmas
            filterBank = [filterBank, getGaussianFilter(s*scaleMultiply)];
        end
        %LoG
        for s=logSigmas
            filterBank = [filterBank, getLOGFilter(s*scaleMultiply)];
        end
        %d/dx, d/dy Gaussians
        for s=dGaussianSigmas
            filterBank = [filterBank, filterDerivativeX(getGaussianFilter(s*scaleMultiply))];
            filterBank = [filterBank, filterDerivativeY(getGaussianFilter(s*scaleMultiply))];
        end
    end
end

function h = getGaussianFilter(sigma)
    h = fspecial('gaussian',ceil(sigma*3*2+1),sigma);
end

function h = getLOGFilter(sigma)
    h = fspecial('log',ceil(sigma*3*2+1),sigma);
end

function hD = filterDerivativeX(h)
    ddx = [-1, 0, 1];
    hD = imfilter(h, ddx); 
end

function hD = filterDerivativeY(h)
    ddy = [-1, 0, 1]';
    hD = imfilter(h, ddy);
end
