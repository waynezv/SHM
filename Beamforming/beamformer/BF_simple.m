% clear all
% close all
% dbstop if error

%% global parameters that may be used in functions in future
global Resource PData TX Trans Receive Recon 

% == input that must match those for data acq ==
Resource.Parameters.speedofSound = 1540;

Trans.frequency = 5.208;        % MHz
Trans.spacingMm = 0.298;        % mm
Trans.numelements = 128;
Receive.samplesPerWave = 4;
TX.steer = [0,0];               % az - 1st, el -2nd, in rad
TX.Origin = [0 0 0];
TX.Delay = zeros(Trans.numelements);

% == intermediate == 
wavelength = Resource.Parameters.speedofSound/Trans.frequency/1000;

Trans.spacing = Trans.spacingMm/wavelength;
Trans.ElementPos = zeros(Trans.numelements, 5);  % x - 1st, y - 2nd, z - 3rd, az - 4th, el - 5th;
Trans.ElementPos(:,1) = -Trans.spacing*(Trans.numelements-1)/2: Trans.spacing: (Trans.numelements-1)/2*Trans.spacing;

PData.PDelta = [0.5 0 0.5];       % wavelenths
PData.startDepth = 5;           % wavelengths
PData.endDepth = 192;           % wavelengths
PData.Origin = [Trans.ElementPos(1,1), 0, P.startDepth];

PData.Size(1) = ceil(PData.endDepth - PData.startDepth)/PData.PDelta(3);    % z - 1st demension, x - 2nd demension, z - 3rd demension 
PData.Size(2) = ceil(Trans.spacing*(Trans.numelements-1)/PData.PDelta(1));
PData.Size(3) = 1;

% == beamforming control ==

Recon.senscutoff = 0.6;
Recon.fnum = 1;

%%

RF = RcvData{1}(:,:,1);
Img = zeros(PData.Size(1), PData.Size(2));


for i=1:PData.Size(1)
    for j=1:PData.Size(2)
        
        % pixel coordinates
        z = PData.Origin(3) + PData.PDelta(3)*(i-1);
        x = PData.Origin(1) + PData.PDelta(1)*(j-1);
        
        % channels contributing to the pixel
        RcvAperture = abs(z/Recon.fnum);
        [~,chStart] = min(abs(Trans.ElementPos(:,1)-(x - RcvAperture/2)));
        [~,chEnd] = min(abs(Trans.ElementPos(:,1)-(x + RcvAperture/2)));
        
        temp = Trans.ElementPos(:,1:3) - repmat(TX.Origin, size(Trans.ElementPos,1),1);
        [~,chOrigin] = min(temp(:,1).^2 + temp(:,2).^2 + temp(:,3).^2);
        
        % pathlengths
        for k = chStart:chEnd
                pathlength_t =  (z-TX.Origin(3))*cos(TX.Steer(1)) + (x-TX.Origin(1))*sin(TX.Steer(1));
                pathlength_r = sqrt((z-Trans.ElementPos(k,3))^2+(x-Trans.ElementPos(k,1)).^2);
            pathlength = (pathlength_t + pathlength_r) + TX.Delay(chOrigin);
            dataIndx = round(pathlength*Receive.samplesPerWave);
                if dataIndx>=1 && dataIndx<=size(RF,1)
                    Img(i,j) = Img(i,j) + RF(dataIndx,k); 
                end
        end
        
    end 
end

%%
% Simple image processing
ImgData1 = abs(Img).^(1/2);
% ImgData1 = 256*(ImgData1-min(ImgData1(:))) / (max(ImgData1(:)-min(ImgData1(:))));  % normalization
figure(1)
imagesc(ImgData1); colormap(gray); axis image;
figure(2)
imagesc(sqrt(ImgData{1}(:,:,1))); colormap(gray); axis image;