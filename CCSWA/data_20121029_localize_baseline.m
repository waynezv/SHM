function [meta, cfg] = data_20121029_localize_baseline(M)

% -------------------------------------------------------------------
% META DATA
% -------------------------------------------------------------------
meta.name    = '20121029_localize_baseline';    % Setup Name
meta.descrip = 'Baseline measurment for localization';   % Description
meta.exptype = 'active';               % Active or passive experiment
meta.fileno  = {M};

% -------------------------------------------------------------------
% DATA FILES
% -------------------------------------------------------------------
if mod(ceil(M/5)-1,2)+1 == 1
    meta.file = {['20121029_localize_baseline_' num2str(M  , '%03i') '.lvm']}; % Data File
else
    meta.file = {['20121029_localize_baseline_' num2str(M  , '%03i') '.lvm']}; % Data File
end

% -------------------------------------------------------------------
% EXPERIMENTAL CONFIGURATION
% -------------------------------------------------------------------
cfg.type = 'plate';
cfg.Rx = [   2.1987   31.1831 ; % 1  - Sensor locations
             4.6693   34.1727 ; % 2
            13.3298   21.4749 ; % 3
            15.2745   10.9300 ; % 4
            18.3249   28.7382 ; % 5
            21.0296    5.7636 ; % 6
            23.4937   10.3590 ; % 7
            24.2559   38.3115 ; % 8
            26.9501   28.4335 ; % 9
            28.1376    7.8686 ; % 10
            29.6309   23.6339 ; % 11
            31.5786   36.9067 ; % 12
            31.6418   21.2292 ; % 13
            33.3560   13.3322 ; % 14
            33.8932    8.9177 ; % 15
            39.4592   36.3576 ; % 16
            45.5611   24.3438 ; % 17
             8.1558    8.7863 ; % 18
             8.0349   38.4063 ; % 19
            13.3162   41.5756 ; % 20
            33.7081   41.9085 ; % 21
            39.5216    5.7762 ; % 22
             2.2029   28.7492 ; % 23
            ]*0.0254;

cfg.bound = [0      0      0      1.2192 ;  % Boundaries
             0      1.2192 0      0      ; 
             1.2192 1.2192 0      1.2192 ;
             0      1.2192 1.2192 1.2192];
cfg.refl   = [1 1 1 1];                     % Reflection coefficients
cfg.pbound = [];                            % Periodic boundaries
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% SAMPLING INFORMATION
% -------------------------------------------------------------------
meta.N    = 10000;                     % Number of samples to extract
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% FILE-EXPERIMENT RELATIONSHIP
% -------------------------------------------------------------------
% SCATTER LOCATION FOR EACH FILE
switch ceil(M/10) 
    case 1
        meta.Sx   = [13.3233   31.1415]*0.0254;   % Known scatter location
    case 2
        meta.Sx   = [23.4795   31.1035]*0.0254;   % Known scatter location
    case 3
        meta.Sx   = [33.2837   31.0668]*0.0254;   % Known scatter location
    case 4
        meta.Sx   = [ 4.6798   23.8715]*0.0254;   % Known scatter location
    case 5
        meta.Sx   = [21.0235   23.7158]*0.0254;   % Known scatter location
    case 6
        meta.Sx   = [26.9534   23.6594]*0.0254;   % Known scatter location
    case 7
        meta.Sx   = [39.4853   23.5400]*0.0254;   % Known scatter location
    case 8
        meta.Sx   = [13.3328   17.0581]*0.0254;   % Known scatter location
    case 9
        meta.Sx   = [23.4891   17.0028]*0.0254;   % Known scatter location
    case 10
        meta.Sx   = [33.3413   16.9492]*0.0254;   % Known scatter location
end

% TRASMISSITERS (Tn) FOR EACH SENSOR
tn = [1:17 1:17];
meta.Tn = {tn(M)};

%  RECIVERS (Rn) FOR EACH SENSOR
if M < 18,       
    rn = 1:8;
    if any(ismember(1:8, tn(M)))
        rn(tn(M)) = 9;
    end
elseif M >= 18
    rn = 10:17;
    if any(ismember(10:17,tn(M)))
        rn(tn(M)-9) = 9;
    end
end
meta.Rn = {rn};

% GET TRANSMITTER AND RECIEVER LOCATIONS FOR THIS FILE
meta.Rx = {cfg.Rx(meta.Rn{1}, :)};
meta.Tx = {cfg.Rx(meta.Tn{1}, :)}; 

% CHANNELS IN FILE
meta.ch = {setdiff((meta.Rn{1} ~= 0).*(1:length(meta.Rn{1})), 0)};


% -------------------------------------------------------------------
% GET DATA
% -------------------------------------------------------------------
% EXTRACT DATA FROM CHOSEN CHANNELS
[x0, s0, meta.Fs, meta.date, meta.time] = readlvm(meta.file{1}, [], meta.ch{1}-1 );
meta.s = {s0};
meta.x = {x0(1:meta.N, :)};  % Truncate to N samples

end