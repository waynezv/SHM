% get_baseline: Returns the meta data and the experimental configurations.
% [meta, cfg] = get_baseline(M)
%   @para: M, files index to extract
%   @return: meta, meta data
%             cfg, configurations
%
% Author: wzhao1#andrew.cmu.edu
% Log   : 05/26/2016 - v1.0 - release: first release

function [meta, cfg] = get_baseline(M)
meta.name = 'baseline data';
meta.descrip = 'Baseline measurement for plate dispersion experiment';
meta.exptype = 'active';
meta.fileno = {M};
meta.N = 10000; % number of samples to extract
meta.Fs = 10e6; % sample rate

cfg.type = 'plate';
cfg.Rx = [  22  4; % 1 - sensor locations, m
            20  7; % 2
            20 17; % .
            21 21; % .
            16 17; % .
            17 13;
            16 10;
            14 14;
            14 19;
            11 15;
            13  9;
            10  5;
             9 13;
             7  9;
             5  4;
             6 15;
             8 19 % 17
          ] * 5e-2;
      
cfg.bound = [0      0      0      1.2 ;  % boundaries
             0      1.2    0      0      ; 
             1.2    1.2    0      1.2 ;
             0      1.2    1.2    1.2];
         
cfg.refl   = [1 1 1 1];                     % reflection coefficients
cfg.pbound = [];                            % periodic boundaries         

% Transmitters
tn = [1:17];
meta.Tn = {tn(M)};
meta.Tx = {cfg.Rx(meta.Tn{1}, :)};

% Recivers
rn = [1:17];
if any(ismember(rn, meta.Tn{1}))
    meta.Rn = {setdiff(rn, meta.Tn{1}, 'stable')};
else
    meta.Rn = {rn};
end
meta.Rx = {cfg.Rx(meta.Rn{1}, :)};

% Channels
% meta.ch = {setdiff((meta.Rn{1} ~= 0).*(1:length(meta.Rn{1})), 0)};
meta.ch = meta.Rn;
end