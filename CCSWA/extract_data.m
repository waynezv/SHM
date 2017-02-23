function [meta, cfg] = extract_data(fnc, M)
%EXTRACT_DATA  Extract data from a particular experiment
%   [meta cfg] = extract_data(fnc, M) extracts experiment information 
%   from the experiments MATLAB function
%
%   INPUTS:
%     FNC: The function associated with an experiment, usually entitled 
%          'data_' followed by the base file name 
%       M: An array of files numbers to extract
%
%   OUTPUTS:
%    META: Structure containing specific experimental information/data
%     CFG: Structure containing related to experiemental setup
%
%   see also: readlvm
%

%   Joel Harley  2014/06/18
%   Last Revision Date: 2014/06/18

%CHANGE LOG
%   

tm = 0;
for m = 1:length(M)
    fprintf('%08i / %08i [Time left: %s]\n', m, length(M), datestr(tm/24/3600*(length(M)-m+1), 'HH:MM:SS')); ts = tic; 
    
    % GET DATA
    [meta0, cfg0] = fnc(M(m));
    
    % CONCATINATE META INFORMATION
    if ~exist('meta' ,'var'), 
        meta = meta0; 
    else
        fields = fieldnames(meta);   
        for i = 1:numel(fields)
            if isfield(meta0, fields{i}) 
                if iscell(meta.(fields{i}))
                    meta.(fields{i}) = [meta.(fields{i}); meta0.(fields{i})];
                end
            end
        end
    end
    
    % CONCATINATE CFG INFORMATION
    if ~exist('cfg' ,'var')
        cfg = cfg0; 
    else
        fields = fieldnames(cfg);   
        for i = 1:numel(fields)
            if isfield(cfg0, fields{i}) 
                if iscell(cfg.(fields{i}))
                    cfg.(fields{i}) = [cfg.(fields{i}); cfg0.(fields{i})];
                end
            end
        end
    end
    
    % REFRESH TIME INFORMATION
    fprintf(repmat('\b', 1, 42)); 
    tm = (toc(ts) + tm)/2;
end

end

