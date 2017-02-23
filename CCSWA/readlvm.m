function [y, s, Fs, d, t] = readlvm(file, N, ch)
%READLVM  Extract data from labview measurement file
%   [ Y, S, FS, D, T ] = COMPILEDATA(FILE, K, CH) collects LabView
%   measurements from FILE and returns the excitation and response 
%   signals, the times of measurement, and the sampling rate. 
%
%       INPUTS:
%        FILE: LabView measurement file with data to extract
%           K: [OPTIONAL] Labview files to extract data from (by number). 
%              If K is a scalar, then data is extract fro all files up to  
%              that number. For example, if FILE = 'file.lvm' and K = 3, 
%              then data will be extracted from 'file_001.lvm', 
%              'file_002.lvm', and 'file_003.lvm'.
%
%              If K is a vector then data is extracted from each value in
%              K. For example, if K = [1 5], then data is extracted from 
%              'file_001.lvm' and 'file_005.lvm'.
%
%              If K is empty (K = []), data is  extracted from FILE without
%              a numerical suffix. By default, K is empty. 
%          CH: [OPTIONAL] A vector of channels to extract data from. By
%              default, data is extracted only from the first channel, 
%              CH = 0. 
%
%       OUTPUTS:
%           Y: An N x C x K matrix of data from files where K is the number
%              of labview files extracted data from and C is the number of
%              channels
%           S: An N x K matrix of excitation signals from files where K is
%              the number of labview files extracted data from
%          FS: An array of length K of sampling rates of data from each
%              file
%           D: A cell array of length K of the date of measurement for each 
%              file
%           T: A cell array of length K of the time of measurement for each 
%              file
%
%   see also: get_data, get_samp, get_time, compilelvm
%

%   Joel Harley  2009/04/10
%   Last Revision Date: 2011/04/28

%CHANGE LOG
%   2009/10/19  Completely rewrote function
%   2010/05/22  Changed output from cellarray to matrix -- easier to use
%   2010/12/27  If N < 1 now, only the raw filename is used
%   2011/04/18  Included channel documentation
%   2011/06/02  Merged with gettime. Removed unnecessary function parameter.
%               Removed directory parameter. 
%   2011/06/22  Changed Y into a 3-dimensional matrix of samples by files
%               by channel for easier use. Made FS an array with the
%               sampling rate for each file. 
%

% SET DEFAULT PARAMETERS
if nargin < 2, N = []; end
if nargin >= 2, MULTI=true; end
if nargin < 3, ch = 0; end

% IF N IS SCALAR, USE ALL NUMBERS UP TO THAT NUMBER 
% IF N IS VECTOR, USE ALL NUMBERS IN THE VECTOR
if length(N) == 1
    Narray = 1:N;
else
    Narray = N;
end

% IF NARRAY IS EMPTY, USE THE FILENAME (NO ADDITIONAL SUFFIX)
if isempty(Narray)
    MULTI=false;
    Narray = 1;
end    

% INITIALIZE CELLS
s1 = cell(length(Narray),1); y1 = cell(length(Narray),1);
dt = cell(length(Narray),1); tm = cell(length(Narray),1); 
Fs0 = cell(length(Narray),1);

% FOR EACH FILE
for k = 1:length(Narray)

    n = Narray(k);
    
    % LOAD DATA
    if MULTI == false
        [s1{n}, y1{n}] = get_data(file, ch);
        [dt{k}, tm{k}] = get_time(file);
        Fs0{k} = get_samp(file);
    else 
        % SEPARATE FILENAME AND EXTENSION
        extloc = regexp(file, '\.[^.]+$');
        filebase = file(1:extloc-1);
        fileext = file(extloc:end);
        
        % ADD A FILE # SUFFIX
        filen = [filebase '_' num2str(n,'%03u') ...
            fileext];

        % EXTRACT DATA
        [s1{k}, y1{k}] = get_data(filen, ch);
        [dt{k}, tm{k}] = get_time(filen);
        Fs0{k} = get_samp(filen);
    end

end

% CONVERT CELLS TO MATRICES
y = cell2mat(permute(y1.', [1 3 2])); s = cell2mat(s1.');
Fs = cell2mat(Fs0);

% KEEP AS CELLS
d = (dt); t = (tm); 

end



function [s, y] = get_data(file, ch)
%GET_DATA  Extract data from a labview measurement file (.lvm)
%   [S Y] = get_data(FILE, CH) extracts data from a labview measurement
%   file from FILE and channel CH. 
%
%   INPUTS:
%    FILE: The .lvm file location to extract data from 
%      CH: [OPTIONAL] A scalar or vector of channels to extract data from.
%      By default, the first channel (0) is taken. 
%
%   OUTPUTS:
%       S: Excitation signal
%       Y: Response signal
%
%   see also: get_time, get_samp, readlvm
%

%   Joel Harley  2009/04/10
%   Last Revision Date: 2011/06/02

%CHANGE LOG
%   2010/11/23  Added an actual description --  been using this function for
%               a while
%   2011/06/02  Updated description and code. Merged with importfile. 
%

    % SET DEFAULT PARAMETERS
    if nargin < 2
        ch = 0;
    end

    % DEFINE CONSTANTS
    DELIMITER = '\t';
    HEADERLINES = 22;

    try
        % IMPORT FILE
        data0 = importdata(file, DELIMITER, HEADERLINES);
        if isfield(data0, 'data')
            data = data0.data;
        else
            data = data0;
        end

        % EXTRACT RESPONSE
        y = data(:,3+ch);
        
        % EXTRACT EXCIATION
        sl = sum(~isnan(data(:,2)));
        s = data(1:sl,2);

    catch ME
        fprintf('ERROR: Unable to access file %s\n', file);
        rethrow(ME)
    end
    
end



function [dt, tm] = get_time(file)
%GET_TIME  Extract measurment time from a labview measurement file (.lvm)
%   [DT TM] = GET_TIME(FILE) extracts the measurement time from FILE (.lvm)
% 
%   INPUTS:
%    FILE: The .lvm file location to extract data from 
%
%   OUTPUTS:
%      DT: Extracted date of measurement
%      TM: Extracted time of measurement
%
%   see also: get_data, get_samp, readlvm
%

%   Joel Harley  2011/04/14
%   Last Revision Date: 2011/06/02

%CHANGE LOG
%   2011/04/14  Wrote function
%   2011/06/02  Renamed IMPORTFILETIME to GET_TIME. Updated description. 


% OPEN FILE
fid = fopen(file);

% EXTRACT DATE INFORMATION
C = textscan(fid, 'Date %f/%f/%f %*[^\n] %*s', 'HeaderLines', 15);
dt = [num2str(C{1}(1)) '/' num2str(C{2}(1)) '/' num2str(C{3}(1))];

% EXTRACT TIME INFORMATION
C = textscan(fid, 'Time %f:%f:%12.10f', 'HeaderLines', 0);
tm = [num2str(C{1}(1)) ':' num2str(C{2}(1)) ':' num2str(C{3}(1))];

% CLOSE FILE
fclose(fid);


end


function Fs = get_samp(file)
%GET_SAMP  Extract sampling rate from a labview measurement file (.lvm)
%   FS = GET_SAMP(FILE) extracts the sampling rate from FILE (.lvm)
%
%   INPUTS:
%    FILE: The .lvm file location to extract data from 
%
%   OUTPUTS:
%      FS: Extracted sampling rate
%
%   see also: get_data, get_time, readlvm
%

%   Last Revision Date: 2011/06/02

%CHANGE LOG
%   2010/11/23  Added a description
%   2011/06/02  Updated desciption and code. Merged with importfile. 


    % DEFINE CONSTANTS
    DELIMITER = '\t';
    HEADERLINES = 19;

    try
        
        % IMPORT FILE
        data0 = importdata(file, DELIMITER, HEADERLINES);
        if isfield(data0, 'data')
            data = data0.data;
        else
            data = data0;
        end

    catch ME 
        fprintf('ERROR: Unable to access file %s\n', file);
        rethrow(ME)
    end

    Fs = 1/data(1);
        
    
end
