function x = velwindow(x, d, v, a, dly)
%VELWINDOW   Applies a velocity window to signals
%   X = VELWINDOW(X,D,V,A,DLY)  applies an exponentially tapered window on
%   a collection of signals. VELWINDOW removes information from each signal 
%   in X that travels slower than it takes for a single sample to travel  
%   distances D with a velocity V. 
%
%   INPUTS:
%       x: A Q x N matrix of Q x 1 signals or a M x 1 cell of signal 
%          matrices
%       d: An N x 1 vector of distance corresponding no each signal or an 
%          M x 1 cell of distance vectors
%       v: A scalar velocity value with units of distance unit per sample
%       a: [OPTIONAL] A scalar exponential attenutation factor. Low values 
%          attenuate at fast rates and high values attenuate at slow rates.
%          By default, a = 100.
%     dly: [OPTIONAL] A scalar values with units of samples. Window is
%          delayed by this amount. By default, dly = 0.
%
%   OUTPUTS:
%       x: A Q x N matrix of Q x 1, windowed signal or a M x 1 cell of
%       windowed signal matrices
%
%

% SET DEFAULT PARAMETERS
if (nargin < 4), a = 100; end
if (nargin < 5), dly = 0; end

% CHECK FOR ERRORS
error(nargchk(3, 5, nargin));
if isscalar(v) ~= 1, error('Error: v should be a scalar value'); end
if iscell(x) ~= iscell(d), error('Error: x and d should have the same structure'); end

% WINDOW ALL SIGNALS
if iscell(x)  % If x has multiple measurement sets
    assert(iscell(d));  % Assert d is also a cell
    M = length(x);      % Number of measurement sets
    for m = 1:M         % Loop over measurement sets
        x{m} = velwindow_sub(x{m}, d{m}, v, a, dly);  % Apply velocity window
    end
else         % If x has one measurement set
    x = velwindow_sub(x, d, v, a, dly);  % Apply velocity window
end

end


function x = velwindow_sub(x, d, v, a, dly)
    Q = size(x, 1);  % Number of frequencies
    N = size(x, 2);  % Number of measurements
    for n = 1:N      % Loop over measurements
        rQ = min(max(1,round( d(n,:)/v + dly )),Q);  % Sample with velcoity v
        rm = rQ:Q;                                   % Samples for tapering
        x(rm,n) = x(rm,n).*exp(-(rm-rm(1)).'./a);    % Taper extraneous arrivals
    end
end
