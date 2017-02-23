function Trans = computeTrans(varargin)
% computeTrans - Assign array specifications for known transducer types
% Arguments:
%    One argument:
%        1st arg - Trans
%        if Trans is a structure, return the full Trans structure. If ID string, return transducer name only.
%    Two arguments:
%        1st arg should be transducer name as string.
%        2nd arg should be parameter desired (eg 'maxHighVoltage').
%        Returns value of parameter.
%
% Use built-in defaults, but allow prior specification of some attributes, such as Trans.frequency.
%
% Vantage version modified to shift all transducer frequencies to
% the nearest value supported with 4X sampling using the 250 MHz Vantage system clock
% 
% Trans.id is a double that holds the logical 32-bit value "0x00FFIIDD", where "FF" is the format 
% version and "IIDD" is the scanhead ID. The format versions known are:
%   0 - ATL/Philips transducers
%   1 - Verasonics vendor 1
%   2 - Verasonics vendor 2
% For example, "Trans.id = hex2dec('00000250');" specifies format version 0
% and the L7-4's ID.  This is specified in the code below as"Trans.id = hex2dec('0250');". 
%
% Trans.impedance is used for estimating transmit output current.
% The output current value is used in checking power dissipation limits for
% transmit devices within the system, but not for the transducer itself.
% It is the user's responsibility to keep transmit power levels below the
% safe operating limits of their transducer.
% Note that a non-realistic default impedance of 20 ohms has been set for MUX probes
%   (which should NEVER be used for push, and perhaps not for any extended transmit bursts), 
%   as a means of severely limiting transmit output levels when long duration profile 5 transmit bursts are used. 
%   In fact, the impedance of an HVMux switch is high enough that the switch will overheat during an
%   extended burst, and the artificially low impedance value is used by software to prevent most operation in profile 5.

% Copyright 2001-2016 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.


% Known transducers and their corresponding IDs and high voltage limits. *** Do not modify these values. ***
KnownTransducers = {'L7-4','0250',96,...
                    'L7-4-2','0250',96,...
                    'L10-5','074C',70,...
                    'L11-4v','1ABB4',75,...
                    'L11-5','0351',96,...
                    'L12-3v','1BBC3',75,...
                    'L12-5 38mm','0755',75,...
                    'L12-5 50mm','0B5B',75,...
                    'L22-8v','18A8A',35,...
                    'L22-14v','2AB18',30,...
                    'L22-14vLF','2AB17',30,...
                    'CL10-5','034D',96,...
                    'CL15-7','035C',96,...
                    'C4-2','20D1',96,...
                    'C5-2','20D9',96,...
                    'C5-2v','1AC52',96,...
                    'C7-4','224E',96,... % no data
                    'C8-4V','228C',96,... % no data
                    'C8-5','22DE',96,... % no data
                    'C9-5ICT','228B',96,...
                    'P3-2','4428',96,... % no data
                    'P4-1','483E',96,...
                    'P4-2','4439',96,...
                    'P4-2v','1AA42',96,...
                    'P5-3','4529',96,... % no data
                    'P6-3','4D3B',96,...
                    'P7-4','462A',96,... 
                    'Adapter Embedded S.T.E','1FFFA',96,... % no data
                    '260 ZIF Backshell Kit','BAD00',96,... % no data
                    '260 ZIF S.T.E. Fixture','1FFFC',96,... % no data
                    '408 ZIF S.T.E. Fixture','1FFFB',96}; % no data
                

% The following "known" probes include ID info only (no specs): C7-4,  C8-4V,  C8-5,  P3-2,  P5-3,  P7-4  
% (Specifications for these probes are planned for a future release).


switch nargin
    
    case 1
        Trans = varargin{1};
        if ~isstruct(Trans)  % if a structure is not provided as input, assume input is ID to translate into string.
            n = find(strcmpi(Trans, KnownTransducers));
            if isempty(n), Trans = 'Unknown';
            else Trans = KnownTransducers{n(1)-1};
            end
            return
        end
        if ~isfield(Trans,'name'), error('computeTrans: Trans.name must be provided in input structure.'); end
        n = find(strcmpi(Trans.name, KnownTransducers));
        if isempty(n), error('computeTrans: Trans.name not recognized as known transducer.'); end
        speedOfSound = 1.540;  % default speed of sound in mm/usec
        verbose = 2;
        if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
            if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
                speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
            end
            if evalin('base','isfield(Resource.Parameters,''verbose'')')
                verbose = evalin('base','Resource.Parameters.verbose');
            end
        end

        % check for user-specified units, and print warning message if not found
        if ~isfield(Trans,'units') || isempty(Trans.units)
            fprintf(2, 'Warning: Trans.units not specified; selecting default units of mm.\n');
            fprintf(2, 'If script requires wavelength units, add an explicit definition of\n');
            fprintf(2, '"Trans.units = ''wavelengths'';" before calling computeTrans.\n');
            Trans.units = 'mm';
        end
        if ~strcmp(Trans.units, 'mm') && ~strcmp(Trans.units, 'wavelengths')
            error('computeTrans: Unrecognized value for Trans.units.  Must be ''mm'' or ''wavelengths''.');
        end

        % if Trans.frequency value has already been specified, we will use
        % it as is.  VSX and update() will confirm the value matches the
        % A/D sample rate constraints, and will exit with an error message
        % to the user if not.  Therefore we do not need to validate the
        % Trans.frequency value here (and could not, since we don't know
        % intended use of 4/3 sampling or interleave, etc.).
        if isfield(Trans,'frequency')
            if isempty(Trans.frequency)
                % if empty, remove it so cases below will assign default frequency
                Trans = rmfield(Trans, 'frequency');
            end
        end
        % also allow user-specified Bandwidth to override the default:
        if isfield(Trans,'Bandwidth')
            if isempty(Trans.Bandwidth)
                % if empty, remove it so cases below will assign default
                % Bandwidth
                Trans = rmfield(Trans, 'Bandwidth');
            end
        end

        Trans.lensCorrection = 0; % specify default value, in case it is not set for a particular transducer;
        
        switch KnownTransducers{n}
            case 'L7-4'
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                % Vantage:  5.208 is closest supported frequency to 5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('0250');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = .250; % width in mm
                Trans.spacingMm = .298;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.887; % in mm units; was 3 wavelengths;
                Trans.impedance = [3 23.9-125i;  3.25 25.4-116i;  3.5 26-106i;  3.75 25.4-98.9i;  4 25.9-89.4i;  4.25 27-79.7i;...
                    4.5 32.8-72.6i;  4.75 39.2-66.2i;  5 46.1-69.6i;  5.25 46.5-72.4i;  5.5 41.9-71.6i;  5.75 43.2-69.8i;...
                    6 42.3-69.8i;  6.25 38.2-71i;  6.5 33.5-66.2i;  6.75 32-59.8i;  7 34.4-54.2i;  7.25 37.4-50.3i;...
                    7.5 42.3-48.2i;  7.75 47.8-47.9i;  8 53-51.3i];
                
                case 'L7-4-2'
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                % Vantage:  5.208 is closest supported frequency to 5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('0250');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = .250; % width in mm
                Trans.spacingMm = .298;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0; % in mm units; was 3 wavelengths;
                Trans.impedance = [3 23.9-125i;  3.25 25.4-116i;  3.5 26-106i;  3.75 25.4-98.9i;  4 25.9-89.4i;  4.25 27-79.7i;...
                    4.5 32.8-72.6i;  4.75 39.2-66.2i;  5 46.1-69.6i;  5.25 46.5-72.4i;  5.5 41.9-71.6i;  5.75 43.2-69.8i;...
                    6 42.3-69.8i;  6.25 38.2-71i;  6.5 33.5-66.2i;  6.75 32-59.8i;  7 34.4-54.2i;  7.25 37.4-50.3i;...
                    7.5 42.3-48.2i;  7.75 47.8-47.9i;  8 53-51.3i];

            case 'L10-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('074C');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(17,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L10-5ApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L11-4v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % nominal frequency in MHz
                % Vantage:  6.25 is closest supported frequency to 6.4286 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.5, 10.5]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('1ABB4');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.270; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.4785; % in mm units; was 5 wavelengths
                Trans.impedance = [ 3 17.4-123i;  3.25 18.6-110i;  3.5 20.5-100i;  3.75 20.4-91i;  4 23.3-83.3i;  4.25 21.1-78.3i;...
                    4.5 21.2-69.8i;  4.75 20.8-65.3i;  5 19.6-57.8i;  5.25 21.5-52.2i;  5.5 20.3-47.6i;  5.75 20.4-40.8i;...
                    6 21.3-37.2i;  6.25 20.4-31.5i;  6.5 22.7-25.8i;  6.75 23.1-22.8i;  7 23.5-17.3i;  7.25 26-13.8i;...
                    7.5 25.7-9.62i;  7.75 28.3-4.22i;  8 31.2-2.05i;  8.25 32.2+1.5i;  8.5 36.8+4.35i;  8.75 37.5+3.61i;...
                    9 38.7+7.27i;  9.25 39.7+7.01i;  9.5 38.5+11.5i;  9.75 42+15.4i;  10 42.9+18.2i;  10.25 47+24i;...
                    10.5 53.9+25.6i;  10.75 60.1+27.2i;  11 70.7+25.1i;  11.25 77.2+16.7i;  11.5 82.2+5.99i;  11.75 77-9.58i;...
                    12 65.4-15.8i;  12.25 54.3-17.7i;  12.5 44-15.6i;  12.75 35.2-10.6i;  13 29.1-3.17i;  13.25 26.2+4.47i;...
                    13.5 25.3+11i;  13.75 25.2+16.5i;  14 25.3+20.8i;  14.25 25+24.8i;  14.5 24.6+29.1i;  14.75 24.6+33.6i;...
                    15 25+38.1i;  15.25 25.6+42.4i;  15.5 26.5+46.5i;  15.75 27.5+50.7i;  16 28.9+54.9i];   

            case 'L11-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('0351');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.235; % width in mm
                Trans.spacingMm = 0.260;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = 50; % using default value 


            case 'L12-3v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 12]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('1BBC3');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .170; % width in mm
                Trans.spacingMm = .200;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 5.0, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(33,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L12-3v_ApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L12-5 38mm'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('0755');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(17,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L12-5_38mmApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L12-5 50mm'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 is closest supported frequency to 8.18 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('0B5B');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = .1703; % width in mm
                Trans.spacingMm = .1953;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,129), ...
                                     'VDASAperture', zeros(17,129,'uint8'));
                for i = 0:128, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                ApertureDataFilename = 'L12-5_50mmApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:129
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L22-8v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [8, 21.5]; end % approx. 90% relative bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('18A8A');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = .0703; % width in mm
                Trans.spacingMm = 0.108;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0; % in mm units; was 12 wavelengths;
                Trans.impedance = 5; % using default value for MUX probe
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 35; end
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 6.5, ... %10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,129), ...
                                     'VDASAperture', zeros(17,129,'uint8'));
                for i = 0:128, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                fid = fopen('L22-8vApertureData.txt');
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:129
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);
                                        
            case 'L22-14v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate), 
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('2AB18');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
                % Note elevation aperture and focus will not be converted
                % to wavelengths regardless of Trans.units
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [14.5, 9.0-21.0i; 15.0, 10.0-22.0i; 15.5, 13.0-22.0i; 16.0, 15.0-22.0i; 16.5, 17.0-21.0i; 17.0, 20.0-21.0i;...
                    17.5, 22.0-21.0i; 18.0, 26.0-16.0i; 18.5, 30.0-16.0i; 19.0, 30.0-11.0i; 19.5, 32.0-7.0i; 20.0, 34.0-4.0i; 20.5, 36.0-4.0i;...
                    21.0, 37.0-3.0i; 21.5, 41.0-5.0i; 22.0, 41.0+1.0i; 22.5, 41.0+7.0i; 23.0, 44.0+10.0i; 23.5, 44.0+12.0i; 24.0, 45.0+16.0i;...
                    24.5, 48.0+20.0i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
                                        
            case 'L22-14vLF'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate), 
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('2AB17');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elevationApertureMm = 3; % Long Focus active elevation aperture in mm
                Trans.elevationFocusMm = 20; % Long Focus nominal elevation focus depth from lens on face of transducer
                % Note elevation aperture and focus will not be converted
                % to wavelengths regardless of Trans.units
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [14.5, 9.0-21.0i; 15.0, 10.0-22.0i; 15.5, 13.0-22.0i; 16.0, 15.0-22.0i; 16.5, 17.0-21.0i; 17.0, 20.0-21.0i;...
                    17.5, 22.0-21.0i; 18.0, 26.0-16.0i; 18.5, 30.0-16.0i; 19.0, 30.0-11.0i; 19.5, 32.0-7.0i; 20.0, 34.0-4.0i; 20.5, 36.0-4.0i;...
                    21.0, 37.0-3.0i; 21.5, 41.0-5.0i; 22.0, 41.0+1.0i; 22.5, 41.0+7.0i; 23.0, 44.0+10.0i; 23.5, 44.0+12.0i; 24.0, 45.0+16.0i;...
                    24.5, 48.0+20.0i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
                
            case 'CL10-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % nominal frequency in MHz (this is the closest allowed frequency)
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.5, 9]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('034D');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.spacingMm = 0.200;   % Spacing between elements in mm.
                Trans.elementWidth = Trans.spacingMm - 0.05; % width in mm
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.1; % in mm units; was 9 wavelengths;
                Trans.impedance = 50; % using default value


            case 'CL15-7'
                if ~isfield(Trans,'frequency'), Trans.frequency = 8.929; end % nominal frequency in MHz
                % Vantage:  8.929 is closest supported frequency to 9.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 9*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('035C');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.16; % width in mm
                Trans.spacingMm = 0.178;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.6899; % in mm units; was 4 wavelengths;
                Trans.impedance = 50; % using default value


            case 'C4-2'
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.976; end % nominal frequency in MHz
                % Vantage:  2.976 is closest supported frequency to 3.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.id = hex2dec('20D1');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end

                Trans.impedance = 50; % using default value

            case 'C5-2'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.125; end % nominal frequency in MHz
                % Vantage:  3.125 is closest supported frequency to 3.2143 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4.5]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.id = hex2dec('20D9');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = [ 1 19.4-255i;  1.25 20.1-186i;  1.5 22.4-136i;  1.75 29-94.8i;  2 43.6-62i;  2.25 64.8-48.9i;...
                    2.5 70.2-54.9i;  2.75 54.3-49.1i;  3 45.4-30.9i;  3.25 42.9-13.3i;  3.5 42.5+1.9i;  3.75 42.1+16.4i;...
                    4 41.5+31.1i;  4.25 41.2+46.3i;  4.5 42.3+62.7i;  4.75 44.8+78.7i;  5 45.2+93i;  5.25 40.9+114i;...
                    5.5 42.6+144i;  5.75 51.2+177i;  6 66.5+214i];

            case 'C5-2v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.7; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2.4, 5]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.id = hex2dec('1AC52');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                radiusMm = 49.57;  % radius in mm.
                spacingMm = .508; % spacing in mm.
                kerf = .048;   % in mm.
                scanangle = 128*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.035; % in mm units; was 3 wavelengths;
                Trans.impedance = [ 1 39.2-216i;  1.25 41-162i;  1.5 41.1-118i;  1.75 47.5-87i;  2 50.8-69.7i;  2.25 44.7-49.5i;...
                    2.5 43.3-27i;  2.75 48.6-8.53i;  3 50.4+4.4i;  3.25 50.7+20.7i;  3.5 56.5+38.2i;  3.75 69.1+51.8i;...
                    4 83.2+53.8i;  4.25 91.9+45.7i;  4.5 81.7+34.6i;  4.75 66+40.1i;  5 56.6+47i;  5.25 42.4+57.2i;...
                    5.5 30.7+75i;  5.75 25.3+95.5i;  6 24.3+115i;  6.25 24.8+132i;  6.5 26+148i;  6.75 27.8+164i;  7 30+179i];

            case 'C9-5ICT'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.id = hex2dec('228B');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 146.677 * (pi/180);    % degrees converted to radians
                radiusMm = 8.511;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .025;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = 50; % using default value

            case 'P4-1'
                % The P4-1 is a 96 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.5, 3.5]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('483E');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 96;
                Trans.Connector = [023 022 021 041 024 042 046 043 045 044 047 018 017 048 013 020 ...
                                   019 014 015 016 049 050 054 051 053 052 009 055 056 011 012 005 ...
                                   006 007 008 010 004 003 002 001 040 039 038 037 033 034 035 036 ...
                                   093 094 095 096 092 091 090 089 128 127 126 125 119 121 122 123 ...
                                   124 117 118 073 074 120 077 076 078 075 079 080 113 114 115 110 ...
                                   109 116 081 112 111 082 085 084 086 083 087 105 088 108 107 106]';
                kerf = .050;   % guess (in mm)
                Trans.spacingMm = 0.2950;   % Spacing between elements in mm.
                Trans.elementWidth = (0.2950-kerf); % width in mm
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:96,1) = Trans.spacingMm*(-((96-1)/2):((96-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.464; % in mm units; was 4 wavelengths;
                Trans.impedance = [ 1 33.9-390i;  1.25 43.1-286i;  1.5 63.9-217i;  1.75 80.1-185i;  2 74.1-159i;...
                    2.25 73.9-126i;  2.5 84-99.2i;  2.75 98-93.3i;  3 89.5-99.4i;  3.25 63.6-89.5i;  3.5 44.2-63.6i;...
                    3.75 32.7-34.5i;  4 25.8-3.09i;  4.25 24.1+29.5i;  4.5 26.3+62.6i;  4.75 33.2+96.9i;  5 45.2+130i];   

            case 'P4-2'
                % The P4-2 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 32-63 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('4439');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = [1:32,97:128]';
                Trans.elementWidth = 0.2950; % width in mm
                Trans.spacingMm = 0.3200;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 3.08; % in mm units; was 5 wavelengths;
                Trans.impedance = [1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
                    4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            case 'P4-2v'
                % The P4-2 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 1-64 are
                %    wired to connector inputs 33-96.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.976; end % nominal frequency in MHz
                % Vantage:  2.976 is closest supported frequency to 3.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('1AA42');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = (33:96)';
                Trans.elementWidth = 0.250; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.587; % in mm units; was 5 wavelengths;
                Trans.impedance = [ 1 33.8-366i;  1.25 49.4-252i;  1.5 64.1-196i;  1.75 75.5-140i;  2 99.5-118i;...
                    2.25 84.2-104i;  2.5 73.9-74.6i;  2.75 73.5-36i;  3 94.5-5.42i;  3.25 115-9.99i;  3.5 111-24.8i;...
                    3.75 74.9-14.8i;  4 66.9+17.5i;  4.25 70.2+22.4i;  4.5 37.7+37.8i;  4.75 26+75.3i;  5 24.8+107i;...
                    5.25 26.1+135i;  5.5 28.5+162i;  5.75 31.7+188i;  6 36.3+215i];   


            case 'P6-3'
                if ~isfield(Trans,'frequency'), Trans.frequency = 4.464; end % nominal frequency in MHz
                % Vantage:  4.464 is closest supported frequency to 4.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3, 6]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('4D3B');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = (.218-.025); % width in mm
                Trans.spacingMm = 0.218;   % Spacing between elements in mm
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:128,1) = Trans.spacingMm*(-((128-1)/2):((128-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.380; % in mm units; was 4 wavelengths;
                Trans.impedance = 48; % Z @ 3.4 Mhz

            case 'P7-4'
                % The P7-4 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 32-63 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.id = hex2dec('462A');
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = [1:32,97:128]';
                Trans.elementWidth = 0.1369; % width in mm
                Trans.spacingMm = 0.1711;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.7; % in mm units; was 5 wavelengths;
                Trans.impedance = 50;%[1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
%                     4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            otherwise
                Trans.id = hex2dec(KnownTransducers{n+1});
                Trans.type = [];
                Trans.frequency = [];
                Trans.spacingMm = [];
                Trans.elementWidth = [];
                Trans.lensCorrection = [];
                Trans.ElementPos = [ 0 0 0 0 ];
                if verbose > 2
                    disp(' ');
                    disp(['computeTrans Status: Data not available for ', Trans.name]);
                    disp('Trans structure must be provided in user script.');
                    disp(' ');
                end
%                 fprintf(2, ['computeTrans: Data not available for ', Trans.name, ';\n']);
%                 fprintf(2, 'Trans structure must be provided in user script.\n');
%                 error(' ');

        end

        % Set a conservative value for maxHighVoltage, if not already defined
        if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

        % Now convert all units as required, based on Trans.units
        scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths

        % regardless of units, always provide spacing and radius in wavelengths
        Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths.
        if Trans.type == 1
            Trans.radius = Trans.radiusMm * scaleToWvl; % convert radiusMm to wavelengths
        end
        if strcmp(Trans.units, 'wavelengths')
            % convert all mm unit variables to wavelengths
            Trans.elementWidth = Trans.elementWidth * scaleToWvl;
            Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
            Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
            Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
            Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
        end
        
    case 2
        % Two inputs provided - 1st input is Trans.name, 2nd is parameter desired (currently only
        % 'maxHighVoltage' allowed).
        Trans = varargin{1};
        if ischar(Trans)
            n = find(strcmpi(Trans, KnownTransducers));
            % if 1st input is not a recognized transducer name, set n to
            % zero to flag as unrecognized transducer
            if isempty(n)
                n = 0; % special case of zero will trigger assignment of default value
            end
            Param = varargin{2};
            switch Param
                case 'maxHighVoltage'
                    if n == 0
                        % unrecognized transducer name, so assign default
                        % maxHighVoltage limit at hw max value
                        Trans = 96;
                    else
                        % return maxHighVoltage from KnownTransducers array
                        Trans = KnownTransducers{n(1)+2};
                    end
                otherwise
                    error('computeTrans: Unrecognized parameter in 2nd argument.');
            end
        else
            error('computeTrans: When called with 2 inputs, 1st input must be transducer name.');
        end
        
    otherwise
        error('computeTrans: computeTrans accepts one or two input arguments.');

end
return
        

