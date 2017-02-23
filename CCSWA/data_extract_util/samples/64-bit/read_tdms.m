% read_tdms: Reads .tdms files using NI provided ReadFile.m.
% data = read_tdms()
%   @para:data_path     
%   @return: data, a cell of raw data from LabVIEW tdms files
%
% Author: wzhao1#andrew.cmu.edu
% Log   : 01/24/2016 - v1.0 - release: first release
%         05/25/2016 - v1.1 - style: some revisions      
%         05/27/2016 - v1.2 - refactor: rewrite the func
function data = read_tdms(data_path)

% +-----------------------------------------------------------------------+
% +                             SET VARIABLES                             +
% +-----------------------------------------------------------------------+
Data_Path_prefix = './data/'; % path storing the raw tdms files
channel_group_name = {'dev0_1-8', 'dev1_9-16'};
dev_ind = 0; % device index

% +-----------------------------------------------------------------------+
% +                           EXTRACT RAW DATA                            +
% +-----------------------------------------------------------------------+
for i = 1:length(channel_group_name)
    Data_Path_temp = fullfile(Data_Path_prefix, channel_group_name{i});
    
    % Get filenames
    filelist_obj  = dir(Data_Path_temp);
    filelist = cell(length(filelist_obj)-2,1); % neglect the first two filenames
                                               % they are the parent paths
    for j = 1:length(filelist)
        filelist{j} = filelist_obj(j+2).name;
        Data_Path = fullfile(Data_Path_temp, filelist{j});
        tmp = regexp(filelist{j}, '\.', 'split');
        file_label = strcat('dev_', num2str(dev_ind), '_file_', tmp{1});
        % Read and save data
        % run './data_extract_util/samples/64-bit/ReadFile.m';
        ReadFile;
    end
    dev_ind = dev_ind + 1;
    
    % Concatenate and reshape data
    data_list = dir(data_path);
    num_files = length(data_list) - 2;
    data_cell = cell(ceil(num_files/2), 1);
    for i = 1:length(data_cell)
        filename1 = data_list(i+2).name; % dev0, channel 8-15
        filename2 = data_list(i+2+17).name; %dev1, channel 0-7
        tmp = regexp(filename1, '_', 'split'); % split name
        ind = str2num(cell2mat(regexp(tmp{4}, '\d{1,2}', 'match'))); % get sensor index
        data_part1 = load(filename2);
        data_part2 = load(filename1);
        data_cell{ind} = [data_part1.data.channel_data data_part2.data.channel_data];
    end  
end

data = data_cell;
end