% RSST Final Analysis User code
% last edit 10/4
% needs rsst.m and RSSTAnalysis.m
% Buzsaki data also requires bload function from hc-2 crcns.com
%% Input
% Data Input
    path = "/Users/Mingda/Desktop/RSST/Data/Buzsaki HC-2/";
    file = "ec013.527.dat";
    plotTitle = "Buzsaki subject 013, trial 527"; % title

    
    %currentFreq = 938; % used for .csv files
    %columnToAnalyze = 2; % used for .csv files

    electrodeNum = 1; % used for .set files

    columnToAnalyze = 4; % used for .dat files
    currentFreq = 20000; % used for .dat files

    datasetNum = 49; % used for .mat files


    % leave start and end as 0 to include the entire dataset
    datasetStart = 0; % starting datapoint in seconds
    datasetEnd = 10; % ending datapoint in seconds


% Plot Settings
    plotAll = 'Off'; % 'On' or 'Off'
    savePlotFile = append(path,file); % leave as '' to not save a file

    % imports time data
    %filePath = append(path,file);
    %timeTable = readtable(filePath);
    %timeTable = timeTable(:,1);
    %currentTimeData = table2array(timeTable);
    currentTimeData = [];

    yMin = 0.3;
    yMax = 160;
    downsample = 'Off';
    downsampleRate = 10;
    listBoxAlgorithm = 'SSWT'; % 'SSWT' or 'CWT'
    switchUnitType = 'Potential'; %'Power' or 'Potential'
    waveParameter = 5; % frequency / stddev(frequency), higher values results in greater frequency resolution
    switchYScale = 'Log'; %'Linear' or 'Log' Y scale
    useMs = 'Off'; % use milliseconds for x-axis

    filterAC = 'Off'; % filters out 60 hz data from AC hum
    removeCOI = 'On'; % adds linear data to front and end to remove the coi from graph for aesthetics. coi imprecision still has an effect.
    editFieldSensitivity = 1; % increases sensitivity non-linearly
    boost = 1; % Boosts sensitivity linearly
    switchSqHz = 'Off'; % increase sensitivity for higher frequencies; Only necessary for analysis of very high frequencies



%% load files
    %converts seconds to datapoint #s
    datasetStart = floor(datasetStart*currentFreq);
    if datasetStart == 0
        datasetStart = 1;
    end
    datasetEnd = floor(datasetEnd*currentFreq);
    % load .csv file
    if contains(file,".csv")
        filePath = append(path,file);
        eegTable = readtable(filePath);
        eegTable = eegTable(:,columnToAnalyze);
        currentData = table2array(eegTable);
    end

    % load .MAT matlab struct format
    % struct contains fields including DatasetName, Fs, and Y
    if contains(file,".mat")
        load(append(path,file));
%         %converts seconds to datapoint #s
%         currentFreq = Database(datasetNum).Fs;
%         %Database(datasetNum).X; unnecessary since we use the sampling rate
%         %to determine time data
%         currentData = Database(datasetNum).Y(:); % EEG data
%         currentData = currentData/1000; % use if data is amplified (i.e. Fried data)
        currentData = FT_data.trial{1,1}(1,1:50000);

    end

    % load .set EEGLAB file
    if contains(file,".set")
        load(append(path,file), '-mat');
        currentData = EEG.chaninterp.data(electrodeNum, :);
        currentData = currentData/10; % use if data is amplified (i.e. Lenartowicz data)
    end

    % load Buzsaki Binary data (HC-2)
    if contains(file,".dat")
        eegArray = bload(append(path,file), [columnToAnalyze datasetEnd*currentFreq]);
        eegArray = eegArray(columnToAnalyze,:);
        currentData = eegArray.';
        currentData = currentData(:)/1000;
        % /1000 to remove amplification of Buzsaki data
    end
    % load Audio data (HC-2)
    if contains(file,".wav")
        [y, Fs] = audioread(append(path,file));
        eegArray = y;
        currentData = eegArray.';
        currentFreq = Fs;
    end
    if datasetEnd == 0
        currentData = currentData(:);
    else
        currentData = currentData(datasetStart:datasetEnd);
    end
    %% set input parameters

    if strcmpi(downsample, 'On')
        currentData = databin(currentData, downsampleRate);
        currentFreq = floor(currentFreq/downsampleRate);
    end

    if isrow(currentData) == true
        currentData = currentData.';
    end
    currentData = rmmissing(currentData);


    if currentFreq <= yMax*2
        disp("Sampling rate is low, higher frequencies will be blank.")
    end
    disp("plotting " + plotTitle);

    if strcmpi(plotAll,'On')
        for i = 1:4
            if i == 1
                listBoxAlgorithm = 'SSWT'; % 'SSWT' or 'CWT'
                switchUnitType = 'Potential'; %'Power';
            elseif i == 2
                listBoxAlgorithm = 'SSWT'; % 'SSWT' or 'CWT'
                switchUnitType = 'Power'; %'Power';
            elseif i == 3
                listBoxAlgorithm = 'CWT'; % 'SSWT' or 'CWT'
                switchUnitType = 'Potential'; %'Power';
            elseif i == 4
                listBoxAlgorithm = 'CWT'; % 'SSWT' or 'CWT'
                switchUnitType = 'Power'; %'Power';
            end
                RSSTAnalysis(plotTitle, currentData, currentFreq,...
                    'algorithm', listBoxAlgorithm,...
                    'unitType', switchUnitType,...
                    'yScale', switchYScale,...
                    'sqHz', switchSqHz,...
                    'removeCOI', removeCOI,...
                    'filterAC', filterAC,...
                    'useMs', useMs,...
                    'waveParameter', waveParameter,...
                    'currentTimeData', currentTimeData,...
                    'sensitivity', editFieldSensitivity,...
                    'yMin', yMin,...
                    'yMax', yMax,...
                    'boost', boost,...
            'savePlotFile',append(path,file));
        end
    else
        RSSTAnalysis(plotTitle, currentData, currentFreq,...
            'algorithm', listBoxAlgorithm,...
            'unitType', switchUnitType,...
            'yScale', switchYScale,...
            'sqHz', switchSqHz,...
            'removeCOI', removeCOI,...
            'filterAC', filterAC,...
            'useMs', useMs,...
            'waveParameter', waveParameter,...
            'currentTimeData', currentTimeData,...
            'sensitivity', editFieldSensitivity,...
            'yMin', yMin,...
            'yMax', yMax,...
            'boost', boost,...
            'savePlotFile',append(path,file));
    end


function [x_int] = databin(x, binSize)
    if isrow(x) == true
        x = x';
    end
    lenX = length(x);

    deleteEnd = mod(lenX, binSize);
    if deleteEnd ~= 0
        x(lenX-deleteEnd+1:lenX) = [];
    end

    reshapeSize = (lenX-deleteEnd)/binSize;
    x_res = reshape(x, [binSize, reshapeSize]);
    x_int = mean(x_res, 1);
end
