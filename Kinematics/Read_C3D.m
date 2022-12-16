function [c3dStruct] = Read_C3D(filename)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename,'r');

%% Creating Header section of C3D reader structure
c3dStruct = Get_c3dHeaders(fileID);

%% Creating Parameter section of C3D reader structure
paramStartBlock = (c3dStruct.Headers.ParameterStartBlock-1)*512;
c3dStruct = Get_C3dParameters(fileID, c3dStruct, paramStartBlock);

%% Retrieving raw 3D point data and Analog data
firstBin = (c3dStruct.Headers.DataStartBlock - 1) * 512;

fseek(fileID, firstBin,'bof');
DataMatrix = fread(fileID, 'float32', 'ieee-le');
fclose(fileID);
%% Creating Data Labels for Data section of C3D reader structure
paramList = [c3dStruct.Parameters.ParameterName{1,:}];

try
    idxH = contains(paramList,'Height');
    idxW = contains(paramList,'Bodymass');
    Heightcell = c3dStruct.Parameters.ParameterName(1:2,idxH);
    Weightcell = c3dStruct.Parameters.ParameterName(1:2,idxW);
   
    Heightcell = strtrim(Heightcell(~cellfun('isempty',Heightcell)));
    HeightLabel = strrep([Heightcell{1}], ' ', ' ');
    HeightLabel = erase(HeightLabel, 'PROCESSING_');
    HeightVal = str2num(strrep(cell2mat(Heightcell{2}), '', ' '));
    Weightcell = strtrim(Weightcell(~cellfun('isempty',Weightcell)));
    WeightLabel = strrep([Weightcell{1}], ' ', ' ');
    WeightLabel = erase(weightLabel, 'PROCESSING_');
    WeightVal = str2num(strrep(cell2mat(Weightcell{2}), '', ' '));
    c3dStruct.Data.(HeightLabel) = HeightVal;
    c3dStruct.Data.(WeightLabel) = WeightVal;
catch
    
end

try
    idxMarkers = contains(paramList,'POINT_LABELS');
    MarkerLabels = c3dStruct.Parameters.ParameterName(:,idxMarkers);
    MarkersList = MarkerLabels(2:end,1);
    MarkerLabels = strtrim(MarkersList(~cellfun('isempty',MarkersList)));
    MarkerLabels = strrep([MarkerLabels{:}], ' ', ' ');
    MarkerLabels = strrep(MarkerLabels(:), '.', ' ');
catch
    MarkerLabels = [];
end

try
    idxAnalogs = contains(paramList,"ANALOG_LABELS");
    AnalogLabels = c3dStruct.Parameters.ParameterName(:,idxAnalogs);
    Analogs = AnalogLabels(2:end,1);
    AnalogLabels = strtrim(Analogs(~cellfun('isempty',Analogs)));
    AnalogLabels = strrep([AnalogLabels{:}], ' ', ' ');
    AnalogLabels = strrep(AnalogLabels(:), '.', ' ');
catch
    AnalogLabels = [];
end

totalFrames = c3dStruct.Headers.LastFrame - c3dStruct.Headers.FirstFrame;

%% Extracting marker and analog data from C3D
[MarkerData,AnalogData] = Extract_C3DData(DataMatrix, c3dStruct.Headers.NumberOfTrajectories, c3dStruct.Headers.AVRatio, c3dStruct.Headers.AnalogChannels);

if size(MarkerLabels,1) < 1
    Markers = 'No Markers present in trial';
else
    loop = 0;
    for i = 1:4:size(MarkerData,1)
        loop = loop + 1;
        try
            Markers.(erase(char(MarkerLabels{loop}),' ')) = MarkerData(i:i+2,:)';
        catch
            Markers.(['unlabeled_Marker', num2str(loop)]) = MarkerData(i:i+2,:)';
        end
    end
end

try
    EMGChan = contains(AnalogLabels(:), 'EMG');
    EMGData = AnalogData(:,EMGChan);
    AnalogData(:,EMGChan) = [];
    AnalogLabels(EMGChan, :) = [];
catch
    EMGData = [];
end

if contains(AnalogLabels(:), 'States')
    removeExtraChans = contains(AnalogLabels(:), 'States');
    AnalogData(:,removeExtraChans) = [];
    AnalogLabels(removeExtraChans, :) = [];
end

try
    for k = 1:size(AnalogLabels,1)
        rawForcePlates.(erase(char(AnalogLabels{k}),{' ', 'r', 'Plate'})) = AnalogData(:,k);
    end
end
%% Calculating COP Data and separating Force Plate data for Data section of C3D reader structure
try
    [ForcePlates, COPChans] = bertec_COP(rawForcePlates, 5, 2048);
    
    COPLabels = {'LeftX', 'LeftY','RightX','RightY','WeightedX','WeightedY'};
    for j = 1:size(COPLabels, 2)
        COP.(COPLabels{j}) = COPChans(:,j);
    end
    
end
%% Saving Data section of C3D reader structure
c3dStruct.Data.Markers = Markers;
c3dStruct.Data.ForcePlates = ForcePlates;
c3dStruct.Data.COP = COP;
c3dStruct.Data.AnalogLabels = AnalogLabels;
c3dStruct.Data.Analogs = AnalogData;
c3dStruct.Data.EMG = EMGData;
end

