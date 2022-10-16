function [c3dStruct] = Read_C3D(filename)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

fileID = fopen(filename,'r');

%% Creating Header section of C3D reader structure
c3dStruct = Get_c3dHeaders(fileID);

%% Creating Parameter section of C3D reader structure
paramStartBlock = (c3dStruct.Headers.ParameterStartBlock-1)*512;
c3dStruct = Get_C3dParameters(fileID, c3dStruct, paramStartBlock);

%% Creating Data Labels for Data section of C3D reader structure
paramList = [c3dStruct.Parameters.ParameterName{1,:}];
idxLabels = contains(paramList,'LABELS');
DataLabels = c3dStruct.Parameters.ParameterName(:,idxLabels);

if size(DataLabels, 2) < 3
    MarkersList = [];
    MarkerLabels = [];
    Analogs = DataLabels(2:end,1);
else
    
    MarkersList = DataLabels(2:end,1);
    MarkerLabels = strtrim(MarkersList(~cellfun('isempty',MarkersList)));
    
    Analogs = DataLabels(2:end,2);
    MarkerLabels = strrep([MarkerLabels{:}], ' ', ' ');
    MarkerLabels = strrep(MarkerLabels(:), '.', ' ');
end

try

AnalogLabels = strtrim(Analogs(~cellfun('isempty',Analogs)));
AnalogLabels = strrep([AnalogLabels{:}], ' ', ' ');
AnalogLabels = strrep(AnalogLabels(:), '.', ' ');

end 

totalFrames = c3dStruct.Headers.LastFrame - c3dStruct.Headers.FirstFrame;

%% Extracting marker and analog data from C3D
[MarkerData,AnalogData] = Extract_C3DData(fileID, c3dStruct.Headers.DataStartBlock, totalFrames, c3dStruct.Headers.NumberOfTrajectories, c3dStruct.Headers.AVRatio, c3dStruct.Headers.AnalogChannels);

if size(MarkerLabels,1) < 1
    Markers = 'No Markers present in trial';
end 

for i = 1:size(MarkerLabels,1)
    try
        Markers.(erase(char(MarkerLabels{i}),' ')) = MarkerData(1:3,:, i)';
    catch
        Markers.(['unlabeled_Marker', num2str(i)]) = MarkerData(1:3,:, i)';
    end
end

EMGChan = contains(AnalogLabels(:), 'EMG');
EMGData = AnalogData(:,EMGChan);

%% Calculating COP Data and separating Force Plate data for Data section of C3D reader structure
[FPData, COPChans] = bertec_COP(AnalogData, 5, 2048);

for k = 1:size(FPData,2)
    ForcePlates.(erase(char(AnalogLabels{k}),{' ', 'r', 'Plate'})) = FPData(:,k);
end

COPLabels = {'LeftX', 'LeftY','RightX','RightY','WeightedX','WeightedY'};
for j = 1:size(COPLabels, 2)
    COP.(COPLabels{j}) = COPChans(:,j);
end

%% Saving Data section of C3D reader structure
c3dStruct.Data.Markers = Markers;
c3dStruct.Data.ForcePlates = ForcePlates;
c3dStruct.Data.COP = COP;
c3dStruct.Data.AnalogLabels = AnalogLabels;
c3dStruct.Data.Analogs = AnalogData;
c3dStruct.Data.EMG = EMGData;
end

