function [rawMarkerData,rawAnalogData] = Extract_C3DData(DataMatrix, MarkerNumber, AVRatio, AnalogNumber)

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

try
points = MarkerNumber * 4 + AnalogNumber * AVRatio;
sortPts = length(DataMatrix) / points;
PointMatrix = reshape(DataMatrix,[],sortPts); 
rawMarkerData = PointMatrix([1:MarkerNumber*4], :);
catch
    rawMarkerData = [];
end

try
    if MarkerNumber < 1
        resAnalogData = PointMatrix([1:end], :);
    else
        resAnalogData = PointMatrix([MarkerNumber*4:end], :);
        resAnalogData(1,:) = [];
    end 
sortAnalog = AVRatio * size(resAnalogData,2);
rawAnalogData = reshape(resAnalogData, [AnalogNumber, sortAnalog])';
catch
    rawAnalogData = [];
end

%%%% Old extraction code replaced by faster matrix operations %%%%
%%% This code won't work without 

% firstBin = (dataStartBin - 1) * 512;
% 
% fseek(fileID, firstBin,'bof');
% DataMatrix = fread(fileID, 'float32', 'ieee-le');
% fclose(fileID);

% rawMarkerData = zeros(4,totalVideoFrames,MarkerNumber);
% AnalogDataPoints = totalVideoFrames * AVRatio;
% rawAnalogData = zeros(AnalogDataPoints, AnalogNumber);
% analogRow = 0;
% 
% for i = 1:totalVideoFrames
%     
%     try
%         rawMarkers = zeros(4,1, MarkerNumber);
%         for j = 1:MarkerNumber
%             xyzCoords = fread(fileID, 4, 'float32', 'ieee-le');
%             rawMarkers(:,:,j) = xyzCoords;
%         end
%         rawMarkerData(:,i,:) = rawMarkers;
%     catch
%     end
%     
%     try
%         
%         for ii = 1:AVRatio
%             analogRow = analogRow + 1;
%             allChannels = zeros(1,15);
%             for jj = 1:AnalogNumber
%                 singleChannel = fread(fileID, 1, 'float32', 'ieee-le');
%                 allChannels(jj) = singleChannel;
%             end
%             rawAnalogData(analogRow,:) = allChannels;
%         end
%     catch
%     end
% end

end

