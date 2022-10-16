function [rawMarkerData,rawAnalogData] = Extract_C3DData(fileID, dataStartBin, totalVideoFrames, MarkerNumber, AVRatio, AnalogNumber)

%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
firstBin = (dataStartBin - 1) * 512;

fseek(fileID, firstBin,'bof');
rawMarkerData = zeros(4,totalVideoFrames,MarkerNumber);
AnalogDataPoints = totalVideoFrames * AVRatio;
rawAnalogData = zeros(AnalogDataPoints, AnalogNumber);
analogRow = 0;

for i = 1:totalVideoFrames
    
    rawMarkers = zeros(4,1, MarkerNumber);
    for j = 1:MarkerNumber
        xyzCoords = fread(fileID, 4, 'float32', 'ieee-le');
        rawMarkers(:,:,j) = xyzCoords;
    end
    rawMarkerData(:,i,:) = rawMarkers;
    
    
    for ii = 1:AVRatio
        analogRow = analogRow + 1;
        allChannels = zeros(1,15);
        for jj = 1:AnalogNumber
            singleChannel = fread(fileID, 1, 'float32', 'ieee-le');
            allChannels(jj) = singleChannel;
        end
        rawAnalogData(analogRow,:) = allChannels;
    end
    
end

end

