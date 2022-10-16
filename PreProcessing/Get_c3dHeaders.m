function [c3dStruct] = Get_c3dHeaders(fileID)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

c3dHeader1 = fread(fileID, 2, 'int8', 'ieee-le');
c3dStruct.Headers.ParameterStartBlock = c3dHeader1(1);

c3dHeader2 = fread(fileID, 5, 'int16', 'ieee-le');
c3dStruct.Headers.NumberOfTrajectories = c3dHeader2(1);
c3dStruct.Headers.FirstFrame = c3dHeader2(3);
c3dStruct.Headers.LastFrame = c3dHeader2(4);
c3dStruct.Headers.MaxInterpolateGap = c3dHeader2(5);

c3dStruct.Headers.ScaleFactor = fread(fileID,1, 'float32', 'ieee-le');
c3dHeader3 = fread(fileID, 2,'int16', 'ieee-le');
c3dStruct.Headers.DataStartBlock = c3dHeader3(1);
c3dStruct.Headers.AVRatio = c3dHeader3(2);

c3dStruct.Headers.AnalogChannels = c3dHeader2(2)/c3dStruct.Headers.AVRatio;

c3dStruct.Headers.VideoFrameRate = fread(fileID, 1,'float32', 'ieee-le');

end

