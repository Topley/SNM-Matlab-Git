function [binPosition,curParamsData] = extractParamData(fileID,param3D, numDim, dataType)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for j = 1:3
    pidx = param3D(j);
    idx = j+1;
    if idx > numDim && pidx == 0
        param3D(j) = 1;
    end
end
flip3DParam = flip(param3D);

curParamsData = {};
for ij = 1:flip3DParam(1)
    for jj = 1:flip3DParam(2)
        
        switch dataType
            case -1
                loopParamData = fread(fileID, flip3DParam(3),'int8', 'ieee-le');
                curParamsData{ij,jj} = string(char(loopParamData'));
            case 1
                loopParamData = fread(fileID, flip3DParam(3), 'int8', 'ieee-le');
                loopParamData = sprintf('%0.6f',loopParamData(1));
                curParamsData{ij,jj} = loopParamData;
            case 2
                loopParamData = fread(fileID, flip3DParam(3), 'int16', 'ieee-le');
                loopParamData = sprintf('%0.6f',loopParamData(1));
                curParamsData{ij,jj} = loopParamData;
            case 4
                loopParamData = fread(fileID, flip3DParam(3), 'float32', 'ieee-le');
                loopParamData = sprintf('%0.6f',loopParamData(1));
                curParamsData{ij,jj} = loopParamData;
        end
        
    end
end
 binPosition = ftell(fileID);
end

