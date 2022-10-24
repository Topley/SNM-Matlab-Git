function [c3dStruct] = Get_C3dParameters(fileID, c3dStruct, paramStartBlock)

fseek(fileID, paramStartBlock,'bof');
c3dParams1 = fread(fileID, 4, 'int8', 'ieee-le');
firstBreak = ftell(fileID);
c3dStruct.Parameters.ParamBlockNumber = c3dParams1(3);
c3dStruct.Parameters.ProcessorType = c3dParams1(4)-83;


blockLoop = firstBreak;
stopLoop = 1;
c3dStruct.Parameters.GroupID = [];
c3dStruct.Parameters.NumberDimensions = {};
c3dStruct.Parameters.DimensionLength = {};
c3dStruct.Parameters.DataType = [];
c3dStruct.Parameters.GroupParamName = {};
c3dStruct.Parameters.NextGroupParam = [];
c3dStruct.Parameters.DescriptionSize = [];
c3dStruct.Parameters.GroupParamDescription = [];
paramDataMatrix = {};
loop = 0;
while stopLoop ~= 0
    
    loop = loop + 1;
%     if loop > 127
%         keyboard
%     end
    fseek(fileID, blockLoop,'bof');
    groupInfo = fread(fileID, 2, 'int8', 'ieee-le');
    
    c3dStruct.Parameters.GroupID = [c3dStruct.Parameters.GroupID;groupInfo(2)];
    stopLoop = groupInfo(2);
    
    binaryParamName = fread(fileID, abs(groupInfo(1)), 'int8', 'ieee-le');
    groupParamName = char(binaryParamName');
    c3dStruct.Parameters.GroupParamName{loop} = string(groupParamName);
    
    nextGroup = fread(fileID, 1, 'int16', 'ieee-le');
    c3dStruct.Parameters.NextGroupParam = [c3dStruct.Parameters.NextGroupParam; nextGroup];
    offsetBlock = ftell(fileID);
    nextBlock = (nextGroup + offsetBlock) - 2;
    
    if stopLoop < 0
        c3dStruct.Parameters.NumberDimensions{loop} = NaN;
        c3dStruct.Parameters.DimensionLength = {c3dStruct.Parameters.DimensionLength;NaN};
        c3dStruct.Parameters.DataType = [c3dStruct.Parameters.DataType;NaN];
        curParamsData = sprintf("Group");
    elseif stopLoop > 0
        dimData = fread(fileID, 2, 'int8', 'ieee-le');
        dataType = dimData(1);
        numDim = dimData(2);
        c3dStruct.Parameters.DataType = [c3dStruct.Parameters.DataType ; double(dataType)];
   
        if numDim > 3
            paramStartDim = fread(fileID, numDim, 'int8', 'ieee-le');
            disp('Dimensions are too high - check Parameter dimension number')
            param3D = paramStartDim(1:3);
        elseif numDim == 3
            paramStartDim = fread(fileID, numDim, 'int8', 'ieee-le');
            param3D = paramStartDim;
        elseif numDim == 2
            paramStartDim = fread(fileID, numDim, 'int8', 'ieee-le');
            param3D = [paramStartDim;0];
        elseif numDim == 1
            paramStartDim = fread(fileID, numDim, 'int8', 'ieee-le');
            param3D = [paramStartDim;0;0];
        elseif numDim == 0
            paramStartDim = fread(fileID, numDim, 'int8', 'ieee-le');
            %numDim = 1;
            param3D = [1,0,0];
        else
        end
        c3dStruct.Parameters.NumberDimensions{loop} = double(numDim);
        c3dStruct.Parameters.DimensionLength = param3D;
        
        holdBlock = ftell(fileID);
        %%%% Function
        [ParamBinPosition,curParamsData] = extractParamData(fileID,param3D, numDim, dataType);
        paramDataMatrix{loop} = curParamsData;
        
     elseif stopLoop == 0
          dimData = fread(fileID, 2, 'int8', 'ieee-le');
        dataType = dimData(1);
        numDim = 0;
        c3dStruct.Parameters.DataType = [c3dStruct.Parameters.DataType ; double(dataType)];
        
          curParamsData = fread(fileID, 1, 'int8', 'ieee-le');
          paramDataMatrix{loop} = curParamsData;
        c3dStruct.Parameters.NumberDimensions{loop} = NaN;
        c3dStruct.Parameters.DimensionLength = {c3dStruct.Parameters.DimensionLength;NaN};
        c3dStruct.Parameters.DataType = [c3dStruct.Parameters.DataType;NaN];
    else
    end
    
 if stopLoop ~= 0
    descriptSize = fread(fileID, 1, 'int8', 'ieee-le');
    c3dStruct.Parameters.DescriptionSize = descriptSize(1);
    groupParamDescript = fread(fileID, descriptSize(1), 'int8', 'ieee-le');
    c3dStruct.Parameters.GroupParamDescription{loop} = char(groupParamDescript');
    blockLoop = nextBlock;
 end 
 
end

for k = 1:size(c3dStruct.Parameters.GroupParamName, 2)
    c3dStruct.Parameters.ParameterName{1,k} = c3dStruct.Parameters.GroupParamName{k};
    
    paramRows = size(paramDataMatrix{k}(:), 1);
    paramCols = size(paramDataMatrix{k}(:), 2);
    
    if paramRows == 1 && paramCols > 1
        paramValue = paramDataMatrix{k}(:);
        c3dStruct.Parameters.ParameterName(2:paramCols+1,k) = paramValue;
    elseif paramRows > 1 && paramCols == 1
        paramValue = paramDataMatrix{k}(:)';
        c3dStruct.Parameters.ParameterName(2:paramRows+1,k) = paramValue;
    else
        c3dStruct.Parameters.ParameterName(2:paramRows+1,k) = paramDataMatrix(k);
        
    end
end

end

