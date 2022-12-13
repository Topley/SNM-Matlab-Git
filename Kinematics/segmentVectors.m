function [segmentVect1] = segmentVectors(trialStruct, jc1, jc2)

%% Hip vector
jc1Size = size(trialStruct(1).(jc1), 2);
jc2Size = size(trialStruct(1).(jc2), 2);

cols = [jc1Size,jc2Size];


if cols(1) == 4 && cols(2) == 4
    cols1 = 2:4;
    cols2 = 2:4;
elseif cols(1) == 4 && cols(2) == 3
    cols1 = 2:4;
    cols2 = 1:3;
elseif cols(1) == 3 && cols(2) == 4
    cols1 = 1:3;
    cols2 = 2:4;
elseif cols(1) == 3 && cols(2) == 3
    cols1 = 1:3;
    cols2 = 1:3;
else
    disp('missing columns')
end


vect1 = trialStruct.(jc1)(:,cols1) - trialStruct.(jc2)(:,cols2);
segmentVect1 = normr(vect1);
%segmentVect1isnan(segmentVect1)

end



