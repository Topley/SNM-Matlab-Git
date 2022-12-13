function [AbdAddAngle, IntExtRotAngle, FlexExtAngle] = EulerAngles(proximalRMat,distalRMat)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

rotatedMat = pagemtimes(proximalRMat, 'transpose', distalRMat, 'none');

ind1 = rotatedMat(1,2,:);

ind2 = rotatedMat(2,2,:);

alphaR = atan2(-ind1,ind2);
%alphaX = permute(alphaR, [3,2,1]);

ind3 = rotatedMat(3,2,:);

b1 = ind1 .^2;

b2 = ind2 .^2;

ind4 = sqrt((b1 + b2));

betaR = atan2(ind3, ind4);
%betaX = permute(betaR(:,:,1:63317), [3,2,1]);

ind5 = rotatedMat(3,1,:);

ind6 = rotatedMat(3,3,:);

gammaR = atan2(-ind5, ind6);
%gammaX = permute(gammaR(:,:,1:63317), [3,2,1]);

radMat = [betaR, gammaR, alphaR];

tRotMat = reshape(pagetranspose(radMat),[], 3);
%tRotMat = squeeze(radMat);

radangles = tRotMat ./ pi;

angles = radangles .* 180;

[b,a] = butter(2,8/1024, 'low');

AbdAddAngle = filtfilt(b, a, angles(:,1));
IntExtRotAngle = filtfilt(b, a, angles(:,2));
FlexExtAngle = filtfilt(b, a, angles(:,3));
%plot(angles(:,1))

end
