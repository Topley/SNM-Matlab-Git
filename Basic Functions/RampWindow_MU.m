function [MUFiring, AVERAGE,firing, TRUEMU] = RampWindow_MU(MUTime, TRUEMU, MUFiring)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
removetrains=[];
loopMuF = {};

try
for i = 1:length(MUFiring)
    loopMuF{i} = MUFiring{i}(MUFiring{i}>MUTime(1).*fsamp & MUFiring{i}<MUTime(2).*fsamp);
end

for jj = 1:length(loopMuF)
    try
        removetrains(jj) = (loopMuF{jj}(end))/fsamp - (loopMuF{jj}(1)/fsamp);
    catch
        removetrains(jj) = 1;
    end
end
TRUEMU(removetrains<4) = [];
loopMuF(removetrains<4) = [];


for k = 1:length(loopMuF)
    splitUnit = round(length(loopMuF{k})./4);
    bigPause(k) = max(diff(loopMuF{k}(splitUnit:end-splitUnit))./fsamp);
end 
loopMuF(bigPause > 1) = [];
TRUEMU(bigPause > 1) = [];
%%% old remove trains arbitrary limits
% removetrains = cellfun(@length, loopMuF);
% TRUEMU(removetrains<10) = [];
% loopMuF(removetrains<10) = [];

%%% This is where the motor units should have their edges cleaned
for i = 1:length(loopMuF) % Create cell array of MU datapoints in seconds
    
    MuF = loopMuF{i};
    preTSAISI = (diff(MuF/fsamp));
    splitUnit = round(length(preTSAISI)./4);
    
    clear r rr;
    try
        %[rr r] = find(preTSAISI(1:splitUnit)>=0.4)
        [rr r] = find(preTSAISI(1:splitUnit)>0.6); % cut-off from Gorasinni 2004
    catch
        [rr r] = find(preTSAISI(1:10)>0.6);% Strict arbitrary cut-off
    end
    if isempty(rr);
        preMUPulses = MuF;
    else
        y = max(r);
        preMUPulses = MuF(y+1:end);
    end
    
    clear h hh
    try
        [hh h] = find(preTSAISI(end-splitUnit:end)>0.6);% cut-off from Gorasinni 2004
    catch
        [hh h] = find(preTSAISI(end-10:end)>0.6);% Strict arbitrary cut-off
    end
    
    if isempty(hh);
        preMUPulses = preMUPulses;
    else
        z = min(h);
        preMUPulses = preMUPulses(1:end-splitUnit+z-2);
    end
    
    MuF = preMUPulses;
    loopMUFiring{i} = MuF/fsamp;
end


TRUEMU = [TRUEMU(~cellfun(@isempty, loopMUFiring)),0];
MUFiring = loopMUFiring(~cellfun(@isempty, loopMUFiring));

%%% Remove spike trains firing for less than 4s (allows for 2s up and 2s down)
removetrains=[];
loopMuF = {};
for i = 1:length(MUFiring)
    loopMuF{i} = MUFiring{i}(MUFiring{i}>MUTime(1) & MUFiring{i}<MUTime(2));
end

for jj = 1:length(loopMuF)
    try
        removetrains(jj) = (loopMuF{jj}(end)) - (loopMuF{jj}(1));
    catch
        removetrains(jj) = 1;
    end
end
TRUEMU(removetrains<4) = [];
MUFiring(removetrains<4) = [];

for k = 1:length(MUFiring)
    bigPause(k) = max(diff(MUFiring{k}(2:end-1)));
end 
MUFiring(bigPause > .5) = [];
TRUEMU(bigPause > .5) = [];
%%% This is where the units are filtered
for i = 1:length(MUFiring)
    MuF = MUFiring{i};
    index = round(MuF*fsamp);
    firing(i,index) = 1;
end

firing(:,length(firing)+fsamp)=0; %adds a second of 0s to the end
AVERAGE = 1./sum(hanning(fsamp*2))*filtfilt(hanning(fsamp*2),1,firing')';%AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);


catch
    disp('Cannot trim Motor unit edges');
end

end

