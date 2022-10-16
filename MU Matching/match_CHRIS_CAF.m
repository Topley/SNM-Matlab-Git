function MATCHED_MU = match_CHRIS_CAF(fullFilename,SIG1,MUTrain1,SIG2,MUTrain2,L,fsamp,DMODE,TH,NS,FHIGH,FLOW,Finterf,THD)
% MATCH Find the same Motor Units across different recordings
%
% MATCHED_MU = match(EMG1,MUTrain1,EMG2,MUTrain2,L,fsamp,DMODE,TH,NS,FHIGH,FLOW,Finterf)
%
% MATCHED_MU provides a metric to track motor unit discharges across
% trials. After user defined filtering and differentiation, the array 
% waveforms are estimated through spike triggered averaging for each 
% of the 2 trials. Two dimensional crosscorelation is used to estimate the
% best fit and lag between trials. The lag is adjusted by moving the second trial window.
% Waveforms are plotted and overlaid. The optimal match for each of MUTrain
% 
% 
% With monopolar settings, a 3D model of the motor unit action potential is
% created **NOTE the current electrode configuration is numbered 
% sequentially as in **
%
% [SIG1]: matrix of multichannel EMG data from first recording
% [MUTrain1]: array of cells corresponding to the motor unit spike times from first recording
% [SIG2]: matrix of multichannel EMG data from second recording
% [MUTrain2]: array of cells corresponding to the motor unit spike times from second recording
% [L]: 1/2 of the length of STA window
% [fsamp]: sampling frequency
% [DMODE]: post-processing differentiation (0=mono, 1=single, 2=double)
% [TH]: threshold for significance of the 2D crosscorelation
% [NS]: minimum number of samples for STA
% [FHIGH]: post-processing high pass filter
% [FLOW]: post-processing low pass filter
% [Finterf]: harmonic line filter 
%
% Current usege
% load 3612C0Ch130510160101.mat; % Raw array data from trial 1 saved as 'EMG'
% EMG1=EMG;
% load 3612C0Ch130510160109.mat; % Raw array data from trial 2 saved as 'EMG'
% EMG2=EMG; 
% load C006_C006_60101_soleus_1_decomposed.mat; % Spike times from trial 1 saved as 'MUFiring'
% MUTrain1 = MUFiring; 
% load C006_C006_60109_soleus_1_decomposed.mat; % Spike times from trial 2 saved as 'MUFiring'
% MUTrain2 = MUFiring;  
% % MATCHED_MU = match(SIG1,MUTrain1,SIG2,MUTrain2,L,fsamp,DMODE,TH,NS,FHIGH,FLOW,Finterf)
% MU_MATCHED = match(SIG1,MUTrain1,SIG2,MUTrain2,25e-3,5120,1,0.85,50,100,500,[60]);
%
% Author: F. Negro (francesco.negro@unibs.it)
% Description and Debugging: C. Thompson (chriskthompson@yahoo.com),
% Debugging: Eduardo Martinez Valdes  mredumartinez@gmail.com
% Debugging: Alessandro Del Vecchio XXXXXXXX@XXXX  

% Revision: 0.1  Date: 2013/07/25  
% Revision: 0.11 Date: 2013/08/05
% Revision: 1.0  Date: 2014/06/15
% Revision 2.0   Date: 2016/01/19   


% error codes
%if (size(EMG1,1)>size(EMG1,2))||(size(EMG2,1)>size(EMG2,2)),
%    error(generatemsgid('InvalidRange'),'EMG signals should have channels in the first dimention.')
%end    

%if size(EMG1,1)~=size(EMG2,1)
%    error(generatemsgid('InvalidRange'),'EMG signals don''t have the same number of channels.')
%end

[fileDir, filename] = fileparts(fullFilename);

try
    pdfdir = [fileDir '\matched_pdf'];
    if exist(pdfdir,'dir') == 0;
        mkdir([fileDir '\matched_pdf']);
    end
end

if L>1
    error(generatemsgid('InvalidRange'),'L should be in seconds.')
end    

%if isempty(intersect(DMODE,[0 1 2])),
%    error(generatemsgid('InvalidRange'),'Wrong DIFF mode.')
%end    

if TH>1 || TH<0,
    error(generatemsgid('InvalidRange'),'TH should be between 0 and 1.')
end    

if FHIGH>FLOW || FHIGH>fsamp/2 || FLOW>fsamp/2,
    error(generatemsgid('InvalidRange'),'Wrong filtering settings.')
end    

L = round(fsamp*L/2); % Conversion to samples

[B,A] = butter(4,[FHIGH,FLOW]*2/fsamp);

Fsamp=fsamp;
sigLength=length(SIG1{2,2});
Frad=round(Fsamp/200);
fInt=round(0:Finterf*sigLength/Fsamp:(sigLength/2));
% Apply line filter to EMG1
if ~isempty(Finterf)
    for R = 1:size(SIG1,1),
        for C = 1:size(SIG1,2),
            if ~isempty(SIG1{R,C}),
                SIG1{R,C} = filtfilt(B,A,SIG1{R,C});
                SIG1{R,C} = real(subtractFreq2(SIG1{R,C},fInt,Frad,0));
            end
        end
    end        
end

%disp('Before Sorting First Recording\n');
%size(MUTrain2,2)
%MUTrain1 = find_doubles(MUTrain1,fsamp,THD);
%disp('After Sorting First Recording\n');
%size(MUTrain2,2)
MUFiring = MUTrain1; % Call spike times


% Spike trigger average EMG1
STA1.data = cell(size(SIG1));
[STA1.data{:}] = deal(zeros(2*L,size(MUFiring,2)));
for MU = 1:size(MUFiring,2)
    index = MUFiring{MU};    
    for R = 1:size(STA1.data,1)
       for C = 1:size(STA1.data,2)
             if ~isempty(SIG1{R,C})
                 p=0;
                 for trigger = 1:length(index)
                        if index(trigger)>L && index(trigger)<length(SIG1{R,C})-L
                                STA1.data{R,C}(:,MU) = STA1.data{R,C}(:,MU) + SIG1{R,C}(index(trigger)-L+1:index(trigger)+L)';
                            p=p+1;         
                        end
                 end
                 STA1.data{R,C}(:,MU) = STA1.data{R,C}(:,MU)/p;
             end
        end
    end
    STA1.TRIGGERS(MU) = p;
    DR = 1./diff(index/fsamp);
    STA1.DR(MU) = mean(DR(DR>4 & DR<50));
end

Fsamp=fsamp;
sigLength=length(SIG2{2,2});
Frad=round(Fsamp/200);
fInt=round(0:Finterf*sigLength/Fsamp:(sigLength/2));

% Apply line filter to EMG2
if ~isempty(Finterf)
    for R = 1:size(SIG2,1),
        for C = 1:size(SIG2,2),
            if ~isempty(SIG2{R,C}),
                SIG2{R,C} = filtfilt(B,A,SIG2{R,C});
                SIG2{R,C} = real(subtractFreq2(SIG2{R,C},fInt,Frad,0));
            end
        end
    end       
end
%disp('Before Sorting Second Recording\n');
%size(MUTrain2,2)
%MUTrain2 = find_doubles(MUTrain2,fsamp,THD);
%disp('After Sorting Second Recording\n');
%size(MUTrain2,2)
MUFiring = MUTrain2;

% Differentiate
%if DMODE>0,
%EMG = diff(EMG,DMODE,1);
%end

% Spike trigger average SIG2
STA2.data = cell(size(SIG2));
[STA2.data{:}] = deal(zeros(2*L,size(MUFiring,2)));
for MU = 1:size(MUFiring,2)
    index = MUFiring{MU};    
    for R = 1:size(STA2.data,1)
    for C = 1:size(STA2.data,2)
       if ~isempty(SIG2{R,C})
            p=0;
            for trigger = 1:length(index)
                if index(trigger)>L && index(trigger)<length(SIG2{R,C})-L
                        STA2.data{R,C}(:,MU) = STA2.data{R,C}(:,MU) + SIG2{R,C}(index(trigger)-L+1:index(trigger)+L)';
                    p=p+1;    
                end   
            end
            STA2.data{R,C}(:,MU) = STA2.data{R,C}(:,MU)/p;
        end
    end
    end
    STA2.TRIGGERS(MU) = p;
    DR = 1./diff(index/fsamp);
    STA2.DR(MU) = mean(DR(DR>4 & DR<50));
end

if isempty(strfind(DMODE,'MONO')),
    STA1.data = spatial_filter(STA1.data,DMODE);
    STA2.data = spatial_filter(STA2.data,DMODE);
end

AVERAGE1.data = cat(3,STA1.data{:});
AVERAGE1.data = shiftdim(AVERAGE1.data,2);
AVERAGE1.TRIGGERS = STA1.TRIGGERS;
AVERAGE1.DR = STA1.DR;
AVERAGE2.data = cat(3,STA2.data{:});
AVERAGE2.data = shiftdim(AVERAGE2.data,2);
AVERAGE2.TRIGGERS = STA2.TRIGGERS;
AVERAGE2.DR = STA2.DR;

%SCALING = max([max(max(max(AVERAGE1.data))),max(max(max(AVERAGE2.data)))]);

SCALING = 500;


cprintf('Black','\nMU1\tMU2\tCORR\t DIST\tMATCH\tNT1\tNT2\n');

% Begin crosscorrelation of spike times
MATCHED_MU =[];
MAXCORRE = [];
MAXDELAY = [];
MAXMUMATCH = [];

NMAT = 0;
NUMAT = 0;
ACUMAT = 0;
ACMAT = 0;

for MU1 = 1:size(AVERAGE1.data,3),
    MAXCORRE(MU1) = eps;

    for MU2 = 1:size(AVERAGE2.data,3),

        CORRE = xcorr2(AVERAGE1.data(:,:,MU1),AVERAGE2.data(:,:,MU2));

        NCORRE = max(max(CORRE,[],2))./sqrt(sum(dot(AVERAGE1.data(:,:,MU1),AVERAGE1.data(:,:,MU1)))*sum(dot(AVERAGE2.data(:,:,MU2),AVERAGE2.data(:,:,MU2))));
        
        [massimo,index] = max(CORRE(:));
        [indr,indc] = ind2sub(size(CORRE),index);
        
        AVERAGE3 = circshift(AVERAGE2.data(:,:,MU2),[0 round(size(AVERAGE2.data(:,:,MU2),2)+indc)]);
        
        DISTANCE = sqrt(sum(sum((AVERAGE1.data(:,:,MU1) - AVERAGE3).^2)))/(norm(AVERAGE1.data(:,:,MU1))+norm(AVERAGE3));
        DISTANCE = 1-DISTANCE;

        %[MU1 MU2 NCORRE DISTANCE]
        
        if NCORRE > TH,
            ACMAT = ACMAT + NCORRE;
            NMAT = NMAT + 1;
        else
            ACUMAT = ACUMAT + NCORRE;
            NUMAT = NUMAT + 1;
        end    
        
        if NCORRE>MAXCORRE(MU1),
            MAXCORRE(MU1) = NCORRE;
            
            MAXDELAY(MU1) = indc;
            MAXMUMATCH(MU1) = MU2;
            MAXDIST(MU1) = DISTANCE;
        end
    end
      
        % Plot 2D motor unit action potentials using show_MUAP function
        figure,
        show_MUAP(filename,STA1,MU1,0,fsamp,4*L,SCALING,'k',MUTrain1,MU1,MU1+100,1);
        hold on,
        axis off;
        %pause%%%%%
        show_MUAP(filename,STA2,MAXMUMATCH(MU1),MAXDELAY(MU1),fsamp,4*L,SCALING,'r',MUTrain2,MU1,MU1+100,2);
        axis off;
        figure(MU1),
        title(num2str([MU1 MAXMUMATCH(MU1) MAXCORRE(MU1) MAXDIST(MU1)]))%%%%%%
        % Print output
        MATCHED_MU(MU1,:) = [MU1 MAXMUMATCH(MU1) MAXCORRE(MU1) MAXDIST(MU1) AVERAGE1.TRIGGERS(MU1) AVERAGE2.TRIGGERS(MAXMUMATCH(MU1)) AVERAGE1.DR(MU1) AVERAGE2.DR(MAXMUMATCH(MU1))];
        saveas(figure(MU1),fullfile(pdfdir,[filename(1:end-4),'_MATCHEDunit_',num2str(MU1),'.pdf']))
  
end   

clc;

for MU = 1:size(MATCHED_MU,1),
        if MATCHED_MU(MU,3)>=TH,
            if MATCHED_MU(MU,5)<NS || MATCHED_MU(MU,6)<NS,
                cprintf('Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'--X--',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
            else
                index = find(MATCHED_MU(:,2)==MATCHED_MU(MU,2));
                if length(index) == 1,    
                    cprintf('Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'--*--',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
                else
                    [minimo,ind] = min(MATCHED_MU(index,4));
                    TEMP = MATCHED_MU(index,:);
                    TEMP = TEMP(ind,:);
                    if MU == TEMP(:,1),
                        cprintf('_Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'--*--',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
                    else
                        cprintf('Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'--*--',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
                    end
                end    
            end
        else
            if MATCHED_MU(MU,5)<NS || MATCHED_MU(MU,6)<NS,
                cprintf('Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'-----',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
            else
                cprintf('Text','%d\t%d\t%1.2f\t%1.2f\t%s\t%d\t%d\n',MATCHED_MU(MU,1),MATCHED_MU(MU,2),MATCHED_MU(MU,3),MATCHED_MU(MU,4),'-----',MATCHED_MU(MU,5),MATCHED_MU(MU,6));
            end
        end
end 
disp(['Average Correlation between MATCHED units = ',num2str(ACMAT/NMAT)]);
disp(['Number of MATCHED units = ',num2str(NMAT)]);
disp(['Average Correlation between UNMATCHED units = ',num2str(ACUMAT/NUMAT)]);
disp(['Number of MATCHED units = ',num2str(NUMAT)]);

save(fullfile(fileDir,[filename(1:end-4),'_MATCHED_MU']),'MATCHED_MU')
end

% Plot 2D motor unit action potentials
function show_MUAP(filename,STA,MU,SHI,fsamp,SCALE1,SCALE2,color,MUFiring,NF1,NF2,NF3)%,filename)
%figure(NF1);
%subplot(1,2,1),
MUAP1 = [];
for R = 1:size(STA.data,1)
    for C = 1:size(STA.data,2)
        hold on,plot([0:length(STA.data{R,C}(:,MU))-1]/fsamp+C*SCALE1/fsamp,circshift(STA.data{R,C}(:,MU),SHI,1)'-R*SCALE2,color);    
        TEMP = circshift(STA.data{R,C}(:,MU),SHI,1);
        for time = 1:length(TEMP)
        MUAP(R,C,time) = TEMP(time);
        end
    end
end
axis off;
title(num2str(MU))
%keyboard
%saveas(gcf,[filename(1:end-4),'_MATCHEDunit_',num2str(MU),'.pdf'])
%save([filename(1:end-4),'_unit_MATCHEDSTA.mat'],'filename','STA','MU','SHI','fsamp','SCALE1','SCALE2','color','MUFiring','NF1','NF2','NF3')

%subplot(2,1,2),
%hold on,plot(MUFiring{MU}(2:end)/fsamp,1./(diff(MUFiring{MU})/fsamp),[color,'*']);xlim([0 max(MUFiring{MU}/fsamp)]);
%ylim([0 50]);
%subplot(2,1,1);

%figure(NF1);
%plot 3D
% 
% subplot(1,2,2);
% MUAP_plane=zeros(5,12,'double');
% [XI,YI] = meshgrid(1:0.05:12,1:0.05:5);
% [Xi,Yi] = meshgrid(1:12,1:5);
% ZI1 = interp2(1:12,1:5,squeeze(MUAP(:,:,end/2))',XI,YI,'cubic');
% mesh(XI,YI,ZI1,'FaceAlpha',0.6); hold on; mesh(Xi,Yi,MUAP_plane,'FaceAlpha',0,'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
% colormap(pink);
% shading interp ;
% axis off;

% figure;
% c = 15
% for time = 1:round(length(MUAP)/5):length(MUAP)
% [XI,YI] = meshgrid(1:0.05:12,1:0.05:5);
% [Xi,Yi] = meshgrid(1:12,1:5);
% ZI1 = interp2(1:12,1:5,squeeze(MUAP(:,:,time))',XI,YI,'cubic');
% ZI1 = ZI1/SCALE2;
% hold on
% image([1:0.05:12]+c,1:0.05:5,ZI1,'CDataMapping','scaled');
% c = c+15;
% %colormap(pink);
% %shading interp ;
% %axis off;
% end
% colormap(jet);

end

% Line filter
function ynew=subtractFreq2(y,Ind,r,draw)
% subtracts the content at the frequencies [Ind-r:Ind+r] from the signal y. With paramter draw>0 it also plots
% the frequency content of original and filtered signal.
% OUTPUT
%   ynew - filtered signal

fy=fft(y);
lenFy=length(fy);
fd=zeros(1,lenFy);  
tInd=[];
% for k=-floor(r/2):floor(r/2)
%     tInd=[tInd Ind+k];
% end

% mFreq=median(fy);
% sFreq=std(fy);
% t2Ind=find(abs(fy) > mFreq+15*sFreq);
% for k=-floor(r/2):floor(r/2)
%     tInd=[tInd t2Ind+k];
% end

WinLen=1000;
for k1=1:WinLen:length(fy)-WinLen;
    mFreq=median(abs(fy(k1+1:k1+WinLen)));
    sFreq=std(abs(fy(k1+1:k1+WinLen)));
    t2Ind=find(abs(fy(k1+1:k1+WinLen)) > mFreq+5*sFreq);
    t2Ind=t2Ind+k1;
    for k=-floor(r/2):floor(r/2)
        tInd=[tInd t2Ind+k];
    end
end


tInd=round(tInd(find(tInd>0 & tInd<=lenFy/2+1)));
for k=tInd
    fd(k)=fy(k);    
end

correct=lenFy-floor(lenFy/2)*2;
fd(lenFy:-1:ceil(lenFy/2)+1)=conj(fd(2:1:ceil(lenFy/2)+1-correct));

ynew=y-ifft(fd(1:lenFy));

if draw>0
fy=fftshift(fy);
fd=fftshift(fd);     
    figure(draw); 
    subplot(2,1,1); plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,real(fy)); hold on; plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,real(fy-fd),'r'); hold off; legend('original signal','filtered signal'); title('real');     
    subplot(2,1,2); plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,imag(fy)); hold on; plot(-ceil(lenFy/2)+1:ceil(lenFy/2)-correct,imag(fy-fd),'r'); hold off; legend('original signal','filtered signal'); title('imag'); 
    %pause;
end
end

% Text colorization
function count = cprintf(style,format,varargin)
% CPRINTF displays styled formatted text in the Command Window
%
% Syntax:
%    count = cprintf(style,format,...)
%
% Description:
%    CPRINTF processes the specified text using the exact same FORMAT
%    arguments accepted by the built-in SPRINTF and FPRINTF functions.
%
%    CPRINTF then displays the text in the Command Window using the
%    specified STYLE argument. The accepted styles are those used for
%    Matlab's syntax highlighting (see: File / Preferences / Colors / 
%    M-file Syntax Highlighting Colors), and also user-defined colors.
%
%    The possible pre-defined STYLE names are:
%
%       'Text'                 - default: black
%       'Keywords'             - default: blue
%       'Comments'             - default: green
%       'Strings'              - default: purple
%       'UnterminatedStrings'  - default: dark red
%       'SystemCommands'       - default: orange
%       'Errors'               - default: light red
%       'Hyperlinks'           - default: underlined blue
%
%       'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White'
%
%    STYLE beginning with '-' or '_' will be underlined. For example:
%          '-Blue' is underlined blue, like 'Hyperlinks';
%          '_Comments' is underlined green etc.
%
%    STYLE beginning with '*' will be bold (R2011b+ only). For example:
%          '*Blue' is bold blue;
%          '*Comments' is bold green etc.
%    Note: Matlab does not currently support both bold and underline,
%          only one of them can be used in a single cprintf command. But of
%          course bold and underline can be mixed by using separate commands.
%
%    STYLE also accepts a regular Matlab RGB vector, that can be underlined
%    and bolded: -[0,1,1] means underlined cyan, '*[1,0,0]' is bold red.
%
%    STYLE is case-insensitive and accepts unique partial strings just
%    like handle property names.
%
%    CPRINTF by itself, without any input parameters, displays a demo
%
% Example:
%    cprintf;   % displays the demo
%    cprintf('text',   'regular black text');
%    cprintf('hyper',  'followed %s','by');
%    cprintf('key',    '%d colored', 4);
%    cprintf('-comment','& underlined');
%    cprintf('err',    'elements\n');
%    cprintf('cyan',   'cyan');
%    cprintf('_green', 'underlined green');
%    cprintf(-[1,0,1], 'underlined magenta');
%    cprintf([1,0.5,0],'and multi-\nline orange\n');
%    cprintf('*blue',  'and *bold* (R2011b+ only)\n');
%    cprintf('string');  % same as fprintf('string') and cprintf('text','string')
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
%    A technical description of the implementation can be found at:
%    <a href="http://undocumentedmatlab.com/blog/cprintf/">http://UndocumentedMatlab.com/blog/cprintf/</a>
%
% Limitations:
%    1. In R2011a and earlier, a single space char is inserted at the
%       beginning of each CPRINTF text segment (this is ok in R2011b+).
%
%    2. In R2011a and earlier, consecutive differently-colored multi-line
%       CPRINTFs sometimes display incorrectly on the bottom line.
%       As far as I could tell this is due to a Matlab bug. Examples:
%         >> cprintf('-str','under\nline'); cprintf('err','red\n'); % hidden 'red', unhidden '_'
%         >> cprintf('str','regu\nlar'); cprintf('err','red\n'); % underline red (not purple) 'lar'
%
%    3. Sometimes, non newline ('\n')-terminated segments display unstyled
%       (black) when the command prompt chevron ('>>') regains focus on the
%       continuation of that line (I can't pinpoint when this happens). 
%       To fix this, simply newline-terminate all command-prompt messages.
%
%    4. In R2011b and later, the above errors appear to be fixed. However,
%       the last character of an underlined segment is not underlined for
%       some unknown reason (add an extra space character to make it look better)
%
%    5. In old Matlab versions (e.g., Matlab 7.1 R14), multi-line styles
%       only affect the first line. Single-line styles work as expected.
%       R14 also appends a single space after underlined segments.
%
%    6. Bold style is only supported on R2011b+, and cannot also be underlined.
%
% Change log:
%    2012-08-09: Graceful degradation support for deployed (compiled) and non-desktop applications; minor bug fixes
%    2012-08-06: Fixes for R2012b; added bold style; accept RGB string (non-numeric) style
%    2011-11-27: Fixes for R2011b
%    2011-08-29: Fix by Danilo (FEX comment) for non-default text colors
%    2011-03-04: Performance improvement
%    2010-06-27: Fix for R2010a/b; fixed edge case reported by Sharron; CPRINTF with no args runs the demo
%    2009-09-28: Fixed edge-case problem reported by Swagat K
%    2009-05-28: corrected nargout behavior sugegsted by Andreas G?b
%    2009-05-13: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%
% See also:
%    sprintf, fprintf

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.08 $  $Date: 2012/10/17 21:41:09 $

  persistent majorVersion minorVersion
  if isempty(majorVersion)
      %v = version; if str2double(v(1:3)) <= 7.1
      %majorVersion = str2double(regexprep(version,'^(\d+).*','$1'));
      %minorVersion = str2double(regexprep(version,'^\d+\.(\d+).*','$1'));
      %[a,b,c,d,versionIdStrs]=regexp(version,'^(\d+)\.(\d+).*');  %#ok unused
      v = sscanf(version, '%d.', 2);
      majorVersion = v(1); %str2double(versionIdStrs{1}{1});
      minorVersion = v(2); %str2double(versionIdStrs{1}{2});
  end

  % The following is for debug use only:
  %global docElement txt el
  if ~exist('el','var') || isempty(el),  el=handle([]);  end  %#ok mlint short-circuit error ("used before defined")
  if nargin<1, showDemo(majorVersion,minorVersion); return;  end
  if isempty(style),  return;  end
  if all(ishandle(style)) && length(style)~=3
      dumpElement(style);
      return;
  end

  % Process the text string
  if nargin<2, format = style; style='text';  end
  %error(nargchk(2, inf, nargin, 'struct'));
  %str = sprintf(format,varargin{:});

  % In compiled mode
  try useDesktop = usejava('desktop'); catch, useDesktop = false; end
  if isdeployed | ~useDesktop %#ok<OR2> - for Matlab 6 compatibility
      % do not display any formatting - use simple fprintf()
      % See: http://undocumentedmatlab.com/blog/bold-color-text-in-the-command-window/#comment-103035
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/1390a26e7ef4aa4d
      % Also see: https://mail.google.com/mail/u/0/?ui=2&shva=1#all/13a6ed3223333b21
      count1 = fprintf(format,varargin{:});
  else
      % Else (Matlab desktop mode)
      % Get the normalized style name and underlining flag
      [underlineFlag, boldFlag, style] = processStyleInfo(style);

      % Set hyperlinking, if so requested
      if underlineFlag
          format = ['<a href="">' format '</a>'];

          % Matlab 7.1 R14 (possibly a few newer versions as well?)
          % have a bug in rendering consecutive hyperlinks
          % This is fixed by appending a single non-linked space
          if majorVersion < 7 || (majorVersion==7 && minorVersion <= 1)
              format(end+1) = ' ';
          end
      end

      % Set bold, if requested and supported (R2011b+)
      if boldFlag
          if (majorVersion > 7 || minorVersion >= 13)
              format = ['<strong>' format '</strong>'];
          else
              boldFlag = 0;
          end
      end

      % Get the current CW position
      cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
      lastPos = cmdWinDoc.getLength;

      % If not beginning of line
      bolFlag = 0;  %#ok
      %if docElement.getEndOffset - docElement.getStartOffset > 1
          % Display a hyperlink element in order to force element separation
          % (otherwise adjacent elements on the same line will be merged)
          if majorVersion<7 || (majorVersion==7 && minorVersion<13)
              if ~underlineFlag
                  fprintf('<a href=""> </a>');  %fprintf('<a href=""> </a>\b');
              elseif format(end)~=10  % if no newline at end
                  fprintf(' ');  %fprintf(' \b');
              end
          end
          %drawnow;
          bolFlag = 1;
      %end

      % Get a handle to the Command Window component
      mde = com.mathworks.mde.desk.MLDesktop.getInstance;
      cw = mde.getClient('Command Window');
      xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);

      % Store the CW background color as a special color pref
      % This way, if the CW bg color changes (via File/Preferences), 
      % it will also affect existing rendered strs
      com.mathworks.services.Prefs.setColorPref('CW_BG_Color',xCmdWndView.getBackground);

      % Display the text in the Command Window
      count1 = fprintf(2,format,varargin{:});

      %awtinvoke(cmdWinDoc,'remove',lastPos,1);   % TODO: find out how to remove the extra '_'
      drawnow;  % this is necessary for the following to work properly (refer to Evgeny Pr in FEX comment 16/1/2011)
      docElement = cmdWinDoc.getParagraphElement(lastPos+1);
      if majorVersion<7 || (majorVersion==7 && minorVersion<13)
          if bolFlag && ~underlineFlag
              % Set the leading hyperlink space character ('_') to the bg color, effectively hiding it
              % Note: old Matlab versions have a bug in hyperlinks that need to be accounted for...
              %disp(' '); dumpElement(docElement)
              setElementStyle(docElement,'CW_BG_Color',1+underlineFlag,majorVersion,minorVersion); %+getUrlsFix(docElement));
              %disp(' '); dumpElement(docElement)
              el(end+1) = handle(docElement);  %#ok used in debug only
          end

          % Fix a problem with some hidden hyperlinks becoming unhidden...
          fixHyperlink(docElement);
          %dumpElement(docElement);
      end

      % Get the Document Element(s) corresponding to the latest fprintf operation
      while docElement.getStartOffset < cmdWinDoc.getLength
          % Set the element style according to the current style
          %disp(' '); dumpElement(docElement)
          specialFlag = underlineFlag | boldFlag;
          setElementStyle(docElement,style,specialFlag,majorVersion,minorVersion);
          %disp(' '); dumpElement(docElement)
          docElement2 = cmdWinDoc.getParagraphElement(docElement.getEndOffset+1);
          if isequal(docElement,docElement2),  break;  end
          docElement = docElement2;
          %disp(' '); dumpElement(docElement)
      end

      % Force a Command-Window repaint
      % Note: this is important in case the rendered str was not '\n'-terminated
      xCmdWndView.repaint;

      % The following is for debug use only:
      el(end+1) = handle(docElement);  %#ok used in debug only
      %elementStart  = docElement.getStartOffset;
      %elementLength = docElement.getEndOffset - elementStart;
      %txt = cmdWinDoc.getText(elementStart,elementLength);
  end

  if nargout
      count = count1;
  end
  return;  % debug breakpoint
end
% Process the requested style information
function [underlineFlag,boldFlag,style] = processStyleInfo(style)
  underlineFlag = 0;
  boldFlag = 0;

  % First, strip out the underline/bold markers
  if ischar(style)
      % Styles containing '-' or '_' should be underlined (using a no-target hyperlink hack)
      %if style(1)=='-'
      underlineIdx = (style=='-') | (style=='_');
      if any(underlineIdx)
          underlineFlag = 1;
          %style = style(2:end);
          style = style(~underlineIdx);
      end

      % Check for bold style (only if not underlined)
      boldIdx = (style=='*');
      if any(boldIdx)
          boldFlag = 1;
          style = style(~boldIdx);
      end
      if underlineFlag && boldFlag
          warning('YMA:cprintf:BoldUnderline','Matlab does not support both bold & underline')
      end

      % Check if the remaining style sting is a numeric vector
      %styleNum = str2num(style); %#ok<ST2NM>  % not good because style='text' is evaled!
      %if ~isempty(styleNum)
      if any(style==' ' | style==',' | style==';')
          style = str2num(style); %#ok<ST2NM>
      end
  end

  % Style = valid matlab RGB vector
  if isnumeric(style) && length(style)==3 && all(style<=1) && all(abs(style)>=0)
      if any(style<0)
          underlineFlag = 1;
          style = abs(style);
      end
      style = getColorStyle(style);

  elseif ~ischar(style)
      error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

  % Style name
  else
      % Try case-insensitive partial/full match with the accepted style names
      validStyles = {'Text','Keywords','Comments','Strings','UnterminatedStrings','SystemCommands','Errors', ...
                     'Black','Cyan','Magenta','Blue','Green','Red','Yellow','White', ...
                     'Hyperlinks'};
      matches = find(strncmpi(style,validStyles,length(style)));

      % No match - error
      if isempty(matches)
          error('YMA:cprintf:InvalidStyle','Invalid style - see help section for a list of valid style values')

      % Too many matches (ambiguous) - error
      elseif length(matches) > 1
          error('YMA:cprintf:AmbigStyle','Ambiguous style name - supply extra characters for uniqueness')

      % Regular text
      elseif matches == 1
          style = 'ColorsText';  % fixed by Danilo, 29/8/2011

      % Highlight preference style name
      elseif matches < 8
          style = ['Colors_M_' validStyles{matches}];

      % Color name
      elseif matches < length(validStyles)
          colors = [0,0,0; 0,1,1; 1,0,1; 0,0,1; 0,1,0; 1,0,0; 1,1,0; 1,1,1];
          requestedColor = colors(matches-7,:);
          style = getColorStyle(requestedColor);

      % Hyperlink
      else
          style = 'Colors_HTML_HTMLLinks';  % CWLink
          underlineFlag = 1;
      end
  end
end
% Convert a Matlab RGB vector into a known style name (e.g., '[255,37,0]')
function styleName = getColorStyle(rgb)
  intColor = int32(rgb*255);
  javaColor = java.awt.Color(intColor(1), intColor(2), intColor(3));
  styleName = sprintf('[%d,%d,%d]',intColor);
  com.mathworks.services.Prefs.setColorPref(styleName,javaColor);
end
% Fix a bug in some Matlab versions, where the number of URL segments
% is larger than the number of style segments in a doc element
function delta = getUrlsFix(docElement)  %#ok currently unused
  tokens = docElement.getAttribute('SyntaxTokens');
  links  = docElement.getAttribute('LinkStartTokens');
  if length(links) > length(tokens(1))
      delta = length(links) > length(tokens(1));
  else
      delta = 0;
  end
end
% fprintf(2,str) causes all previous '_'s in the line to become red - fix this
function fixHyperlink(docElement)
  try
      tokens = docElement.getAttribute('SyntaxTokens');
      urls   = docElement.getAttribute('HtmlLink');
      urls   = urls(2);
      links  = docElement.getAttribute('LinkStartTokens');
      offsets = tokens(1);
      styles  = tokens(2);
      doc = docElement.getDocument;

      % Loop over all segments in this docElement
      for idx = 1 : length(offsets)-1
          % If this is a hyperlink with no URL target and starts with ' ' and is collored as an error (red)...
          if strcmp(styles(idx).char,'Colors_M_Errors')
              character = char(doc.getText(offsets(idx)+docElement.getStartOffset,1));
              if strcmp(character,' ')
                  if isempty(urls(idx)) && links(idx)==0
                      % Revert the style color to the CW background color (i.e., hide it!)
                      styles(idx) = java.lang.String('CW_BG_Color');
                  end
              end
          end
      end
  catch
      % never mind...
  end
end
% Set an element to a particular style (color)
function setElementStyle(docElement,style,specialFlag, majorVersion,minorVersion)
  %global tokens links urls urlTargets  % for debug only
  global oldStyles
  if nargin<3,  specialFlag=0;  end
  % Set the last Element token to the requested style:
  % Colors:
  tokens = docElement.getAttribute('SyntaxTokens');
  try
      styles = tokens(2);
      oldStyles{end+1} = styles.cell;

      % Correct edge case problem
      extraInd = double(majorVersion>7 || (majorVersion==7 && minorVersion>=13));  % =0 for R2011a-, =1 for R2011b+
      %{
      if ~strcmp('CWLink',char(styles(end-hyperlinkFlag))) && ...
          strcmp('CWLink',char(styles(end-hyperlinkFlag-1)))
         extraInd = 0;%1;
      end
      hyperlinkFlag = ~isempty(strmatch('CWLink',tokens(2)));
      hyperlinkFlag = 0 + any(cellfun(@(c)(~isempty(c)&&strcmp(c,'CWLink')),tokens(2).cell));
      %}

      styles(end-extraInd) = java.lang.String('');
      styles(end-extraInd-specialFlag) = java.lang.String(style);  %#ok apparently unused but in reality used by Java
      if extraInd
          styles(end-specialFlag) = java.lang.String(style);
      end

      oldStyles{end} = [oldStyles{end} styles.cell];
  catch
      % never mind for now
  end
  
  % Underlines (hyperlinks):
  %{
  links = docElement.getAttribute('LinkStartTokens');
  if isempty(links)
      %docElement.addAttribute('LinkStartTokens',repmat(int32(-1),length(tokens(2)),1));
  else
      %TODO: remove hyperlink by setting the value to -1
  end
  %}

  % Correct empty URLs to be un-hyperlinkable (only underlined)
  urls = docElement.getAttribute('HtmlLink');
  if ~isempty(urls)
      urlTargets = urls(2);
      for urlIdx = 1 : length(urlTargets)
          try
              if urlTargets(urlIdx).length < 1
                  urlTargets(urlIdx) = [];  % '' => []
              end
          catch
              % never mind...
              a=1;  %#ok used for debug breakpoint...
          end
      end
  end
  
  % Bold: (currently unused because we cannot modify this immutable int32 numeric array)
  %{
  try
      %hasBold = docElement.isDefined('BoldStartTokens');
      bolds = docElement.getAttribute('BoldStartTokens');
      if ~isempty(bolds)
          %docElement.addAttribute('BoldStartTokens',repmat(int32(1),length(bolds),1));
      end
  catch
      % never mind - ignore...
      a=1;  %#ok used for debug breakpoint...
  end
  %}
  
  return;  % debug breakpoint
end
% Display information about element(s)
function dumpElement(docElements)
  %return;
  numElements = length(docElements);
  cmdWinDoc = docElements(1).getDocument;
  for elementIdx = 1 : numElements
      if numElements > 1,  fprintf('Element #%d:\n',elementIdx);  end
      docElement = docElements(elementIdx);
      if ~isjava(docElement),  docElement = docElement.java;  end
      %docElement.dump(java.lang.System.out,1)
      disp(' ');
      disp(docElement)
      tokens = docElement.getAttribute('SyntaxTokens');
      if isempty(tokens),  continue;  end
      links = docElement.getAttribute('LinkStartTokens');
      urls  = docElement.getAttribute('HtmlLink');
      try bolds = docElement.getAttribute('BoldStartTokens'); catch, bolds = []; end
      txt = {};
      tokenLengths = tokens(1);
      for tokenIdx = 1 : length(tokenLengths)-1
          tokenLength = diff(tokenLengths(tokenIdx+[0,1]));
          if (tokenLength < 0)
              tokenLength = docElement.getEndOffset - docElement.getStartOffset - tokenLengths(tokenIdx);
          end
          txt{tokenIdx} = cmdWinDoc.getText(docElement.getStartOffset+tokenLengths(tokenIdx),tokenLength).char;  %#ok
      end
      lastTokenStartOffset = docElement.getStartOffset + tokenLengths(end);
      txt{end+1} = cmdWinDoc.getText(lastTokenStartOffset, docElement.getEndOffset-lastTokenStartOffset).char;  %#ok
      %cmdWinDoc.uiinspect
      %docElement.uiinspect
      txt = strrep(txt',sprintf('\n'),'\n');
      try
          data = [tokens(2).cell m2c(tokens(1)) m2c(links) m2c(urls(1)) cell(urls(2)) m2c(bolds) txt];
          if elementIdx==1
              disp('    SyntaxTokens(2,1) - LinkStartTokens - HtmlLink(1,2) - BoldStartTokens - txt');
              disp('    ==============================================================================');
          end
      catch
          try
              data = [tokens(2).cell m2c(tokens(1)) m2c(links) txt];
          catch
              disp([tokens(2).cell m2c(tokens(1)) txt]);
              try
                  data = [m2c(links) m2c(urls(1)) cell(urls(2))];
              catch
                  % Mtlab 7.1 only has urls(1)...
                  data = [m2c(links) urls.cell];
              end
          end
      end
      disp(data)
  end
end
% Utility function to convert matrix => cell
function cells = m2c(data)
  %datasize = size(data);  cells = mat2cell(data,ones(1,datasize(1)),ones(1,datasize(2)));
  cells = num2cell(data);
end
% Display the help and demo
function showDemo(majorVersion,minorVersion)
  fprintf('cprintf displays formatted text in the Command Window.\n\n');
  fprintf('Syntax: count = cprintf(style,format,...);  click <a href="matlab:help cprintf">here</a> for details.\n\n');
  url = 'http://UndocumentedMatlab.com/blog/cprintf/';
  fprintf(['Technical description: <a href="' url '">' url '</a>\n\n']);
  fprintf('Demo:\n\n');
  boldFlag = majorVersion>7 || (majorVersion==7 && minorVersion>=13);
  s = ['cprintf(''text'',    ''regular black text'');' 10 ...
       'cprintf(''hyper'',   ''followed %s'',''by'');' 10 ...
       'cprintf(''key'',     ''%d colored'',' num2str(4+boldFlag) ');' 10 ...
       'cprintf(''-comment'',''& underlined'');' 10 ...
       'cprintf(''err'',     ''elements:\n'');' 10 ...
       'cprintf(''cyan'',    ''cyan'');' 10 ...
       'cprintf(''_green'',  ''underlined green'');' 10 ...
       'cprintf(-[1,0,1],  ''underlined magenta'');' 10 ...
       'cprintf([1,0.5,0], ''and multi-\nline orange\n'');' 10];
   if boldFlag
       % In R2011b+ the internal bug that causes the need for an extra space
       % is apparently fixed, so we must insert the sparator spaces manually...
       % On the other hand, 2011b enables *bold* format
       s = [s 'cprintf(''*blue'',   ''and *bold* (R2011b+ only)\n'');' 10];
       s = strrep(s, ''')',' '')');
       s = strrep(s, ''',5)',' '',5)');
       s = strrep(s, '\n ','\n');
   end
   disp(s);
   eval(s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
% - Fix: Remove leading space char (hidden underline '_')
% - Fix: Find workaround for multi-line quirks/limitations
% - Fix: Non-\n-terminated segments are displayed as black
% - Fix: Check whether the hyperlink fix for 7.1 is also needed on 7.2 etc.
% - Enh: Add font support

function xcor=fXcorr(Ind1,Ind2,Lim)
% cross correlation function
xcor=zeros(2*Lim+1,1);
for k=-Lim:Lim
    xcor(k+Lim+1)=length(intersect(Ind1,Ind2+k));
end
end
    
function MUFiring_NEW = find_doubles(MUFiring,fsamp,TH)
        % sort the sources and take the ones with lowest CoV
LW1 = 100e-3;
LW2 = 0.5e-3;
p=1;
DOPPI=[];
B = 1;
BUONA = [];

for MU1 = 1:size(MUFiring,2)-1,
    COV = std(diff(MUFiring{MU1}))/mean(diff(MUFiring{MU1}))*100;
    good = [];
    good = [MU1,COV];
    if isempty(intersect(MU1,DOPPI)),
        for MU2 = (MU1 + 1):size(MUFiring,2),
            if isempty(intersect(MU2,DOPPI)),
                Firing1 = MUFiring{MU1};
                Firing2 = MUFiring{MU2};

                [minimo,indmin] = min([length(Firing1),length(Firing2)]);
                
                %size(Firing1)
                %size(Firing2)
                if indmin==1,
                    Firing2 = Firing2(Firing2>min(Firing1) & Firing2<max(Firing1));
                else
                    Firing1 = Firing1(Firing1>min(Firing2) & Firing1<max(Firing2));
                end    
                %size(Firing1)
                %size(Firing2)
                
                CORR = fXcorr(Firing1,Firing2,round(LW1*fsamp));
                [massimo,indmax] = max(CORR);
                if indmax > round(LW2*fsamp)+1 && indmax<length(CORR)-round(LW2*fsamp),
                    SENS = sum(CORR(indmax-round(LW2*fsamp):indmax+round(LW2*fsamp)))/max([length(Firing1),length(Firing2)])*100;
                   % figure(100),plot(CORR);title([MU1,MU2,SENS]);pause(0.5);
                
                    if SENS>TH, % Soglia per riconoscere le sorgenti uguali, 30 % cross correlation
                        DOPPI(p) = MU2;  % salva doppioni
                        COV = std(diff(MUFiring{MU2}))/mean(diff(MUFiring{MU2}))*100;
                        good = [good;[MU2,COV]];
                        p=p+1;
                    end
                end
            end
        end
    [minimo,index] = min(good(:,2));
    %good(index,:)
    BUONA{B} = MUFiring{good(index,1)};    

    B = B+1;
    end
end   
MUFiring_NEW = BUONA;
end

