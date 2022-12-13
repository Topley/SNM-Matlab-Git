function [ForcePlates, COP_vars, COP_var_HP] = bertec_COP(raw_forces_moments, low_pass, Fs, varargin)

% For this function to work, the variable raw_forces_moments needs to be a
% 12 column matrix, with columns 1-12 corresponding to left
% forces (x,y,z) and moments (x,y,z) and right forces and moments (in that
% order). 

% Outputs are the corresponding forces and moments after low-pass filtering
% as well as the left, right, and net COP signals. 

% Optional argument to provide the same filtered COP variables after applying
% a high-pass filter. Specify 'high_pass' followed by the cutoff frequency.


% % % to test
% raw_forces_moments=vicon_data_upsampled;
% low_pass=5;
% Fs=2048;
% varargin='high';

high_pass=[];
if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'high_pass'
        high_pass = varargin{i+1};
      otherwise
        errordlg('unknown argument')
    end
  end
else
end


lFx_raw = raw_forces_moments.LeftFx; % left plate data
lFy_raw = raw_forces_moments.LeftFy;
lFz_raw = raw_forces_moments.LeftFz;
lMx_raw = raw_forces_moments.LeftMx;
lMy_raw = raw_forces_moments.LeftMy;
lMz_raw = raw_forces_moments.LeftMz;

rFx_raw=raw_forces_moments.RightFx; % right plate data
rFy_raw=raw_forces_moments.RightFy;
rFz_raw=raw_forces_moments.RightFz;
rMx_raw=raw_forces_moments.RightMx;
rMy_raw=raw_forces_moments.RightMy;
rMz_raw=raw_forces_moments.RightMz;

% filter raw force and moment data
cut=low_pass/(Fs*0.5);
[B,A]=butter(5,cut);

lFx=filtfilt(B,A,lFx_raw);
lFy=filtfilt(B,A,lFy_raw);
lFz=filtfilt(B,A,lFz_raw);
lMx=filtfilt(B,A,lMx_raw);
lMy=filtfilt(B,A,lMy_raw);
lMz=filtfilt(B,A,lMz_raw);

rFx=filtfilt(B,A,rFx_raw);
rFy=filtfilt(B,A,rFy_raw);
rFz=filtfilt(B,A,rFz_raw);
rMx=filtfilt(B,A,rMx_raw);
rMy=filtfilt(B,A,rMy_raw);
rMz=filtfilt(B,A,rMz_raw);

% calculate COP - LEFT plate
L_COPx_biased=(((lMy)./lFz)*100)-50; %ML
L_COPy_biased=(((lMx)./lFz)*100); %AP
% debias left COPx and COPy
L_COPx=L_COPx_biased-mean(L_COPx_biased);
L_COPy=L_COPy_biased-mean(L_COPy_biased);

% calculate COP - RIGHT plate
R_COPx_biased=(((rMy)./rFz)*100)+50; 
R_COPy_biased=(((rMx)./rFz)*100); 
% debias right COPx and COPy
R_COPx=R_COPx_biased-mean(R_COPx_biased);
R_COPy=R_COPy_biased-mean(R_COPy_biased); 

% Weighted COP (COP net)
COPx_weighted_biased=[];
COPy_weighted_biased=[];
for iii=1:length(R_COPx_biased)
    total_Fz=rFz(iii)+lFz(iii);
    R_COP_weight=rFz(iii)/total_Fz;
    L_COP_weight=lFz(iii)/total_Fz;

    COPx_weighted_biased(iii,1)=(L_COPx_biased(iii)*L_COP_weight)+(R_COPx_biased(iii)*R_COP_weight);
    COPy_weighted_biased(iii,1)=(L_COPy_biased(iii)*L_COP_weight)+(R_COPy_biased(iii)*R_COP_weight);
end

% debias COP net
COPx_weighted=COPx_weighted_biased-mean(COPx_weighted_biased(1:100)); %ML
COPy_weighted=COPy_weighted_biased-mean(COPy_weighted_biased(1:100)); %AP

filt_forces_moments=[lFx,lFy,lFz,lMx,lMy,lMz,rFx,rFy,rFz,rMx,rMy,rMz];

ForcePlates.LeftFx = lFx; 
ForcePlates.LeftFy = lFy; 
ForcePlates.LeftFz = lFz; 
ForcePlates.LeftMx = lMx; 
ForcePlates.LeftMy = lMy; 
ForcePlates.LeftMz = lMz; 
ForcePlates.RightFx = rFx; 
ForcePlates.RightFy = rFy; 
ForcePlates.RightFz = rFz; 
ForcePlates.RightMx = rMx; 
ForcePlates.RightMy = rMy;
ForcePlates.RightMz = rMz; 

COP_vars=[L_COPx_biased,L_COPy_biased,R_COPx_biased,R_COPy_biased,COPx_weighted_biased,COPy_weighted_biased];

if ~isempty(high_pass)
    highcut = high_pass/(Fs*0.5);
    [B2,A2] = butter(2,highcut,'high');
    
    R_COPx_HP = filtfilt(B2,A2,R_COPx);
    R_COPy_HP = filtfilt(B2,A2,R_COPy);
    
    L_COPx_HP = filtfilt(B2,A2,L_COPx);
    L_COPy_HP = filtfilt(B2,A2,L_COPy);
    
    COPx_weighted_HP = filtfilt(B2,A2,COPx_weighted);
    COPy_weighted_HP = filtfilt(B2,A2,COPy_weighted);
    
    COP_var_HP=[R_COPx_HP,R_COPy_HP,L_COPx_HP,L_COPy_HP,COPx_weighted_HP,COPy_weighted_HP];
end