function [TQwin,AnalogWinTime, cutTorque, cutTorqueFeedback,varargout] = RampWindow_Analog(Window, Torque, TorqueFeedback, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fsamp = 2048;
TQwin = [(Window(1)*fsamp)+1 (Window(2)*fsamp)];

try
    cutTorque = Torque(TQwin(1):TQwin(2));
    cutTorqueFeedback = TorqueFeedback(TQwin(1):TQwin(2));
    AnalogWinTime = TQwin(1)/fsamp:1/fsamp:Window(2);
    
    if size(varargin) > 0 
        try
            dfTorque = (cutTorqueFeedback >= 0);
            
            pfTQfeedback = cutTorqueFeedback;
            dfTQfeedback = cutTorqueFeedback;
            
            dfTQfeedback(~dfTorque) = nan;
            pfTQfeedback(dfTorque) = nan;
            varargout = [pfTQfeedback,dfTQfeedback];
        catch
        end
    end
    
catch
    disp('Error with function RampWindow_Analog Variables');
    
end

end

