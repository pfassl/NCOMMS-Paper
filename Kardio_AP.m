function [RMP, AP_Amp, APD20, APD50, APD90] = Kardio_AP(varargin)

%[RMP, AP_Amp, APD20, APD50, APD90] = Kardio_AP(varargin)
%
%This function analyzes action potentials elicited in primary
%cardiomyocytes. Specifically, it returns the resting membrane potential
%(RMP), action potential amplitude (AP_Amp, [mV]) and the AP durations 
%(widths, [ms]) at 20%, 50% and 90% repolarzation (APD20, APD50, APD90).
%
%---Usage---
%
%The function can be called with either one or two input arguments. If only
%one input argument is passed, the input has to be the path of the source
%file as a string. E.g. 
%
%[RMP, AP_Amp, APD20, APD50, APD90] = 
%Kardio_AP('/mnt/raw/pascal/data/Kardio_Patch/151210/151210_001_AP_05Hz.mat')
%
%It is also possible the add 'plot' as the second input argument to the
%function which will cause the script to plot the analysis results.
%
%The output values are stored as a vector where each value corresponds to
%the calculated RMP, AP_Amp and APDs obtained from each sweep. RMP and
%AP_Amp are given in mV. APDs are given in ms.
%
%The function only accepts .mat files.
%
%(c) Pascal Geschwill 14.12.2015

%% Input Handling
switch nargin
    case 0
        error('Kardio_AP:NoInput','Input expected.')
    case 1
        if ischar(varargin{1}) && isempty(strfind(varargin{1},'2hz'))
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = false;
            TwoHzFlag = false;
        elseif ischar(varargin{1}) && ~isempty(strfind(varargin{1},'2hz'))
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = false;
            TwoHzFlag = true;
        else
            error('Kardio_AP:NoStringInput','String input expected.')
        end
    case 2
        if ischar(varargin{1}) && strcmp(varargin{2},'plot') && isempty(strfind(varargin{1},'2hz'))
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = true;
            TwoHzFlag = false;
        elseif ischar(varargin{1}) && strcmp(varargin{2},'plot') && ~isempty(strfind(varargin{1},'2hz'))
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = true;
            TwoHzFlag = true;
        else
            error('Kardio_AP:NoStringInput','Unexpected Input.')
        end
    otherwise
        error('Kardio_AP:TooManyInputs','Too many imput arguments.')
end

%% Analysis of AP Widths

dt = 1/20000; %Sampling frequency; 1 datapoint corresponds to 0.05ms

NoOfSweeps = size(Data,2);

if TwoHzFlag == false
    
    %Initialize Matrices, AP Duration at 20%, 50% and 90% Repolarization
    APD20_mat = zeros(NoOfSweeps,1);
    APD50_mat = zeros(NoOfSweeps,1);
    APD90_mat = zeros(NoOfSweeps,1);
    RMP_mat = zeros(NoOfSweeps,1); %Resting Membrane Potential
    AP_Amp_mat = zeros(NoOfSweeps,1); %Action Potential Amplitudes
    
else
    
    %Initialize Matrices, AP Duration at 20%, 50% and 90% Repolarization
    APD20_mat = zeros(2*NoOfSweeps,1);
    APD50_mat = zeros(2*NoOfSweeps,1);
    APD90_mat = zeros(2*NoOfSweeps,1);
    RMP_mat = zeros(2*NoOfSweeps,1); %Resting Membrane Potential
    AP_Amp_mat = zeros(2*NoOfSweeps,1); %Action Potential Amplitudes
    
end

for cS = 1:NoOfSweeps
    
    RMP_mat(cS) = mean(Data(1500:2000,cS)); %Calculate resting potential
    
    AP_Amp_mat(cS) = max(Data(1500:4000,cS)) - RMP_mat(cS); %Calculate AP Amplitude
    Repol90 = AP_Amp_mat(cS) - AP_Amp_mat(cS)*0.9 + RMP_mat(cS);  %Use the Amplitude of the AP to calculate the percent repolarization
    Repol50 = AP_Amp_mat(cS) - AP_Amp_mat(cS)*0.5 + RMP_mat(cS);
    Repol20 = AP_Amp_mat(cS) - AP_Amp_mat(cS)*0.2 + RMP_mat(cS);
    
    %Use smoothed data for estimating the width of the AP because of
    %oversampling
    SmoothData = smooth(Data(:,cS),10);
       
    %Find the intersections of the Repolarization Value with the trace to
    %get the width of the AP
    [x90, ~] = intersections([0 4000], [Repol90 Repol90], 1:4000, SmoothData(1:4000));
    [x50, ~] = intersections([0 4000], [Repol50 Repol50], 1:4000, SmoothData(1:4000));
    [x20, ~] = intersections([0 4000], [Repol20 Repol20], 1:4000, SmoothData(1:4000));

%     [x90, ~] = intersections([0 4000], [Repol90 Repol90], 1:4000, Data(1:4000,cS));
%     [x50, ~] = intersections([0 4000], [Repol50 Repol50], 1:4000, Data(1:4000,cS));
%     [x20, ~] = intersections([0 4000], [Repol20 Repol20], 1:4000, Data(1:4000,cS));
    
    APD90_mat(cS) = (x90(end) - x90(1))*dt*1000;    %Calculate width of AP by taking the difference of the first and the last intersection
    APD50_mat(cS) = (x50(end) - x50(1))*dt*1000;
    APD20_mat(cS) = (x20(end) - x20(1))*dt*1000;
    
%     %Use max of dv/dt as starting point for APD measurement
%     Data_derived = diff(SmoothData(1000:4000));
%     [~, DiffMaxIdx1] = max(Data_derived);
%     
%     APD90_mat(cS) = (x90(end) - (DiffMaxIdx1+1000))*dt*1000;    %Calculate width of AP by taking the difference of the last intersection and the max dv/dt
%     APD50_mat(cS) = (x50(end) - DiffMaxIdx)*dt*1000;    %Possibly not as accurate as the differente of the first and last intersection
%     APD20_mat(cS) = (x20(end) - DiffMaxIdx)*dt*1000;
    
    if cS >= 15 && TwoHzFlag == true
        
        for cS2 = 16:2*NoOfSweeps
            
            RMP_mat(cS2) = mean(Data(11500:12000,cS2-15)); %Calculate resting potential
            
            AP_Amp_mat(cS2) = max(Data(11500:14000,cS2-15)) - RMP_mat(cS2); %Calculate AP Amplitude
            Repol90 = AP_Amp_mat(cS2) - AP_Amp_mat(cS2)*0.9 + RMP_mat(cS2);  %Use the Amplitude of the AP to calculate the precent repolarization
            Repol50 = AP_Amp_mat(cS2) - AP_Amp_mat(cS2)*0.5 + RMP_mat(cS2);
            Repol20 = AP_Amp_mat(cS2) - AP_Amp_mat(cS2)*0.2 + RMP_mat(cS2);
            
            %Use smoothed data for estimating the width of the AP because of
            %oversampling
            SmoothData = smooth(Data(:,cS2-15),10);
            
            %Find the intersections of the Repolarization Value with the trace to
            %get the width of the AP
            [x90_2, ~] = intersections([10000 14000], [Repol90 Repol90], 10001:14000, SmoothData(10001:14000));
            [x50_2, ~] = intersections([10000 14000], [Repol50 Repol50], 10001:14000, SmoothData(10001:14000));
            [x20_2, ~] = intersections([10000 14000], [Repol20 Repol20], 10001:14000, SmoothData(10001:14000));
            
            %     [x90, ~] = intersections([0 4000], [Repol90 Repol90], 1:4000, Data(1:4000,cS));
            %     [x50, ~] = intersections([0 4000], [Repol50 Repol50], 1:4000, Data(1:4000,cS));
            %     [x20, ~] = intersections([0 4000], [Repol20 Repol20], 1:4000, Data(1:4000,cS));
            
            APD90_mat(cS2) = (x90_2(end) - x90_2(1))*dt*1000;    %Calculate width of AP by taking the difference of the first and the last intersection
            APD50_mat(cS2) = (x50_2(end) - x50_2(1))*dt*1000;
            APD20_mat(cS2) = (x20_2(end) - x20_2(1))*dt*1000;
            
%             %Use max of dv/dt as starting point for APD measurement
%             Data_derived = diff(SmoothData(11000:14000));
%             [~, DiffMaxIdx2] = max(Data_derived);
%             
%             APD90_mat(cS2) = (x90_2(end) - (DiffMaxIdx2+11000))*dt*1000;    %Calculate width of AP by taking the difference of the last intersection and the max dv/dt
%             APD50_mat(cS2) = (x50_2(end) - DiffMaxIdx)*dt*1000;    %Possibly not as accurate as the differente of the first and last intersection
%             APD20_mat(cS2) = (x20_2(end) - DiffMaxIdx)*dt*1000;
            
        end
        
    end
    
    if PlotFlag == 1
        
        hold on
        figure(1)
        
        subplot(2,3,[1,2,4,5])
        plot(1000*dt-cS*0.03:dt:4000*dt-cS*0.03, Data(1000:4000,cS) - (cS-1) * 30, 'color', [0.5 0.5 0.5])      
        hold on
        plot(1000*dt-cS*0.03:dt:4000*dt-cS*0.03, SmoothData(1000:4000) - (cS-1) * 30, 'g')
        plot([1500*dt-cS*0.03 2000*dt-cS*0.03], [RMP_mat(cS)- (cS-1) * 30 RMP_mat(cS)- (cS-1) * 30],'b', 'linewidth', 2)
        plot([x90(1)*dt-cS*0.03 x90(end)*dt-cS*0.03], [Repol90 - (cS-1) * 30 Repol90 - (cS-1) * 30],'r', 'linewidth', 2)
        plot([x50(1)*dt-cS*0.03 x50(end)*dt-cS*0.03], [Repol50 - (cS-1) * 30 Repol50 - (cS-1) * 30],'r', 'linewidth', 2)
        plot([x20(1)*dt-cS*0.03 x20(end)*dt-cS*0.03], [Repol20 - (cS-1) * 30 Repol20 - (cS-1) * 30],'r', 'linewidth', 2)

%         This is for dvdt APD measurement
%         plot([(DiffMaxIdx1+1000)*dt-cS*0.03 (DiffMaxIdx1+1000)*dt-cS*0.03], [-90 - (cS-1) * 30 50 - (cS-1) * 30],'k')
%         plot([(DiffMaxIdx1+1000)*dt-cS*0.03 x90(end)*dt-cS*0.03], [Repol90 - (cS-1) * 30 Repol90 - (cS-1) * 30],'r', 'linewidth', 2)

%         if cS >= 15 && TwoHzFlag == true
% %             
%             hold on
%             plot(5000*dt-(cS2-15)*0.03:dt:8000*dt-(cS2-15)*0.03, Data(11000:14000,cS2-15) - (cS2-15-1) * 30, 'color', [0.5 0.5 0.5])
%             hold on
%             plot(5000*dt-(cS2-15)*0.03:dt:8000*dt-(cS2-15)*0.03, SmoothData(11000:14000) - (cS2-15-1) * 30, 'g')
%             plot([5500*dt-(cS2-15)*0.03 6000*dt-(cS2-15)*0.03], [RMP_mat(cS2)- (cS2-15-1) * 30 RMP_mat(cS2)- (cS2-15-1) * 30],'b', 'linewidth', 2)
%             %         plot([x90(1)*dt-cS*0.03 x90(end)*dt-cS*0.03], [Repol90 - (cS-1) * 30 Repol90 - (cS-1) * 30],'r', 'linewidth', 2)
% %             plot([x50(1)*dt-(cS2-15)*0.03 x50(end)*dt-(cS2-15)*0.03], [Repol50 - (cS2-15-1) * 30 Repol50 - (cS2-15-1) * 30],'r', 'linewidth', 2)
% %             plot([x20(1)*dt-(cS2-15)*0.03 x20(end)*dt-(cS2-15)*0.03], [Repol20 - (cS2-15-1) * 30 Repol20 - (cS2-15-1) * 30],'r', 'linewidth', 2)
%             plot([(DiffMaxIdx2+5000)*dt-(cS2-15)*0.03 (DiffMaxIdx2+5000)*dt-(cS2-15)*0.03], [-90 - (cS2-15-1) * 30 50 - (cS2-15-1) * 30],'k')
%             plot([(DiffMaxIdx2+5000)*dt-(cS2-15)*0.03 x90(end)*dt-(cS2-15)*0.03], [Repol90 - (cS2-15-1) * 30 Repol90 - (cS2-15-1) * 30],'r', 'linewidth', 2)  
% %             
%         end
            
        xlabel('Time [s]')
        ylabel('Membrane Potential [mV]') 
        text(-0.3, 0, varargin{1})
        box off
        
        if cS == NoOfSweeps
            
        subplot(2,3,[3 6])
        bar([median(APD20_mat) median(APD50_mat) median(APD90_mat)])
        ylim([0 70])
        set(gca, 'XTick', 1:3, 'XTickLabel', {'APD20', 'APD50', 'APD90'});
        text(0.75, max(APD90_mat)-2, ['Median RMP = ', num2str(round(median(RMP_mat), 2)),'mV'])
        text(0.75, max(APD90_mat)-6, ['Median AP Amplitude = ', num2str(round(median(AP_Amp_mat), 2)),'mV'])
        ylabel('Median Action Potential Duration [ms]')
        box off
             
        end
        
    end
    
%     keyboard
%     clf
    
end

%% Output

RMP = RMP_mat;
AP_Amp = AP_Amp_mat;
APD20 = APD20_mat;
APD50 = APD50_mat;
APD90 = APD90_mat;

end