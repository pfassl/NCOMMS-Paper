function [Rs, Rm, Cap, Ih, I_HCN] = Kardio_ResCap(varargin)
    
%[Rs, Rm, Cap, Ih] = Kardio_ResCap(varargin)
%
%This function reads a recorded data file from a cardiomyocyte and
%calculates the access resistance (Rs), the membrane resistance (Rm), the
%capacitance (Cap), the holding current (Ih) as well as the HCN (I_HCN) 
%current of the recorded cell. 
%The algorithm uses the capacitive current which is elicited by the *first*
%+5mV pulse of the double rectangular test pulse that is applied at the
%beginning of each sweep.
%
%---Usage---
%
%The function can be called with either one or two input arguments. If only
%one input argument is passed, the input has to be the path of the source
%file as a string. E.g. 
%
%[Rs, Rm, Cap, Ih. I_HCN] =
%Kardio_ResCap('/mnt/raw/pascal/data/Kardio_Patch/151105/151105_001.mat').
%
%It is also possible the add 'plot' as the second input argument to the
%function which will cause the script to plot the analysis results.
%
%The output values are stored as a vector where each value corresponds to
%the calculated Rs, Rm, Cap and Ih obtained from each sweep. Rs and Rm are 
%given in MOhms, Cap is given in pF and Ih/I_HCN are given in pA.
%
%The function only accepts .mat files.
%
%(c) Pascal Geschwill 13.11.2015

%% Input Handling
switch nargin
    case 0
        error('Kardio_ResCap:NoInput','Input expected.')
    case 1
        if ischar(varargin{1})
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = false;
        else
            error('Kardio_ResCap:NoStringInput','String input expected.')
        end
    case 2
        if ischar(varargin{1}) && strcmp(varargin{2},'plot')
            Data = Kardio_getdata(varargin{1}); %Only reads DC trace
            PlotFlag = true;
        else
            error('Kardio_ResCap:NoStringInput','Unexpected Input.')
        end
    otherwise
        error('Kardio_ResCap:TooManyInputs','Too many imput arguments.')
end

%% Calculation of Rs, Rm, Cm, Ih and I_HCN

dt = 1/20000; %Sampling frequency; 1 datapoint corresponds to 0.05ms

%Definition of exponential fit function for calculating the membrane time
%constant from the first capacitive transient (+5mV step)
%Adapted from: http://www.mathworks.com/matlabcentral/fileexchange/41774-synapticconductance/content/SynapticConductance/iv_curve_ABF.m
fitwindow = 0:dt:dt*(400-1);
f2 = fittype('a*exp(-b*x) + c','dependent',{'y'},'independent',{'x'},'coefficients',{'a', 'b', 'c'});
fo = fitoptions('method','NonlinearLeastSquares','normalize', 'off','startpoint',[1 200 1]);

%Initialize output matrices
Rs_mat      = zeros(15,1);
Rm_mat      = zeros(15,1);
Cm_mat      = zeros(15,1);
tau_mat     = zeros(15,1);
Ih_mat      = zeros(15,1);
I_HCN_mat   = zeros(15,1);

for cS = 1:15 %Number of sweeps
    
    Ih_mat(cS) = mean(Data(1500:2000,cS)); %Calculate holding current
    I_HCN_mat(cS) = mean(Data(69000:70000,cS)) - Ih_mat(cS); %100ms at the end of the hyperpolarizing step, corrected for holding current
    Rs_mat(cS) = 1000 * (5 / (max(Data(1:4000,cS)) - Ih_mat(cS))); %5mV Step *1000 to get MOhm value
    Rm_mat(cS) = 1000 * (5 / (mean(Data(3500:4000,cS)) - Ih_mat(cS))) - Rs_mat(cS); %5mV Step *1000 to get MOhm value
    
    [~, i] = max(Data(1:4000,cS)); %Index of the max of the capacitive transient
    CapFit = fit(fitwindow', Data(i:i+400-1,cS), f2, fo);
    
    tau_mat(cS) = 1000/CapFit.b;
%     Cm_mat(cS) = tau_mat(cS)/Rm_mat(cS) * 1000; %Inverse of the CapFit.b value gives time constant in seconds, therefore multiply by 1000 to get ms.
    %Then again multiply by 1000 to get value in pF
 
    Cm_mat(cS) = 1000*tau_mat(cS) * (1/Rm_mat(cS) + 1/Rs_mat(cS)); %Implemented from Gentet et al. Biophys J. 2000 Jul; 79(1): 314?320. 
    
%     Cm_mat(cS) = 1000*tau_mat(cS) * (1/Rs_mat(cS)); %Implemented from Gentet et al. Biophys J. 2000 Jul; 79(1): 314?320. 

    if PlotFlag == true
        
        figure(1)
        
%         subplot(2,3,1)
%         hold on
%         plot(0-cS*0.01:dt:dt*(4000-1)-cS*0.01, Data(1:4000,cS))
%         plot([1500*dt 2000*dt],[mean(Data(1500:2000,cS)) mean(Data(1500:2000,cS))], 'r', 'linewidth', 3)
%         plot([3500*dt 4000*dt],[mean(Data(3500:4000,cS)) mean(Data(3500:4000,cS))], 'r','linewidth', 3)
%         xlabel('Time [s]')
%         zlabel('Membrane Current [pA]')
%         title('Test pulses')
        
        subplot(2,2,1)
        hold on
        plot(0-cS*0.01:dt:dt*(4000-1)-cS*0.01, Data(1:4000,cS) - (cS-1) * 40)
        plot([1500*dt-cS*0.01 2000*dt-cS*0.01],[mean(Data(1500:2000,cS))-(cS-1)*40 mean(Data(1500:2000,cS))-(cS-1)*40], 'r', 'linewidth', 3)
        plot([3500*dt-cS*0.01 4000*dt-cS*0.01],[mean(Data(3500:4000,cS))-(cS-1)*40 mean(Data(3500:4000,cS))-(cS-1)*40], 'r','linewidth', 3)        
        xlabel('Time [s]')
        ylabel('Membrane Current [pA]')
        title('Test pulses, sorted')
        
        subplot(2,2,2)
        hold on 
        plot(CapFit, fitwindow, Data(i:i+400-1,cS))
        xlabel('Time [s]')
        ylabel('Membrane Current [pA]')
        legend(gca, 'off')
        title('Exponential fits')
        
        subplot(2,2,3)
        hold on
        plot(0:dt:dt*(100000-1), Data(:,cS))
        plot([69000*dt 70000*dt], [mean(Data(69000:70000,cS)) mean(Data(69000:70000,cS))],'r', 'linewidth', 3)
        ylim([I_HCN_mat(1)+Ih_mat(1)-100 500])
        xlabel('Time [s]')
        ylabel('Membrane Current [pA]')
        title('HCN Currents')
        
        subplot(2,2,4)
        hold on
        plot(-140:10:0, I_HCN_mat./Cm_mat,'-o')
        xlabel('Membrane Potential [mV]')
        ylabel('Current Density [pA/pF]')
        title('Current Density')
        
    end
    
    
end

%% Output
Rs = Rs_mat;
Rm = Rm_mat;
Cap = Cm_mat;
Ih = Ih_mat;
I_HCN = I_HCN_mat;

end
    