function [DC, AC, Voltage] = Kardio_getdata(varargin)

%[DC, AC, Voltage] = Kardio_getdata(varargin)
%
%This function will plot data acquired from patch-clamp recordings in
%primary cardiomyocytes.
%
%---Usage---
%
%The function can be called with either 1 or 0 input arguments. Calling it
%with no inputs opens a prompt to select the file which should be read
%out. When calling the function with one input argument, the path and the
%filename should be parsed.
%The function will provide three outputs: The unfiltered raw trace (DC), 
%the AC filtered trace (AC) and the voltage trace (Voltage).
%The function only accepts .mat files.
%
%(c) Pascal Geschwill 13.11.2015

%% Input Handling
switch nargin
    case 0
        [FileName,PathName] = uigetfile;  %Choose the file
        Data = load([PathName,FileName]);
    case 1
        if ischar(varargin{1}) && strcmp(varargin{1}(end-3:end),'.mat')           
            FileName = varargin{1}; %Filename is passed as argument  
            Data = load(FileName);
        else
            error('Kardio_getdata:NoStringInput','String input expected or wrong filetype.')
        end
    otherwise
        error('Kardio_getdata:TooManyInputs','Too many imput arguments')
end

%% Readout
DataFieldname = fieldnames(Data); %Data is imported as an overly convoluted struct

%Read data
DC = Data.(DataFieldname{1}).values(:,1,:);
AC = Data.(DataFieldname{1}).values(:,2,:);
Voltage = Data.(DataFieldname{1}).values(:,3,:);

%Reformat data
DC = reshape(DC, length(DC),size(DC,3));
AC = reshape(AC, length(AC),size(AC,3));
Voltage = reshape(Voltage, length(Voltage),size(Voltage,3));

end