% nanoRFE_vna6GHz_plot.m
%{
Description: Import *.s1p and *.s2p files exported from NanoRFE VNA6000

The format for .s2p files are:
*.s2p Files Each record contains 1 stimulus value and 4 S-parameters 
(total of 9 values):
[1=Stim,2=Re(S11),3=Im(S11),4=Re(S21),5=Im(S21),...
 6=Re(S12),7=Im(S12),8=Re(S22),9=Im(S22)]

%}
% Author: Walter V. Unglaub
% Date Last Modified: 12/16/2024
%==========================================================================
clear all;
clc;
clf; 
% close all;
% warning off;
tic;
%==========================================================================
LW=1; FS=14; MS=10; AFS=FS; TFS=FS; azim=45; pol=45; fyr=0.05; publish=1;
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',AFS);
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultTextFontSize',TFS);
kilo=1e-3; Mega=1e-6; Giga=1e-9; Tera=1e-12;
milli=1e3; micro=1e6; nano=1e9; pico=1e12; femto=1e15; % unit prefixes
%==========================================================================
convert_read = 1; % '1'=convert, then read,'2'=read converted file

if convert_read==1
% Copy the files over with a new name.
fileName = 'port2.s2p';

% % First create the folder B, if necessary.
% outputFolder = fullfile(pwd, 'data')
% if ~exist(outputFolder, 'dir')
%   mkdir(outputFolder);
% end
% Prepare the input filename.
inputFullFileName = fullfile(pwd, fileName);
outputBaseFileName = sprintf('%s.txt', fileName(1:end-4)); % prepare output filename
outputFullFileName = fullfile(outputBaseFileName);
copyfile(inputFullFileName, outputFullFileName); % copy/rename s2p file
end

% Read in data
d = table2array(readtable('port2.txt'));

f   = d(:,1);            % [Hz] frequency
s11 = d(:,2)+1i*d(:,3);  % S11 complex signal
s21 = d(:,4)+1i*d(:,5);  % S21 complex signal
s12 = d(:,6)+1i*d(:,7);  % S12 complex signal
s22 = d(:,8)+1i*d(:,9);  % S22 complex signal

gain  = 20*log10(abs(s21));
phase = phase(s21);


ttl=['NanoRFE VNA6000 (50kHz-6GHz)'];

%{
figure(1); clf;
loglog(f,abs(s11),'-k','linewidth',LW); hold on; grid on;
xlabel(['Frequency, {\itf} (Hz)'],'Fontsize',FS);
ylabel(['|S11|'],'Fontsize',FS);
title(ttl,'Fontweight','normal','Fontsize',FS);
%}

% LOGMAG
figure(2); clf;
semilogx(f*Mega,gain,'-k','linewidth',LW); hold on; grid on;
xlim([min(f),max(f)]*Mega); ylim([-80,5]);
xlabel(['Frequency, {\itf} (MHz)'],'Fontsize',FS);
ylabel(['|S21| (dB)'],'Fontsize',FS);
title(ttl,'Fontweight','normal','Fontsize',FS);

%{
% PHASE
figure(3); clf;
semilogx(f,phase,'-k','linewidth',LW); hold on; grid on;
xlabel(['Frequency, {\itf} (Hz)'],'Fontsize',FS);
ylabel(['\phi(S21) (rads)'],'Fontsize',FS);
title(ttl,'Fontweight','normal','Fontsize',FS);
%}

%{
figure(3); 
loglog(f,abs(s12),'-k','linewidth',LW); hold on; grid on;
xlabel(['Frequency, {\itf} (Hz)'],'Fontsize',FS);
ylabel(['|S12|'],'Fontsize',FS);
title(ttl,'Fontweight','normal','Fontsize',FS);
figure(4); 
loglog(f,abs(s22),'-k','linewidth',LW); hold on; grid on;
xlabel(['Frequency, {\itf} (Hz)'],'Fontsize',FS);
ylabel(['|S22|'],'Fontsize',FS);
title(ttl,'Fontweight','normal','Fontsize',FS);
%}



























