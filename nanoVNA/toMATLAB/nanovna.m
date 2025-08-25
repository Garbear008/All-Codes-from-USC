%    MATLAB scripts to read data from NanoVNA vector network analyzer.
%    It allows to establish connection with NanoVNA, send commands,
%    obtain S11 and S21 parameters, save it to S2P file and show
%    Logmag and Smith chart with no needs to use external applications.
%
%    https://github.com/qrp73/NanoVNA-MATLAB
%
%    Copyright (C) 2019 alex_m
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
% function nanovna()
    clear; close all; clc;

    model=1; %'1'=10kHz-1.5GHz,'2'=10kHz-3GHz,'3'=50kHz-6GHz

    if model==1
        hCom = nanoOpen('COM5'); % DeepElec NanoVNA-F (10kHz-1.5GHz)
    elseif model==2
        hCom = nanoOpen('COM6'); % Seesii NanoVNA V2 (10kHz-3GHz)
    elseif model==3
        hCom = nanoOpen('COM5'); % 6GHz NanoRFE VNA6000-A (50kHz-6GHz)
    end
    
    cleaner = onCleanup(@() nanoClose(hCom));

    fmin = 50; % [kHz] start frequency
    fmax = 300; % [MHz] stop frequency
    
    % nanoCommand(hCom, ['sweep start ',num2str(fmin),'000']);
    % nanoCommand(hCom, ['sweep stop ',num2str(fmax),'000000']);

    nanoCommand(hCom, 'sweep start 50000');
    nanoCommand(hCom, 'sweep stop 300000000');

    netwk = nanoGetData(hCom); % extract data
    
    % save network data
    %rfwrite(netwk.Parameters, netwk.Frequencies, 'data.s2p');
    %netwk = sparameters('data.s2p');

    % plot S11 LOGMAG
    fig = figure('Name','LOGMAG', 'NumberTitle','off');
    %set(fig,'Position',[10 10 320 240])
    rfplot(netwk,1,1, '-r')
    ylim([-90 10])
    hold on
    rfplot(netwk,2,1, '-b')
    hold off
    
    % plot S11 SMITH CHART
    fig = figure('Name','SMITH', 'NumberTitle','off');
    %set(fig,'Position',[10 10 320 240])
    smith(netwk,1,1);

    % plot S11 PHASE
    fig = figure('Name','PHASE', 'NumberTitle','off');
    %set(fig,'Position',[10 10 320 240])
    rfplot(netwk,1,1,'angle','-r')
    
    save netwk netwk;
    
    %===plot TDR step response===
    
    % fit to a rational function object
    tdrfit = rationalfit(netwk.Frequencies, rfparam(netwk, 1,1));
    
    % Parameters for a step signal
    Ts = 3e-12;
    N = 32768;
    Trise = 1e-10;
    
    % Calculate the step response for TDR and plot it
    [tdr,t1] = stepresp(tdrfit,Ts,N,Trise);
    tdrz = arrayfun(@(x) 50 * (1+x)/(1-x), tdr);
    
    fig = figure('Name','TDR Step Response', 'NumberTitle','off');
    plot(t1*1e9, tdrz, 'LineWidth',2);
    ylabel('TDR Step Response [Ohm]');
    xlabel('Delay [ns]');
    xlim([0 Ts*N*1e9])
    ylim([10 60])
    grid on;
    grid minor;
    
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







