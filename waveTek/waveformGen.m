clear all; clc;

resourceList = visadevlist
% Parameters
amp = 0.1;                 % Vpp
offset = 0;                 % DC offset
sRate = 1e6;                % 1 MSa/s
name = 'Wave';              % Arbitrary waveform name

% --- Generate waveform
arb = sin(2*pi*5e3*(0:1/sRate:1e-3));  % 5 kHz sine, 1 ms
IPAddress='169.254.5.22';

arbTo33500(arb,IPAddress,amp,sRate,name)
