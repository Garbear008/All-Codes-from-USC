clear; clc;
%check if the oscilloscope is embedded as a VISA object
resourceList = visadevlist
% use existing visadev object
osc = visadev("USB0::0x0699::0x0367::C103716::0::INSTR");  

% configure waveform source and format
writeline(osc, 'DATA:SOURCE CH1');      %use Channel 1
writeline(osc, 'DATA:ENC RPB');         % RPBinary = signed, fastest
writeline(osc, 'DATA:WIDTH 1');         %1 byte per point

%read waveform preamble
writeline(osc, 'WFMPRE:XINCR?');  
xIncr = str2double(readline(osc));
writeline(osc, 'WFMPRE:YOFF?');  
yOff  = str2double(readline(osc));
writeline(osc, 'WFMPRE:YMULT?');  
yMult = str2double(readline(osc));
writeline(osc, 'WFMPRE:YZERO?');  
yZero = str2double(readline(osc));

%get waveform data
writeline(osc, 'CURVE?');

%get binary waveform block from waveform data
rawData = readbinblock(osc, 'uint8');
readline(osc);  %read and discard the LF terminator

%convert binary data to voltage
voltage = (double(rawData) - yOff) * yMult + yZero;
time = (0:length(voltage)-1) * xIncr;

% plot waveform
figure;
plot(time, voltage);
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Tektronix TDS 2012B - Channel 1');
grid on;

