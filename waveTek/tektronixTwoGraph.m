clear all; clc;
set(0,'defaultAxesFontSize', 16);

resourceList = visadevlist
%get oscilloscope VISA object
osc = visadev("USB0::0x0699::0x0367::C103716::0::INSTR");  

%read CH1
[t1, v1] = readChannel(osc, 'CH1');

%read CH2
[t2, v2] = readChannel(osc, 'CH2');

%plot
figure;
plot(t1, v1, 'b', t2, v2, 'r');
legend('CH1', 'CH2');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Waveforms from CH1 and CH2');
grid on;

%plot
figure;
plot(tm, v1m, 'b', tm, v2m, 'r');
legend('CH1', 'CH2');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Multiple Waveforms from CH1 and CH2');
grid on;

%function to read channels
function [t, v] = readChannel(osc, chan)
    writeline(osc, ['DATA:SOURCE ' chan]);
    writeline(osc, 'DATA:ENC RPB');      % Fast binary
    writeline(osc, 'DATA:WIDTH 1');      % 1 byte per point

    writeline(osc, 'WFMPRE:XINCR?'); xIncr = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YOFF?');  yOff  = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YMULT?'); yMult = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YZERO?'); yZero = str2double(readline(osc));
    writeline(osc, 'HOR:RECORDLENGTH?');
    maxPoints = str2double(readline(osc));
    writeline(osc, sprintf('DATA:STOP %d', maxPoints));

    writeline(osc, 'CURVE?');
    raw = readbinblock(osc, 'uint8');
    readline(osc); % Clear LF

    v = (double(raw) - yOff) * yMult + yZero;
    t = (0:length(v)-1) * xIncr;
end
