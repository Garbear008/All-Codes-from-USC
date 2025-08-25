clear; clc; close all;
set(0,'defaultAxesFontSize', 16);

resourceList = visadevlist

%parameters
fStart = 1e3; %start of sweep frequency
fStop  = 1e6; %end of sweep frequency
nPoints = 10; %amount of data points from oscilloscope
amplitude = 0.5; %amplitude of wave
settleTime = 2; %amount of time to pause to display accurate figures on oscillscope (usually set at 1)
numOfPeriodsInDisplay = 5; %amount of periods in display for oscilloscope
sampleRate = 1e9; %sample rate of oscilloscope

%waveform generator and oscilloscope
wavegen = visadev("TCPIP0::169.254.5.22::inst0::INSTR:"); 
%for dell laptop wavegen use "USB0::0x0957::0x2C07::MY52800530::INSTR"
%otherwise use ethernet "TCPIP0::169.254.5.22::inst0::INSTR:"

scope   = visadev("USB0::0x0699::0x0367::C103716::0::INSTR");

%frequency vector
frequencies = logspace(log10(fStart), log10(fStop), nPoints);

%set a 1 by nPoints matrix for gain
gain_dB = zeros(1, nPoints);
phase=zeros(1,nPoints);

%configure waveform generator
writeline(wavegen, '*RST');
writeline(wavegen, 'FUNC SIN');
writeline(wavegen, sprintf('VOLT %f', amplitude));
writeline(wavegen, 'VOLT:OFFSET 0');
writeline(wavegen, 'OUTPut ON');

%turn channels on oscilloscope and allow probing
writeline(scope, '*RST');
writeline(scope, 'SELect:CH1 ON');
writeline(scope, 'SELect:CH2 ON');
writeline(scope, 'CH1:PROBE 1');
writeline(scope, 'CH2:PROBE 1');

%trigger setup
writeline(scope, 'TRIGger:MAIn:EDGE:SOURce CH1');
writeline(scope, 'TRIGger:MAIn:EDGE:SLOPe RISing');
writeline(scope, 'TRIGger:MAIn:LEVel 0');

%set up measurements
writeline(scope, 'MEASUrement:MEAS1:TYPe PK2pk');
writeline(scope, 'MEASUrement:MEAS1:SOURCE CH1');

writeline(scope, 'MEASUrement:MEAS2:TYPe PK2pk');
writeline(scope, 'MEASUrement:MEAS2:SOURCE CH2');
pause(1);

for k = 1:nPoints
    f = frequencies(k);
    fprintf("Sweeping: %.2f Hz (%.1f%% done)\n", f, 100*k/nPoints);

    % Update freq and scope timescale
    writeline(wavegen, sprintf('FREQ %e', f));
    % Estimate period
    period = 1 / f;
    horScale = (numOfPeriodsInDisplay * period) / 10;  % 10 divisions on screen
    writeline(scope, sprintf('HORizontal:MAIn:SCAle %e', horScale));
    
    % 2. Vertical scale
    v_div = amplitude / 3;  
    writeline(scope, sprintf('CH1:SCAle %e', v_div));
    writeline(scope, sprintf('CH2:SCAle %e', v_div));
    
    % 3. Center vertically
    writeline(scope, 'CH1:POSition 0');
    writeline(scope, 'CH2:POSition 0');
    
    if (k==1)
        pause(1);
    end
    pause(settleTime);

    % Read input
    writeline(scope, 'MEASUrement:MEAS1:VALue?'); v_in = str2double(readline(scope));
    writeline(scope, 'MEASUrement:MEAS2:VALue?'); v_out = str2double(readline(scope));

    writeline(scope, sprintf('CH2:SCAle %e', v_out/4));
    pause(1);
    writeline(scope, 'MEASUrement:MEAS2:VALue?'); v_out = str2double(readline(scope));
    gain_dB(k) = 20 * log10(v_out / v_in);

    clear; clc; close all;
% Connect to TDS2012B
scope = visadev("USB0::0x0699::0x0367::C103716::0::INSTR"); % Change ID to match your scope
fopen(scope);

% Get CH1
% Force CH1 ASCII and read raw string
writeline(scope, 'DATA:SOURCE CH1');
writeline(scope, 'DATA:ENCdg ASCII'); % full name, not ASC
writeline(scope, 'DATA:WIDTH 1');
rawX = query(scope, 'CURVE?');
x = str2double(split(strip(rawX), ','));

% Get CH2
writeline(scope, 'DATA:SOURCE CH2');
writeline(scope, 'DATA:ENCdg ASC');  % force ASCII again!
writeline(scope, 'DATA:WIDTH 1');
writeline(scope, 'DATA:START 1');
writeline(scope, 'DATA:STOP 2500');
rawY = query(scope, 'CURVE?');
y = str2double(split(strip(rawY), ','));

clear scope;

% Remove offsets
x = x - mean(x);
y = y - mean(y);

% Peak-to-peak in Y (A)
A = max(y) - min(y);

% Vertical intercept at X=0 (B)
% We'll interpolate to find where X is closest to 0
[~, idx0] = min(abs(x)); 
B = abs(y(idx0));

% Calculate phase in degrees
phase(k) = asin(B / A).*(180/pi);
    fprintf("Freq: %.0f Hz ; Vin: %.3f V ; Vout: %.3f V; Phase: %.0f degrees\n",f, v_in, v_out, phase(k));
end

%plot
figure;
frequencies_interp = logspace(log10(frequencies(1)), log10(frequencies(end)), length(frequencies)*10);
gain_dB_interp = interp1(frequencies, gain_dB, frequencies_interp, 'pchip');     % 'pchip' gives smooth curve

semilogx(frequencies_interp, gain_dB_interp, 'b-', 'LineWidth', 2); 
grid on;
xlabel('Frequency (Hz)'); 
ylabel('Magnitude (dB)');
ylim([-30,0]);
title('Transfer Function - Gain');

figure;
phase_interp = interp1(frequencies, phase, frequencies_interp, 'pchip');     % 'pchip' gives smooth curve

semilogx(frequencies_interp, phase_interp, 'b-', 'LineWidth', 2); 
grid on;
xlabel('Frequency (Hz)'); 
ylabel('Phase (degrees)');
ylim([-180,180]);
title('Transfer Function - Phase');

%clear the instruments
clear wavegen scope

%function to read channels
function [t, v] = readChannel(osc, chan)
    writeline(osc, ['DATA:SOURCE ' chan]);
    writeline(osc, 'DATA:ENC RPB');      % Fast binary
    writeline(osc, 'DATA:WIDTH 1');      % 1 byte per point

    writeline(osc, 'WFMPRE:XINCR?'); xIncr = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YOFF?');  yOff  = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YMULT?'); yMult = str2double(readline(osc));
    writeline(osc, 'WFMPRE:YZERO?'); yZero = str2double(readline(osc));

    writeline(osc, 'CURVE?');
    raw = readbinblock(osc, 'uint8');
    readline(osc); % Clear LF

    v = (double(raw) - yOff) * yMult + yZero;
    t = (0:length(v)-1) * xIncr;
end

function signal_filtered=fastFourier(signal,f, time)
%==========================================================================
% Time series setup
dt=time(2)-time(1); %compute sampling time from defined number of samples
% time=tmin:dt:tmax; %time array (s)
nt=length(time); %number of time samples

%computing FFT parameters
fs=1/dt; %sampling frequency (Hz)
fN=fs/2; %Nyquist frequency (Hz)
freq=linspace(-fN,fN,nt); %frequency array (must have same length as time)

%==========================================================================
% Fourier transform (with shift)

%compute shifted Fast Fourier Transform (shifted FFT)
FFT=fftshift(fftn(ifftshift(signal)))/sqrt(nt); %take FFT and shift
%create a custom top-hat filter
filter=zeros(1,nt);
f_low=0.9*f; %lower frequency cutoff threshold (high-pass)
f_high=1.1*f; %higher frequency cutoff threshold (low-pass)
%use both cutoff thresholds to create custom top-hat bandpass filter
for i=1:nt
    if abs(freq(i))>f_low && abs(freq(i))<f_high
        filter(i)=1;
    end
end

%apply filter to FFT of signal (convolution), then compute inverse FFT
FFT_filtered=FFT.*filter; %apply filter
%perform inverse FFT operations (IFFT) to extract cleaned-up signal (check this!!)
signal_filtered=fftshift(ifftn(ifftshift(FFT_filtered)))*sqrt(nt);
end