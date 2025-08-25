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
x = str2double(split(strip(rawX), ','))

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
phi_rad = asin(B / A);
phi_deg = phi_rad * (180/pi);

fprintf('Phase magnitude from Lissajous: %.2f degrees\n', phi_deg);
