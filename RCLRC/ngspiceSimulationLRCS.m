%DOWNLOAD ngspice
%AND ALSO:
%Copy the path of the folder containing the ngspice exe
%Search in Windows: Environment Variables
%Click "Edit the system environment variables"
%Click the "Environment Variables..." button
%Under "System variables", find and select Path
%Click Edit...
%Click New and paste the path:
%Then restart matlab and this code will work now

clear all; clc;
clf; close all;
warning off;
tic

%make/edit file names
netlist = 'LRCseries.cir';
raw_file = 'LRCseries.csv';

%edit theoretical values for the theoretical circuit
R=4.3;
C=14.177*10^(-9);
L=28.351*10^(-6);

%edit SPICE netlist for simulation (also edit nominal values)
R1='4.3';
C1='14.7n';
L1='27.7u';
fid = fopen(netlist,'w','native','UTF-8');
fprintf(fid, '* RC Low-Pass Filter AC Sweep\n');
fprintf(fid, 'Vin in GND AC 1\n');
fprintf(fid, 'L1 in aa %s\n',L1);
fprintf(fid, 'C1 aa out %s\n',C1);
fprintf(fid, 'R1 out GND %s\n',R1);
fprintf(fid, '.ac dec 100 1 1e6\n');
fprintf(fid, '.control\n');
fprintf(fid, 'run\n');
fprintf(fid, 'wrdata %s v(out) \n', raw_file );  % CSV-like output
fprintf(fid, 'quit\n');
fprintf(fid, '.endc\n');
fprintf(fid, '.end\n');
fclose(fid);

%=======================================================================
%run ngspice
[status, result] = system(sprintf('ngspice -b %s', netlist));

%lets compiler know if there is error
if status ~= 0
    disp('Ngspice failed!');
    disp(result);  % show error message
else
    disp('Ngspice ran successfully!');
end

%to import the csv data as a matrix
%the matrix has 3 columns starting at row 1 and column 1
%col 1 has the frequency
%col 2 has the real part of Vout
%col 3 has the imaginary part of Vout

M=dlmread(raw_file,'',0,0);

% Figure parameters =====================================================%
def_size = 1; %switch to return to default figure size
LW  = 1; %linewidth for 2D plots
FS	= 14;  %fontsize for axis labels
FSN = 14;  %fontsize for tick mark labels
FST = 16;  %fontsize for figure title
frcx = 0.45; %fraction of screen length along x-direction
frcy = 0.55; %fraction of screen length along y-direction
rpy = 0.1; %range percentage for y-axis limits (shortcut)
fullscreen = get(0,'ScreenSize'); %get dimensions of monitor screen
pxi=(1-frcx)*fullscreen(3)/2; %lower left corner x-coordinate
pyi=(1-frcy)*fullscreen(4)/2; %lower left corner y-coordinate
pxs=frcx*fullscreen(3); %custom base dimension of figure
pys=frcy*fullscreen(4); %custom height dimension of figure

%compute relevant transfer function characteristics
f=M(:,1);% Frequency in Hz
Vcomplex=M(:,2)+1i*M(:,3);% Construct complex v(out)
Vout=20*log10(abs(Vcomplex));% Gain in dB
Vphase=angle(Vcomplex)*(180/pi);% Phase in degrees

%transfer functions for LRCS circuits
amp_function =R ./ sqrt(R.^2 + (2*pi*f*L - 1./(2*pi*f*C)).^2);
phase_function =(180/pi).*( -atan( (2*pi*f*L./R) - (1./(2*pi*f*R*C)) ));

%plot the data
    clf;
    if def_size==1
    fig1=figure(1);
    else
	fig1=figure('Position',[pxi,pyi,pxs,pys]);
end
%First plot has both theoretical and nominal of amplitude function
    subplot(2,2,1);
    semilogx(f,20*log10(amp_function), '-b', 'Linewidth',LW);
    hold on;
    semilogx(f,Vout,'--w','Linewidth',LW);
    axis([1,1e6,-100,0]);
    legend('Theoretical', 'Nominal', 'Location', 'northoutside');
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title(['Discrepancy Between Two Amplitude Functions'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

%Second plot has the first plot's difference vs. frequency
    subplot(2,2,2);
    semilogx(f,Vout-20*log10(amp_function), '-g', 'Linewidth',LW);
    axis([1,1e6,-1,1]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title(['Difference of Amplitude as Function'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

%Third plot has theoretical and nominal of phase function
    subplot(2,2,3);
    axis([1,1e6,-100,100]);
    semilogx(f,phase_function,'-b','Linewidth',LW);
    hold on;    
    semilogx(f,Vphase,'--w','Linewidth',LW);
    legend('Theoretical', 'Nominal', 'Location', 'northoutside');
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Phase (\deg)','Fontsize',FS);
    title(['Discrepancy Between Two Phase Functions'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

%Last plot has third plot's difference vs frequency
    subplot(2,2,4);
    semilogx(f,Vphase-phase_function,'-g','Linewidth',LW);
    axis([1,1e6,-100,100]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Phase (\deg)','Fontsize',FS);
    title(['Difference of Phase as Functions'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;
    drawnow;

toc
