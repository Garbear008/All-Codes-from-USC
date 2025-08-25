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
netlist = 'RC1R.cir';
raw_file = 'RC1R.csv';

%edit theoretical values and transfer function for RC
R=1000;
C=10^(-6);

%create a SPICE netlist for simulation
R1='1k';
C1='1u';
fid = fopen(netlist,'w','native','UTF-8');
fprintf(fid, '* RC Low-Pass Filter AC Sweep\n');
fprintf(fid, 'Vin in GND AC 1\n');
fprintf(fid, 'R1 in out %s\n',R1);
fprintf(fid, 'C1 out GND %s\n',C1);
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
transfer_function= 1 ./ sqrt(1 + ((2*pi.*f) .* R .* C).^2);
Phase= angle(1 ./ (1 + 1i * 2*pi.*f * R * C)) * (180 / pi);

%plot the data
    clf;
    if def_size==1;
    fig1=figure(1);
    else
	fig1=figure('Position',[pxi,pyi,pxs,pys]);
end
    subplot(2,1,1);
    semilogx(f,Vout,'-b','Linewidth',LW);
    axis([1,1e6,-100,0]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title(['AC Sweep Analysis'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

    subplot(2,1,2);
    semilogx(f,Vphase,'-g','Linewidth',LW);
    axis([1,1e6,-100,0]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Phase (\deg)','Fontsize',FS);
    set(gca,'Fontsize',FSN);
    grid on;

    drawnow;

toc
