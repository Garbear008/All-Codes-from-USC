clear all; clc;
clf; close all;
warning off;
tic
%this is code to graph AC signal analysis for filters
%csv must only contain the frequency,v(phase), and v(output)

%get the csv file
csvfile='.csv';

%to import the csv data as a matrix
%the matrix has 3 columns starting at row 2 and column 1
%row 1 has the labels

M=dlmread(csvfile,'',1,0);

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
Vphase=M(:,2);
Vout=M(:,3);
f=M(:,1);

%plot the data
    clf;
    if def_size==1
    fig1=figure(1);
    else
	fig1=figure('Position',[pxi,pyi,pxs,pys]);
end
    subplot(2,1,1);
    semilogx(f,Vout,'-b','Linewidth',LW);
    axis([1,1e6,-100,0]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    set(gca,'Fontsize',FSN);
    grid on;
    subplot(2,1,2);
    semilogx(f,Vphase,'-g','Linewidth',LW);
    axis([1,1e6,-100,0]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Phase (\deg)','Fontsize',FS);
    set(gca,'Fontsize',FSN);
    grid on;
    title(['AC Sweep Analysis'],'Fontsize',FST);
    drawnow;
toc
