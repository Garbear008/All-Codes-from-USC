%TO USE kicad-cli and ngspice
%Copy the path from the kicad_cli.exe
%Search in Windows: Environment Variables
%Click "Edit the system environment variables"
%Click the "Environment Variables..." button
%Under "System variables", find and select Path
%Click Edit...
%Click New and paste the path:
%repeat for ngspice.exe
%Restart Matlab

clear; clc; %clr all variables, clr command window
clf; close all; % clr fig, close all fig windows
warning off;
tic % start timer
set(0,'defaultAxesFontSize', 16);

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


% === USER CONFIGURABLE VARIABLES ===
schematicFile = 'OpAmp.kicad_sch'; %ensure your schematic has Vout and Vin labeled
spiceFile = 'InvertingAmp.cir';
raw_file = 'opAmp.csv';

analysisType = 'tran'; % 'tran' or 'ac'
tran_params = '1m 5'; % <timestep> <tstop> [tstart [tmaxstep]]
ac_params = 'dec 100 1 1e6'; % <dec|oct|lin> <numpoints> <fstart (Hz)> <fstop (Hz)>

spiceIncludeFile = 'C:\USCProjects\KiCad\9.0\symbols\kicad-symbols-master\Simulation_SPICE.sp'; %Simulation_SPICE.sp necessary to simulate opAmp

% The default voltage source value string used if none is found in the schematic (it will be mentioned in output)
defaultVSourceValue = 'SIN(0 1 5 0 0 0) AC 1'; %'SIN(<offset voltage> <amplitude> <Freq (Hz)> <Time delay (s)> <Theta (damping factor> <Phi (phase, degrees)>) [AC <mag>]'
%Any VSIN will use this ^^^^^



exportSpiceFromKicad(schematicFile, spiceFile, raw_file, analysisType, tran_params, ac_params, spiceIncludeFile, defaultVSourceValue);

% Run with ngspice
[status, result] = system(sprintf('ngspice -b "%s"', spiceFile));
disp(result);

M=dlmread(raw_file,'',0,0);
if strcmp(analysisType, 'tran')
    time = M(:, 1);
    vout = M(:, 2);
    vin = M(:, 4);
end

%compute relevant transfer function characteristics for ac
if strcmp(analysisType, 'ac')
    f=M(:,1);% Frequency in Hz
    Vcomplex=M(:,2)+1i*M(:,3);% Construct complex v(out)
    Vout=20*log10(abs(Vcomplex));% Gain in dB
    Vphase=angle(Vcomplex)*(180/pi);% Phase in degrees
end

%tran plot
if strcmp(analysisType, 'tran')
    if def_size==1
    fig1=figure(1);
        else
	    fig1=figure('Position',[pxi,pyi,pxs,pys]);
    end
    plot(time, vin, 'r--', 'LineWidth', LW); 
    hold on;
    plot(time, vout, 'b', 'LineWidth', LW);
    grid on;
    xlabel('Time (s)', 'FontSize', FS);
    ylabel('Voltage (V)', 'FontSize', FS);
    title('Transient Response of Operational Amplifier', 'FontSize', FST);
    legend('V_{in}', 'V_{out}');
    set(gca, 'FontSize', FSN);
%ac plots
elseif strcmp(analysisType, 'ac')
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
    axis([1,1e6,-100,100]);
    xlabel('Frequency, f (Hz)','Fontsize',FS);
    ylabel('Phase (\deg)','Fontsize',FS);
    set(gca,'Fontsize',FSN);
    grid on;
    drawnow;
else
    error('Unsupported analysis type');
end

toc % stop timer


% Functions
function exportSpiceFromKicad(schematicFile, spiceFile, raw_file, analysisType, tran_params, ac_params, spiceIncludeFile, defaultVSourceValue)
    netlistFile = 'temp_netlist.net';

    cmd = sprintf('kicad-cli sch export netlist "%s" -o "%s"', schematicFile, netlistFile);
    [status, output] = system(cmd);
    if status ~= 0
        error("KiCad netlist export failed:\n%s", output);
    end

    exportSpiceFromKicadNetfile(netlistFile, spiceFile, raw_file, analysisType, tran_params, ac_params, spiceIncludeFile, defaultVSourceValue);
end


function exportSpiceFromKicadNetfile(netFile, spiceFile, raw_file, analysisType, tran_params, ac_params, spiceIncludeFile, defaultVSourceValue)
    txt = fileread(netFile);

    compPattern = '\(comp\s+\(ref\s+"(?<ref>[^"]+)"\).*?\(value\s+"(?<value>[^"]+)"\).*?(?:\(property\s+"Spice_Model"\s+"(?<spicemodel>[^"]+)"\))?';
    comps = regexp(txt, compPattern, 'names');

    netBlockStrings = extractNetBlocks(txt);
    netBlocks = cell(1, length(netBlockStrings));

    for i = 1:length(netBlockStrings)
        block = netBlockStrings{i};
        nameMatch = regexp(block, '\(name\s+"([^"]+)"\)', 'tokens', 'once');
        if isempty(nameMatch), continue; end
        netName = nameMatch{1};
        netBlocks{i} = {netName, block};
    end

    % Map nets to safe SPICE node names
    netMap = containers.Map();
    for i = 1:length(netBlocks)
        if isempty(netBlocks{i}), continue; end
        name = netBlocks{i}{1};
        if strcmpi(name, 'GND') || strcmpi(name, '0')
            netMap(name) = '0';
        else
            safeName = regexprep(name, '[^a-zA-Z0-9_]', '_');
            netMap(name) = safeName;
        end
    end

    % Map ref:pin to net name
    refPinToNode = containers.Map();
    for i = 1:length(netBlocks)
        if isempty(netBlocks{i}), continue; end
        name = netBlocks{i}{1};
        body = netBlocks{i}{2};
        nodeBlocks = regexp(body, '\(node\s+(?:\([^)]+\)\s*)+\)', 'match');

        for k = 1:length(nodeBlocks)
            nodeText = nodeBlocks{k};
            refMatch = regexp(nodeText, '\(ref\s+"([^"]+)"\)', 'tokens', 'once');
            pinMatch = regexp(nodeText, '\(pin\s+"([^"]+)"\)', 'tokens', 'once');
            if ~isempty(refMatch) && ~isempty(pinMatch)
                ref = refMatch{1};
                pin = pinMatch{1};
                key = sprintf('%s:%s', ref, pin);
                refPinToNode(key) = netMap(name);
            end
        end
    end

    % Build SPICE netlist lines
    spiceLines = {};
    hasOpamp = false;
    for i = 1:length(comps)
        ref = comps(i).ref;
        value = comps(i).value;
        pins = getPinsForComponent(refPinToNode, ref);

        switch upper(ref(1))
            case {'R','C','L'}
                if length(pins) < 2, continue; end
                line = sprintf('%s %s %s %s', ref, pins{1}, pins{2}, value);

            case 'V'
                if length(pins) < 2, continue; end
                val = '';
                if isfield(comps(i), 'spicemodel') && ~isempty(comps(i).spicemodel)
                    val = comps(i).spicemodel;
                elseif ~isempty(value)
                    val = value;
                else
                    val = 'DC 0';
                end
                % If the value is exactly "VSIN" (case-insensitive), replace with default voltage source string
                if strcmpi(val, 'VSIN')
                    val = defaultVSourceValue;
                    fprintf('%s value defaulted to %s\n', ref, val);
                end
                line = sprintf('%s %s %s %s', ref, pins{1}, pins{2}, val);

            case 'U' % Op-amp
                if length(pins) < 5
                    warning('Op-amp %s has fewer than 5 pins. Skipping.', ref);
                    continue;
                end
                line = sprintf('X%s %s %s %s %s %s kicad_builtin_opamp', ...
                    ref, pins{1}, pins{2}, pins{3}, pins{4}, pins{5});
                hasOpamp = true;

            otherwise
                warning('Unsupported component: %s', ref);
                continue;
        end

        spiceLines{end+1} = line; %#ok<AGROW>
    end

    % Write SPICE netlist file
    fid = fopen(spiceFile, 'w', 'native', 'UTF-8');
    fprintf(fid, '* SPICE netlist generated from KiCad .net\n');
    fprintf(fid, '.title KiCad schematic\n');
    fprintf(fid, '.include "%s"\n', spiceIncludeFile);

    for i = 1:length(spiceLines)
        fprintf(fid, '%s\n', spiceLines{i});
    end

    if strcmp(analysisType, 'tran')
        fprintf(fid, '.tran %s\n', tran_params);
    elseif strcmp(analysisType, 'ac')
        fprintf(fid, '.ac %s\n', ac_params);
    else
        error('Unsupported analysis type');
    end

    fprintf(fid, '.control\n');
    fprintf(fid, 'run\n');
    fprintf(fid, 'wrdata %s v(OUT) v(IN)\n', raw_file);
    fprintf(fid, 'quit\n');
    fprintf(fid, '.endc\n');
    fprintf(fid, '.end\n');
    fclose(fid);

    fprintf('SPICE netlist written to: %s\n', spiceFile);
end


function netBlocks = extractNetBlocks(txt)
    startIdxs = strfind(txt, '(net ');
    netBlocks = cell(1, length(startIdxs));
    for i = 1:length(startIdxs)
        startPos = startIdxs(i);
        pos = startPos;
        depth = 0;
        while pos <= length(txt)
            if txt(pos) == '('
                depth = depth + 1;
            elseif txt(pos) == ')'
                depth = depth - 1;
                if depth == 0
                    break
                end
            end
            pos = pos + 1;
        end
        netBlocks{i} = txt(startPos:pos);
    end
end


function pins = getPinsForComponent(map, ref)
    pins = {};
    pinNum = 1;
    while true
        key = sprintf('%s:%d', ref, pinNum);
        if isKey(map, key)
            pins{end+1} = map(key); %#ok<AGROW>
            pinNum = pinNum + 1;
        else
            break;
        end
    end
    if isempty(pins)
        pins = repmat({'0'}, 1, 2);
    end
end
