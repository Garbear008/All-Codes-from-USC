%clear command window
clear; clc;
clf; close all;
warning off;
tic

schematicFile = 'butterworthnominal.kicad_sch';
spiceFile = 'butterworth.cir';
raw_file = 'butterworth.csv';

function exportSpiceFromKicad(schematicFile, spiceFile, raw_file)
    % Temporary netlist file
    netlistFile = 'temp_netlist.net';

    % Fix paths for Windows
    schematicFile = strrep(schematicFile, '\', '/');
    netlistFile = strrep(netlistFile, '\', '/');
    spiceFile = strrep(spiceFile, '\', '/');

    % Step 1: Export KiCad netlist (XML)
    cmd = sprintf('kicad-cli sch export netlist "%s" -o "%s"', schematicFile, netlistFile);
    [status, output] = system(cmd);
    if status ~= 0
        error("KiCad netlist export failed:\n%s", output);
    end
    % Step 2: Export SPICE netlist from the temporary netlist file
    exportSpiceFromKicadNetfile(netlistFile, spiceFile, raw_file);
end

function exportSpiceFromKicadNetfile(netFile, spiceFile, raw_file)
    % Read entire netlist text
    txt = fileread(netFile);

    % Extract components
    compPattern = '\(comp\s+\(ref\s+"(?<ref>[^"]+)"\).*?\(value\s+"(?<value>[^"]+)"\)';
    comps = regexp(txt, compPattern, 'names');

    % Extract full balanced (net ...) blocks
    netBlockStrings = extractNetBlocks(txt);
    netBlocks = cell(1, length(netBlockStrings));

    for i = 1:length(netBlockStrings)
        block = netBlockStrings{i};
        
        % Extract net name
        nameMatch = regexp(block, '\(name\s+"([^"]+)"\)', 'tokens', 'once');
        if isempty(nameMatch)
            continue
        end
        netName = nameMatch{1};
        
        netBlocks{i} = {netName, block};
    end

    % Build net-to-node mapping
    netMap = containers.Map();
    nodeCounter = 1;
    for i = 1:length(netBlocks)
        if isempty(netBlocks{i}), continue; end
        name = netBlocks{i}{1};
        if strcmpi(name, 'GND') || strcmpi(name, '0')
            netMap(name) = '0';
        else
            netMap(name) = num2str(nodeCounter);
            nodeCounter = nodeCounter + 1;
        end
    end

    % Build ref+pin to node map
    refPinToNode = containers.Map();
    for i = 1:length(netBlocks)
        if isempty(netBlocks{i}), continue; end
        name = netBlocks{i}{1};
        body = netBlocks{i}{2};

        % Find all (node ...) inside the net block
        nodeBlocks = regexp(body, '\(node\s+(?:\([^)]+\)\s*)+\)', 'match');

        fprintf('Net "%s": found %d node(s)\n', name, length(nodeBlocks));

        for k = 1:length(nodeBlocks)
            nodeText = nodeBlocks{k};
            % Extract ref and pin
            refMatch = regexp(nodeText, '\(ref\s+"([^"]+)"\)', 'tokens', 'once');
            pinMatch = regexp(nodeText, '\(pin\s+"([^"]+)"\)', 'tokens', 'once');
            if ~isempty(refMatch) && ~isempty(pinMatch)
                ref = refMatch{1};
                pin = pinMatch{1};
                key = sprintf('%s:%s', ref, pin);
                refPinToNode(key) = netMap(name);
                fprintf('Mapping %s:%s to node %s\n', ref, pin, netMap(name));
            end
        end
    end

    disp('Mapped pins:');
    disp(keys(refPinToNode));

    % Create SPICE netlist
    spiceLines = {};
    for i = 1:length(comps)
        ref = comps(i).ref;
        value = comps(i).value;
        pins = getPinsForComponent(refPinToNode, ref);
        if length(pins) < 2
            warning('Component %s has <2 pins. Skipping.', ref);
            continue;
        end

        % SPICE line
        switch upper(ref(1))
            case {'R', 'C', 'L'}
                line = sprintf('%s %s %s %s', ref, pins{1}, pins{2}, value);
            case 'V'
                line = sprintf('%s %s %s AC 2', ref, pins{1}, pins{2});
            otherwise
                warning('Unsupported component: %s', ref);
                continue;
        end

        spiceLines{end+1} = line; %#ok<AGROW>
    end

    % Write to file
    fid = fopen(spiceFile, 'w','native','UTF-8');
    fprintf(fid, '* SPICE netlist generated from KiCad .net\n');
    fprintf(fid, '.title Auto-generated netlist\n');
    for i = 1:length(spiceLines)
        fprintf(fid, '%s\n', spiceLines{i});
    end
    fprintf(fid, '.ac dec 100 1e5 1e9\n');
    fprintf(fid, '.control\n');
    fprintf(fid, 'run\n');
    fprintf(fid, 'wrdata %s v(5) \n', raw_file);  % CSV-like output
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
    % Return node numbers for pins "1" and "2"
    pins = {};
    for pin = ["1", "2"]
        key = sprintf('%s:%s', ref, pin);
        if isKey(map, key)
            pins{end+1} = map(key); %#ok<AGROW>
        else
            pins{end+1} = '999'; % fallback node
        end
    end
end

nodeTest = '(node (ref "C1") (pin "2") (pintype "passive"))';
match = regexp(nodeTest, '\(node\s+(?:\([^)]+\)\s*)+\)', 'match');
disp(match);
    
exportSpiceFromKicad(schematicFile, spiceFile, raw_file);

% Optional: run with ngspice
[status, result] = system(sprintf('ngspice -b "%s"', spiceFile));
disp(result);

%to import the csv data as a matrix
%the matrix has 3 columns starting at row 1 and column 1
%col 1 has the frequency
%col 2 has the real part of Vout
%col 3 has the imaginary part of Vout
M=dlmread(raw_file,'',0,0);

% Figure parameters =====================================================%
def_size = 1; %switch to return to default figure size
LW  = 1; %linewidth for 2D plots
FS	= 12;  %fontsize for axis labels
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
f=10^(-6).*M(:,1);% Frequency in Hz
Vcomplex=M(:,2)+1i*M(:,3);% Construct complex v(out)
Vout=20*log10(abs(Vcomplex));% Gain in dB
Vphase=(angle(Vcomplex));% Phase in degrees
Vphase = (unwrap(Vphase))*(180/pi);

%transfer Functions for butterworth
cutoff_frequency = 80*10^6;
n=5;
amp_function = 1 ./ sqrt(1 + (f./(80)).^(2*n));
x = f / cutoff_frequency;                % Normalized frequency

% Exact real-valued phase expression (in degrees)
phase_function = -(180/pi) * ( ...
    atan((x - 0.9511) ./ 0.3090) + ...
    atan((x - 0.5878) ./ 0.8090) + ...
    atan((x) ./ 1.0000) + ...
    atan((x + 0.5878) ./ 0.8090) + ...
    atan((x + 0.9511) ./ 0.3090) ...
);


%plot the data
    clf;
    if def_size==1
    fig1=figure(1);
    else
	fig1=figure('Position',[pxi,pyi,pxs,pys]);
end

%First plot has both theoretical and nominal of amplitude function
    subplot(2,2,1);
    semilogx(f,20.*log10(amp_function), '-r', 'Linewidth',LW);
    grid on;
    hold on;
    semilogx(f,Vout,'--k','Linewidth',LW);
    axis([1,1e3,-100,0]);
    legend('Theoretical', 'Nominal', 'Location', 'northoutside');
    xlabel('Frequency, {\itf} (MHz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title('Discrepancy Between Amplitudes','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

%Second plot has the first plot's difference vs. frequency
    subplot(2,2,2);
    semilogx(f,abs(Vout-20.*log10(amp_function)), '-g', 'Linewidth',LW);
    axis([1e6,1e9,0,5]);
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title('Difference of Amplitude as Function','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;

%Third plot has theoretical and nominal of phase function
    subplot(2,2,3);
    semilogx(f,phase_function,'-r','Linewidth',LW);
    axis([1e6,1e9,-450,0]);
    grid on;
    hold on;
    semilogx(f,Vphase,'--k','Linewidth',LW);
    legend('Theoretical', 'Nominal', 'Location', 'northoutside');
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel(['Phase (' char(176) ')'],'Fontsize',FS);
    title('Discrepancy Between Phases','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;
    

%Last plot has third plot's difference vs frequency
    subplot(2,2,4);
    semilogx(f,abs(Vphase-phase_function),'-g','Linewidth',LW);
    axis([1e6,1e9,0,25]);
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel(['Phase (' char(176) ')'],'Fontsize',FS);
    title('Difference of Phase as Functions','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;
    
    drawnow;

toc