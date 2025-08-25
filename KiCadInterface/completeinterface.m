%DOWNLOAD ngspice
%AND ALSO:
%Copy the path from folder containing the ngspice exe
%Search in Windows: Environment Variables
%Click "Edit the system environment variables"
%Click the "Environment Variables..." button
%Under "System variables", find and select Path
%Click Edit...
%Click New and paste the path: 

%clear command window
clear; clc;
clf; close all;
warning off;
tic 

%define the original schematic path
original = 'LRCparallel.kicad_sch';

%define a new folder (path) and the new schematic filename
new_folder = '\LRCvariant\';
new_file = fullfile(new_folder, 'LRCvariant.kicad_sch');
spiceFile = 'example.cir';
raw_file = 'example.csv';

%define component values to update
component_values = struct('R1', '1k', 'C1', '100n', 'L1', '10m', 'L2', '100n', 'C2', '100m');

%===========================================================================
%make sure the destination folder exists
if ~exist(new_folder, 'dir')
    mkdir(new_folder);  %create it if not
end

%copy original schematic to new file
copyfile(original, new_file);

%read in the copied schematic file
fid = fopen(new_file, 'r', 'n', 'UTF-8');  %open file with UTF-8 encoding
lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = lines{1};  %extract cell array of lines

%prepare to process lines and store updated ones
new_lines = {};  %to store modified lines
inside_component = false;

%loop through all lines in the schematic
i = 1;
while i <= length(lines)
    raw_line = lines{i};
    line = strtrim(raw_line);  %removes leading/trailing spaces

    %check if the line is the beginning of a component (symbol)
    if startsWith(line, '(lib_id "Device:')
        inside_component = true;
        current_ref = '';
        j = i + 1;

        %look for the reference property in the next few lines
        while j <= length(lines) && ~contains(lines{j}, '(property "Reference"')
            j = j + 1;
        end
        if j <= length(lines)
            ref_line = lines{j};
            tokens = regexp(ref_line, '"([^"]+)"', 'tokens');
            if length(tokens) >= 2
                current_ref = tokens{2}{1};  %extract component reference (e.g., R1)
            end
        end

        %look for the Value property after Reference
        while j <= length(lines) && ~contains(lines{j}, '(property "Value"')
            j = j + 1;
        end

        %to replace matched component ref
        if isfield(component_values, current_ref) && j <= length(lines)
            new_val = component_values.(current_ref);
            %preserve indentation from original line
            indent = regexp(lines{j}, '^\s*', 'match', 'once');
            lines{j} = sprintf('%s(property "Value" \"%s\"', indent, new_val);
        end
    end

    %append current line to output (whether changed or not)
    new_lines{end+1, 1} = lines{i};
    i = i + 1;
end

%write updated lines back into the new schematic
fid = fopen(new_file, 'w', 'n', 'UTF-8');
for i = 1:length(new_lines)
    fprintf(fid, '%s\n', new_lines{i});
end
fclose(fid);
fprintf('New schematic created at: %s\n', new_file);


nodeTest = '(node (ref "C1") (pin "2") (pintype "passive"))';
match = regexp(nodeTest, '\(node\s+(?:\([^)]+\)\s*)+\)', 'match');
disp(match);
    
exportSpiceFromKicad(new_file, spiceFile, raw_file);


%run Ngspice based on our exported netlist
run_ngspice_command = sprintf('ngspice -b %s', spiceFile);
[status2, result] = system(run_ngspice_command);
if status2 ==0
    fprintf('Simulation successful! \n');
else
    error('Failed to simulate netlist: %s', result);
end

%to import the csv data as a matrix
%the matrix has 3 columns starting at row 1 and column 1
%col 1 has the frequency
%col 2 has the real part of Vout
%col 3 has the imaginary part of Vout

M=dlmread(raw_file,'',0,0);

%compute relevant transfer function characteristics
f=M(:,1);% Frequency in Hz
Vcomplex=M(:,2)+1i*M(:,3);% Construct complex v(out)
Vout=20*log10(abs(Vcomplex));% Gain in dB
Vphase=angle(Vcomplex)*(180/pi);% Phase in degrees

%Theoretical values of the circuit (to be graphed for comparison)
R=1559.5;
C=5.85*10^(-9);
L=0.068*10^(-3);

%Transfer Functions
amp_function = 1 ./ sqrt(1 + (R .* (2*pi*f*C - 1./(2*pi*f*L))).^2);
phase_function=180/pi.*(-atan( 2*pi*f*R*C - R./(2*pi*f*L)));
cutoff = -3;

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
    hold on;
    semilogx(f,Vout,'--w','Linewidth',LW);
    hold on;
    semilogx(f,cutoff * ones(size(f)), '--b', 'Linewidth',LW);
    axis([1,1e6,-100,0]);
    legend('New', 'Original', 'Location', 'northoutside');
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title(['Discrepancy Between Two Amplitude Functions'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    set(gca, 'XTick', 10.^(0:6)); %1 to 1MEG Hz
    setFigureLightMode();
    grid on;

%Second plot has the first plot's difference vs. frequency
    subplot(2,2,2);
    semilogx(f,Vout-20.*log10(amp_function), '-g', 'Linewidth',LW);
    axis([1,1e6,-2,2]);
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel('Gain (dB)','Fontsize',FS);
    title('Difference of Amplitude as Function','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    setFigureLightMode();
    axis tight;
    set(gca, 'XTick', 10.^(0:6)); %1 to 1MEG Hz
    grid on;

%Third plot has theoretical and nominal of phase function
    subplot(2,2,3);
    semilogx(f,phase_function,'-r','Linewidth',LW);
    axis([1,1e6,-100,100]);
    hold on;
    semilogx(f,Vphase,'--w','Linewidth',LW);
    legend('New', 'Original', 'Location', 'northoutside');
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel(['Phase (' char(176), ')'],'Fontsize',FS);
    title('Discrepancy Between Two Phase Functions','Fontsize',FST);
    set(gca,'Fontsize',FSN);
    setFigureLightMode();
    set(gca, 'XTick', 10.^(0:6)); %1 to 1MEG Hz
    grid on;

%Last plot has third plot's difference vs frequency
    subplot(2,2,4);
    semilogx(f,Vphase-phase_function,'-g','Linewidth',LW);
    axis([1,1e6,-100,100]);
    xlabel('Frequency, {\itf} (Hz)','Fontsize',FS);
    ylabel(['Phase (' char(176), ')'],'Fontsize',FS);
    title(['Difference of Phase as Functions'],'Fontsize',FST);
    set(gca,'Fontsize',FSN);
    grid on;
    setFigureLightMode(); 
    axis tight;
    set(gca, 'XTick', 10.^(0:6)); %1 to 1MEG Hz
    drawnow;

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
                line = sprintf('%s %s %s AC 1', ref, pins{1}, pins{2});
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
    fprintf(fid, '.ac dec 100 1 1e6\n');
    fprintf(fid, '.control\n');
    fprintf(fid, 'run\n');
    fprintf(fid, 'wrdata %s v(3) \n', raw_file);  % CSV-like output
    fprintf(fid, 'quit\n');
    fprintf(fid, '.endc\n');
    fprintf(fid, '.end\n')
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


function setFigureLightMode(fig)
    % If no figure is passed, use current figure
    if nargin < 1
        fig = gcf;
    end

    % Set figure background to white
    set(fig, 'Color', 'white');

    % Apply settings to all axes in the figure
    ax = findall(fig, 'Type', 'axes');
    for i = 1:length(ax)
        set(ax(i), 'Color', 'white', ...               % Axes background
                   'XColor', 'black', 'YColor', 'black', 'ZColor', 'black', ...  % Axes lines
                   'GridColor', [0.6 0.6 0.6], ...      % Light gray grid
                   'MinorGridColor', [0.8 0.8 0.8], ...
                   'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on', ...
                   'FontSize', 12, ...
                   'FontWeight', 'normal');
    end

    % Set legend and colorbar background/text if they exist
    legends = findall(fig, 'Type', 'legend');
    for l = legends'
        set(l, 'TextColor', 'black', 'Color', 'white', 'EdgeColor', 'black');
    end

    colorbars = findall(fig, 'Type', 'ColorBar');
    for cb = colorbars'
        set(cb, 'Color', 'black', 'FontSize', 12);
    end
    titles = findall(fig, 'Type', 'text', 'Tag', '');  % Title, xlabel, ylabel, etc.
    for t = titles'
         set(t, 'Color', 'black');
    end
end
toc