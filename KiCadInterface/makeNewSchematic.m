
%clear command window
clear; clc;
clf; close all;
warning off;
tic 

%define the original schematic path
original = 'butterworth.kicad_sch';

%define a new folder (path) and the new schematic filename
new_folder = 'C:\USCProjects\KiCad\9.0\projects\Passive Butterworth\butterworth\';
new_file = fullfile(new_folder, 'butterworthnominal.kicad_sch');

%define component values to update
component_values = struct('Rsource1', '50', 'Rload1','50','C1', '62p', 'L1', '56n', 'L2', '180n', 'C2', '62p', 'L3','56n');

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

toc