function prt = bv_prt_write(prt, varargin)
%function prt = bv_prt_write(prt, varargin)
%Writes a Brainvoager protocol file (.prt) based on the input structure prt
%which must contain the following fields:
%.name (the full name of the protocol file to be created)
%.NrOfConditions
%.Condition{:} (cell array of conditions, each one containing the fields
%               name, ntpts, estart, eend, color)
%SEE ALSO bv_read_prt.m
%Jens.Schwarzbach@form.unitn.it 20071121

if ~isfield(prt, 'FileVersion'), prt.FileVersion = 2; end
% if ~isfield(prt, 'ResolutionOfTime'), prt.ResolutionOfTime = 'volumes'; end
if ~isfield(prt, 'ResolutionOfTime'), prt.ResolutionOfTime = 'msec'; end
if ~isfield(prt, 'Experiment'), prt.Experiment = '<none>'; end
if ~isfield(prt, 'BackgroundColor'), prt.BackgroundColor = [0 , 0, 0]; end
if ~isfield(prt, 'TextColor'), prt.TextColor = [255, 255, 255]; end
if ~isfield(prt, 'TimeCourseColor'), prt.TimeCourseColor = [255, 255, 255]; end
if ~isfield(prt, 'TimeCourseThick'), prt.TimeCourseThick = 2; end
if ~isfield(prt, 'ReferenceFuncColor'), prt.ReferenceFuncColor = [0, 255, 0]; end
if ~isfield(prt, 'ReferenceFuncThick'), prt.ReferenceFuncThick = 1; end


cfg = [];
if ~isempty(varargin)
    cfg = varargin{1};
end
if ~isfield(cfg, 'CreateEmptyCondition'), cfg.CreateEmptyCondition = 0; end
if ~isfield(cfg, 'AlwaysOverwrite'), cfg.AlwaysOverwrite = 0; end



%CHECK FOR EXISTENCE OF OUTPUT FILE
if ~cfg.AlwaysOverwrite
    if exist(prt.name, 'file') == 2
        msg = sprintf('File %s already exists. Overwrite?', prt.name);
        ButtonName=questdlg(msg, 'Warning', 'Yes','No','No');
        if strcmp(ButtonName, 'No')
            prt = [];
            return
        end
    end
end


%CREATE AN EXTRA CONDITION FROM EMPTY ENTRIES
if cfg.CreateEmptyCondition
    switch strtok(prt.ResolutionOfTime)
        case 'msec'
            tres = 1000*cfg.TR;
        case 'volumes'
            tres = 1;
    end
    bl.vec = ones(cfg.VolumesInRun*tres, 1);
    
    %CREATE A CONDITION FROM ALL EMPTY TIMEPOINTS
    estart = [];
    eend = [];
    for c = 1: prt.NrOfConditions
        estart = [estart(:); prt.Condition{c}.estart(:)];
        eend = [eend(:); prt.Condition{c}.eend(:)];
        for t = 1:length(estart)
            bl.vec(estart(t):eend(t)) = 0;
        end
    end
    cases = find(bl.vec == 1);
    find(diff(cases) ~= 1)
end

%WRITE PROTOCOL FILE
fid = fopen(prt.name, 'wt');
if fid == -1
    msg = sprintf('The file %s cannot be opened.', prt.name);
    errordlg(msg, 'FileOpen Error', 1);
    return
end

fprintf(1, 'WRITING %s ...', prt.name)

%[dummy, prt.FileVersion] = fscanf(fid, '%s %d', 1);
fprintf(fid, '\n');
fprintf(fid, 'FileVersion:         %d\n\n', prt.FileVersion);
fprintf(fid, 'ResolutionOfTime: %s\n\n', prt.ResolutionOfTime);
fprintf(fid, 'Experiment: %s\n\n', prt.Experiment);
fprintf(fid, 'BackgroundColor:     %3d %3d %3d\n', prt.BackgroundColor);
fprintf(fid, 'TextColor:           %3d %3d %3d\n\n', prt.TextColor);
fprintf(fid, 'TimeCourseColor:     %3d %3d %3d\n', prt.TimeCourseColor);
fprintf(fid, 'TimeCourseThick:         %d\n', prt.TimeCourseThick);
fprintf(fid, 'ReferenceFuncColor:  %3d %3d %3d\n', prt.ReferenceFuncColor);
fprintf(fid, 'ReferenceFuncThick:         %d\n\n', prt.ReferenceFuncThick);
fprintf(fid, 'NrOfConditions:         %d\n\n', prt.NrOfConditions);

%LOOP THROUGH CONDITIONS
for c = 1:prt.NrOfConditions
    fprintf(fid, '%s\n', prt.Condition{c}.name);
    fprintf(fid, '%d\n', prt.Condition{c}.ntpts);
    if prt.Condition{c}.ntpts == 0
        %EMPTY CONDITIONS PRODUCE A WARNING BUT ARE NEVERTHELESS WRITTEN TO
        %FILE BECAUSE THEY MAY BE NEEDED FOR GLM
        fprintf(1, '\nWARNING! Condition ''%s'' in file ''%s'' is empty.\n', prt.Condition{c}.name, prt.name);
    end
    for t = 1:prt.Condition{c}.ntpts
        fprintf(fid, '%8d %8d\n', round(prt.Condition{c}.estart(t)), round(prt.Condition{c}.eend(t)));
    end
    fprintf(fid, 'Color:  %3d %3d %3d\n\n', prt.Condition{c}.color);
    
    
end
fclose( fid);
fprintf(1, 'DONE\n');
