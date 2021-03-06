
%This is the program of the experiment in which we are trying to decode
%motion direction in cortical bind patients.
%e.g motionDots(1,5); %subject 1 and run 5;
%By Ruyuan Zhang last updated 08/28/2015

%%


%% left-right dirction

function motionDots(subj_number, run_number)
%e.g motionDots(1,5); %subject 1 and run 5;

%check inputs
if nargin==0 
    subj_number=0;
    run_number=0;
    DEBUG = 1;
elseif nargin==1 
    run_number=0;
    DEBUG = 1;
elseif nargin==2 
    DEBUG = 0;
else
    error('Inappropriate number of input arguments')
end;

clc;close all;

warning('off','MATLAB:dispatcher:InexactMatch');
ListenChar(2);
 
%% ---------All parameters should be here---------------
% KEY Monitor Parameters
Params.General.View_Distance           = 119;          %cm
Params.General.Resolution              = [1400 1050]; %
Params.General.MonitorWidth            = 42.8;        %cm
Params.General.Linearize 	           = 1;           %whether monitor is linearized. if=1, program will look for "MyGammaTable"
Params.General.Break_if_mismatch       = 1;       % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)
Params.General.Scale_factor            = Params.General.MonitorWidth*60/(Params.General.Resolution(1)*Params.General.View_Distance *tand(1));     % most important parameter - how many acrmin is one screen pixel?

% experiment parameters
Params.General.Experiment              ='MotionDots';
Params.General.Sub                     = subj_number;
Params.General.Run                     = run_number;  %
Params.General.N_run                   = 1;  % how many runs in this program
Params.General.Dur_stimuli             = 8; %secs, duration for a motion stimuli
Params.General.Dur_blank               = 8; %secs, duration for blank after stimuli;
Params.General.Spatial_envelope        = 0;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
Params.General.Trials_condition        = 4; %how many trials /direction/retinal location in a run
Params.General.N_direction             = 2; %number of directions, left/right here;
Params.General.N_location              = 2; %number of retinal locatin we want to test;
Params.General.EyeTrack                = 0; %1, use eye tracking; 0, no eye tracking

% motion dots parameters
Params.General.Background              = [128 128 128];%rgb value
Params.General.Dot_color               = [0 0 0];%rgb value
Params.General.Dot_background          = [128 128 128];%rgb value
Params.General.Dot_size                = 0.4;    % deg, radius of each dot
Params.General.Env_radius              = 3;      % eg, radius of whole stimuli;
Params.General.Stimuli_ecc             = [6 3];  %[horizantol(deg),vertical(deg)];eccentricity of the stimuli

Params.General.Dots_density            = 1.2;    % number of dots in per deg^2
Params.General.Contrast                = 100;    % 100%
Params.General.Speed                   = 5;      %deg/sec
Params.General.F_kill                  = 0; %fraction of dots to kill each frame (limited lifetime)

%do some simple computation
Params.General.White                   = 255;
Params.General.Dot_size                = round(Params.General.Dot_size*2*60/Params.General.Scale_factor); %dot diameter in pixel
Params.General.N_total_dots            = round(Params.General.Dots_density*pi*Params.General.Env_radius^2);  % number of dots
Params.General.H_ecc_stim              = round(Params.General.Stimuli_ecc(1)*60/Params.General.Scale_factor);
Params.General.V_ecc_stim              = round(Params.General.Stimuli_ecc(2)*60/Params.General.Scale_factor);


%% ---------------------------------------
savefile              =1;
if subj_number == 0,savefile = 0 ; end;
tme                   =clock;
filename = strcat([Params.General.Experiment '_Sub' int2str(Params.General.Sub) '_Run' int2str(Params.General.Run)],'.mat');

%% -----info on prt file;
%prt for four conditions
prt.name = [Params.General.Experiment '_Sub' int2str(Params.General.Sub) '_Run' int2str(Params.General.Run) '.prt'];
prt.NrOfConditions = Params.General.N_direction * Params.General.N_location;
prt.Condition{:}.ntpts=Params.General.Trials_condition;
prt.Condition{1}.name='DirLeft_FieldLeft'; %direction_field;
prt.Condition{2}.name='DirLeft_FieldRight'; %direction_field;
prt.Condition{3}.name='DirRight_FieldLeft'; %direction_field;
prt.Condition{4}.name='DirRight_FieldRight'; %direction_field;
% generate the colors (jet is a 64x3 color table blue-green-yellow-red):
colors = 255*jet; close(gcf);
prt.colors = round(linspace(1,length(colors),prt.NrOfConditions+1));%%one more color for eyeTracking data
prt.colors = round(colors(prt.colors,:));


%prt for individual trials
prtSingleTrial.name = [Params.General.Experiment '_Sub' int2str(Params.General.Sub) '_Run' int2str(Params.General.Run) '_SingleTrial.prt'];
prtSingleTrial.NrOfConditions = Params.General.N_direction * Params.General.N_location * Params.General.Trials_condition;
prtSingleTrial.Condition{:}.ntpts=1; % for signal trial
prtSingleTrial.Condition{1}.name='DirLeft_FieldLeft_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{2}.name='DirLeft_FieldLeft_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{3}.name='DirLeft_FieldLeft_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{4}.name='DirLeft_FieldLeft_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{5}.name='DirLeft_FieldRight_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{6}.name='DirLeft_FieldRight_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{7}.name='DirLeft_FieldRight_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{8}.name='DirLeft_FieldRight_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{9}.name='DirRight_FieldLeft_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{10}.name='DirRight_FieldLeft_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{11}.name='DirRight_FieldLeft_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{12}.name='DirRight_FieldLeft_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{13}.name='DirRight_FieldRight_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{14}.name='DirRight_FieldRight_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{15}.name='DirRight_FieldRight_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{16}.name='DirRight_FieldRight_trial4'; %direction_field_trial4;

% generate the colors (jet is a 64x3 color table blue-green-yellow-red):
colors = 255*jet; close(gcf);
prtSingleTrial.colors = round(linspace(1,length(colors),prtSingleTrial.NrOfConditions+1)); %one more color for eyeTracking data
prtSingleTrial.colors = round(colors(prtSingleTrial.colors,:));



%%
%create data structure template for each run and we use this template
%to do randomize trial order in each run
trials_run       = Params.General.Trials_condition * Params.General.N_direction * Params.General.N_location;%how many trials/run.
direction_list   = sort(rem(1:trials_run,Params.General.N_direction)+1);
location_list    = rem(1:trials_run,Params.General.N_location)+1;


try
    %open Screen windows,
    Screen('Preference', 'SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens); %present in max screeen numbers
    [w,screen_rect]=Screen('OpenWindow',screenNumber,0,[],[],2);
    Params.General.IFI=Screen('GetFlipInterval', w);
    Params.General.Frame_rate=1/Params.General.IFI;	% frames per second

    
    %Now I know frame_rate, do some simple computations
    Params.General.MvLength_stimuli       = round(Params.General.Frame_rate*Params.General.Dur_stimuli);
    Params.General.MvLength_blank         = round(Params.General.Frame_rate*Params.General.Dur_blank);
    Params.General.PixStep                = Params.General.Speed*60/(Params.General.Frame_rate * Params.General.Scale_factor);        %moving speed, how many pixel/frame
    %housekeeping stuff
    stimulus_radius  = round(60* Params.General.Env_radius/Params.General.Scale_factor);%how many pixels for stimuli radius
    
    if Params.General.Linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, Params.General.Background);Screen('Flip', w);
    Screen('TextSize',w,30);Screen('TextFont',w,'Charcoal');
    sr_hor = round(screen_rect(3)/2); sr_ver = round(screen_rect(4)/2);
    
    
    
    % make the spatial envelope
    [x,y]=meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle=((stimulus_radius)^2-(x.^2+y.^2));
    for i=1:bps;
        for j =1:bps;
            if circle(i,j) < 0; circle(i,j) = 0;
            else
                circle(i,j) = 1;
            end;
        end;
    end;
    if Params.General.Spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/2)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif Params.General.Spatial_envelope == 2 %raised cosine evelope
        R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);
    end
    %make a square mask to exclude some dots pixels outside the square
    %aperture.
    square_mask = zeros(bps+2*Params.General.Dot_size, bps+2*Params.General.Dot_size,4);
    for color=1:3 %rgb
        square_mask(:,:,color)=ones(bps+2*Params.General.Dot_size,bps+2*Params.General.Dot_size)*Params.General.Background(color);
    end
    square_mask(Params.General.Dot_size+1:end-Params.General.Dot_size,Params.General.Dot_size+1:end-Params.General.Dot_size,4)=1;
    square_mask(:,:,4)=(1-square_mask(:,:,4))*255;
    maskSquare=Screen('MakeTexture',w, square_mask);
    
    
    %create cosine mask to blend image below
    circleCOS = zeros(bps, bps,4);
    for color=1:3 %rgb
        circleCOS(:,:,color)=ones(bps,bps)*Params.General.Background(color);
    end
    circleCOS(:,:,4)=(1-circle)*255;
    maskCOS=Screen('MakeTexture',w, circleCOS);
    
    %  create stimulus rectangles
    movie_rect = [0 0 bps bps];
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    screen_rect_middle = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
    %screen_patch = screen_rect_middle + [Params.General.H_ecc_stim,V_ecc_stim,Params.General.H_ecc_stim,V_ecc_stim];
    
    exit_index = 0;
    
    
    %% 1st set up eyetracker
    % eyetracking
    if ~DEBUG && Params.General.EyeTrack
        HideCursor;
        elSettings = SetupEyeLink(Params,w,screen_rect);
        EyelinkDoDriftCorrect(elSettings.el);
        gaze_failed = []; %we initial a empty matrix to have gaze failed data
    elseif length(Screen('Screens'))==1 % can't see stim being created if only one screen
        %sca;
    end
    
    % some code that can be evaluated quickly in the main script
    get_time = @(startTime,since) (since-startTime)*1000;
    check_eyes = [...
        'if ~DEBUG&&Params.General.EyeTrack '...
        '    elSettings = confirmFixation(elSettings); '... % checks if eye is in bounds
        '    if ~gaze_ok && elSettings.endstatus(end) == 1 '... % if eye was out of bounds but is back in
        '        gaze_failed(end,2) = get_time(startTime,GetSecs); '... % end the prt condition timer
        '        gaze_ok = 1; '... % acknowledge that the eye is back in bounds
        '    elseif gaze_ok && elSettings.endstatus(end) == 3 '... % gaze was in bounds but is out now
        '        gaze_failed(end+1,1) = get_time(startTime,GetSecs); '... % start the prt condition timer
        '        gaze_ok = 0; '... % acknowledge that the eye is out of bounds
        '    end; '...
        'end'
        ];
    gaze_ok = 1; %index whether gaze is out of bound or not;
    
    
    
    %% MAIN LOOP
    HideCursor;
    tic; skip_first = 1;trial = 1;
    
    Screen('FillRect',w, Params.General.Background);
    Screen('DrawText',w,'Press any key to start the experiment',sr_hor-300,sr_ver-80,Params.General.Dot_color);Screen('Flip',w);
    KbWait(-1);
    Screen('Flip', w);
    FlushEvents('keyDown');
    
    
    
    % IS THIS DRAW BOX WORKING? start recoding from beginning of the
    % experiment
    if ~DEBUG && Params.General.EyeTrack
        Eyelink('StartRecording');
        WaitSecs(.050);
        Eyelink('Command', 'draw_box %d %d %d %d 15', elSettings.bounds(1),elSettings.bounds(2),...
            elSettings.bounds(3),elSettings.bounds(4));
    end
    
    
    %% wait for trigger
    disp('Waiting for trigger from scanner or the tilde button (~)');
    KbName('UnifyKeyNames');
    trigger = KbName('`~');
    [~,~,key] = KbCheck(-1);
    while ~key(trigger);
        WaitSecs(.001);
        Screen('DrawText',w,'Waiting for trigger....',sr_hor-300,sr_ver-80,Params.General.Dot_color);Screen('Flip',w);
        [~,~,key] = KbCheck(-1);
    end
    startTime=GetSecs;
    
    
    %% Main loop
    for run=1:Params.General.N_run
        
        if exit_index == 1;
            ShowCursor
            break;
        end
        
        %% eyetracking part
        if ~DEBUG && Params.General.EyeTrack; Eyelink('Message', 'TRIALID %d', trial); end
        if ~DEBUG && Params.General.EyeTrack; Eyelink('Message','SCANNER SYNCTIME'); end
%         while GetSecs - start < wait_beginning % wait for rest time to be up
%          eval(check_eyes) % check fixation while we wait
%         end
        %%
        
        
        %randomize trial sequence in this run
        trial_order = randperm(trials_run);
        temp_data   = [direction_list' location_list'];
        temp_data   = temp_data(trial_order,:);%[direction, location] in each trial
        
        
        Screen('FillRect',w, Params.General.Background);
        Screen('FillOval', w,Params.General.Dot_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
        Screen('DrawText',w,'Keep your fixation on the center spot...',sr_hor-300,sr_ver-120,Params.General.Dot_color);
        vbl=Screen('Flip', w); %end of the blank before this run

        
        % check whether eyetrack is still running
        %eyetrackCounter = 1;
        
        Screen('FillRect',w, Params.General.Background);
        vbl=Screen('Flip', w,vbl+Params.General.Dur_blank);
        %preblanck = vbl-startTime
        prt.data = zeros(trials_run,5); %direction, location, condition, start time, end time
        trialEndTime = 0;  %set a default value 
        for trial=1:trials_run
            t1=GetSecs;
            
            %determine direction and retinal location of the stimulus in
            %this trial
            if temp_data(trial,1)==1 %left
                direction = -1;
            elseif temp_data(trial,1)==2 %right
                direction = 1;
            end
            prt.data(trial,1)=temp_data(trial,1); %save direction
            
            switch temp_data(trial,2)
                case 1 % left
                    screen_patch=screen_rect_middle+[-Params.General.H_ecc_stim Params.General.V_ecc_stim -Params.General.H_ecc_stim Params.General.V_ecc_stim];
                case 2 % right
                    screen_patch=screen_rect_middle+[Params.General.H_ecc_stim Params.General.V_ecc_stim Params.General.H_ecc_stim Params.General.V_ecc_stim];
            end
            prt.data(trial,2)= temp_data(trial,2);%save location
            prt.data(trial,3)= (temp_data(trial,1)-1)*Params.General.N_location+temp_data(trial,2);%save condition (1-8);
            
            %------------------------make the movie----------------------------
            xy = diag([bps, bps]) * rand( 2, round( Params.General.N_total_dots) );
            dxdy = repmat([direction*Params.General.PixStep;0],1,Params.General.N_total_dots);                 % change in x and y per frame (pixels),let's fix to move to right here.
            
            
            %set color of dots
            luminance = Params.General.Dot_color'*ones( 1, size(xy,2) );
            color_arg = [luminance; Params.General.Contrast/100*255*ones( 1, size(xy,2))];
            
            
            %make the first image frame
            Screen('FillRect',w, Params.General.Dot_background);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);%blend function to ensure making smooth dots
            Screen('DrawDots', w, xy , Params.General.Dot_size, color_arg, [], 1);
            
            
            screen_patch_dots=repmat([screen_patch(1);screen_patch(2)],1,Params.General.N_total_dots);
            
            clear frame;
            
            
            
            %----------play the movie----------------
            %  initiate trial
            Screen('FillRect',w, Params.General.Background);
            Screen('Flip', w);
            mm = 31;
            for i=0:4
                nn = mm-i*4;
                Screen('FrameOval', w,Params.General.Dot_color,[sr_hor-nn, sr_ver-nn, sr_hor+nn, sr_ver+nn],2,2)
                Screen('Flip', w);
                WaitSecs(0.05);
            end
            Screen('FillOval', w,Params.General.Dot_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            Screen('Flip', w);
            WaitSecs(0.2);
           
            % play the movie
            
            priorityLevel=MaxPriority(w);Priority(priorityLevel);
            %testblank = GetSecs-trialEndTime
            trialStartTime = GetSecs;
            prt.data(trial,4) = (GetSecs-startTime)*1000;%starting time;
            for frame = 1:Params.General.MvLength_stimuli
                
                xy=xy+dxdy;
                % check to see which dots have gone beyond the borders of the annuli
                dot_out = find(xy(1,:) > bps |xy(1,:) <0 | rand(1,Params.General.N_total_dots) < Params.General.F_kill);	% dots to reposition
                nout = length(dot_out);
                if nout
                    % choose new coordinates
                    xy(1,dot_out) = xy(1,dot_out)-direction*bps; %horizontal location for new dots,we just consider horizontal for left and right
                    xy(2,dot_out) = bps * rand( 1, round(nout) ); %vertical location for new dots
                    
                end;
                xy_present = xy+screen_patch_dots;
                Screen('FillRect',w, Params.General.Dot_background,screen_patch);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);%blend function to ensure making smooth dots
                Screen('DrawDots', w, xy_present, Params.General.Dot_size, color_arg, [], 1);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture',w,maskSquare,movie_rect+[0 0 2*Params.General.Dot_size 2*Params.General.Dot_size],screen_patch+[-Params.General.Dot_size -Params.General.Dot_size Params.General.Dot_size Params.General.Dot_size],[],[]);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[]);
                Screen('FillOval', w,Params.General.Dot_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
                vbl = Screen('Flip',w);
                
                % run eyetracker
                if ~DEBUG && Params.General.EyeTrack
                   eval(check_eyes);
                end
            end
               
             
            %test timing
            trialEndTime= GetSecs;
            preDelay = trialStartTime-t1;
            
            
            prt.data(trial,5)=(trialEndTime-startTime)*1000;%end time;
            
            t =  trialEndTime - trialStartTime
            %t = ( trialEndTime - trialStartTime) / Params.General.MvLength_stimuli
            
            
            Screen('FillRect',w, Params.General.Background);
            Screen('FillOval', w,Params.General.Dot_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            Screen('Flip',w);
            
            Priority(0);
            FlushEvents('keyDown');
            
            Screen('FillRect',w, Params.General.Background);
            Screen('FillOval', w,Params.General.Dot_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            postDelay=GetSecs-trialEndTime;
            vbl=Screen('Flip', w, trialEndTime+Params.General.Dur_blank-preDelay-postDelay);
                        
            FlushEvents('keyDown');
            % Close movies
            %Screen(movie_play, 'Close');
            clear frame;
            
        end
    end
    ShowCursor;
    if ~DEBUG && Params.General.EyeTrack
        Eyelink('Message','!V IAREA RECTANGLE %d %d %d %d %d %s',1,elSettings.bounds(1),...
            elSettings.bounds(2),elSettings.bounds(3),elSettings.bounds(4),'fixation');
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        Eyelink('ReceiveFile');
        % this was edited by alex 4/8/15)
        
%         %movefile(edf,sub_folder,'f')
%         mkdir('\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\ouput')
%         if exist(['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\ouput\' num2str(Params.General.Sub) '_R1.edf'],'file') > 1;
%             copyfile(['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\ouput\' num2str(Params.General.Sub) '_R1.edf'], ['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\output' Params.General.Experiment num2str(Params.General.Sub) '_R1.edf']);
%         end
%         if exist(['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\ouput\' num2str(Params.General.Sub) '_R2.edf'],'file') > 1;
%             copyfile(['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\ouput\' num2str(Params.General.Sub) '_R2.edf'], ['\\files\cantlonData\NEUROPSYCHLAB\Brad_ExpData\MRI\Retinotopy\output' Params.General.Experiment num2str(Params.General.Sub) '_R2.edf']);
%         end
    end
    
    
    
    
    if savefile~=0
        %compute prt for four conditions
        for i =1:prt.NrOfConditions
            prt.Condition{i}.ntpts = 4;% 4 trials/condition
            prt.Condition{i}.estart = prt.data(prt.data(:,3)==i,4);
            prt.Condition{i}.eend = prt.data(prt.data(:,3)==i,5);
            prt.Condition{i}.color=prt.colors(i,:);
        end
        %add gaze condition
        if ~DEBUG && Params.General.EyeTrack
            prt.NrOfConditions=prt.NrOfConditions+1;
            prt.Condition{prt.NrOfConditions}.name = 'gaze';
            prt.Condition{prt.NrOfConditions}.ntpts = size(gaze_failed,1);
            prt.Condition{prt.NrOfConditions}.estart = gaze_failed(:,1);
            prt.Condition{prt.NrOfConditions}.eend = gaze_failed(:,2);
            prt.Condition{prt.NrOfConditions}.color=prt.colors(end,:);
        end
        bv_prt_write (prt);
        
        %compute prt for all individual trials
        for i = 1: prtSingleTrial.NrOfConditions
            prtSingleTrial.Condition{i}.ntpts=1; % for signal trial
            prtSingleTrial.Condition{i}.color=prtSingleTrial.colors(i,:);
            condition_temp=fix((i-1)/4)+1;
            rem_temp = rem(i-1,4)+1;
            estart_temp = prt.data(prt.data(:,3)==condition_temp,4);
            eend_temp   = prt.data(prt.data(:,3)==condition_temp,5);
            prtSingleTrial.Condition{i}.estart =estart_temp(rem_temp);
            prtSingleTrial.Condition{i}.eend = eend_temp(rem_temp);
            
        end
        prtSingleTrial.data=prt.data;
        if ~DEBUG && Params.General.EyeTrack
            %add gaze condition
            prtSingleTrial.NrOfConditions=prtSingleTrial.NrOfConditions+1;
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.name = 'gaze';
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.ntpts = size(gaze_failed,1);
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.estart = gaze_failed(:,1);
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.eend = gaze_failed(:,2);
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.color=prtSingleTrial.colors(end,:);
        end
        bv_prt_write (prtSingleTrial);
        save(filename);
    end
        
      
    
    
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    sca;
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    TotalTime = GetSecs-startTime
catch
    ListenChar(1);
    ddd = lasterror;
    ddd.message
    ddd.stack(1,1).line
    psychrethrow(lasterror);
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    sca;
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch..



