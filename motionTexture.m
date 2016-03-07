
function motionTexture(subj_number, run_number)
% function for direction identification of moving texture
% this function only present one moving stimuli in either the blind field
% or the good field
%e.g motionTexture(1,5); %subject 1 and run 5;

%check inputs
if ~exist('subj_number','var')||isempty('subj_number')
    subj_number = 99;
    debug = 1; 
end

if ~exist('run_number','var')||isempty('run_number')
    run_number = 1;
end


clc;close all;

warning('off','MATLAB:dispatcher:InexactMatch');
ListenChar(2);

%% -------------parameters you may want to change ------------------
Params.General.EyeTrack                = 0; %1, use eye tracking; 0, no eye tracking
Params.General.Env_radius              = 3;      % deg, radius of the envelope of stimuli;
Params.General.Stimuli_ecc             = [6 3];  %[horizantol(deg),vertical(deg)];eccentricity of the stimuli
...horizontol should always be a positive value; vertical ecc can be negative (upper visual field) and positive (lower visual field)
Params.General.SFRange                 = [0 1];   %cycles/deg, cut for spatial frequency band for texture stimuli



%% ---------All parameters should be here---------------
% KEY Monitor Parameters
Params.General.View_Distance           = 119;          %cm
Params.General.Resolution              = [1400 1050]; %width, height
Params.General.MonitorWidth            = 42.8;        %cm
Params.General.Linearize 	           = 1;           %whether monitor is linearized. if=1, program will look for "MyGammaTable"
Params.General.Break_if_mismatch       = 1;       % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)
Params.General.Scale_factor            = Params.General.MonitorWidth*60/(Params.General.Resolution(1)*Params.General.View_Distance *tand(1));     % most important parameter - how many acrmin is one screen pixel?


% experiment parameters
Params.General.Experiment              ='TextureMotion';
Params.General.Sub                     = subj_number;
Params.General.Run                     = run_number;  %
Params.General.N_run                   = 1;  % how many runs in this program
Params.General.Dur_stimuli             = 8;  %secs, duration for a motion stimuli
Params.General.Dur_blank               = 8;  %secs, duration for blank after stimuli;
Params.General.Spatial_envelope        = 0;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
Params.General.Trials_condition        = 4;  %how many trials /direction/retinal location in a run
Params.General.N_direction             = 2;  %number of directions, left/right here;
Params.General.N_location              = 2;  %number of retinal locatin we want to test;

% motion texure parameters
Params.General.Background              = [128 128 128];%rgb value
Params.General.Character_color         = [0 0 0];%rgb value


Params.General.Dots_density            = 8;    % number of dots in per deg^2
Params.General.Contrast                = 100;  % 100%
Params.General.Speed                   = 5;    %deg/sec
Params.General.F_kill                  = 0; %fraction of dots to kill each frame (limited lifetime)

%do some simple computation
Params.General.White                   = 255;
Params.General.Amplitude               = 128;
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
    Screen('TextSize',w,30);Screen('TextFont',w,'Arial');
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
    if ~debug && Params.General.EyeTrack
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
        'if ~debug&&Params.General.EyeTrack '...
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
    Screen('DrawText',w,'Press any key to start the experiment',sr_hor-300,sr_ver-80,Params.General.Character_color);Screen('Flip',w);
    KbWait(-1);
    Screen('Flip', w);
    FlushEvents('keyDown');
    
    
    
    % IS THIS DRAW BOX WORKING? start recoding from beginning of the
    % experiment
    if ~debug && Params.General.EyeTrack
        Eyelink('StartRecording');
        WaitSecs(.050);
        Eyelink('Command', 'draw_box %d %d %d %d 15', elSettings.bounds(1),elSettings.bounds(2),...
            elSettings.bounds(3),elSettings.bounds(4));
    end
    
    
    
    %% ------------------------make the movie----------------------------
    % we premake the movie before experiment starts
    % make a wide image patch with dots.
    bigPatch_X = ceil( bps + Params.General.PixStep * Params.General.MvLength_stimuli);
    big_im     = rand(bps,bigPatch_X);
    %big_im(big_im<0.5)=0;
    %big_im(big_im>=0.5)=1;

    %big_im = Params.General.Amplitude + Params.General.Amplitude * big_im;

    % filtered image
    samplingRate =  1*60/Params.General.Scale_factor;
    F = fftshift(fft2(big_im));
    cutBound = (Params.General.SFRange/samplingRate)'*[bps, bigPatch_X];
    [m,n]=meshgrid(-bigPatch_X/2:bigPatch_X/2-1,-stimulus_radius:stimulus_radius);

    if Params.General.SFRange(1) == 0 && Params.General.SFRange(2) ~= 0%low pass
        c = sqrt((m/cutBound(2,2)).^2+(n/cutBound(2,1)).^2)<=1;

    elseif Params.General.SFRange(2) == inf && Params.General.SFRange(1) ~= 0 %high pass
        c=sqrt((m/cutBound(1,2)).^2+(n/cutBound(1,1)).^2)>=1;
    else %band pass
        c1=sqrt((m/cutBound(1,2)).^2+(n/cutBound(1,1)).^2)>=1;
        c2=sqrt((m/cutBound(2,2)).^2+(n/cutBound(2,1)).^2)<=1;
        c=c1.*c2;
    end
    F=F.*c;
    big_im = real(ifft2(ifftshift(F)));
    big_im = (big_im-min(big_im(:)))/(max(big_im(:))-min(abs(big_im(:)))); %normalize to (0,1) range
    big_im = Params.General.Amplitude + Params.General.Amplitude * (2*big_im-1);
    %big_im = im2uint8(big_im);

    bigPatch = Screen('MakeTexture', w, big_im); %make texture of big matrix


    
    
    %% wait for trigger
    disp('Waiting for trigger from scanner or the tilde button (~)');
    KbName('UnifyKeyNames');
    trigger = KbName('`~');
    [~,~,key] = KbCheck(-1);
    while ~key(trigger);
        WaitSecs(.001);
        Screen('DrawText',w,'Waiting for trigger....',sr_hor-300,sr_ver-80,Params.General.Character_color);Screen('Flip',w);
        [~,~,key] = KbCheck(-1);
    end
    startTime=GetSecs;
    
    
    %% Main loop
    for run=1:Params.General.N_run
        
        if exit_index == 1;
            ShowCursor;
            break;
        end
        
        %% eyetracking part
        if ~debug && Params.General.EyeTrack; Eyelink('Message', 'TRIALID %d', trial); end
        if ~debug && Params.General.EyeTrack; Eyelink('Message','SCANNER SYNCTIME'); end
%         while GetSecs - start < wait_beginning % wait for rest time to be up
%          eval(check_eyes) % check fixation while we wait
%         end
        %%
        
        
        %randomize trial sequence in this run
        trial_order = randperm(trials_run);
        temp_data   = [direction_list' location_list'];
        temp_data   = temp_data(trial_order,:);%[direction, location] in each trial
        
        
        Screen('FillRect',w, Params.General.Background);
        Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
        Screen('DrawText',w,'Keep your fixation on the center spot...',sr_hor-300,sr_ver-120,Params.General.Character_color);
        vbl=Screen('Flip', w); %end of the blank before this run

        
       
        Screen('FillRect',w, Params.General.Background);
        Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
        vbl=Screen('Flip', w,vbl+Params.General.Dur_blank - 1);% 
        prt.data = zeros(trials_run,5); %direction, location, condition, start time, end time
        trialEndTime = 0;  %set a default value 
        for trial=1:trials_run
            t1=GetSecs;
            
            %determine direction and retinal location of the stimulus in
            %this trial
            if temp_data(trial,1)==1 %left
                rotation_angle = 0;
            elseif temp_data(trial,1)==2 %right
                rotation_angle = 180;
            end
            prt.data(trial,1)=temp_data(trial,1); %save direction
            
            switch temp_data(trial,2)
                case 1 % left
                   screen_patch=screen_rect_middle+[-Params.General.H_ecc_stim Params.General.V_ecc_stim -Params.General.H_ecc_stim Params.General.V_ecc_stim];
                case 2 % right               
                   screen_patch=screen_rect_middle+[Params.General.H_ecc_stim Params.General.V_ecc_stim Params.General.H_ecc_stim Params.General.V_ecc_stim];
            end
            prt.data(trial,2)= temp_data(trial,2);%save location
            prt.data(trial,3)= (temp_data(trial,1)-1)*Params.General.N_location+temp_data(trial,2);%save condition (1-4);
            
            % play the animation of foveal
            Screen('FillRect',w, Params.General.Background);
            Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            Screen('Flip', w);
            mm = 31;
            for i=0:4
                nn = mm-i*4;
                Screen('FrameOval', w,Params.General.Character_color,[sr_hor-nn, sr_ver-nn, sr_hor+nn, sr_ver+nn],2,2)
                Screen('Flip', w);
                WaitSecs(0.05);
            end
            Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            Screen('Flip', w);
            WaitSecs(0.2);
            
            % play the movie
            priorityLevel=MaxPriority(w);Priority(priorityLevel);
            %testblank = GetSecs-trialEndTime
            trialStartTime = GetSecs;
            prt.data(trial,4) = (GetSecs-startTime)*1000;%starting time;
            for frame = 1:Params.General.MvLength_stimuli
                
                offset=(frame-1)*Params.General.PixStep;  %how many pixel moved in this frame
                %because of left motion, source rect should move to from left
                %to right
                %offset
                srcrect=[offset 0 bps+offset bps];
                
                
                %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture',w,bigPatch,srcrect,screen_patch,rotation_angle);
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                Screen('DrawTexture',w,maskCOS,movie_rect,screen_patch,[],[]);
                
                Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
                if frame == 1
                    vbl=Screen('Flip',w,startTime+(trial-1)*16+8);%control the timing for first frame
                    fprintf('Stimulus onset time is %.10f ms\n',(vbl-startTime)*1000);
                else
                    Screen('Flip',w);
                end
                % run eyetracker
                if ~debug && Params.General.EyeTrack
                   eval(check_eyes);
                end
            end
                
             
            %test timing
            trialEndTime= GetSecs;
            
            prt.data(trial,4) = (vbl-startTime)*1000;%starting time 
            prt.data(trial,5)=(trialEndTime-startTime)*1000;%end time;
            
            %t =  trialEndTime - trialStartTime
            if trial>1
                trial
                t = (prt.data(trial,4) -  prt.data(trial-1,5))/1000;
                fprintf('the previous blank is %.10f second\n',t);
            end
            
            
            Screen('FillRect',w, Params.General.Background);
            Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            vbl=Screen('Flip',w);
            %fprintf('blank onset time is %.10f ms',(vbl-startTime)*1000);
            
            Priority(0);
            FlushEvents('keyDown');
            
            Screen('FillRect',w, Params.General.Background);
            Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            postDelay=GetSecs-trialEndTime;
            if trial == trials_run
                vbl=Screen('Flip', w, vbl+Params.General.Dur_blank);% we don't need to compensate delay in the last trial
            else
                vbl=Screen('Flip', w, vbl+Params.General.Dur_blank - 1); %we need to compensate delayed time at the begnining of the trial
            end
                        
            FlushEvents('keyDown');
            % Close movies
            %Screen(movie_play, 'Close');
            clear frame;
            
        end
    end
    ShowCursor;
    TotalTime = GetSecs-startTime
    
    if ~debug && Params.General.EyeTrack
        Eyelink('Message','!V IAREA RECTANGLE %d %d %d %d %d %s',1,elSettings.bounds(1),...
            elSettings.bounds(2),elSettings.bounds(3),elSettings.bounds(4),'fixation');
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        Eyelink('ReceiveFile');

    end
    
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    sca;
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
    
    
    
    if savefile
        save(filename);
        %compute prt for four conditions
        for i =1:prt.NrOfConditions
            prt.Condition{i}.ntpts = 4;% 4 trials/condition
            prt.Condition{i}.estart = prt.data(prt.data(:,3)==i,4);
            prt.Condition{i}.eend = prt.data(prt.data(:,3)==i,5);
            prt.Condition{i}.color=prt.colors(i,:);
        end
        %add gaze condition
        if ~debug && Params.General.EyeTrack
            prt.NrOfConditions=prt.NrOfConditions+1;
            prt.Condition{prt.NrOfConditions}.name = 'gaze';
            prt.Condition{prt.NrOfConditions}.ntpts = size(gaze_failed,1);
            if prt.Condition{prt.NrOfConditions}.ntpts==0
                prt.Condition{prt.NrOfConditions}.estart = 0;
                prt.Condition{prt.NrOfConditions}.eend = 0;
            else
                prt.Condition{prt.NrOfConditions}.estart = gaze_failed(:,1);
                prt.Condition{prt.NrOfConditions}.eend = gaze_failed(:,2);
            end
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
        if ~debug && Params.General.EyeTrack
            %add gaze condition
            prtSingleTrial.NrOfConditions=prtSingleTrial.NrOfConditions+1;
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.name = 'gaze';
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.ntpts = size(gaze_failed,1);
            if prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.ntpts == 0
                  prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.estart = 0;
                  prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.eend = 0;
            else
                  prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.estart = gaze_failed(:,1);
                  prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.eend = gaze_failed(:,2); 
            end
            prtSingleTrial.Condition{prtSingleTrial.NrOfConditions}.color=prtSingleTrial.colors(end,:);
        end
        bv_prt_write (prtSingleTrial);
        save(filename);
    end

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



