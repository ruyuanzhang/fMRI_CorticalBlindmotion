%% program for decoding plaid motion on cortical blind patients
%






%%
%This is the program of the experiment in which we are trying to decode
%motion direction in cortical bind patients. 

%By Ruyuan Zhang 07/08/2015

%%


%% left-right dirction

clc;clear all;close all;

warning('off','MATLAB:dispatcher:InexactMatch');
ListenChar(2);

%% ---------All parameters should be here---------------
% KEY Monitor Parameters
p.view_Distance          = 58;
p.scale_factor           = 2;     % most important parameter - how many acrmin is one screen pixel? for SONY monitor at [1024 640] resolution, and 30.4in vieving distance, use scale_factor=2
p.frame_rate             = 120;        %
p.resolution             = [1280 720]; %
p.linearize 	         = 1;         % whether monitor is linearized. if=1, program will look for "MyGammaTable"
p.break_if_mismatch      = 1;       % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)

% experiment parameters
p.subject_initials       ='test';
p.n_run                  = 10;  % number of total runs
p.dur_stimuli            = 8; %secs, duration for a motion stimuli
p.dur_blank              = 1; %secs, duration for blank after stimuli;
p.spatial_envelope       = 2;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
p.trials_condition       = 1; %how many trials /direction/retinal location in a run
p.n_direction            = 2; %number of directions, left/right here;
p.n_location             = 4; %number of retinal locatin we want to test;

% motion dots parameters
p.background             = 126;%rgb value
p.env_radius             = 5;      % eg, radius of whole stimuli; 
p.stimuli_ecc            = [7 7];  %[horizantol(deg),vertical(deg)];eccentricity of the stimuli

p.contrast               = 100;    % 100% for plaid so 50% for each component
p.speed                  = 5;      %deg/sec
p.SF                     = 0.5 ;     %cycle/deg,spatial frequency

%do some simple computation
p.white                  = 255;
p.TF                     = p.SF *p.speed;
p.mvLength_stimuli       = round(p.frame_rate*p.dur_stimuli);
p.mvLength_blank         = round(p.frame_rate*p.dur_blank);
p.H_ecc_stim             = round(p.stimuli_ecc(1)*60/p.scale_factor);
p.V_ecc_stim             = round(p.stimuli_ecc(2)*60/p.scale_factor);            
p.orientation            =0;%we can set orientation to vertical but rotate when present it.


%% ---------------------------------------
savefile                 =1;
tme                      =clock;
% filename = strcat(subject_initials,'_block',num2str(block),'_',num2str(speed),'.mat');
% IsExist = exist(filename,'file');
% if IsExist && block ~= 0
%     ListenChar(1);
%     error('data file name exists')
% end
% if block==0
%     savefile=0;
% end



%%
%create data structure template for each run and we we use this template
%to do randomize trial order in each run
trials_run       = p.trials_condition * p.n_direction * p.n_location;%how many trials/run.
direction_list   = sort(rem(1:trials_run,p.n_direction)+1);
location_list    = rem(1:trials_run,p.n_location)+1;


%%

try
    %open Screen windows,
    Screen('Preference', 'SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens); %present in max screeen numbers
    w=Screen('OpenWindow',screenNumber,0,[],[],2);
    w
    screen_rect = Screen('Rect',w);
    
    if p.linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, p.background);Screen('Flip', w);
    Screen('TextSize',w,30);Screen('TextFont',w,'Charcoal');
    sr_hor = round(screen_rect(3)/2); sr_ver = round(screen_rect(4)/2);
    
    
    % MAIN LOOP
    %HideCursor;
    tic; skip_first = 1;trial = 1;
    %Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    %Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
    w
    Screen('DrawText',w,['Press any key to start the experiment'],sr_hor-300,sr_ver-80,0);Screen('Flip',w);
    KbWait(-1);
    
    GetChar;
    
    %Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    %Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
    Screen('Flip', w);
    FlushEvents('keyDown');
    
    %housekeeping stuff
    stimulus_radius  = round(60* p.env_radius/p.scale_factor);%how many pixels for stimuli radius
    %housekeeping stuff for drift grating    
    TFstep = (2*pi*p.TF)/p.frame_rate;
    f=(p.SF*p.scale_factor/60)*2*pi;
    p.orientation=p.orientation*pi/180;
    a=cos(p.orientation)*f; b=sin(p.orientation)*f;
    q1 = 1; q2 = 1;
    amplitude = p.background*p.contrast/100;

    
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
    if p.spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/2)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif p.spatial_envelope == 2 %raised cosine evelope
        R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);
    end

    
    %  create stimulus rectangles
    movie_rect = [0 0 bps bps];
    
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    screen_rect_middle = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
    %screen_patch = screen_rect_middle + [p.H_ecc_stim,V_ecc_stim,p.H_ecc_stim,V_ecc_stim];
    
    %% wait for trigger
    disp('Waiting for trigger from scanner or the tilde button (~)');
    KbName('UnifyKeyNames');
    trigger = KbName('`~');
    [~,~,key] = KbCheck(-1);
    while ~key(trigger);
        WaitSecs(.001);
        Screen('DrawText',w,'Waiting for trigger....',sr_hor-300,sr_ver-80,0);Screen('Flip',w);
        [~,~,key] = KbCheck(-1);
    end
    startTime=GetSecs;
    
    
    exit_index = 0;
    HideCursor;
    
    for run=1:p.n_run
        
        if exit_index == 1;
            ShowCursor;
            break;
        end
        
        %randomize trial sequence in this run
        trial_order = randperm(trials_run);
        temp_data   = [direction_list' location_list'];
        temp_data   = temp_data(trial_order,:);%[direction, location] in each trial
      
        
        for trial=1:trials_run
            
            %determine direction and retinal location of the stimulus in
            %this trial
            if temp_data(trial,1)==1 %left
                rotation_angle = 135;
            elseif temp_data(trial,1)==2 %right
                rotation_angle = 315;
            end
            
            switch temp_data(trial,2)
                case 1 % upright
                    screen_patch=screen_rect_middle+[p.H_ecc_stim -p.V_ecc_stim p.H_ecc_stim -p.V_ecc_stim];
                case 2 % downright
                    screen_patch=screen_rect_middle+[p.H_ecc_stim p.V_ecc_stim p.H_ecc_stim p.V_ecc_stim];
                case 3 % downleft
                    screen_patch=screen_rect_middle+[-p.H_ecc_stim p.V_ecc_stim -p.H_ecc_stim p.V_ecc_stim];
                case 4 % upleft
                    screen_patch=screen_rect_middle+[-p.H_ecc_stim -p.V_ecc_stim -p.H_ecc_stim -p.V_ecc_stim];
            end
            
            
            %------------------------make the movie----------------------------
            motion_step(1) = rand*2*pi;
            for i=2:p.mvLength_stimuli;
                motion_step(i) = motion_step(i-1)+TFstep;
            end
            
            for i = 1:p.mvLength_stimuli;
                grating1=round(((sin(a*x+b*y+motion_step(i)).*circle*amplitude)+p.background));
                grating2=grating1';
                movie_temp=(grating1+grating2)/2;
                movie_play{i} = Screen('MakeTexture',w, movie_temp);
            end
            
            
            %----------play the movie----------------
            %  initiate trial
            t1 = GetSecs;
            Screen('FillRect',w, p.background);
            Screen('Flip', w);
            mm = 19;
            for i=0:4
                nn = mm-i*4;
                Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn, sr_hor+nn, sr_ver+nn],2,2)
                Screen('Flip', w);
                WaitSecs(0.05);
            end
            Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn, sr_hor+nn, sr_ver+nn],2,2)
            Screen('Flip', w);
            WaitSecs(0.36);
            
            Screen('FillRect',w, p.background);
            Screen('Flip', w);
            WaitSecs(0.3);
            
            % play the movie
            priorityLevel=MaxPriority(w);Priority(priorityLevel);
            blah=GetSecs;
                      
           
            for frame = 1:p.mvLength_stimuli
                Screen('DrawTexture', w, movie_play{frame}, movie_rect, screen_patch, rotation_angle);
                Screen('Flip',w);            
            end
            %test timing
            t=GetSecs-blah
            
            Screen('FillRect',w, p.background);
            Screen('Flip',w);
            Priority(0);
            FlushEvents('keyDown');
         
            Screen('FillRect',w, p.background);
            vbl=Screen('Flip', w);
            
            
            
            
                 
            FlushEvents('keyDown');
            
            % Close movies
            for frame = 1:p.mvLength_stimuli
                Screen('Close',  movie_play{frame});         
            end
            clear grating1 grating2 movie_temp;
            
        end
         
    end
    ShowCursor;
    if savefile~=0;
        save(filename);
    end
    
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
catch
    ListenChar(1);
    ddd = lasterror;
    ddd.message
    ddd.stack(1,1).line
    psychrethrow(lasterror);
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    Screen('CloseAll');
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch..
time = toc/60;


%%%%finally, do some computation of staircase,it should output thresholds
%%%%and initial values for next training blocks

