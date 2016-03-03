function motion_localizer_texture
% motion localizer using low spatial frequency
% Left-static-right-static cycles for MRI experiment, MT detection, 
% no input and output

% history
% 02/28/2016 Ruyuan Zhang created this function


clc;close all;

warning('off','MATLAB:dispatcher:InexactMatch');
%ListenChar(2);

%% -------------parameters you may want to change ------------------
Params.General.Env_radius              = 6;      % deg, radius of the envelope of stimuli;
Params.General.Stimuli_ecc             = [0 0];  %[horizantol(deg),vertical(deg)];eccentricity of the stimuli
...horizontol should always be a positive value; vertical ecc can be negative (upper visual field) and positive (lower visual field)
Params.General.SFRange                 = [0 1];   %cycles/deg, cut for spatial frequency band for texture stimuli


%% ---------All parameters should be here---------------
% General experiment parameters
Params.General.Experiment               = 'Motion_Localizer_Texture'; % name of the experiment.


% KEY Monitor Parameters
Params.General.Resolution              = [1400 1050]; %width, height
Params.General.MonitorWidth            = 42.8;        %cm
Params.General.Linearize 	           = 1;           %whether monitor is linearized. if=1, program will look for "MyGammaTable"
Params.General.Scale_factor             = 2 ; % scale factors

% experiment parameters
Params.General.N_run                   = 1;  % how many runs in this program
Params.General.Dur_stimuli             = 8;  %secs, duration for a motion stimuli
Params.General.Dur_blank               = 8;  %secs, duration for blank after stimuli;
Params.General.Spatial_envelope        = 0;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
Params.General.nTrials                 = 16; % how many trials we want, each trial has a left/static/right/static loop.

% motion texure parameters
Params.General.Background              = [127 127 127];%rgb value
Params.General.Character_color         = [0 0 0];%rgb value
Params.General.Contrast                = 100;    % 100%
Params.General.Speed                    = 6; %deg/sec;

%do some simple computation
Params.General.White                   = 254;
Params.General.Amplitude               = 127;
Params.General.H_ecc_stim              = round(Params.General.Stimuli_ecc(1)*60/Params.General.Scale_factor);
Params.General.V_ecc_stim              = round(Params.General.Stimuli_ecc(2)*60/Params.General.Scale_factor);

direction_list                          =  rem(2:Params.General.nTrials+1,2);   %directions, 0,left,1,right;

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
    %housekeeping stuff
    stimulus_radius  = round(60* Params.General.Env_radius/Params.General.Scale_factor);%how many pixels for stimuli radius    
    if Params.General.Linearize
        screen_clut = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
        Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    end
    Screen('FillRect',w, Params.General.Background);Screen('Flip', w);
    Screen('TextSize',w,30);Screen('TextFont',w,'Arial');
    sr_hor = round(screen_rect(3)/2); sr_ver = round(screen_rect(4)/2);
    
    
    %Now I know frame_rate, do some simple computations
    Params.General.MvLength_stimuli       = round(Params.General.Frame_rate*Params.General.Dur_stimuli);
    Params.General.MvLength_blank         = round(Params.General.Frame_rate*Params.General.Dur_blank);
    Params.General.PixStep                = Params.General.Speed*60/(Params.General.Frame_rate * Params.General.Scale_factor);        %moving speed, how many pixel/frame
       
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
    screen_rect = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
      
    exit_index = 0;
    
    
    %% premake the motion movie
    %------------------------make the movie----------------------------
    % make a wide image patch with dots.
    bigPatch_X = ceil( bps + Params.General.PixStep * Params.General.MvLength_stimuli);
    big_im     = rand(bps,bigPatch_X);
    
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

    
    
    %% MAIN LOOP
    HideCursor;
    tic; trial = 1;
        
    Screen('FillRect',w, Params.General.Background);
    Screen('DrawText',w,'Press any key to start the experiment',sr_hor-300,sr_ver-80,Params.General.Character_color);
    Screen('DrawText',w,'Please keep your fixation on center spot throught the experiment',sr_hor-300,sr_ver-40,Params.General.Character_color);Screen('Flip',w);
    KbWait(-1);
    Screen('Flip', w);
    FlushEvents('keyDown');
    
    %% wait for trigger
    disp('Waiting for trigger from scanner or the tilde button (~)');
    KbName('UnifyKeyNames');
    trigger = KbName('`~');
    [~,~,key] = KbCheck(-1);
    while ~key(trigger);
        WaitSecs(.001);
        Screen('DrawText',w,'Waiting for trigger...',sr_hor-300,sr_ver-80,Params.General.Character_color);Screen('Flip',w);
        [~,~,key] = KbCheck(-1);
    end
    startTime = GetSecs;
    
    % we put the first blank out of cycle;
    Screen('FillRect',w, Params.General.Background);
    Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
    vbl=Screen('Flip', w);
    
    % first blank at begining
    Screen('FillRect',w, Params.General.Background);
    Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
    vbl_blank=Screen('Flip', w,vbl+Params.General.Dur_blank); % 0.6295 is estimated time for preDelay below
    
    
    %% Main loop    
    for trial=1:Params.General.nTrials
           
        %% play the movie to the left
        priorityLevel=MaxPriority(w);Priority(priorityLevel);
        trialStartTime = GetSecs;
        preDelay =  trialStartTime-vbl_blank
        prt.data(trial,4) = (GetSecs-startTime)*1000;%starting time;
        
        if trial ==1
            fprintf('This blank lasts %.10f\n', GetSecs-vbl); 
        else
            fprintf('This blank lasts %.10f\n', GetSecs-trialEndTime); 
        end
   
        
        for frame = 1:Params.General.MvLength_stimuli
            offset=(frame-1)*Params.General.PixStep;  %how many pixel moved in this frame
            %because of left motion, source rect should move to from left
            %to right
            %offset
            srcrect=[offset 0 bps+offset bps];

            %Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawTexture',w,bigPatch,srcrect,screen_rect,direction_list(trial)*180);
            Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('DrawTexture',w,maskCOS,movie_rect,screen_rect,[],[]);
            Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
            vbl=Screen('Flip',w);                               
        end

        %% 
        trialEndTime = GetSecs;

        prt.data(trial,5)=(trialEndTime-startTime)*1000;%end time;
        
        fprintf('This trial lasts %.10f seconds \n',trialEndTime - trialStartTime);

            trial
%            t = (prt.data(trial,4) -  prt.data(trial-1,5))/1000 


        % blank after left motion    
        Screen('FillRect',w, Params.General.Background);
        Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);
        blankStartTime=Screen('Flip',w);
        Priority(0);
        Screen('FillRect',w, Params.General.Background);
        Screen('FillOval', w,Params.General.Character_color,[sr_hor-15, sr_ver-15, sr_hor+15, sr_ver+15]);        
        vbl_blank=Screen('Flip', w, trialEndTime+Params.General.Dur_blank-preDelay); % give some room for preDelay below

        FlushEvents('keyDown');
        % Close movies
        %Screen(movie_play, 'Close');
        clear frame;
    end
    ShowCursor;
    TotalTime = GetSecs-startTime;
    fprintf('The whole run took %.10f\n', TotalTime)
    
    
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    sca;
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    ListenChar(1);
       
    save(filename);
    
    %% set up prt
    % some info for prt file
    % for motion localizer prt. We use this to define motion selective
    % voxels
    prt.name = [Params.General.Experiment '_localizerPRT'];
    prt.NrofCodnitions                      = 2;%left/right/blank;
    prt.Condition{:}.ntpts                  = 8;
    prt.Condition{1}.name                   ='Left';
    prt.Condition{2}.name                   ='Right';
    colors = 255*jet;close(gcf);
    prt.colors = round(linspace(1,length(colors),prt.NrOfConditions));
    prt.colors = round(colors(prt.colors,:));
    
    
    % for single trials prt
    % we need single trial prt in order to perform signal trial decoding
    % analysis and know the best achievable decoding performance.
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
    prtSingleTrial.Condition{16}.name='DirRight_FieldRight_trial4'; %direction_field_trial4


    
    %compute prt for three conditions (left/right/blank);
    for i =1:prt.NrOfConditions
        prt.Condition{i}.ntpts = 8;% 8 trials/condition
        prt.Condition{i}.estart = prt.data(prt.data(:,3)==i,4);
        prt.Condition{i}.eend = prt.data(prt.data(:,3)==i,5);
        prt.Condition{i}.color=prt.colors(i,:);
    end
    %add gaze condition
    if ~DEBUG && Params.General.EyeTrack
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
    if ~DEBUG && Params.General.EyeTrack
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

catch
    ListenChar(1);
    ddd = lasterror;
    ddd.message
    psychrethrow(lasterror);
    screen_clut = [0:255; 0:255; 0:255]'./(255);
    Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);
    sca;
    Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
    Priority(0);
end %try..catch..



