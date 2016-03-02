function elSettings = SetupEyeLink(Params,w,rect)

elSettings.gazeDepend = 0; % enforce fixation (0=don't wait for fixation to continue)
elSettings.flipToWarn = 0; % don't change the screen when checking fixation
elSettings.printInfo = 0; % don't print out a bunch of info during eye checking
elSettings.bounds = [rect(3)/2-200,rect(4)/2-200,rect(3)/2+200,rect(4)/2+200]; % allowed bounds for eyes
elSettings.el = EyelinkInitDefaults(w); % initialize eyelink (el) defaults

EyelinkInit(0); % initializes connection
Eyelink('OpenFile',sprintf('%s%d_R%d.edf',Params.General.Experiment(1:3),Params.General.Sub,Params.General.Run)); % name the edf file

Eyelink('Command', 'add_file_preamble_text ''%s''',Params.General.Experiment);

% give experimenter a few seconds to skip the custom calibration (and
% continue on to drift correct in RunMapper), but you can hit any keyboard
% button to definitely do it
tic; wait = 10;
while toc < wait
    Screen('DrawText',w,sprintf('Click to skip custom calibration: %d',ceil(wait-toc)),50,50);
    Screen('Flip',w);
    [~,~,buttons] = GetMouse(w);
    if any(buttons) % click the mouse, skip the calibration
        return
    end
    if KbCheck % hit the keyboard, do the calibration
        break
    end
end

Eyelink('Command', 'screen_pixel_coords = %ld %ld %ld %ld', 0, 0, rect(3)-1, rect(4)-1); % sets screen size
Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, rect(3)-1, rect(4)-1); % sets display size

Eyelink('Command', 'calibration_type = HV5');
Eyelink('Command', 'generate_default_targets = NO');

% instructions: move the mouse around and talk to subject until you find a
% good place. When you have a good place for a calibration spot, click the
% mouse.
samples = 5;
B = 25;
S = 5;
xc = []; yc = [];
found = 0;
for i = 1:samples % scroll through each calibration dot location
    while length(xc) < i % continue this code while a new point hasn't been set
        [x,y,buttons] = GetMouse(w);
        if any(buttons) % if the mouse was clicked, save the location
            xc(i) = x; yc(i) = y;
            found = found + 1;
            while any(buttons); [~,~,buttons] = GetMouse(w); end % wait for button release
        end
        % draw current location of mouse
        Screen('DrawText',w,sprintf('%d / %d',i,samples),50,50);
        Screen('FillRect',w,[128 128 128]);
        Screen('FillOval',w,[0 0 0],[x-B, y-B, x+B, y+B]);
        Screen('FillOval',w,[128 128 128],[x-S, y-S, x+S, y+S]);
        Screen('Flip',w);
    end
end

% sca;

Eyelink('Command','calibration_samples = %d',samples+1);
Eyelink('Command','calibration_sequence = 0,1,2,3,4,5');
Eyelink('Command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',...
    xc(1),yc(1),  xc(2),yc(2),  xc(3),yc(3),  xc(4),yc(4),  xc(5),yc(5) );
Eyelink('Command','validation_samples = %d',samples);
Eyelink('Command','validation_sequence = 0,1,2,3,4,5');
Eyelink('Command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',...
    xc(1),yc(1),  xc(2),yc(2),  xc(3),yc(3),  xc(4),yc(4),  xc(5),yc(5) );

Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'); % what info is available in the edf
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'); % what info is available on-line
Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA,HTARGET,GAZERES,STATUS,INPUT');

if Eyelink('IsConnected')~=1; % make sure we're still connected.
    error('EyeLink no longer connected\n');
end

EyelinkUpdateDefaults(elSettings.el); % update the defaults with our newly edited settings
EyelinkDoTrackerSetup(elSettings.el); % do calibration and validation
