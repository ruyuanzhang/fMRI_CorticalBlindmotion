%% Confirmation of Gaze
% The purpose of this function is to ensure that the participant is fixated
% on the center of the screen during behavior/MRI testing using psychtoolbox.
% The major strength of this function is that it will give the participant
% a warning and pause the program when he or she diverts gaze. The major
% deficiency of this function comes from its generalizability: It cannot
% have any decision-making power with regards to how to continue stimulus
% presentation in the orginal program. Some experiments may require the
% current stimulus to be thrown out, while others may simply carry on the
% presentation. How this function is integrated into the full experiment
% will have influence (but not total control) over how to respond to
% diverted gaze during a task.
% 
% Confirmation, or a confirmed gaze, is when gaze coordinates are within
% the allowed bounding box for a gaze-dependent program. Momentary
% missing data due to blinking or some error does NOT qualify as
% non-confirmation. Confirmation is checked if and only if data is not
% missing. In certain situations based on the parameters you can set, the
% function will either try for re-confirmation (checking again for fixation
% after finding the gaze outside the bounding box), or will simply leave the
% function non-confirmed (that is, there is data showing where the gaze is,
% but that gaze is not inside the allowed fixation bounding box).
%
% ENDSTATUS gives information about the gaze. The function always returns
% one of the following:
% (1) 'confirmed' - gaze in bounding box on first try
% (2) 'diverted but confirmed' - gaze left the bounding box but was
% re-confirmed
% (3) 'outside bounds' - gaze is outside bounds (only for use in
% non-gaze-dependent mode)
% (4) 'elSettings not found' - there was no input argument elSettings
% (5) 'timeout' - gaze was diverted but was not reconfirmed for too long,
% the duration of which can be set in the parameters
% (6) 'keyboard' - gaze confirmation was aborted with keyboard press ('q')
% (7) 'missing data' - couldn't find an eye or missing data so there was no
% confirmation attempt
% (8) 'no new float sample' - additional data not collected from eyetracker
% 
% In certain cases during fixation confirmation, the screen changes to
% warn the participant of a non-confirmed gaze. However, it would be nice
% put back the old screen once we have re-confirmation. This will require
% duplicating the pre-confirmation check PTB Screen function in an IF
% statement following the call to the current function, given certain
% values of ENDSTATUS.
%   % elSettings = confirmationFixation(elSettings);
%   % switch elSettings.endstatus(end)
%   %     case {2, 3, 5, 6} %these are the endstatuses that alter the visible screen
%   %         Screen(...);
%   %         Screen('Flip',...);
%   % end
% 
% Function compiled by EG Gaffin-Cahn 10/4/12, partially adapted from demos
% and examples from the Eyelink Toolbox in Psychtoolbox.
%
% TO DO:
%   -specific button to exit gaze confirmation - STILL NOT WORKING?
%   -Standalone mode - start doing that
%   -drift correction from timeout is having trouble re-integrating into
%   the function and re-checking gaze - STILL NOT WORKING?
%
% Last update: 7/16/13 1:00PM

function elSettings = confirmFixation(elSettings)
%% Standalone mode
% Should only enter this IF statement if this function is run in
% standalone, rather than being called by a psychophysical behavioral task.
if ~Eyelink('IsConnected');
    standalone = 1;
    %connect to EyeLink...
    
%% Set the parameters from elSettings
% ELSE will be invoked when this function is being called from another
% function where the eyelink setup has already been accomplished, including
% calibration and setting of the parameters.
else
    try %takes the elSettings structure argument and ensures everything is there
        Eyelink('Message', 'Attempting to confirm fixation');
        
        if ~isfield(elSettings,'el')
            disp('You will need to pass elSettings.el in your argument.')
            disp('''el'' contains the settings and parameters for the eye tracker')
            return
        end
        el = elSettings.el;
        
        [width, height] = Screen('WindowSize', elSettings.el.window);
        if isfield(elSettings,'xc') 
            width = elSettings.xc*2;
        end
        if isfield(elSettings,'yc') 
            height = elSettings.yc*2;
        end

        if ~isfield(elSettings,'bounds') %allowed boundaries for eye fixation
            halfside = 25;
            elSettings.bounds = [width/2-halfside, height/2-halfside, width/2+halfside, height/2+halfside];
        end
        bounds = elSettings.bounds;
        
        %field in elSettings structure        %default value                      %description of setting
        %parameters
        if ~isfield(elSettings,'thick');      elSettings.thick=2;                 end %thickness of the boundary box
        if ~isfield(elSettings,'saveGaze');   elSettings.saveGaze=0;              end %1 saves gaze in gazeLoc (BEWARE: very large array)
        if ~isfield(elSettings,'gazeDepend'); elSettings.gazeDepend=1;            end %1 prevents continuing until gaze is in bounds
        if ~isfield(elSettings,'drawGaze');   elSettings.drawGaze=0;              end %1 draws a circle around current gaze if outside bounds
        if ~isfield(elSettings,'flipToWarn'); elSettings.flipToWarn=1;            end %1 will go through with function flipToFixation
        if ~isfield(elSettings,'timeout');    elSettings.timeout=1;               end %1 timeout will kill confirmation, otherwise, drift correction
        if ~isfield(elSettings,'expiry');     elSettings.expiry = 10;             end %time until timeout when cannot confirm gaze
        if ~isfield(elSettings,'bkgdcolor');  elSettings.bkgdcolor=[128 128 128]; end %color of the background during warnMSG
        if ~isfield(elSettings,'centerMSG');  elSettings.centerMSG=1;             end %centers warnMSG on gaze if fixation is broken
        if ~isfield(elSettings,'printInfo');  elSettings.printInfo=1;             end %prints out endstatus
        %other
        if ~isfield(elSettings,'eye_used');   elSettings.eye_used=-1;             end %which eye being used, default -1 will be fixed later on
        if ~isfield(elSettings,'rect');       elSettings.rect=[0 0 width height]; end %size of screen
        if ~isfield(elSettings,'xc');         elSettings.xc=width/2;              end %y center of screen
        if ~isfield(elSettings,'yc');         elSettings.yc=height/2;             end %x center of screen
        if ~isfield(elSettings,'gazeloc');    elSettings.gazeloc=[];              end %TIMEx2 array of gaze locations
        if ~isfield(elSettings,'endstatus');  elSettings.endstatus = [];          end %makes sure endstatus exists
        
        %eye_used is too bulky a variable name when appended to the
        %elSettings structure. To convert this to a matlab-friendly index,
        %we add +1 when refering using eye_used
        eye_used = elSettings.eye_used;
        elSettings.warnMSG = 'Please focus on the center of the screen';
        quitCheck  = KbName('q');
        standalone = 0; %1 for standalone mode
        gazeOut    = 1; %True when confirming fixation so that certain
                        %commands will not be called over again when gaze 
                        %leaves the bounding box and the program waits in
                        %the while loop for the gaze to re-center.
    catch %#ok<CTCH>
        disp('No elSettings argument provided..')
        Eyelink('Message', '4: elSettings not found');
        elSettings.endstatus(end+1) = 4;
        return
    end
end
    
%% Gaze check prep
% This section prepares for gaze fixation confirmation by getting the
% newest float sample, the eye being tracked, the location of the gaze,
% assuming it is not missing.
WaitSecs(.002); %this is the ONLY WaitSecs if neither gazeDepend nor standalone
if Eyelink('NewFloatSampleAvailable') > 0; %checks if new eyelink data is available
    evt = Eyelink('NewestFloatSample'); %returns special eyelink structure from newest dataset
    
    %If the eye being tracked has not been set yet. Sometimes this fails to
    %catch on for whatever reason. If it passes 3 seconds without finding
    %which eye, it times out and displays the problem. Maybe we should
    %include going back to the Camera Setup screen?
    timer = GetSecs;
    while eye_used == -1
        eye_used = Eyelink('EyeAvailable'); %find out which eye is being tracked
        if eye_used == el.BINOCULAR %if both are being tracked, only record RIGHT
            eye_used = el.RIGHT_EYE;
        end
        if GetSecs-timer > 3;
            Eyelink('Message','Cannot find which eye to use');
            disp('Cannot find which eye to use after 3 seconds...');
            break
        end
    end
    
    x = evt.gx(eye_used+1); %gets eye x location
    y = evt.gy(eye_used+1); %gets eye y location
    
    % checks to see if data has been collected and eye is found
    if x~=el.MISSING_DATA && y~=el.MISSING_DATA && evt.pa(eye_used+1)>0
        eyedriftStart = GetSecs;
        endstatus = 1;
        
        %% Gaze check
        % This is the WHILE loops that checks if the gaze location is
        % within the bounding box allowed. The syntax for bounds is
        % [-x, -y, +x, +y], ie. [top left, top right, bottom
        % right, bottom left]
        while x < bounds(1) || y < bounds(2) || x > bounds(3) || y > bounds(4)
            
            %breaks the while loop if not in gaze-contingent mode
            if ~elSettings.gazeDepend
                endstatus = 3;
                Eyelink('Message','3: outside bounds');
                break;
            end
            
            % only enters IF gazeOut the first time into the function, it does
            % not repeat during the encompassing while loop.
            if gazeOut && elSettings.flipToWarn
                gazeOut = 0;
                elSettings.x = x; elSettings.y = y;
                flipToFixation(elSettings);
            end
            
            %re-measures eye location to check if gaze is in bounding box.
            WaitSecs(.01);
            if Eyelink('NewFloatSampleAvailable') > 0;
                evt = Eyelink('NewestFloatSample');
                x = evt.gx(eye_used+1);
                y = evt.gy(eye_used+1);
            end

            %set the eye locations to -1 so that we can log missing data
            if x == el.MISSING_DATA || y == el.MISSING_DATA;
                disp('missing data now...it was fine before')
                x = -1; y = -1;
            end
            
            %Keeps track of eye locations during non-confirmed gaze. If,
            %during the wait for the participant's eyes to move back to the
            %bounding box, the experimenter decides to continue anyway, he
            %or she will have the option to press (any) key to exit the
            %necessity for confirmation. Only saves the gaze locations if
            %in saveGaze mode.
            if elSettings.saveGaze
                tempEyeLoc = round([x,y]);
                elSettings.gazeloc = vertcat(elSettings.gazeloc, tempEyeLoc);
                WaitSecs(0.1);
            end
            endstatus = 2;
            
            %This allows exit from confirmationFixation even if gaze has
            %not been re-confirmed. If eyes are outside the bounding box
            %for x consecutive seconds, something, may be wrong.
            if GetSecs - eyedriftStart > elSettings.expiry;
                if elSettings.timeout %allows function to timeout
                    endstatus = 5;
                    Eyelink('Message','5: timeout');
                    break;
                else
                    EyelinkDoDriftCorrection(elSettings.el);
                end
                %Put bounding box back after drift correction and reset the
                %drift clock, otherwise it will do drift correction on the
                %very next gaze location check
                elSettings.x = x; elSettings.y = y;
                if elSettings.flipToWarn
                    flipToFixation(elSettings);
                end
                eyedriftStart = GetSecs;
                x = evt.gx(eye_used+1); %gets eye x location
                y = evt.gy(eye_used+1); %gets eye y location
            end
            
            %This allows exit from confirmationFixation even if gaze has
            %not been re-confirmed. Checks for the quit keyboard button ('q')
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(quitCheck)
                endstatus = 6;
                Eyelink('Message','6: keyboard');
                break;
            end;
        end
        
%% End status record
% This section records the endstatus of the gaze, assuming there was no
% timeout,  keyboard-triggered manual stop, or gaze is not confirmed and we
% are not in gazeDepend mode.
        
    else
        endstatus = 7;
        Eyelink('Message','7: missing data');
    end
else
    x = -1; y = -1;
    endstatus = 8;
    Eyelink('Message','8: no new float sample');
end

if endstatus == 1
        Eyelink('Message', '1: confirmed');
elseif endstatus == 2
        Eyelink('Message', '2: diverted but confirmed');
end

if standalone
    WaitSecs(0.05);
    Eyelink('StopRecording');
    Eyelink('CloseFile'); 
    Eyelink('Shutdown');
    sca;
end

%% Final saving
% This section prints out the endstatus and logs the output information
% from the gaze confirmation.

% Prints out the endstatus if it is different than the previously
% registered endstatus.
if elSettings.printInfo && length(elSettings.endstatus) > 1 && endstatus ~= elSettings.endstatus(end)
    fprintf('\nEnd status code: %d\n',endstatus);
end

%keeps track of eye locations
if elSettings.saveGaze
    tempEyeLoc = [round(x) round(y)];
    elSettings.gazeloc = vertcat(elSettings.gazeloc, tempEyeLoc);
end
elSettings.endstatus(end+1) = endstatus; %log the final gaze status

%% Show bounding box
% When gazeDepend is on and the participant's gaze leave the allowed box,
% this function is called to show the boundary box and give a warning
% message to keep fixated on the center of the screen.
function flipToFixation(elSettings)

window     = elSettings.el.window;
xc         = elSettings.xc;
yc         = elSettings.yc;
x          = elSettings.x;
y          = elSettings.y;
bounds     = elSettings.bounds;
warnMSG    = elSettings.warnMSG;
centerMSG  = elSettings.centerMSG;
thick      = elSettings.thick;
bkgdcolor  = elSettings.bkgdcolor;
rect       = elSettings.rect;
gazeCircle = [x-10 y-10 x+10 y+10]; %size of the gaze location circle in gazeDepend mode
% pointList  = [...
%     bounds(1) bounds(2);...
%     bounds(1) bounds(4);...
%     bounds(3) bounds(4);...
%     bounds(3) bounds(2)];

%Puts the warning message in the center of the screen, or where the
%non-confirmed gaze is, depending on the setting called centerMSG
if centerMSG %put the message on the center of the screen
    MSGxpos = xc-100;
    MSGypos = yc-100;
else %put the message at the current gaze location
    MSGxpos = x-100;
    MSGypos = y;
end

%Draws a black circle at current gaze if the gaze leaves the bounding box
%and we are in drawGaze mode
if elSettings.drawGaze
    Screen('FillOval',window,[0 0 0],gazeCircle);
end
    
%Draws the fixation box
Screen('FillRect',window,bkgdcolor,rect);
Screen('FrameRect',window,[220 220 220],bounds,thick);
Screen('FillPoly',window,[0 0 0],...
    [xc-6,yc+2; xc-2,yc+2; xc-2,yc+6;...
     xc+2,yc+6; xc+2,yc+2; xc+6,yc+2;...
     xc+6,yc-2; xc+2,yc-2; xc+2,yc-6;...
     xc-2,yc-6; xc-2,yc-2; xc-6,yc-2]);
Screen('DrawText',window,warnMSG,MSGxpos,MSGypos,[0 0 0]);
Screen('Flip',window);