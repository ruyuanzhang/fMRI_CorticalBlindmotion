function motion_localizer
% radial dot motion localizer
% moving-static cycles for MRI experiment, MT detection
% uses DrawDots.mex which is a modification of the DOTS.mex package from
% Shadlen

% written by Keith Schneider, 1999 September 24,26
% ks cleaned up 10/12/99
% ks modification for alternating directions (per Huk, Dougherty & Heeger)
% ks modified for PTB 3
% edited by DT

% -----------------------------
% set the experiment parameters
% -----------------------------

p.mon_width  	= 42.8;	% horizontal dimension of viewable screen (cm)
p.v_dist 		= 79;	% viewing distance (cm) (31in, measured by frank
p.norm_speed 	= 8;	% normal dot speed (deg/sec)
p.ndots 		= 1000;	% number of dots
p.max_d 		= 15;	% maximum radius of  annulus (deg)
p.min_d 		= 0.75;	% minumum
p.alt_time		= 1;	% time between inward/outward direction alternation (sec)
p.dot_d			= 0.1;	% dot diameter (deg)
p.move_dur		= 16;	% duration of moving dots (sec)
p.static_dur	= 16;	% duration of static dots (sec)
p.ncycles		= 8;	% number of stimulus cycles
p.fkill			= 0.05;	% faction of dots to terminate per frame, to produce limited lifetime

%esc_key = KbName('escape'); % press this key to abort

% ---------------
% open the screen
% ---------------

screens=Screen('Screens');
screenNumber=max(screens);
[w, rect] = Screen('OpenWindow', screenNumber, 0);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

xc = rect(3)/2; % x center
yc = rect(4)/2; % y center

HideCursor;	% Hide the mouse cursor
bg = [0 0 0];	% background color
dot_color = [255 255 255];  % dot color`
black_index = [0 0 0];
white_index = [255 255 255];

Screen('TextSize', w, 24);
Priority(MaxPriority(w));

Screen('FillRect', w, bg);

ifi=Screen('GetFlipInterval', w);
p.frame_rate=1/ifi;	% frames per second

% -------------
% set up arrays
% -------------

pix_per_deg = pi * rect(3) / atan(p.mon_width/p.v_dist/2) / 360;	% pixels per degree
ds_to_pf = pix_per_deg / p.frame_rate;	% conversion factor from deg/sec to pixels/frame
pfs = p.norm_speed * ds_to_pf;	% normal dot speed (pixels/frame)
dot_d = pix_per_deg*p.dot_d;

% -----------------------------------
% set up dot positions and velocities
% -----------------------------------

rmax = p.max_d * pix_per_deg;	% maximum radius of annulus (pixels from center)
rmin = p.min_d * pix_per_deg; 	% minimum
r = (rmax - rmin) * sqrt(rand(p.ndots,1)) + rmin;	% r
t = 2*pi*rand(p.ndots,1);					% theta
cs = [cos(t), sin(t)];
xy = [r r] .* cs;	% dot positions in Cartesian pixel coordinates relative to center

mdir = 2 * floor(rand+0.5) * ones(p.ndots,1) - 1;	% motion direction (in or out) for each dot
dr = pfs * mdir;	% change in radius (in pixels) per frame
dxdy = [dr dr] .* cs;	% change in x and y (in pixels) per frame


% --------------------
% start experiment now: draw fixation point and text and wait for key press to begin
% --------------------

%basic set up
KbName('UnifyKeyNames'); % make KbName function work across platforms
commandwindow;
exitScript = KbName('q'); % to quit the program early, press Q
trigger = KbName('`~');

% wait for trigger
Screen('TextSize',w,40);
Screen('DrawText',w,'Waiting for trigger...',50,50,255);
Screen('Flip',w);
[~,~,keyCode] = KbCheck;
while ~keyCode(trigger)
    [~,~,keyCode] = KbCheck;
end

%after trigger is received...
Screen('Flip',w); %flip screen to remove text

txt = 'Fixate.  No task.';
normBoundsRect = Screen('TextBounds', w, txt);
txtloc = [xc - normBoundsRect(3)/2, yc + normBoundsRect(4)/2];
[newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),white_index);

% draw fixation point
Screen(w, 'FillRect', black_index, [xc-3 yc-3 xc+3 yc+3]);
Screen(w, 'FillRect', white_index, [xc-2 yc-2 xc+2 yc+2]);
Screen(w, 'FillRect', black_index, [xc-1 yc-1 xc+1 yc+1]);

Screen('Flip', w);
pause(16);

switch_mov = cumsum(repmat([p.move_dur p.static_dur],1,p.ncycles));
qmov = 1;
thetime = 0;
imov = 1;
nextalt = p.alt_time;
start_time = GetSecs;

while thetime <= p.ncycles*(p.move_dur + p.static_dur)

    Screen('DrawDots', w, xy', dot_d, dot_color, [xc yc], 1);
    % draw fixation point
    Screen('FillRect', w, black_index, [xc-3 yc-3 xc+3 yc+3]);
    Screen('FillRect', w, white_index, [xc-2 yc-2 xc+2 yc+2]);
    Screen('FillRect', w, black_index, [xc-1 yc-1 xc+1 yc+1]);
    Screen('Flip', w);

    if thetime >= switch_mov(imov)
        imov = imov + 1;
        qmov = 1 - qmov;
    end
    if thetime >= nextalt
        nextalt = nextalt + p.alt_time;
        dxdy = -dxdy;
        dr = -dr;
    end

    if qmov

        % -------------------------------------------------
        % update the dot positions for the next video frame
        % -------------------------------------------------

        xdots(1:2:2*p.ndots-1,1:2) = xy;	% old dot positions to be erased
        xy = xy + dxdy;						% move dots
        r = r + dr;							% update polar coordinates too

        % check to see which dots have gone beyond the borders of the annuli

        L = find(r > rmax | r < rmin | rand(p.ndots,1) < p.fkill);
        nL = length(L);

        if nL

            % choose new r and angle

            r(L) = (rmax-rmin) * sqrt(rand(nL,1)) + rmin;
            t(L) = 2*pi*(rand(nL,1));

            % now convert the polar coordinates to Cartesian

            cs(L,:) = [cos(t(L)), sin(t(L))];
            xy(L,:) = [r(L) r(L)] .* cs(L,:);

            % compute the new cartesian velocities

            dxdy(L,:) = [dr(L) dr(L)] .* cs(L,:);
        end
    end

    [keyIsDown,secs,keyCode] = KbCheck;
%    if keyCode(esc_key)
%        break
%    end


    thetime = GetSecs - start_time;
end

Priority(0);
Screen('Flip',w); %flip screen to remove text
pause(16);
ShowCursor;
Screen('Close', w);
