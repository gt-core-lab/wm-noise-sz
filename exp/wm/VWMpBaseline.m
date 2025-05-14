%% VISUAL WORKING MEMORY PRECISION TASK (WITH ORIENTATION) - BASELINE RUN & VWM INSTRUCTIONS/PRACTICE
%
% This version should be run prior to the spatial delayed response WM task that requires subject to report
% orientation of probe without memory burden
%
% Megan Ichinose & Woon Park 
% April 20th, 2016
%
% Updated May 27th, 2016 - Megan Ichinose
% Update includes USING GRIFFIN POWERMATE DIAL, not arrow keys
% 
% =========================================================================
%               INPUT PARAMETERS                 
% =========================================================================

warning('OFF'); clc; sca;
Screen('CloseAll');
warning('off','MATLAB:dispatcher:InexactMatch')
oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'SkipSyncTests', 1);
PsychDefaultSetup(2);

%Specify display parameters (added by WP 05/10/16)
    
    scale_factor     = 1.7;        % most important parameter - how many acrmin is one screen pixel? 
    frame_rate       = 120;        % 
    resolution       = [1024 768]; % [1024 640];
    linearize 	     = 1;          % whether monitor is linearized. if=1, program will look for "MyGammaTable"
    
    %Fixation parameters (added by MI 05/19/16)
    V_ecc_fix = 0;
    H_ecc_fix = 0;
    V_ecc_fix = V_ecc_fix*60/scale_factor;
    H_ecc_fix = H_ecc_fix*60/scale_factor;
    l = 7; % fixation cross line length 
    ww = 4; % fixation cross line width 

%Specify some task parameters

    %WM set size, number of trials/set, blocks, etc. 
    set_size        = 1;      %Set sizes of 1,2,4
    num_set         = 24;     %Number of trials per set size  
    block_length    = 24;     %Equal number of trials in each set per block
    num_blocks      = 1;      %Number of blocks

    T_num = set_size*num_set;  %Total number of trials, also = block_length*num_blocks

    %Specify text 
    txtsz           = 25;      %Size of all presented text
    txtfont         = 'Arial'; %Text font
    
    %Specify hypotenuse for circle array (Model stim will be presented at
    %length of Hypot from center fix (edited by WP 05/10/16)
    Hypot_deg       = 4;       %Degree (Distance between center/fixation & stimuli)
    Hypot           = round(Hypot_deg*60/scale_factor); %pixels (original: 150)
    
    %Specify timing of stimuli
    PresDur         = 2;       %Stimulation duration (s)
    delay           = 1;       %Delay time (s)
    
    %Stimuli parameters (edited by WP 05/10/16) 
    LineCol = 255; %black=0, white=255
%     LineWidth = 10; %Pixels of line width
%     LineHeight = 55; %Pixels of line hieght
    SF                  = 1;    % cycle/deg
    radius              = 1;    % deg
    sigma               = 0.75; % SD for gaussian envelope 
    rAngle              = 0;    % reference angle (0)
    grating_contrast    = 1;    % 1 for 100%
    stimulus_radius     = round(60*radius/scale_factor);
    f                   = (SF*scale_factor/60)*2*pi;
    rAngle              = rAngle*pi/180;
    a                   = cos(rAngle)*f; 
    b                   = sin(rAngle)*f;
    
    %Specify output directory name
    dirout = '/Data/';
    
% =========================================================================
%               WRITE SUBJECT DATA FILE              
% =========================================================================


ID = input('Enter subject name or ID number:', 's');
SN = str2num(ID);
if ~exist('ID','var')
    ID=999;
end

%Warn if duplicate subID
FileName =  strcat(ID, '_VWMpBaseline.csv');
if exist(strcat('Data/',FileName),'file')
    resp=input(['The file ' FileName ' already exists./nDo you want to overwrite it? [Type ok for overwrite]'], 's');
    if ~strcmp(resp,'ok') %abort experiment if overwriting not confirmed
        disp('experiment aborted')
        return
    end
end

%Create datafile
FID = fopen(FileName, 'w');
if FID == -1
    fprintf(1,'File Not Open');
end
fprintf(FID, 'Trial, T1Loc, ProbeAngle, Theta, RespAngle, RespError, RT\n');
fclose(FID);


% =========================================================================
%               OPEN & INITIALIZE SCREEN 
% =========================================================================

screenNumber = max(Screen('Screens')); %Set screen number

numBuffers = 2; %Double buffering (use >2 for development, 2 for experiment)
screenColor = 126; %Light gray (edited by WP 05/10/16; previously: [155 155 155])
[w, wRect] = Screen('OpenWindow', screenNumber, screenColor, [], [], numBuffers);

%Create this for fixation
sr_hor = round(wRect(3)/2); sr_ver = round(wRect(4)/2);

Screen('Flip', w);
[CX,CY] = RectCenter(wRect); %Get center coordinates of window

if linearize % added by WP (05/10/16)
    load('mygammatable.mat');
    Screen('LoadNormalizedGammaTable', screenNumber, mygammatable);
end

Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); %Set blend function for screen
maxPriorityLevel = MaxPriority(w); %Query max priority level
Priority(maxPriorityLevel); %Set priority level
timing=Screen('GetFlipInterval',w); % get the flip rate of current monitor.

rand('twister',sum(100*clock)); %Reset random generator
HideCursor;

% =========================================================================
%               CREATE STIMULI & SET KEY RESPONSES
% =========================================================================

% edited by WP 05/10/16

% % Create line stimulus for memory array
% rectMat=repmat(LineCol,LineWidth,LineHeight); %Line params from above
% MemLine = Screen('MakeTexture', w, rectMat); %Creates line

% Make spatial envelope

[x,y] = meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
bps = (stimulus_radius)*2+1;
circle = ((stimulus_radius)^2-(x.^2+y.^2));

for i = 1:bps
    for j = 1:bps
        if circle(i,j) < 0
            circle(i,j) = 0;
        else
            circle(i,j) = 1;
        end;
    end;
end;

R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);

% Specify keys used
KeyboardIndex = GetKeyboardIndices; %Query keyboard device number
Larrow = KbName('LeftArrow'); %Used for angle rotation
Rarrow = KbName('RightArrow'); %Used for angle rotation
SubmitResp = KbName('UpArrow'); %Submit angle response

 
% =========================================================================
%               CONDITION INITIALIZATION (Create design matrix)
% =========================================================================

%Set variables for output
Order = 1:T_num;
T1_Loc = repmat(1:8, num_set/8, 1); %Target location (1-8)
ProbeAngle = repmat(9999,T_num,1); %Same as target angle
Theta = repmat(9999,T_num,1); %Raw response angle value
RespAngle = repmat(9999,T_num,1); %Converted response angle that is closest to target
RespError = repmat(9999,T_num,1); 
RT = zeros(1, T_num);

%Randomize location trials within blocks

Order = reshape(Order,block_length,num_blocks);

for i=0 : num_blocks-1
    
    N = Order(:,i+1);
    Order(:,i+1) = N(randperm(end),:);
    
end

Order = reshape(Order,T_num,1);


% =========================================================================
%               SETUP POSSIBLE DISPLAY LOCATIONS
% =========================================================================

Block = 1; %Initialize starting block number
num_locs = 8; % Number of possible target locations in the circular memory array
nsearch = num_locs*2;
searchradius = Hypot;

angles = (1:16) * 2*pi/nsearch;
[xPos, yPos] = pol2cart(angles, searchradius); %convert
    %now back to Cartesian
    xPos = xPos + CX;
    yPos = yPos + CY;
    xPos = xPos';
    yPos = yPos';
% PresLocs = [xPos-(LineWidth/2), yPos-(LineHeight/2), xPos+(LineWidth/2), yPos+(LineHeight/2)];
PresLocs = [xPos-round(bps/2), yPos-round(bps/2), xPos+round(bps/2), yPos+round(bps/2)]; % edited by WP 05/10/16

Loc1 = PresLocs(1,:);
Loc2 = PresLocs(3,:);
Loc3 = PresLocs(5,:);
Loc4 = PresLocs(7,:);
Loc5 = PresLocs(9,:);
Loc6 = PresLocs(11,:);
Loc7 = PresLocs(13,:);
Loc8 = PresLocs(15,:);

DispLocs = [Loc1;Loc2;Loc3;Loc4;Loc5;Loc6;Loc7;Loc8];
% DistCenter = [CX-(LineWidth/2), CY-(LineHeight/2), CX+(LineWidth/2), CY+(LineHeight/2)];
DistCenter = [CX-round(bps/2), CY-round(bps/2), CX+round(bps/2), CY+round(bps/2)]; % edited by WP 05/10/16

% =========================================================================
%               BEGIN EXPERIMENT WITH INSTRUCTIONS
% =========================================================================

while KbCheck; end %Call KbCheck at the beginning since it's slow to load the 1st time

% Draw Instructions 
Screen(w,'FillRect',screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
Instructions = ['In this 1st task you will see a small patch with lines in it appear somewhere on the screen.\n\n'...
                'Another patch with lines will then appear at the center of the screen.\n\n.'...
                'Use the DIAL to turn the angle of this center patch to match\n'...
                'the EXACT angle of the lines in the patch that is still on the screen.\n\n'...
                'Once you make the lines in your patch match the existing one as PRECISELY as you can,\n'...
                'press the DIAL to submit your response.'];
DrawFormattedText(w, Instructions, 'center', 'center', [255 255 255]);
Screen('Flip', w);
KbWait(KeyboardIndex,2);

% Open Powermate
handle = PsychPowerMate('Open');


%% ========================================================================
%                BEGIN TRIALS
% =========================================================================

% Fixation to begin
Screen('FillRect', w, screenColor);
Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
        Screen('Flip', w);
        WaitSecs(2);

% added by WP 05/10/16
background = screenColor;
amplitude = background*grating_contrast;

for Trial = 1 : T_num

% Make stimulus (added by WP 05/10/16)
% grating = round(((sin(a*x+b*y+rand*2*pi).*circle*amplitude)+background));
dice = round(rand(1));
if dice == 1
    grating = round(((sin(a*x+b*y+pi/2).*circle*amplitude)+background));
elseif dice == 0
    grating = round(((sin(a*x+b*y-pi/2).*circle*amplitude)+background));
end
MemLine = Screen('MakeTexture',w,grating); % grating texture 

% -------------- PREPARE THE MEMORY ARRAY  -------------------------------

% This section of code will set locations and rotations angles for each
% target in the memory array. The rotation is set to be randomly sampled
% from 0:180, with the parameter than no two items are within 5 degrees from each
% other. This section also selects the target to be probed after the delay,
% and randomly selects an orientation for that probe. 

PossLocs = 1:8; %Possible circle locations
PossAngles = -90:89; %Possible angle rotations (180)

%This resets possible angle rotations for each memory target with each
%trial
T1_Angle = randsample(PossAngles,1);

T1L = T1_Loc(Order(Trial));%get the location coords for T1
Screen('DrawTexture', w, MemLine, [], DispLocs(T1L,:), T1_Angle);

% present memory array with fixation
Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
Screen('Flip', w);
WaitSecs(0.5);


% --------------------- TEST PROBE  -----------------------------------

% Randomize theta angle & get start time
th = datasample(-90:89,1); % Initial rotation angle (degrees)
rot_spd = 2;

%Present Memory Probe
Screen('DrawTexture', w, MemLine, [], DistCenter, th);
Screen('DrawTexture', w, MemLine, [], DispLocs(T1L,:), T1_Angle);
Screen('Flip', w);

bdown = 0;
[~, dialStart] = PsychPowerMate('Get', handle); %Initiate dial position
Begin_Time = GetSecs; 

%Get Response (Rotation)
while(~any(bdown))    
    [mx,my,bdown] = GetMouse;
    
    [button, dialPos] = PsychPowerMate('Get', handle); 
    if dialPos > dialStart %Clockwise
        th = th + rot_spd; dialStart = dialPos;
    end
    if dialPos < dialStart %Counterclockwise
        th = th - rot_spd; dialStart = dialPos;
    end
    Screen('DrawTexture', w, MemLine, [], DistCenter, th);
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1L,:), T1_Angle);
    Screen('Flip', w);
    if button == 1
        theta = th; %Use for later
        End_Time = GetSecs;
        RT(Trial) = (End_Time-Begin_Time)*1000;
        break
    end   
    WaitSecs(0.001);  

end


%Change degree format for probe & response angles from [-90:90) to [-180:180)

%Compute which was the probed (correct) angle
Probe_Angle = wrapTo180(T1_Angle); %This won't actually change it, but keep in for logical step
Resp_Angle = wrapTo180(theta); %This WILL change response (since theta can be >180 in neg and pos directions)

AltResp = 180-abs(Resp_Angle); %Find alternative ('flipped by 180') response
if abs(abs(Probe_Angle)-abs(Resp_Angle)) > abs(abs(Probe_Angle)-abs(AltResp)) %Choose the closest angle
    Resp_Angle = AltResp;
end

%Compute response error
Resp_Error = abs(Probe_Angle)-abs(Resp_Angle);

%Save some trial data
ProbeAngle(Trial) = Probe_Angle; 
Theta(Trial) = th; 
RespAngle(Trial) = Resp_Angle;
RespError(Trial) = Resp_Error;

% ------------ FINISH UP TRIAL & GET DATA --------------------------------

% Clear the screen to finish trial
Screen(w, 'FillRect', screenColor);

% Present ITI/Fixation Display before next trial
    mm = 19;
    for i = 0:4
        nn = mm-i*4;
        Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
        Screen('Flip', w);
        WaitSecs(0.05);
    end
    Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
    Screen('Flip', w);
    WaitSecs(0.7);

% Save trial data
FID = fopen(FileName, 'a');
Save = [Trial; T1_Loc(Order((Trial))); ProbeAngle(Trial); Theta(Trial); RespAngle(Trial); RespError(Trial); RT(Trial)];
fprintf(FID, '%7d, %7d, %7d, %7d, %7d, %3.4f, %3.4f\r\n', Save); 
fclose(FID);

% close screen (added by WP 05/10/16)
Screen('Close', MemLine);
clear grating;

end %Closes the for Trial=1:T_num loop, ending experiment

movefile(FileName,strcat(pwd,dirout)); 
%movefile(matfname,strcat(pwd,dirout)); 

% ------------------------------------------------------------------------
%                     BEGIN INSTRUCTIONS FOR WM TASK 
% ------------------------------------------------------------------------

Screen(w, 'FillRect', screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
MemInstructions = ['Good Job!\n\n'...
                   'In the next, longer task, you will be doing something a little different...\n'...
                   'Patches will first appear somewhere on the screen for a brief time & disappear.\n'...
                   'You must pay close attention & REMEMBER the angles of the lines\n'...
                   'in the patches presented at each location.\n'...
                   'After a very short delay, a patch will then reappear at\n'...
                   'one of the prior locations.\n\n'...
                   'Use the DIAL to change the angle\n'...
                   'of this patch to EXACTLY MATCH the angle of the lines\n'...
                   'in the patch that was previously shown in THAT SPECIFIC location.'];
DrawFormattedText(w, MemInstructions, 'center', 'center', [255 255 255]);
Screen('Flip', w);
KbWait(KeyboardIndex,2);

Screen(w, 'FillRect', screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
MemInstructions2 = ['You are not being timed, so try to be AS PRECISE as possible\n'...
                    'when matching the angle of the patch.\n\n'...
                   'Once you have matched the angle as precisely as possible,\n'...
                   'press the DIAL to submit your response.\n\n'...
                   'You will complete 450 trials total, with periodic breaks.\n'...
                   'See the following slide for an example of what the trials will look like...'];
DrawFormattedText(w, MemInstructions2, 'center', 'center', [255 255 255]);
Screen('Flip', w);
KbWait(KeyboardIndex,2);

WMInstructions = imread('WMinstructions.jpg','jpg');
Screen('PutImage', w, WMInstructions, wRect); % put image on screen
Screen('Flip',w); % now visible on screen
KbWait(KeyboardIndex,2);


%---- BEGIN PRACTICE WM TRIALS (3) ----------------------------------

% added by WP 05/10/16
background = screenColor;
amplitude = background*grating_contrast;

% Fixation to begin
Screen('FillRect', w, screenColor);
Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
Screen(w, 'Flip');
WaitSecs(2);

for Trial = 1 : 3

SetSize=[1,2,4];
% Make stimulus (added by WP 05/10/16)
% grating = round(((sin(a*x+b*y+rand*2*pi).*circle*amplitude)+background));
dice = round(rand(1));
if dice == 1
    grating = round(((sin(a*x+b*y+pi/2).*circle*amplitude)+background));
elseif dice == 0
    grating = round(((sin(a*x+b*y-pi/2).*circle*amplitude)+background));
end
MemLine = Screen('MakeTexture',w,grating); % grating texture 

% -------------- PREPARE THE MEMORY ARRAY  -------------------------------

% This section of code will set locations and rotations angles for each
% target in the memory array. The rotation is set to be randomly sampled
% from 0:180, with the parameter than no two items are within 5 degrees from each
% other. This section also selects the target to be probed after the delay,
% and randomly selects an orientation for that probe. 

PossLocs = 1:8; %Possible circle locations
PossAngles = -90:89; %Possible angle rotations (180)

%This resets possible angle rotations for each memory target with each
%trial (no angle within 5 degrees of another)
    T1_Angle = randsample(PossAngles,1); 
        PossAngles(find(PossAngles>T1_Angle-6 & PossAngles<T1_Angle+6))=[]; 
    T2_Angle = randsample(PossAngles,1);
        PossAngles(find(PossAngles>T2_Angle-6 & PossAngles<T2_Angle+6))=[];
    T3_Angle = randsample(PossAngles,1);
        PossAngles(find(PossAngles>T3_Angle-6 & PossAngles<T3_Angle+6))=[];
    T4_Angle = randsample(PossAngles,1);
        
if SetSize(Trial) == 1 %set size 1
    T_Locs = randsample(PossLocs,1);
    T1_Loc = T_Locs(1); 
    T2_Loc = 0; T3_Loc = 0; T4_Loc = 0; %Fill unused spots with blanks    
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1_Loc,:), T1_Angle);
    ProbeLoc = T1_Loc; %One of memory target locations (1:8)
    T2_Angle = 9999; T3_Angle = 9999; T4_Angle = 9999; %Fill unused spots with blanks
    
elseif SetSize(Trial) == 2 %set size 2
    T_Locs = randsample(PossLocs,2);
    T1_Loc = T_Locs(1);
    T2_Loc = T_Locs(2);
    T3_Loc = 0; T4_Loc = 0; %Fill unused spots with blanks
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1_Loc,:), T1_Angle); 
    Screen('DrawTexture', w, MemLine, [], DispLocs(T2_Loc,:), T2_Angle); 
    ProbeLoc = randsample(T_Locs,1);
    T3_Angle = 9999; T4_Angle = 9999; %Fill unused spots with blanks
    
elseif SetSize(Trial) == 4 %set size 4
    T_Locs = randsample(PossLocs,4);
    T1_Loc = T_Locs(1);
    T2_Loc = T_Locs(2);
    T3_Loc = T_Locs(3);
    T4_Loc = T_Locs(4);
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1_Loc,:), T1_Angle); 
    Screen('DrawTexture', w, MemLine, [], DispLocs(T2_Loc,:), T2_Angle); 
    Screen('DrawTexture', w, MemLine, [], DispLocs(T3_Loc,:), T3_Angle); 
    Screen('DrawTexture', w, MemLine, [], DispLocs(T4_Loc,:), T4_Angle);
    ProbeLoc = randsample(T_Locs,1);
    
end

% present memory array with fixation
Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
Screen('Flip', w);
WaitSecs(PresDur); %Present memory array for specified presentation time


% --------------------- DELAY PERIOD  ---------------------------------


Screen(w, 'FillRect', screenColor);
Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
Screen('Flip', w);
WaitSecs(delay); %Present fixation for specified delay time


% --------------------- TEST PROBE  -----------------------------------


% Randomize theta angle & get start time
th = datasample(-90:89,1); % Initial rotation angle (degrees)
rot_spd = 2;

bdown = 0;
[~, dialStart] = PsychPowerMate('Get', handle); %Initiate dial position
Begin_Time = GetSecs; 

%Get Response (Rotation)
while(~any(bdown))    
    [mx,my,bdown] = GetMouse;
    
    [button, dialPos] = PsychPowerMate('Get', handle); 
    if dialPos > dialStart %Clockwise
        th = th + rot_spd; dialStart = dialPos;
    end
    if dialPos < dialStart %Counterclockwise
        th = th - rot_spd; dialStart = dialPos;
    end
    Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
    Screen('DrawTexture', w, MemLine, [], DispLocs(ProbeLoc,:),th);
    Screen('Flip', w);
    
    if button == 1
        theta = th; %Use for later
        End_Time = GetSecs;
        RT(Trial) = (End_Time-Begin_Time)*1000;
        break
    end   
    WaitSecs(0.0001);  

end

% Clear the screen to finish trial & present fixation for ITI
Screen(w, 'FillRect', screenColor);
    mm = 19;
    for i = 0:4
        nn = mm-i*4;
        Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
        Screen('Flip', w);
        WaitSecs(0.05);
    end
    Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
    Screen('Flip', w);
    WaitSecs(0.7);

end

%--------------------------------------------------------------------
%	Ending experiment
%--------------------------------------------------------------------
 
%matfile = [T1_Loc(Order)' ProbeAngle Theta RespAngle RespError RT'];
%matfname = strcat(ID,'_VWMdat.mat');
%save(matfname,'matfile');

%Move data file to data directory and create a .mat file with data

%movefile(FileName,strcat(pwd,dirout)); 
%movefile(matfname,strcat(pwd,dirout)); 

% Close Powermate
PsychPowerMate('Close', handle);

ListenChar(0);

Screen(w, 'FillRect', screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
EndText = ('Finished with the practice! Any questions?');
DrawFormattedText(w, EndText, 'center', 'center', [255 255 255]);
Screen('Flip', w);
WaitSecs(4);

% added by WP 05/10/16
screen_clut = [0:255; 0:255; 0:255]'./(255);
Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);

Screen('CloseAll');
Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);
sca;
ShowCursor;






