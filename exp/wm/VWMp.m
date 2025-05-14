%% VISUAL WORKING MEMORY PRECISION TASK (WITH ORIENTATION)
%
% A spatial delayed response WM task that requires subject to report
% orientation of probe 
%
% Megan Ichinose & Woon Park 
% April 20th, 2016
% 
% Updated May 27th, 2016 - Megan Ichinose
% Update includes USING GRIFFIN POWERMATE DIAL, not arrow keys
% Update also includes not sampling rotations of gratings within 8 degrees of each other 
% 
% =========================================================================
%               INPUT PARAMETERS                 
% =========================================================================

warning('off','MATLAB:dispatcher:InexactMatch'); clc; sca;
Screen('CloseAll');
Screen('Preference', 'SuppressAllWarnings', 1);
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
    set_size        = 3;       %Set sizes of 1,2,4
    num_set         = 150;     %Number of trials per set size  
    block_length    = 45;      %Equal number of trials in each set per block
    num_blocks      = 10;      %Number of blocks

    T_num = set_size*num_set;  %Total number of trials, also = block_length*num_blocks

    %Specify timing of stimuli
    PresDur         = 2;       %Stimulation duration (s)
    delay           = 1;       %Delay time (s)
     
    %Specify text 
    txtsz           = 25;      %Size of all presented text
    txtfont         = 'Arial'; %Text font
    
    %Specify hypotenuse for circle array (edited by WP 05/10/16)
    Hypot_deg       = 4;       %Degree (Distance between center/fixation & stimuli)
    Hypot           = round(Hypot_deg*60/scale_factor); %pixels (original: 150)
    
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
FileName =  strcat(ID, '_VWMp.csv');
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
fprintf(FID, 'Trial, Order, SetSize, T1_Loc, T2_Loc, T3_Loc, T4_Loc, T1_Angle, T2_Angle, T3_Angle, T4_Angle, ProbeLoc, ProbeAngle, Theta, RespAngle, RespError, RT\n');
fclose(FID);


% =========================================================================
%               OPEN & INITIALIZE SCREEN 
% =========================================================================


screenNumber = max(Screen('Screens')); %Set screen number

screenColor = 126; %Light gray (edited by WP 05/10/16; previously: [135 135 135])
[w, wRect] = Screen('OpenWindow', screenNumber, screenColor);

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
timing=Screen('GetFlipInterval',w);

rand('twister',sum(100*clock)); %Reset random generator
%HideCursor;


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
SetSize = repmat([1;2;4],num_set,1);
T1_Loc = zeros(T_num,1);
T2_Loc = zeros(T_num,1);
T3_Loc = zeros(T_num,1);
T4_Loc = zeros(T_num,1);
Probe_Loc = zeros(T_num,1);
Probe_Angle = repmat(9999,T_num,1);
Theta = repmat(9999,T_num,1);
Resp_Angle = repmat(9999,T_num,1);
Resp_Error = repmat(9999,T_num,1);
RT = zeros(T_num,1);

%Randomize set size trials within blocks

Order = reshape(Order,block_length,num_blocks);

for i=0 : num_blocks-1
    
    N = Order(:,i+1);
    Order(:,i+1) = N(randperm(end),:);
    
end

Order = reshape(Order,T_num,1); %Now trials are randomized within blocks


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

% % Create stimulus rectangles 
% movie_rect = [0,0,bps,bps];
% scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
% scr_top = fix(screen_rect(4)/2)-round(bps/2);
% screen_rect_middle = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
% screen_patch = screen_rect_middle;  %may change this if eccentricity is not 0;

Loc1 = PresLocs(1,:);
Loc2 = PresLocs(3,:);
Loc3 = PresLocs(5,:);
Loc4 = PresLocs(7,:);
Loc5 = PresLocs(9,:);
Loc6 = PresLocs(11,:);
Loc7 = PresLocs(13,:);
Loc8 = PresLocs(15,:);

DispLocs = [Loc1;Loc2;Loc3;Loc4;Loc5;Loc6;Loc7;Loc8];


% =========================================================================
%               BEGIN EXPERIMENT WITH INSTRUCTIONS
% =========================================================================

while KbCheck; end %Call KbCheck at the beginning since it's slow to load the 1st time
%ListenChar(2);

% Open Powermate Dial
handle = PsychPowerMate('Open');

% Draw Prompt 
Screen(w, 'FillRect', screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
Screen('DrawText', w, 'Press the spacebar to begin the experiment.', CX-200, CY, [255 255 255]);
Screen(w, 'Flip');
% GetChar;
KbWait(KeyboardIndex,2);


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
% from 0:180, with the parameter than no two items are within 8 degrees from each
% other. This section also selects the target to be probed after the delay,
% and randomly selects an orientation for that probe. 

PossLocs = 1:8; %Possible circle locations
PossAngles = -90:89; %Possible angle rotations (180)

%This resets possible angle rotations for each memory target with each
%trial (no angle within 8 degree buffer of another)
buff=8; %Set angle buffer
    T1_Angle = randsample(PossAngles,1); 
       PossAngles(find(PossAngles>T1_Angle-8 & PossAngles<T1_Angle+8))=[];
       if T1_Angle <= -90+buff 
           PossAngles(find(PossAngles>abs(T1_Angle)))=[];
       elseif T1_Angle >= 90-buff
           PossAngles(find(PossAngles<-T1_Angle))=[];
       end
    T2_Angle = randsample(PossAngles,1);
       PossAngles(find(PossAngles>T2_Angle-8 & PossAngles<T2_Angle+8))=[];
       if T2_Angle <= -90+buff 
           PossAngles(find(PossAngles>abs(T2_Angle)))=[];
       elseif T2_Angle >= 90-buff
           PossAngles(find(PossAngles<-T2_Angle))=[];
       end
    T3_Angle = randsample(PossAngles,1);
       PossAngles(find(PossAngles>T3_Angle-8 & PossAngles<T3_Angle+8))=[];
       if T3_Angle <= -90+buff 
           PossAngles(find(PossAngles>abs(T3_Angle)))=[];
       elseif T3_Angle >= 90-buff
           PossAngles(find(PossAngles<-T3_Angle))=[];
       end
    T4_Angle = randsample(PossAngles,1);
        
if SetSize(Order(Trial)) == 1 %set size 1
    T_Locs = randsample(PossLocs,1);
    T1_Loc = T_Locs(1); 
    T2_Loc = 0; T3_Loc = 0; T4_Loc = 0; %Fill unused spots with blanks    
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1_Loc,:), T1_Angle);
    ProbeLoc = T1_Loc; %One of memory target locations (1:8)
    T2_Angle = 9999; T3_Angle = 9999; T4_Angle = 9999; %Fill unused spots with blanks
    
elseif SetSize(Order(Trial)) == 2 %set size 2
    T_Locs = randsample(PossLocs,2);
    T1_Loc = T_Locs(1);
    T2_Loc = T_Locs(2);
    T3_Loc = 0; T4_Loc = 0; %Fill unused spots with blanks
    Screen('DrawTexture', w, MemLine, [], DispLocs(T1_Loc,:), T1_Angle); 
    Screen('DrawTexture', w, MemLine, [], DispLocs(T2_Loc,:), T2_Angle); 
    ProbeLoc = randsample(T_Locs,1);
    T3_Angle = 9999; T4_Angle = 9999; %Fill unused spots with blanks
    
elseif SetSize(Order(Trial)) == 4 %set size 4
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

% ListenChar(2)
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
    if dialPos > dialStart % Rotate clockwise
        th = th + rot_spd; dialStart = dialPos;
    end
    if dialPos < dialStart % Rotate counterclockwise
        th = th - rot_spd; dialStart = dialPos;
    end
    Screen('FrameOval', w,60,[sr_hor-3, sr_ver-3+V_ecc_fix, sr_hor+3, sr_ver+3+V_ecc_fix],2,2)
    Screen('DrawTexture', w, MemLine, [], DispLocs(ProbeLoc,:),th);
    Screen('Flip', w);
    
    if button == 1
        End_Time = GetSecs;
        RT(Trial) = (End_Time-Begin_Time)*1000;
        theta = th; %Use for later
        break
    end   
    WaitSecs(0.0001);  

end


%Change degree format for probe & response angles from [-90:90) to [-180:180)

%Compute which was the probed (correct) angle 
AngleOpts = [T1_Angle,T2_Angle,T3_Angle,T4_Angle];
ProbeAngle = AngleOpts(find(T_Locs==ProbeLoc));

ProbeAngle = wrapTo180(ProbeAngle); %This won't actually change it, but keep in for logical step
RespAngle = wrapTo180(theta); %This WILL change response (since theta can be >180 in neg and pos directions)

AltResp = 180-abs(RespAngle); %Find alternative ('flipped by 180') response
if abs(abs(ProbeAngle)-abs(RespAngle)) > abs(abs(ProbeAngle)-abs(AltResp)) %Choose the closest angle
    RespAngle = AltResp;
end

%Compute response error
RespError = abs(ProbeAngle)-abs(RespAngle);

%Save some trial data for later .mat file
Probe_Loc(Trial) = ProbeLoc;
Probe_Angle(Trial) = ProbeAngle;
Theta(Trial) = th;
Resp_Angle(Trial) = RespAngle;
Resp_Error(Trial) = RespError;

% ------------ FINISH UP TRIAL & GET DATA --------------------------------

% Save trial data
FID = fopen(FileName, 'a');
Save = [Trial; Order(Trial); SetSize(Order(Trial)); T1_Loc; T2_Loc; T3_Loc; T4_Loc; T1_Angle; T2_Angle; T3_Angle; T4_Angle; Probe_Loc(Trial); Probe_Angle(Trial); Theta(Trial); Resp_Angle(Trial); Resp_Error(Trial); RT(Trial)];
fprintf(FID, '%7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %7d, %3.4f, %3.4f\n', Save); 
fclose(FID);


% ---------------- ADD BREAKS FOR BLOCKS --------------------------------


if mod(Trial, block_length) == 0;
    if Trial ~= T_num
    Block = Block+1; %Increase block by 1
    FlushEvents('KeyDown');
    WaitSecs(0.5);
    Screen(w, 'FillRect', screenColor);
    Screen('TextSize', w, txtsz);
    Screen('TextFont', w, txtfont);
    BreakText = sprintf('You have %i more to go!\n\nPress the space bar when you are ready to continue', (T_num-(Trial)));
    DrawFormattedText(w, BreakText, 'center', 'center', [255 255 255]);
    Screen('Flip', w);
    KbWait(KeyboardIndex,2);
    end
end

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

% close screen (added by WP 05/10/16)
Screen('Close', MemLine);
clear grating;

end %Closes the for Trial=1:T_num loop, ending experiment

%--------------------------------------------------------------------
%	Ending experiment
%--------------------------------------------------------------------

% Close Powermate
PsychPowerMate('Close', handle);

matfile = [SetSize(Order) Probe_Loc Probe_Angle Theta Resp_Angle Resp_Error RT];
matfname = strcat(ID,'_VWMdat.mat');
save(matfname,'matfile');

%Move data file to data directory
movefile(FileName,strcat(pwd,dirout));
movefile(matfname,strcat(pwd,dirout));

ListenChar(0);

Screen(w, 'FillRect', screenColor);
Screen('TextSize', w, txtsz);
Screen('TextFont', w, txtfont);
EndText = ('You have completed the experiment! Please see experimenter');
DrawFormattedText(w, EndText, 'center', 'center', [255 255 255]);
Screen('Flip', w);
WaitSecs(5);

% added by WP 05/10/16
screen_clut = [0:255; 0:255; 0:255]'./(255);
Screen('LoadNormalizedGammaTable',screenNumber,screen_clut);

Screen('CloseAll');
sca;
ShowCursor;






