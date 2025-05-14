% This code measures Threshold versus Contrast (TvC) using FAST
% 2 AFC; gabor stimuli are sandwiched by two noise frames
% Participants judge whether the gabor is clockwise or counterclockwise
% tilted

% updated by WP 04/19/2016

clc;
clear all;
close all;

% change this before each experiment --------
subject_initials      = 'test';
% -------------------------------------------


ListenChar(2);
warning('off','MATLAB:dispatcher:InexactMatch')
KbName('UnifyKeyNames');
keys=KbName({'LeftArrow' 'RightArrow'}); % for response
space=KbName('space');


% KEY Monitor Parameters
scale_factor          = 1.7;        % most important parameter - how many acrmin is one screen pixel? for SONY monitor at [1024 640] resolution, and 30.4in viewing distance, use scale_factor=2
frame_rate            = 120;        % 
resolution            = [1024 768]; % [1024 640];
linearize 	          = 1;          % whether monitor is linearized. if=1, program will look for "MyGammaTable"
break_if_mismatch     = 1;          % if 1, this will quit the program if monitor resolution does not match above settings (use 1 for the ACTUAL EXP, 0 for testing on a different monitor, eg laptop)

% trials 
trialsPerCondition    = 40; % 40 for Vandy; 50 gives 800 total; 30 gives 480 total;
trialsperblock        = 80; % added by WP 09/04/2014

% break duration (added by WP 09/04/2014)
breakdur              = 30; 

% control monitor black out (added by WP 09/04/2014)
blackout              = 0; 

% keyboard settings (added by WP 09/04/2014)
keysetting = 1; % 1: home&exp / 2: office
keyindice = GetKeyboardIndices;
keyboardnum = keyindice(1);
if keysetting == 1
    keyboardnum = keyindice(1);
elseif keysetting == 2
    keyboardnum = keyindice(2);
end
ApplyKbFilter;

% check filename availability
data_file_name = strcat(subject_initials,'|PTM|data.mat');
IsExist = exist(data_file_name,'file');
if IsExist ~= 0
    error('data file name exists')
    ListenChar(1);
end

% create a data folder if it doesn't exist (added by WP 09/04/2014)
dirname = strcat('@data');
if isdir(dirname) == 0 
    mkdir(dirname);
end


% stimulus and time parameters
spatial_envelope      = 2;      % 0 = disk, 1 = Gabor patch, 2 = raised cosine envelope
background            = 126;    % background intensity, in gray scale untis
SF                    = 1;      % cycle/deg, 
radius                = 1;      % deg (3)
eccentricity          = 0;      % degree from target center to fixation
rAngle                = 0;      % reference angle, 0
dAngle                = 45;     % angle for rotation
sigma                 = 0.75;   % standard deviation for gaussian envelope
noise_factor          = 4;      % noise pixel size
dur                   = 0.0167; % sec, gabor and each noise frame display time
isi                   = 0.3;   % sec, interval between response and next trial

% if you change noise level, make sure to change it accordingly in
% Estimate_TvC.m and computeOverallLieli.m
% Ne =[0 0.003 0.006 0.0123 0.0249 0.043 0.0744 0.11];  
% Ne =[0 0.003 0.006 0.0123 0.0249 0.060 0.14 0.21];
% Ne = [0 0.003 0.0061 0.0124 0.0251 0.051 0.103 0.21]; % ROC ASD 
% Ne = [0 0.021 0.041 0.083 0.124 0.165 0.248 0.33]; % RUYUAN's PNAS PAPER
% Ne = [0 0.01 0.018 0.032 0.057 0.1 0.184 0.33]; % undergrad pilot
Ne = [0 0.01 0.0166 0.0276 0.0458 0.0761 0.1264 0.21];

stimulus_radius = round(60*radius/scale_factor);
f = (SF*scale_factor/60)*2*pi;
rAngle = rAngle*pi/180;
a = cos(rAngle)*f; 
b = sin(rAngle)*f;
white = background*2;
dur_frame  = round(frame_rate*dur);

V_ecc_fix = 0;
H_ecc_fix = 0;
V_ecc_fix = V_ecc_fix*60/scale_factor;
H_ecc_fix = H_ecc_fix*60/scale_factor;

% FAST structure parameters

upper_bound           = 1; % upper contrast value
lower_bound           = 0.001; % This is very important since psychometric function usually defined as (0,inf)
n_alternatives        = 2; % 2 or 4.;%2,2 AFC; 4, 4 AFC; It matters just when we update fast
pthreshold            = [0.7071 0.7937];%two performance level
total_trials          = trialsPerCondition*length(Ne)*length(pthreshold);

fastSS_79 = fastFull( n_alternatives,...
        'func_TvC_79',...
        'psyWeibull',...
        { [0.5, 2], [0.1, 2], [0.000001, 0.5], [0.0001, 0.5],[0.0001, 20]} ); %fast structure for mask
fastSS_70 = fastFull( n_alternatives,...
        'func_TvC_70',...
        'psyWeibull',...
        { [0.5, 2], [0.1, 2], [0.000001, 0.5], [0.0001, 0.5],[0.0001, 20]} ); %fast structure for mask


% data structure
data = zeros(total_trials, 8); % trial number,accuracy,,rotation,noise level,contrast,rs_key,correct,RT
data(:,1) = (1:total_trials)';     
data(:,2) = [ones(total_trials/2,1);ones(total_trials/2,1)*2]; % 1,70%;2,79%
[data(:,2) index] = Shuffle(data(:,2));
data(:,3) = rem(randperm(total_trials),2)*2-1; % -1,counter-clockwise rotated;1 clockwise rotated 
noise_list = repmat((1:length(Ne))',total_trials/length(Ne),1); % 8 noise level
noise_list = noise_list(index);

l = 7; % fixation cross line length 
ww = 4; % fixation cross line width

%---------------------------------------

try
    % open Screen windows
     Screen('Preference', 'SkipSyncTests', 1);
    oldVisualDebugLevel = Screen('Preference', 'VisualDebugLevel', 3);
    oldSupressAllWarnings = Screen('Preference', 'SuppressAllWarnings', 1);
    screens=Screen('Screens');
    screenNumber=max(screens);
%     screenNumber = 1;
    w=Screen('OpenWindow',screenNumber,0,[],[],2);
    
    screen_rect = Screen('Rect',w);
    if linearize
        load('mygammatable.mat');
        Screen('LoadNormalizedGammaTable', screenNumber, mygammatable);
    end
    Screen('FillRect',w, background);
    Screen('Flip', w);
    Screen('FillRect',w, background);
    Screen('TextSize',w,30);

    sr_hor = round(screen_rect(3)/2); sr_ver = round(screen_rect(4)/2);
    
    % black out the control screen (added by WP 09/02/2014)
    if blackout 
        if size(screens,2) == 2
            w2=Screen('OpenWindow',0,0,[],[],2);
            Screen('FillRect',w2, 0); 
            Screen('Flip', w2);
        end
    end
    
    
    % make spatial envelope
    [x,y] = meshgrid(-stimulus_radius:stimulus_radius,-stimulus_radius:stimulus_radius);
    bps = (stimulus_radius)*2+1;
    circle = ((stimulus_radius)^2-(x.^2+y.^2));
    
    for i=1:bps;
        for j =1:bps;
            if circle(i,j) < 0; circle(i,j) = 0;
            else
                circle(i,j) = 1;
            end;
        end;
    end;
    noise_circle=circle;
    if spatial_envelope == 1
        circle = (exp(-(((x)/(sqrt(2)*Gaussian_stdev/2)).^2)-((y/(sqrt(2)*Gaussian_stdev/2)).^2)).*circle);
    elseif spatial_envelope == 2
        R = (sqrt(x.^2 + y.^2) + eps).*circle;R = R/max(max(R));
        cos2D = (cos(R*pi)+1)/2;circle = (cos2D.*circle);
    end
    
    % create stimulus rectangles
    movie_rect = [0,0,bps,bps];
    scr_left_middle = fix(screen_rect(3)/2)-round(bps/2);
    scr_top = fix(screen_rect(4)/2)-round(bps/2);
    screen_rect_middle = movie_rect + [scr_left_middle, scr_top, scr_left_middle, scr_top];
    screen_patch = screen_rect_middle;  %may change this if eccentricy is not 0;
    
    
    
    % set up welcome interface
    tic; 

    Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);

    HideCursor;
    Screen('DrawText',w,[int2str(total_trials),'  trials. Press 0 key to start the experiment'],sr_hor-300,sr_ver-80,0);Screen('Flip',w);
    
    KbWait(keyboardnum,3);
    
    Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
    Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
    Screen('Flip', w);
    FlushEvents('keyDown');
    %
  
    trial = 1;
    while trial <= total_trials
        
 
        t1 = GetSecs;
        
        current_noise = Ne(noise_list(trial)); % choose a noise level
        
        if data(trial,2)==1
            current_con = fastChooseYp( fastSS_70, current_noise,pthreshold(1));
%             current_con = fastChooseY(fastSS_70, current_noise);
        elseif data(trial,2)==2
            current_con = fastChooseYp( fastSS_79, current_noise,pthreshold(2));
%             current_con = fastChooseY(fastSS_79, current_noise);
        end

        if current_con > upper_bound  % 100 for contrast upper bound
            current_con = upper_bound;
            
        end
        if current_con < lower_bound  % 0.1 for contrast lower bound
            current_con = lower_bound;
            
        end

        
        amplitude = background*current_con;
  

        % make stimulus
        % first grating
        grating = round(((sin(a*x+b*y+rand*2*pi).*circle*amplitude)+background));
        frame = zeros(bps,bps,3);
        for i = 1:ceil(dur_frame/3)
            if frame_rate == 360
                for j=1:3
                    if ((i-1)*3+j)>dur_frame
                        %frame(:,:,j) = ones(bps)*background;
                        switch j
                            case 1
                                frame(:,:,3) = ones(bps)*background;
                            case 2
                                frame(:,:,1) = ones(bps)*background;
                            case 3
                                frame(:,:,2) = ones(bps)*background;
                        end
                    else
                        switch j
                            case 1
                                frame(:,:,3) = grating;
                            case 2
                                frame(:,:,1) = grating;
                            case 3
                                frame(:,:,2) = grating;
                        end
                    end
                end
                movie_signal{i} = Screen('MakeTexture',w,frame);
            else
               movie_signal{i} = Screen('MakeTexture',w,grating); 
            end
        end       
        % create noise
        for j=1:2
            NN = round(bps/noise_factor)+1; %# of pixel
            img = current_noise*randn(NN);
            while 1
                out = abs(img(:))>1; nout=sum(out);
                if nout == 0, break; end
                img(out) = current_noise*randn(nout,1);
            end
            img = Expand(img,noise_factor);
            img = img(1:size(noise_circle,1),1:size(noise_circle,2));%match size of img and circle mask
            img = img.*noise_circle*background+background;
            noiseImg{j} = img;
        end
        
        % make noise movie
        frame = zeros(bps,bps,3);
        for i = 1:ceil(dur_frame/3)
            if frame_rate == 360
                for j=1:3
                    if ((i-1)*3+j)>dur_frame
                        %frame(:,:,j) = ones(bps)*background;
                        switch j
                            case 1
                                frame(:,:,3) = ones(bps)*background;
                            case 2
                                frame(:,:,1) = ones(bps)*background;
                            case 3
                                frame(:,:,2) = ones(bps)*background;
                        end
                    else
                        switch j
                            case 1
                                frame(:,:,3) = noiseImg{1};
                            case 2
                                frame(:,:,1) = noiseImg{1};
                            case 3
                                frame(:,:,2) = noiseImg{1};
                        end
                    end
                end
                movie_noise1{i} = Screen('MakeTexture',w,frame);
            else
                movie_noise1{i} = Screen('MakeTexture',w,noiseImg{1});
            end
        end
        frame = zeros(bps,bps,3);
        for i = 1:ceil(dur_frame/3)
            if frame_rate == 360
                for j=1:3
                    if ((i-1)*3+j)>dur_frame
                        %frame(:,:,j) = ones(bps)*background;
                        switch j
                            case 1
                                frame(:,:,3) = ones(bps)*background;
                            case 2
                                frame(:,:,1) = ones(bps)*background;
                            case 3
                                frame(:,:,2) = ones(bps)*background;
                        end
                    else
                        switch j
                            case 1
                                frame(:,:,3) = noiseImg{2};
                            case 2
                                frame(:,:,1) = noiseImg{2};
                            case 3
                                frame(:,:,2) = noiseImg{2};
                        end
                    end
                end
                movie_noise2{i} = Screen('MakeTexture',w,frame);
            else
               movie_noise2{i} = Screen('MakeTexture',w,noiseImg{2}); 
            end
        end
        
        %check noise image by saving a noise image
%         imwrite(noiseImg{1}./255,'noiseImg1.png','png');
%         imwrite(noiseImg{2}./255,'noiseImg2.png','png');
%         imwrite(grating./255,'grating.png','png');

        
        %  initiate trial
        t2 = GetSecs - t1;
        
        FlushEvents('keyDown');
        WaitSecs(1-t2);

        Screen('FillRect',w, background);
        Screen('Flip', w);
        mm = 19;
        for i = 0:4
            nn = mm-i*4;
            Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
            Screen('Flip', w);
            WaitSecs(0.05);
        end
        Screen('FrameOval', w,60,[sr_hor-nn, sr_ver-nn+V_ecc_fix, sr_hor+nn, sr_ver+nn+V_ecc_fix],2,2)
        Screen('Flip', w);
        WaitSecs(0.36);
        Screen('FillRect',w, background);
        Screen('Flip', w);
        WaitSecs(0.3);
        

        % play the movie
        priorityLevel=MaxPriority(w);Priority(priorityLevel);
        blah=GetSecs;
        Screen('Flip',w);
        
        %show first noise
        for i = 1:ceil(dur_frame/3);
            Screen('DrawTexture', w, movie_noise1{i}, movie_rect, screen_patch);
            Screen('Flip',w);
        end;
        
        %show signal
        for i = 1:ceil(dur_frame/3);
            Screen('DrawTexture', w, movie_signal{i}, movie_rect, screen_patch, dAngle*data(trial,3));
            Screen('Flip',w);
        end
        
        %show second noise
        for i = 1:ceil(dur_frame/3);
            Screen('DrawTexture', w, movie_noise2{i}, movie_rect, screen_patch);
            Screen('Flip',w);
        end;
        
        Screen('FillRect',w, background);
        Screen('Flip',w);
        ttt1 = GetSecs;

        Priority(0);
        
        % get the response && Update QUEST
        answer = data(trial,3);
        [rs, endSecs, keyCode] = get_response(answer,keyboardnum);
        
        Screen('FillRect',w, background);
        vbl=Screen('Flip', w);
        ampS = 1/4;
        
        if rs == 1
            Snd('Play',sin((1:1000)/4)*ampS);
            Snd('Wait');
        end

        
        %record data and update FAST
        rs_key = find(keyCode == 1);
        % Avoid dimension mismatch error (pressing
        % 2 keys at once) - ADDED BY MCI 7/14/16
        rs_key = rs_key(1);
        
        data(trial,4:7) = [current_noise,current_con, rs_key,rs];
        data(trial,8) = endSecs - ttt1;
        
        if data(trial,2)==1
            [fastSS_70, resample] = fastUpdate( fastSS_70, [current_noise, current_con, rs] );
        elseif data(trial,2)==2
            [fastSS_79, resample] = fastUpdate( fastSS_79, [current_noise, current_con, rs] );
        end

        FlushEvents('keyDown');
        
        % close movie
        for i = 1:ceil(dur_frame/3);
            % if time_gauss(i) ~= amplitude;
            Screen(movie_signal{i}, 'Close');
            Screen(movie_noise1{i}, 'Close');
            Screen(movie_noise2{i}, 'Close');
            % end
        end
        clear noiseImg img
        
        % set a break
        if rem(trial,trialsperblock)==0 && trial ~=total_trials %% modified by WP 09/04/2014

            %break countdown (added by WP 09/04/2014)
            JustBreak(w,breakdur,keyboardnum);
            Snd('Play',sin((1:1000)/4));
            WaitSecs(.5)
            Snd('Play',sin((1:1000)/4));
            WaitSecs(.5)
            Snd('Play',sin((1:1000)/4));
           
            Screen(w,'DrawLine',0,sr_hor-l+H_ecc_fix,sr_ver+V_ecc_fix,sr_hor+l+H_ecc_fix,sr_ver+V_ecc_fix,ww);
            Screen(w,'DrawLine',0,sr_hor+H_ecc_fix,sr_ver-l+V_ecc_fix,sr_hor+H_ecc_fix,sr_ver+l+V_ecc_fix,ww);
            Screen('DrawText',w,[int2str(trial),'  trials. Press 0 key to continue'],sr_hor-300,sr_ver-80,0);Screen('Flip',w);
            KbWait(keyboardnum,3);

        end
        
    cd(dirname);
        save(data_file_name);
    cd ..

    trial = trial+1;
    end
    
    ShowCursor;
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

time = toc/60

