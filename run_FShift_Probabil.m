function [] = run_FShift_Probabil(sub,flag_training, flag_isolum, flag_block)
% run_FShift_PerIrr(sub,flag_training, flag_isolum, flag_block)
%   runs experiment SSVEP_FShift_Probabil
%       sub:            participant number
%       flag_training:  1 = do training
%       flag_isolum:    1 = do isoluminance adjustment
%       flag_block:     1 = start with block 1
%           e.g. run_FShift_Probabil(1,1, 0, 1)
% 
% current version includes two central relevant and one central irrelevant RDK with hue changes as target to be discriminated
% target within each trial


% Christopher Gundlach, Maria Dotzer,  Leipzig, 2023,2021, 2020

if nargin < 4
    help run_FShift_PerIrr
    return
end

%% parameters
% sub = 1; flag_training = 1; flag_block = 1; flag_isolum = 1;
% design
p.sub                   = sub;                  % subject number
p.flag_block            = flag_block;           % block number to start
p.flag_training         = flag_training;        % do training

p.ITI                   = [1000 1000];          % inter trial interval in ms
p.targ_respwin          = [200 1000];           % time window for responses in ms

% screen
p.scr_num               = 1;                    % screen number
p.scr_res               = [1920 1080];          % resolution
p.scr_refrate           = 480;                  % refresh rate in Hz (e.g. 85)
p.scr_color             = [0.05 0.05 0.05 1];      % default: [0.05 0.05 0.05 1]; ; color of screen [R G B Alpha]
p.scr_imgmultipl        = 4;

% some isoluminace parameters
p.isol.TrlAdj           = 5;                    % number of trials used for isoluminance adjustment
p.isol.MaxStd           = 10;                   % standard deviation tolerated
p.isol.run              = false;                % isoluminance run?
p.isol.override         = [];                   % manually set colors for RDK1 to RDKXs e.g. []
p.isol.override         = [0.0862745098039216 0.0345098039215686 0 1;0 0.0627450980392157 0.156862745098039 1;0 0.0823529411764706 0 1 ];
p.isol.override         = [0.0862745098039216 0.0345098039215686 0 1;0 0.0627450980392157 0.156862745098039 1;0 0.0823529411764706 0 1 ].*2; p.isol.override(p.isol.override>1)=1;
p.isol.override         = [0.0862745098039216 0.0345098039215686 0 1;0 0.0627450980392157 0.156862745098039 1;0.02 0.0823529411764706 0.02 1 ].*2; p.isol.override(p.isol.override>1)=1;
% % color 1 lum 68 chroma 70
% p.isol.override         = [0.151573078382719	0.492713810961937	0.169273961794666 1; ... % green --> increases are hard
%                             0.768808896052382	0.235396366987767	0.0676584621046665 1; ...% orange
%                             0.700395021369019	0.219797886668636	0.690819078346578 1 ]; % purple/magenta
% color 1 lum 68 chroma 60 [currently preferred #2]
p.isol.override         = [0.179103922145369	0.479134122907170	0.193900675407810 1; ... % green --> increases are hard
                            0.715591409278279	0.256239038997824	0.0950458443741471 1; ...% orange
                            0.648338268946474	0.245821973632775	0.639442314393116 1 ]; % purple/magenta
p.chroma_target = [-90 83];
p.isol.bckgr            = [0.38 0.38 0.38 1];          % isoluminant to background or different color?

% % color 2 lum 68 chroma 70
% p.isol.override         = [0.204088651658460	0.477390994465875	0.0695047784619124 1; ... % green
%                             0.768808896052382	0.235396366987767	0.0676584621046665 1; ...% orange
%                             0.700395021369019	0.219797886668636	0.690819078346578 1 ]; % purple/magenta
% % color 2 lum 68 chroma 60
% p.isol.override         = [0.243399878157271	0.456929342428427	0.119526503676649 1; ... % green
%                             0.768808896052382	0.235396366987767	0.0676584621046665 1; ...% orange
%                             0.700395021369019	0.219797886668636	0.690819078346578 1 ]; % purple/magenta


% % color 3 lum 68 chroma 70
% p.isol.override         = [0.110803064969462	0.505084699244808	0.240138958629960 1; ... % green
%                             0.768808896052382	0.235396366987767	0.0676584621046665 1; ...% orange
%                             0.700395021369019	0.219797886668636	0.690819078346578 1 ]; % purple/magenta

%%%%
% color version 2 green red/orange blue/purple lum 48 chroma 50 [currently preferred #1]
p.isol.override         = [0.0794654838295763	0.212979965009343	0.0694822333587645 1; ... % green
                            0.354923915689341	0.0960941504677729	0.0525769382198423 1; ...% orange/red
                            0.189288097847939	0.139241714541135	0.439109578562397 1 ]; % blue/magenta
p.chroma_target = [-90 100];
p.isol.bckgr            = [0.168 0.168 0.168 1];          % isoluminant to background or different color?
%%%%

% p.isol.bckgr            = [0.38 0.38 0.38 1];          % isoluminant to background or different color?
% p.isol.bckgr            = p.scr_color;          % isoluminant to background or different color?


% stimplan and design
p.stim.condition        = [1 2 3 4 5 6];    % attend RDK 1 valid, attend RDK 2 valid; attend RDK 1 invalid, attend RDK 2 invalid; neutral neutral
p.stim.RDK2attend       = [1 2 1 2 3 3];    % defines which RDK to attend in which condition
p.stim.targetRDK        = [1 2 2 1 1 2];    % defines which RDK will feature the target
p.stim.con_repeats_e    = [180 180 80 80 80 80];    % number of repetitions for each condition | experiment
                            % [240 240 100 100 100 100] would amount to 480/680 = 0.7059 per cent valid cues
                            % [180 180 80 80 80 80] would amount to 480/680 = 0.692 per cent valid cues
p.stim.con_repeats_t    = [5 5 2 2 2 2];    % number of repetitions for each condition | training
p.trial_irreg_prob      = [0.1];            % multiplier for trial being regular (events appear in p.targ_win_reg)

p.precue_time_n         = [1500 2000];          % pre cue time range of normal trials in ms
p.postcue_time          = [2750 2750];          % post cue time
p.ITI                   = [1000 1000];          % inter trial interval in ms

p.targ_duration         = [0.200];              % time for target luminance change in s
p.targ_win_reg          = [1500 2500];          % presentation time window for regular trials
p.targ_win_ctrl         = [300 1500];           % presentation time window for control trials trials
p.targ_respwin          = [200 1200];           % response window in ms

% (sum(p.stim.con_repeats_e)+sum(p.stim.con_repeats_e*p.trial_irreg_prob))*(mean(p.precue_time_n)+mean(p.postcue_time)+mean(p.ITI))/(1000*60)
% [240 240 100 100 100 100] = 88 min
% [180 180 80 80 80 80] = 69 min; 748 trials; 680 regular; 68 catch trials

p.stim.blocknum         = 17;                   % 44 trials per block and a length of around 4 minutes


% introduce RDK structure
RDK.RDK(1).size         = [360 360];                    % width and height of RDK in pixel; only even values [360 = 11.22°]
RDK.RDK(1).centershift  = [0 0];                        % position of RDK center; x and y deviation from center in pixel
RDK.RDK(1).col          = [1 0.4 0 1; p.scr_color(1:3) 0];% "on" and "off" color
RDK.RDK(1).freq         = 18;                            % flicker frequency, frequency of a full "on"-"off"-cycle
RDK.RDK(1).mov_freq     = 120;                          % Defines how frequently the dot position is updated; 0 will adjust the update-frequency to your flicker frequency (i.e. dot position will be updated with every "on"-and every "off"-frame); 120 will update the position for every frame for 120Hz or for every 1. quadrant for 480Hz 
RDK.RDK(1).num          = 120;                           % number of dots
RDK.RDK(1).mov_speed    = 1;                            % movement speed in pixel
RDK.RDK(1).mov_dir      = [0 1; 0 -1; -1 0; 1 0];       % movement direction  [0 1; 0 -1; -1 0; 1 0] = up, down, left, right
RDK.RDK(1).dot_size     = 10;                           % size of dots
RDK.RDK(1).shape        = 0;                            % 1 = square RDK; 0 = ellipse/circle RDK;
RDK.RDK(1).chromatarget = p.chroma_target;              % percent changes in chroma defined as events

RDK.RDK(2:3) = deal(RDK.RDK(1));
%[RDK.RDK(2:3).col] = deal([0 0.4 1 1; p.scr_color(1:3) 0], [0 1 0 1; p.scr_color(1:3) 0]);
[RDK.RDK.col] = deal([p.isol.override(1,:); p.scr_color(1:3) 0], [p.isol.override(2,:); p.scr_color(1:3) 0], [p.isol.override(3,:); p.scr_color(1:3) 0]);
[RDK.RDK(2:3).freq] = deal(21, 24);
[RDK.RDK(2:3).chromatarget] = deal(p.chroma_target, p.chroma_target);

%plot_colorwheel([1 0.4 0; 0 0.4 1; 0 1 0],'ColorSpace','propixxrgb','LAB_L',50,'NumSegments',60,'AlphaColWheel',1,'LumBackground',100)
 
RDK.event.type          = 'colorchange';       % event type color change
RDK.event.duration      = p.targ_duration;      % time of color onset
RDK.event.coherence     = 0.8;                    % ration of dots coherently changing color

% fixation cross
p.crs.color             = [0.8 0.8 0.8 1];      % color of fixation cross
p.crs.size              = 12;                   % size of fixation
p.crs.width             = 2;                    % width of fixation cross
p.crs.cutout            = 0;                    % 1 = no dots close to fixation cross

% trigger
p.trig.rec_start        = 253;                  % trigger to start recording
p.trig.rec_stop         = 254;                  % trigger to stop recording
p.trig.tr_start         = 77;                   % trial start; main experiment
p.trig.tr_stop          = 88;                   % trial end; main experiment
p.trig.tr_con_type      = [1 2 3 4 5 6]*10;     % condition type
p.trig.tr_timing_type   = [0 1];                % trial timing type: regular or catch trial
p.trig.event_side       = [100 200];            % left right
p.trig.event_cuevalid   = [10 10 20 20 30 30];  % valid invalid neutral
p.trig.event_type       = [1 2];                % hue increase hue decrease
p.trig.button           = 150;                   % button press


% possible condition triggers:
% {[1 101 201 111 121 211 221]; [2 102 202 112 122 212 222]; [3 103 203 113 123 213 223]; ...
% [4 104 204 114 124 214 224]; [5 105 205 115 125 215 225]; [6 106 206 116 126 216 226]}

% logfiles
p.log.path              = '/home/stimulationspc/matlab/User/christopher/stim_ssvep_fshift_probabil/logfiles';
p.log.exp_name          = 'SSVEP_FShift_Probabil';
p.log.add               = '_a';


%% check for logfile being present
filecheck=dir(sprintf('%sVP%02.0f_timing*',p.log.path,p.sub));
if ~isempty(filecheck)
    reply = input(sprintf('\nVP%02.0f existiert bereits. Datei überschreiben? [j/n]... ',p.sub),'s');
    if strcmp(reply,'j')
        p.filename = sprintf('VP%02.0f_timing',p.sub);
    else
        [temp name_ind]=max(cellfun(@(x) numel(x), {filecheck.name}));
        p.filename = sprintf('%s%s',filecheck(name_ind).name(1:end-4),p.log.add);
    end
else
    p.filename = sprintf('VP%02.0f_timing',p.sub);
end

t.isol = {};
% routine to check for older isoluminance adjustments
for i_file = 1:numel(filecheck)
    t.in = load(fullfile(filecheck(i_file).folder,filecheck(i_file).name));
    t.datenum{i_file} = filecheck(i_file).datenum;
    t.isol{i_file} = t.in.p.isol;
    
end



%% Screen init
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.size/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

%% keyboard and ports setup ???
% keyboard setup
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('UpArrow') KbName('DownArrow') KbName('1!') KbName('2@') KbName('3#')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.INCREASE, key.DECREASE] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5));

%% start experiment
% initialize randomization of stimulation frequencies and RDK colors
% inititalize RDKs [RDK1 and RDK2 task relevant at center;  RDK3 RDK4 RDK5 not and in periphery]
% rand('state',p.sub)
rng(p.sub,'v4')

[RDK.RDK.col_init] = deal(RDK.RDK.col);


% randomize frequencies
[RDK.RDK.freq] = deal(RDK.RDK(randperm(numel(RDK.RDK))).freq);

% randomize colors? yes
t.colridx = randperm(numel(RDK.RDK));
[RDK.RDK.col] = deal(RDK.RDK(t.colridx).col_init);
[RDK.RDK.chromatarget] = deal(RDK.RDK(t.colridx).chromatarget);

% initialize blank variables
timing = []; button_presses = []; resp = []; randmat = [];

%% initial training
% if p.flag_training
%     fprintf(1,'\nTraing starten mit q')
%     inp.prompt_check = 0;
%     while inp.prompt_check == 0             % loop to check for correct input
%         [key.keyisdown,key.secs,key.keycode] = KbCheck;
%         if key.keycode(key.SECRET)==1
%             flag_trainend = 0; inp.prompt_check = 1;
%         end
%         Screen('Flip', ps.window, 0);
%     end
%     
%     
%     i_bl = 1;
%     flag_trainend = 0;
%     while flag_trainend == 0 % do training until ended
%         rand('state',p.sub*i_bl) % determine randstate
%         randmat.training{i_bl} = rand_FShift_Probabil(p, RDK,  1);
%         [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
%             pres_FShift_Probabil(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
%         save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
%         pres_feedback(resp.training{i_bl},p,ps, key,RDK)
%                
%         % loop for training to be repeated
%         fprintf(1,'\nTraing wiederholen? (j/n)')
%         inp.prompt_check = 0;
%         while inp.prompt_check == 0             % loop to check for correct input
%             [key.keyisdown,key.secs,key.keycode] = KbCheck; 
%             if key.keycode(key.YES)==1
%                 i_bl = i_bl + 1; flag_trainend = 0; inp.prompt_check = 1;
%             elseif key.keycode(key.NO)==1
%                 flag_trainend = 1; inp.prompt_check = 1;
%             end
%             Screen('Flip', ps.window, 0);
%         end  
%         
%     end
% end

%% then isoluminance adjustment
% do the heterochromatic flicker photometry
if flag_isolum == 1
%     
%     PsychDefaultSetup(2);
%     Datapixx('Open');
%     Datapixx('SetPropixxDlpSequenceProgram', 0);
%     Datapixx('RegWrRd');
     
    p.isol.init_cols = cell2mat(cellfun(@(x) x(1,1:3)', {RDK.RDK.col}, 'UniformOutput', false))';
    
    % start isoluminance script only RGB output (no alpha)
    [Col2Use] = PRPX_IsolCol_480_Lab(...
        [p.isol.bckgr(1:3); p.isol.init_cols],...
        p.isol.TrlAdj,...
        p.isol.MaxStd,...
        0,...
        RDK.RDK(1).size.*2);
    
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = [Col2Use(1+i_RDK,:) 1];
    end
    % index function execution
    p.isol.run = sprintf('originally run: %s',datestr(now));
    p.isol.coladj = [Col2Use(2:end,:) ones(size(Col2Use,1)-1,1)];
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nadjusted colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
    
%     Screen('CloseAll')
%     Datapixx('SetPropixxDlpSequenceProgram', 0);
%     Datapixx('RegWrRd');
%     Datapixx('close');
else
    % select colors differently
    fprintf(1,'\nKeine Isoluminanzeinstellung. Wie soll verfahren werden?')
    % specify options
    % option1: use default values
    isol.opt(1).available = true;
    t.cols = cell2mat({RDK.RDK(:).col}');
    isol.opt(1).colors = t.cols(1:2:end,:);
    isol.opt(1).text = sprintf('default: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(1).colors(:,1:3)'));
    % option2: use isoluminance values of previously saved dataset
    if ~isempty(t.isol) % file loaded?
        [t.t t.idx] = max(cell2mat(t.datenum));
        isol.opt(2).available = true;
        isol.opt(2).colors = t.isol{t.idx}.coladj(1:end,:);
        isol.opt(2).text = sprintf('aus gespeicherter Datei: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(2).colors(:,1:3)'));
    else
        isol.opt(2).available = false;
        isol.opt(2).colors = [];
        isol.opt(2).text = [];
    end
    % option3: use manual override
    if ~isempty(p.isol.override)
        isol.opt(3).available = true;
        isol.opt(3).colors = p.isol.override(t.colridx,:);
        isol.opt(3).text = sprintf('manuell definiert in p.isol override: %s',sprintf('[%1.2f %1.2f %1.2f] ',isol.opt(3).colors(:,1:3)'));
    else
        isol.opt(3).available = false;
        isol.opt(3).colors = [];
        isol.opt(3).text = [];
    end
    % check for buttons
    IsoButtons = Buttons(6:8);
    isol.prompt.idx = find([isol.opt(:).available]);
    t.prompt = [];
    for i_prompt = 1:numel(isol.prompt.idx)
        t.prompt = [t.prompt sprintf('\n(%1.0f) %s',i_prompt,isol.opt(isol.prompt.idx(i_prompt)).text)];
    end
    
    % display options
    fprintf('%s',t.prompt)
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if any(key.keycode)
            inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    Col2Use = isol.opt(isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)).colors;
    % for troubleshooting/testing
    % Col2Use =p.isol.override(t.colridx,:);


    % use selected colors
    for i_RDK = 1:numel(RDK.RDK)
        RDK.RDK(i_RDK).col(1,:) = Col2Use(i_RDK,:);
    end
    % index function execution
    switch isol.prompt.idx(key.keycode(IsoButtons(1:numel(isol.prompt.idx)))==1)
        case 1
            p.isol.run = sprintf('default at %s',datestr(now));
        case 2
            p.isol.run = sprintf('reloaded at %s from %s',datestr(now),datestr(t.datenum{t.idx}));
        case 3
            p.isol.run = sprintf('override at %s',datestr(now));
    end
    p.isol.coladj = Col2Use;
%     save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    
    fprintf('\nselected colors:\n')
    for i_col = 1:size(p.isol.coladj,1)
        fprintf('RDK%1.0f [%1.4f %1.4f %1.4f]\n', i_col,p.isol.coladj(i_col,:))
    end
end

%% redo initialization
ps.input = struct('ScrNum',p.scr_num,'RefRate',p.scr_refrate,'PRPXres',p.scr_res,'BckGrCol',p.scr_color,'PRPXmode',2);
[~, ps.screensize, ps.xCenter, ps.yCenter, ps.window, ps.framerate, ps.RespDev, ps.keymap] = PTExpInit_GLSL(ps.input,1);

% some initial calculations
% fixation cross
ps.center = [ps.xCenter ps.yCenter];
p.crs.half = p.crs.size/2;
p.crs.bars = [-p.crs.half p.crs.half 0 0; 0 0 -p.crs.half p.crs.half];

% shift into 4 quadrants (running with 480 Hz)
ps.shift = [-ps.xCenter/2, -ps.yCenter/2; ps.xCenter/2, -ps.yCenter/2;... % shifts to four quadrants: upper left, upper right, lower left, lower right
    -ps.xCenter/2, ps.yCenter/2; ps.xCenter/2, ps.yCenter/2];

p.crs.lines = [];
for i_quad=1:p.scr_imgmultipl
    p.crs.lines = cat(2, p.crs.lines, [p.crs.bars(1,:)+ps.shift(i_quad,1); p.crs.bars(2,:)+ps.shift(i_quad,2)]); %array with start and end points for the fixation cross lines, for all four quadrants
end

% keyboard setup
KbName('UnifyKeyNames')
Buttons = [KbName('ESCAPE') KbName('Q') KbName('SPACE') KbName('UpArrow') KbName('DownArrow') KbName('j') KbName('n')];
RestrictKeysForKbCheck(Buttons);
key.keymap=false(1,256);
key.keymap(Buttons) = true;
key.keymap_ind = find(key.keymap);
[key.ESC, key.SECRET, key.SPACE, key.INCREASE, key.DECREASE, key.YES, key.NO] = deal(...
    Buttons(1),Buttons(2),Buttons(3),Buttons(4),Buttons(5), Buttons(6), Buttons(7));


%% do training again?
% loop for training to be repeated
fprintf(1,'\nTraing starten (j/n)')
inp.prompt_check = 0;
while inp.prompt_check == 0             % loop to check for correct input
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    if key.keycode(key.YES)==1
        flag_trainend = 0; inp.prompt_check = 1;
    elseif key.keycode(key.NO)==1
        flag_trainend = 1; inp.prompt_check = 1;
    end
    Screen('Flip', ps.window, 0);
end

if ~exist('i_bl'); i_bl = 1; end
while flag_trainend == 0 % do training until ended
    rng(p.sub*i_bl,'v4') % determine randstate
    randmat.training{i_bl} = rand_FShift_Probabil(p, RDK,  1);
    [timing.training{i_bl},button_presses.training{i_bl},resp.training{i_bl}] = ...
        pres_FShift_Probabil(p, ps, key, RDK, randmat.training{i_bl}, i_bl,1);
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
    pres_feedback(resp.training{i_bl},p,ps, key,RDK)
    
    % loop for training to be repeated
    fprintf(1,'\nTraing wiederholen? (j/n)')
    inp.prompt_check = 0;
    while inp.prompt_check == 0             % loop to check for correct input
        [key.keyisdown,key.secs,key.keycode] = KbCheck;
        if key.keycode(key.YES)==1
            i_bl = i_bl + 1; flag_trainend = 0; inp.prompt_check = 1;
        elseif key.keycode(key.NO)==1
            flag_trainend = 1; inp.prompt_check = 1;
        end
        Screen('Flip', ps.window, 0);
    end
    
end


%% present each block
% randomization
% rand('state',p.sub);                         % determine randstate
rng(p.sub,'v4')
randmat.experiment = rand_FShift_Probabil(p, RDK,  0);    % randomization
for i_bl = p.flag_block:p.stim.blocknum
    % start experiment
    [timing.experiment{i_bl},button_presses.experiment{i_bl},resp.experiment{i_bl}] = ...
        pres_FShift_Probabil(p, ps, key, RDK, randmat.experiment, i_bl,0);
    % save logfiles
    save(sprintf('%s%s',p.log.path,p.filename),'timing','button_presses','resp','randmat','p', 'RDK')
          
    pres_feedback(resp.experiment{i_bl},p,ps, key, RDK)    
end

fprintf(1,'\n\nENDE\n')

%Close everything
Datapixx('SetPropixxDlpSequenceProgram', 0);
Datapixx('RegWrRd');
Datapixx('close');
ppdev_mex('Close', 1);
ListenChar(0);
sca;


end

