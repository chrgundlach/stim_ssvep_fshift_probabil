function [timing,key,resp] = pres_FShift_Probabil(p, ps, key, RDK, conmat, blocknum, flag_training)
% presents experiment SSVEP_FShiftBase
%   p               = parameters
%   ps              = screen parameters
%   RDK             = RDK parameters
%   blocknum        = number of block
%   flag_training   = flag for trainig (1) or experiment (0)

%% adaptations for training
if flag_training == 1
    blocknum_present = blocknum;
    blocknum = 1;
end

%% prepare datapixx output
bufferAddress = 8e6;
% % samplesPerTrigger = 1;
triggersPerRefresh = 1; % send one per refresh (i.e. p.scr_refrate at 480 Hz)


%% initialize some RDK settings
% define input for RDK init function
RDKin.scr = ps; RDKin.scr.refrate = p.scr_refrate;
RDKin.Propixx = p.scr_imgmultipl;
RDKin.RDK = RDK;
RDKin.crs = p.crs;

%% loop for each trial
trialindex = find(conmat.mats.block==blocknum);

% send start trigger
if flag_training~=1
    PRPX_sendRecTrigger('start')
end

% Wait for release of all keys on keyboard, then sync us to retrace:
KbWait(ps.RespDev,1); % wait for releasing keys on indicated response device

% create keayboard queue
KbQueueCreate(ps.RespDev)

if flag_training~=1
    fprintf(1,'\nexperiment block %1.0f - Praesentation - Trial:', blocknum)
else
    fprintf(1,'\ntraining block %1.0f - Praesentation - Trial:', blocknum_present)
end
ttt=WaitSecs(0.3);

% loop across trials
for i_tr = 1:numel(trialindex)
    fprintf('%1.0f',mod(i_tr,10))
    inittime=GetSecs;
    %% initialize trial structure, RDK, cross, logs
    RDKin.trial = struct( ...
        'duration',conmat.trials(trialindex(i_tr)).post_cue_times+conmat.trials(trialindex(i_tr)).pre_cue_times,...
        'frames',conmat.trials(trialindex(i_tr)).post_cue_frames+conmat.trials(trialindex(i_tr)).pre_cue_frames,...
        'cue',conmat.trials(trialindex(i_tr)).pre_cue_frames+1);
    RDKin.trial.event = struct( ...
        'onset',conmat.trials(trialindex(i_tr)).event_onset_frames,...
        'color',conmat.trials(trialindex(i_tr)).event_col, ...
        'RDK',conmat.trials(trialindex(i_tr)).event_pos);
    RDKin.RDK.RDK = RDK.RDK;
    [colmat,dotmat,dotsize,rdkidx,frames] = RDK_init_FShift_Probabil(RDKin.scr,RDKin.Propixx,RDKin.RDK,RDKin.trial,RDKin.crs);
    
  
    % initialize fixation cross
    colmat_cr = repmat(p.crs.color' ,[1 size(p.crs.lines,2), size(colmat,3)]);
    if conmat.trials(trialindex(i_tr)).cue_direction ==3
        % cue both colors
        t.idx = randsample(2,2);
        colmat_cr(:,[1:4:end 2:4:end],(conmat.trials(trialindex(i_tr)).pre_cue_frames/p.scr_imgmultipl)+1:end)=...
            repmat(RDK.RDK(t.idx(1)).col(1,:)', ...
            [1,size(p.crs.lines,2)/2,conmat.trials(trialindex(i_tr)).post_cue_frames/p.scr_imgmultipl]);
        colmat_cr(:,[3:4:end 4:4:end],(conmat.trials(trialindex(i_tr)).pre_cue_frames/p.scr_imgmultipl)+1:end)=...
            repmat(RDK.RDK(t.idx(2)).col(1,:)', ...
            [1,size(p.crs.lines,2)/2,conmat.trials(trialindex(i_tr)).post_cue_frames/p.scr_imgmultipl]);
    else
        % from cue onwards: color of to be attended RDK
        colmat_cr(:,:,(conmat.trials(trialindex(i_tr)).pre_cue_frames/p.scr_imgmultipl)+1:end)=...
            repmat(RDK.RDK(conmat.trials(trialindex(i_tr)).cue_direction).col(1,:)', ...
            [1,size(p.crs.lines,2),conmat.trials(trialindex(i_tr)).post_cue_frames/p.scr_imgmultipl]);
    end
    
    % preallocate timing
    timing(i_tr) = struct('VBLTimestamp',NaN(1,frames.flips),'StimulusOnsetTime',NaN(1,frames.flips),...
        'FlipTimestamp',NaN(1,frames.flips),'Missed',NaN(1,frames.flips));
    
    %% set up responses and bookkeeping
    %setup key presses
    key.presses{i_tr}=nan(size(colmat,3),sum(key.keymap));
    key.presses_t{i_tr}=nan(size(colmat,3),sum(key.keymap));
    
    resp(i_tr).trialnumber              = trialindex(i_tr);
    resp(i_tr).blocknumber              = conmat.trials(trialindex(i_tr)).blocknum;
    resp(i_tr).condition                = conmat.trials(trialindex(i_tr)).condition;
    resp(i_tr).trial_timing_type        = conmat.trials(trialindex(i_tr)).trial_timing_type;
    resp(i_tr).cue_validity             = conmat.trials(trialindex(i_tr)).cue_validity; 
    resp(i_tr).cue_validity_label       = conmat.trials(trialindex(i_tr)).cue_validity_label; 
    resp(i_tr).cue                      = conmat.trials(trialindex(i_tr)).cue_direction; % attended RDK
    resp(i_tr).color_attended           = RDK.RDK(resp(i_tr).cue).col(1,:);
    resp(i_tr).freq_attended            = RDK.RDK(resp(i_tr).cue).freq;
    resp(i_tr).cue_onset_fr             = conmat.trials(trialindex(i_tr)).pre_cue_frames + 1;
    resp(i_tr).cue_onset_t_est          = (conmat.trials(trialindex(i_tr)).pre_cue_frames + 1)/p.scr_refrate*1000;
    resp(i_tr).cue_onset_t_meas         = nan; % measured onset time for cue
    resp(i_tr).pre_cue_frames           = conmat.trials(trialindex(i_tr)).pre_cue_frames;
    resp(i_tr).pre_cue_times            = conmat.trials(trialindex(i_tr)).pre_cue_times;
    resp(i_tr).post_cue_times           = conmat.trials(trialindex(i_tr)).post_cue_frames;
    resp(i_tr).post_cue_frames          = conmat.trials(trialindex(i_tr)).post_cue_times;
    resp(i_tr).eventpos                 = conmat.trials(trialindex(i_tr)).event_pos; % target at which RDK?
    resp(i_tr).event_type               = conmat.trials(trialindex(i_tr)).event_type; % 1 = targethue decrease; 2 = increase
    resp(i_tr).eventRDK_color           = cell2mat(cellfun(@(x) x(1,:),{RDK.RDK(resp(i_tr).eventpos(resp(i_tr).eventpos>0)).col},...
        'UniformOutput',false)');
    resp(i_tr).event_color              = conmat.trials(trialindex(i_tr)).event_col(1,:); 
    resp(i_tr).eventRDK_freq            = [RDK.RDK(resp(i_tr).eventpos(resp(i_tr).eventpos>0)).freq]';
    resp(i_tr).event_onset_frames       = conmat.trials(trialindex(i_tr)).event_onset_frames;
    resp(i_tr).event_onset_times        = conmat.trials(trialindex(i_tr)).event_onset_times;
    
    %% set up datapixx trigger vector
    % prepare datapixx scheduler
    % trigger signal needs to encompass first frame after stimulation has ended to send trial stop trigger
    if size(colmat_cr,3)*p.scr_imgmultipl>...
            conmat.trials(trialindex(i_tr)).pre_cue_frames+conmat.trials(trialindex(i_tr)).post_cue_frames+1
        doutWave = zeros(1,size(colmat_cr,3)*p.scr_imgmultipl); % each entry corresponds to a trigger
    else
        doutWave = zeros(1,size(colmat_cr,3)*p.scr_imgmultipl+p.scr_imgmultipl); % each entry corresponds to a trigger
    end
    
    % write trigger numbers into doutwave
    % trial start trigger
    doutWave(1) = p.trig.tr_start;
    % trial end trigger (presented at first frame after stimulation)
    doutWave(size(colmat_cr,3)*p.scr_imgmultipl+1) = p.trig.tr_stop;
    % condition trigger
    resp(i_tr).triggernum = ...
        p.trig.tr_con_type(resp(i_tr).condition); % condition 1 2 3 4 5 6
    switch resp(i_tr).trial_timing_type
        case 'regular'
            resp(i_tr).triggernum = resp(i_tr).triggernum + p.trig.tr_timing_type(1);
        case 'catch'
            resp(i_tr).triggernum = resp(i_tr).triggernum + p.trig.tr_timing_type(2);
    end
    doutWave(resp(i_tr).cue_onset_fr) = resp(i_tr).triggernum;
    % event trigger
    t.trigger = p.trig.event_type(resp(i_tr).event_type)+ ...
        p.trig.event_side(resp(i_tr).eventpos) + ...
        p.trig.event_cuevalid(resp(i_tr).condition);
    doutWave(resp(i_tr).event_onset_frames) = t.trigger;
    
    doutWave = [doutWave;zeros(triggersPerRefresh-1,numel(doutWave))]; doutWave=doutWave(:);
    samplesPerFlip = triggersPerRefresh * p.scr_imgmultipl;
    % figure; plot(doutWave)        
    
    % draw fixation cross
    Screen('DrawLines', ps.window, p.crs.lines, p.crs.width, p.crs.color, ps.center, 0);
    
    % send 0 before again to reset everything
    Datapixx('SetDoutValues', 0);
    Datapixx('RegWrRd');
    
    % write outsignal
    Datapixx('WriteDoutBuffer', doutWave, bufferAddress);
    % disp(Datapixx('GetDoutStatus'));
    Datapixx('SetDoutSchedule', 0, [samplesPerFlip, 2], numel(doutWave), bufferAddress); % 0 - scheduleOnset delay, [samplesPerFlip, 2] - sheduleRate in samples/video frame, framesPerTrial - maxSheduleFrames
    Datapixx('StartDoutSchedule');
    
    %% keyboard
    % start listening to keyboard
    KbQueueStart(ps.RespDev);
    KbQueueFlush(ps.RespDev);
    
    % flip to get everything synced
    Screen('Flip', ps.window, 0);
    
    %% loop across frames
    for i_fl = 1:frames.flips
        %% Drawing
        % RDK
        Screen('DrawDots', ps.window, dotmat(:,:,i_fl), dotsize(:,i_fl), colmat(:,:,i_fl), ps.center, 0, 0);
        % fixation cross
        Screen('DrawLines', ps.window, p.crs.lines, p.crs.width, colmat_cr(:,:,i_fl), ps.center, 0);
        
        %% start trigger schedule and start listening to response device
        if i_fl == 1 % send the trigger with the start of the 1st flip
            Datapixx('RegWrVideoSync');
        end

        %trs
        if i_fl > conmat.trials(i_tr).event_onset_frames/4 & i_fl < conmat.trials(i_tr).event_onset_frames/4+3
%             fprintf('event!')
%             unique(colmat(:,rdkidx(:,i_fl)==2,i_fl))
%             t.tt = colmat(:,rdkidx(:,i_fl)==1,i_fl);
            t.tt2 = conmat.trials(i_tr).event_type;
        end
        
        % Flip
        [timing(i_tr).VBLTimestamp(i_fl), timing(i_tr).StimulusOnsetTime(i_fl), timing(i_tr).FlipTimestamp(i_fl), timing(i_tr).Missed(i_fl)] = Screen('Flip', ps.window, 0);
        
        % send trigger/save timing/ reset timing
        if i_fl == 1
            % start trigger
            starttime=GetSecs;
            KbEventFlush(ps.RespDev); % flush keyboard
        end
        
        %% check for button presses
        [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
        key.presses{i_tr}(i_fl,:)=key.firstPress(key.keymap)>1;
        key.presses_t{i_tr}(i_fl,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl,:);
        if any(key.firstPress(key.keymap)>1)
            lptwrite(1,find(key.firstPress(key.keymap),1,'first'),500);
        end
    end
    %% ITI
    % draw fixation cross again
    Screen('DrawLines', ps.window, p.crs.lines, p.crs.width, p.crs.color, ps.center, 0);
    
    % flip
    Screen('Flip', ps.window, 0);
    
    % get time
    crttime = GetSecs;

       
    % add waiting period to check for late button presses
    ttt=WaitSecs(p.targ_respwin(2)/1000- ...
        ((conmat.trials(trialindex(i_tr)).post_cue_times + conmat.trials(trialindex(i_tr)).pre_cue_times) ...
        - conmat.trials(trialindex(i_tr)).event_onset_times) - ...
        p.targ_duration);

    
    % check for button presses
    [key.pressed, key.firstPress]=KbQueueCheck(ps.RespDev);
    key.presses{i_tr}(i_fl+1,:)=key.firstPress(key.keymap)>1;
    key.presses_t{i_tr}(i_fl+1,:)=(key.firstPress(key.keymap)-starttime).*key.presses{i_tr}(i_fl+1,:);
    
    
    
    
    %%%%%%%%
    % do behavioral calculation
    key.presses{i_tr}(1,:)=[];
    key.presses_t{i_tr}(1,:)=[];
    
    % get frame and timing of button press onsets
    resp(i_tr).button_presses_fr=nan(max(sum(key.presses{i_tr})),size(key.presses{i_tr},2));
    resp(i_tr).button_presses_t=resp(i_tr).button_presses_fr;
    for i_bt = 1:size(key.presses{i_tr},2)
        try
            resp(i_tr).button_presses_fr(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                find(key.presses{i_tr}(:,i_bt));
            resp(i_tr).button_presses_t(1:sum(key.presses{i_tr}(:,i_bt)),i_bt)=...
                key.presses_t{i_tr}(find(key.presses{i_tr}(:,i_bt)),i_bt)*1000; % in ms
        catch
            resp(i_tr).button_presses_fr(:,i_bt)=nan;
        end
    end
    
    % check for hits {'hit','miss','error'}
    resp(i_tr).event_response_type = {}; %{'hit','miss','error'}
    resp(i_tr).event_response_RT = []; %reaction time or nan
    % define correct responses
    if resp(i_tr).event_type == 1 % 1 = targethue decrease; 2 = increase
        button_index = key.keymap_ind == key.DECREASE;
    elseif resp(i_tr).event_type == 2 % 1 = targethue decrease; 2 = increase
        button_index = key.keymap_ind == key.INCREASE;
    end
    if ~isempty(resp(i_tr).button_presses_t) % any button presses?
        % target-button presses in time window?
        hitindex = (resp(i_tr).button_presses_t(:,button_index) > (resp(i_tr).event_onset_times*1000+p.targ_respwin(1)) & ...
            resp(i_tr).button_presses_t(:,button_index) < (resp(i_tr).event_onset_times*1000+p.targ_respwin(2)))& button_index;
        % button press of any other button in response window?
        errorindex = (resp(i_tr).button_presses_t > (resp(i_tr).event_onset_times*1000+p.targ_respwin(1)) & ...
            resp(i_tr).button_presses_t < (resp(i_tr).event_onset_times*1000+p.targ_respwin(2))) & ~button_index;
        if any(hitindex) % any hits?
            resp(i_tr).event_response_type = 'hit'; % save hits (button presses in correct time window)
            resp(i_tr).event_response_RT = resp(i_tr).button_presses_t(hitindex) - resp(i_tr).event_onset_times*1000;
        else % no hit
            if any(errorindex) % but an error?
                resp(i_tr).event_response_type = 'error';
                resp(i_tr).event_response_RT = resp(i_tr).button_presses_t(errorindex) - resp(i_tr).event_onset_times*1000;
            else % no, than it's a miss!
                resp(i_tr).event_response_type = 'FA';
            end
        end
    else % no button presses at all, then it is a miss!
        resp(i_tr).event_response_type = 'miss';
    end
    
    
    % wait
    ttt=WaitSecs(max(p.ITI)/1000-((GetSecs - crttime)));
    crttime2 = GetSecs; % for troubleshooting

    
    
end

% stop response recording
KbQueueRelease(ps.RespDev);

% send start trigger
if flag_training~=1
    PRPX_sendRecTrigger('stop')
end



end