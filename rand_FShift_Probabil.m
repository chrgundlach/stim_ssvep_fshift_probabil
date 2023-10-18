function [ conmat ] = rand_FShift_Probabil(p,RDK,flag_training)
%rand_FShiftBase randomizes experimental conditions
% move onset only works for constant frequency for all RDKs (i.e. 120)




% set trial number etc
if flag_training~=0
    conmat.totaltrials = ceil(sum(p.stim.con_repeats_t)*(1+p.trial_irreg_prob));
    conmat.totalblocks = 1;
else
    conmat.totaltrials = round(sum(p.stim.con_repeats_e)*(1+p.trial_irreg_prob));
    conmat.totalblocks = p.stim.blocknum;
end
conmat.trialsperblock = conmat.totaltrials/conmat.totalblocks;

% matrix with onset times of on framesfor RDKs
t.onframesonset = nan(numel(RDK.RDK),p.scr_refrate*max(p.postcue_time/1000));
t.onframesonset_times = t.onframesonset; % onset times in s
for i_rdk = 1:numel(RDK.RDK)
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2));
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end

% move
t.movonset_frames=nan(1,p.scr_refrate*max(p.postcue_time/1000));
t.movonset_times=nan(1,p.scr_refrate*max(p.postcue_time/1000));
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;


%% start randomization
% randomize condition
if flag_training==0
    % distribute conditions
    conmat.mats.condition = [cell2mat(arrayfun(@(x,y) repmat(x,1,y), p.stim.condition, p.stim.con_repeats_e, 'UniformOutput', false)) ...
        cell2mat(arrayfun(@(x,y) repmat(x,1,y), p.stim.condition, p.stim.con_repeats_e*p.trial_irreg_prob, 'UniformOutput', false))];
    % define trial timing type
    conmat.mats.trial_timing_type = [repmat({'regular'},1,sum(p.stim.con_repeats_e)) ...
        repmat({'catch'},1,sum(p.stim.con_repeats_e*p.trial_irreg_prob))];
else % training
    % distribute conditions
    t.mat1 = cell2mat(arrayfun(@(x,y) repmat(x,1,y), p.stim.condition, p.stim.con_repeats_t, 'UniformOutput', false));
    t.mat2 = cell2mat(arrayfun(@(x,y) repmat(x,1,y), ...
        p.stim.condition, ceil(p.stim.con_repeats_t*p.trial_irreg_prob*(1/prod(p.stim.con_repeats_t*p.trial_irreg_prob))), ...
        'UniformOutput', false));
    conmat.mats.condition = [t.mat1 t.mat2(randsample(numel(t.mat2),ceil(sum(p.stim.con_repeats_t*p.trial_irreg_prob))))];
    % define trial timing type
    conmat.mats.trial_timing_type = [repmat({'regular'},1,sum(p.stim.con_repeats_t)) ...
        repmat({'catch'},1,ceil(sum(p.stim.con_repeats_t*p.trial_irreg_prob)))];
end

% determine attended RDK
conmat.mats.cue_direction = nan(1,conmat.totaltrials);
conmat.mats.cue_direction_label = cell(1,conmat.totaltrials);
conmat.mats.cue_validity = nan(1,conmat.totaltrials);
conmat.mats.cue_validity_label = cell(1,conmat.totaltrials);
t.cue_validity = [1 1 2 2 3 3]; t.cue_validity_label = {'valid';'invalid';'neutral'}; t.cue_cue_direction_label = {'RDK1','RDK2','both'};
for i_con = 1:numel(p.stim.condition)
    % determine cue direction [1 = left, 2 = right, 3 = both]
    conmat.mats.cue_direction(conmat.mats.condition==p.stim.condition(i_con)) = p.stim.RDK2attend(i_con);
    conmat.mats.cue_direction_label(conmat.mats.condition==p.stim.condition(i_con)) = t.cue_cue_direction_label(p.stim.RDK2attend(i_con));
    % randomize cue conditions (valid, invalid, neutral) [1,2,3]
    conmat.mats.cue_validity(conmat.mats.condition==p.stim.condition(i_con)) = t.cue_validity(i_con);
    conmat.mats.cue_validity_label(conmat.mats.condition==p.stim.condition(i_con)) = t.cue_validity_label(t.cue_validity(i_con));
end

% determine target RDKs
conmat.mats.event_pos=nan(1,conmat.totaltrials);
conmat.mats.event_pos_label = cell(1,conmat.totaltrials);
for i_con = 1:numel(p.stim.condition)
    % index all trials
    conmat.mats.event_pos(conmat.mats.condition == p.stim.condition(i_con)) = p.stim.targetRDK(i_con);
    conmat.mats.event_pos_label(conmat.mats.condition == p.stim.condition(i_con)) = {sprintf('RDK%1.0f',p.stim.targetRDK(i_con))};
end


% randomize post-cue target timing (catch trials and regular trials)
conmat.mats.event_onset_frames = nan(1,conmat.totaltrials);
% all RDKs on in control window
t.poss_frames_ctrl = find((p.targ_win_ctrl(1)/1000)<t.movonset_times & ...
    t.movonset_times<(p.targ_win_ctrl(2)/1000) & ...
    all(~isnan(t.onframesonset),1));

% all RDKs on in regular window
t.poss_frames_reg = find((p.targ_win_reg(1)/1000)<t.movonset_times & ...
    t.movonset_times<(p.targ_win_reg(2)/1000) & ...
    all(~isnan(t.onframesonset),1));

% randomize for all control trials
for i_con = 1:numel(p.stim.condition)
    % index condition and control trials
    t.idx = conmat.mats.condition == p.stim.condition(i_con) & strcmp(conmat.mats.trial_timing_type, 'catch');
    % randomly assign frames
    conmat.mats.event_onset_frames(t.idx) = [...
        repmat(t.poss_frames_ctrl,1,floor(sum(t.idx)/numel(t.poss_frames_ctrl))), ...
        t.poss_frames_ctrl(randsample(numel(t.poss_frames_ctrl),mod(sum(t.idx),numel(t.poss_frames_ctrl))))];
end

% randomize for all regular trials
for i_con = 1:numel(p.stim.condition)
    % index condition and control trials
    t.idx = conmat.mats.condition == p.stim.condition(i_con) & strcmp(conmat.mats.trial_timing_type, 'regular');
    % randomly assign frames
    conmat.mats.event_onset_frames(t.idx) = [...
        repmat(t.poss_frames_reg,1,floor(sum(t.idx)/numel(t.poss_frames_reg))), ...
        t.poss_frames_reg(randsample(numel(t.poss_frames_reg),mod(sum(t.idx),numel(t.poss_frames_reg))))];
end
conmat.mats.event_onset_times = conmat.mats.event_onset_frames./p.scr_refrate;
% % graphical check
% figure; subplot(2,1,1);histogram(conmat.mats.event_onset_frames(:),50);subplot(2,1,2);histogram(conmat.mats.event_onset_times(:),50)
% figure; subplot(2,1,1);histogram(diff(conmat.mats.event_onset_frames),50);subplot(2,1,2);histogram(conmat.mats.event_onset_times(:),50)
% 
% for i_tr = 1:100
% test(i_tr,:,:) = conmat.mats.event_onset_times;
% end
% figure; subplot(2,1,1); histogram(test(:)); subplot(2,1,2); histogram(diff(test,1,2))


% randomize pre-cue times
% all possible pre_cue_frames
t.allframes = (p.precue_time_n(1)/1000)*p.scr_refrate:(p.precue_time_n(2)/1000)*p.scr_refrate;
t.allframes = t.allframes(mod(t.allframes,p.scr_imgmultipl)==0); % only frames that are integers of frames per flip (i.e. 4)
if conmat.totaltrials<numel(t.allframes)
    conmat.mats.pre_cue_frames = t.allframes(randsample(1:numel(t.allframes),conmat.totaltrials));
else
    conmat.mats.pre_cue_frames = [repmat(t.allframes,1,floor(conmat.totaltrials/numel(t.allframes))) ...
        t.allframes(round(linspace(1,numel(t.allframes),mod(conmat.totaltrials,numel(t.allframes)))))];
end
conmat.mats.pre_cue_frames = conmat.mats.pre_cue_frames(randperm(numel(conmat.mats.pre_cue_frames)));
conmat.mats.pre_cue_times = conmat.mats.pre_cue_frames./p.scr_refrate;

% add pre-cue frames to events
conmat.mats.event_onset_times = conmat.mats.event_onset_times+conmat.mats.pre_cue_times;
conmat.mats.event_onset_frames = conmat.mats.event_onset_frames + conmat.mats.pre_cue_frames;

% randomize whether it's a chroma increase or decrease as target
% 1 = hue decrease; 2 = hue increase
conmat.mats.event_type = nan(1,conmat.totaltrials);
% randomize for all control trials
t.type = [1 2];
for i_con = 1:numel(p.stim.condition)
    % index condition and control trials
    t.idx = conmat.mats.condition == p.stim.condition(i_con) & strcmp(conmat.mats.trial_timing_type, 'catch');
    % randomly assign frames
    t.mat  = [...
        repmat(t.type,1,floor(sum(t.idx)/numel(t.type))), ...
        t.type(randsample(numel(t.type),mod(sum(t.idx),numel(t.type))))];
    conmat.mats.event_type(t.idx) = t.mat(randperm(numel(t.mat)));
end

% randomize for all regular trials
t.type = [1 2];
for i_con = 1:numel(p.stim.condition)
    % index condition and control trials
    t.idx = conmat.mats.condition == p.stim.condition(i_con) & strcmp(conmat.mats.trial_timing_type, 'regular');
    % randomly assign frames
    t.mat  = [...
        repmat(t.type,1,floor(sum(t.idx)/numel(t.type))), ...
        t.type(randsample(numel(t.type),mod(sum(t.idx),numel(t.type))))];
    conmat.mats.event_type(t.idx) = t.mat(randperm(numel(t.mat)));
end

% define RDK color
% cols2plot_lab = xyz2lab(rgb2xyz_custom(cols2plot,[0.665 0.321], [0.172 0.726], [0.163 0.039], [0.3127 0.3290]));
% c = makecform('lab2lch');
% lchpict = applycform(labpict,c);

% define RDK chroma increases and decreases
for i_RDK = 1:numel(RDK.RDK)
    t.xyzcol = rgb2xyz_custom(RDK.RDK(i_RDK).col(1,1:3),[0.665 0.321], [0.172 0.726], [0.163 0.039], [0.3127 0.3290]);
    t.lchcol = colorspace('LCH<-XYZ',t.xyzcol);
    t.targetcols_lch = [ t.lchcol+[0 t.lchcol(2)*(RDK.RDK(i_RDK).chromatarget(1)/100) 0]; ...
        t.lchcol+[0 t.lchcol(2)*(RDK.RDK(i_RDK).chromatarget(2)/100) 0]];
%     t.targetcols_lch = [ t.lchcol+[0 t.lchcol(2)*(RDK.RDK(i_RDK).chromatarget(1)/100) 0]; ...
%         t.lchcol+[0 t.lchcol(2)*(RDK.RDK(i_RDK).chromatarget(2)/100) 0]];

    t.targetcols_xyz = colorspace('XYZ<-LCH',t.targetcols_lch);
    t.targetcols_rgb = xyz2rgb_custom(t.targetcols_xyz,[0.665 0.321], [0.172 0.726], [0.163 0.039], [0.3127 0.3290]);
    t.targetcols_rgb_clipped = t.targetcols_rgb;
    t.targetcols_rgb_clipped(t.targetcols_rgb_clipped<0)=0;
    RDK.RDK(i_RDK).eventcol = {[t.targetcols_rgb_clipped(1,:) 1; RDK.RDK(i_RDK).col(2,:)]; ...
        [t.targetcols_rgb_clipped(2,:) 1; RDK.RDK(i_RDK).col(2,:)]};

    % bookkeeping
    RDK.RDK(i_RDK).col_lch = t.lchcol;
    RDK.RDK(i_RDK).eventcol_lch = t.targetcols_lch;

%     % retest
%     t.lchcol_t = colorspace('LCH<-XYZ', ...
%         rgb2xyz_custom(t.targetcols_rgb_clipped, ...
%         [0.665 0.321], [0.172 0.726], [0.163 0.039], [0.3127 0.3290]));
end

% define RGB target color
conmat.mats.event_col = cell(1,conmat.totaltrials);
for i_tr = 1:conmat.totaltrials
    conmat.mats.event_col{i_tr} = RDK.RDK(conmat.mats.event_pos(i_tr)).eventcol{conmat.mats.event_type(i_tr)};
end


%% randomize all information across experiment
t.tidx = randperm(conmat.totaltrials);
conmat.mats.condition = conmat.mats.condition(:,t.tidx);
conmat.mats.trial_timing_type = conmat.mats.trial_timing_type(t.tidx);
conmat.mats.cue_direction = conmat.mats.cue_direction(:,t.tidx);
conmat.mats.cue_direction_label = conmat.mats.cue_direction_label(:,t.tidx);
conmat.mats.cue_validity = conmat.mats.cue_validity(:,t.tidx);
conmat.mats.cue_validity_label = conmat.mats.cue_validity_label(:,t.tidx);
conmat.mats.event_onset_frames = conmat.mats.event_onset_frames(:,t.tidx);
conmat.mats.event_onset_times = conmat.mats.event_onset_times(:,t.tidx);
conmat.mats.event_pos = conmat.mats.event_pos(:,t.tidx);
conmat.mats.event_pos_label = conmat.mats.event_pos_label(:,t.tidx);
conmat.mats.event_type = conmat.mats.event_type(:,t.tidx);
conmat.mats.event_col = conmat.mats.event_col(:,t.tidx);
conmat.mats.pre_cue_frames = conmat.mats.pre_cue_frames(:,t.tidx);
conmat.mats.pre_cue_times = conmat.mats.pre_cue_times(:,t.tidx);

conmat.mats.block = repmat(1:conmat.totalblocks,conmat.trialsperblock,1);
conmat.mats.block = conmat.mats.block(:)';

conmat.RDK = RDK;

%% write all information into trial structure
% create frame mat, onset time for events

for i_tr = 1:conmat.totaltrials
    % trialnumber
    conmat.trials(i_tr).trialnum = i_tr;
    
    % block number
    conmat.trials(i_tr).blocknum = conmat.mats.block(i_tr);
    
    % condition
    conmat.trials(i_tr).condition = conmat.mats.condition(i_tr);
    
    % trial timing type
    conmat.trials(i_tr).trial_timing_type = conmat.mats.trial_timing_type{i_tr};

    % cue direction
    conmat.trials(i_tr).cue_direction = conmat.mats.cue_direction(i_tr);
    conmat.trials(i_tr).cue_direction_label = conmat.mats.cue_direction_label{i_tr};

    % cue validity
    conmat.trials(i_tr).cue_validity = conmat.mats.cue_validity(i_tr);
    conmat.trials(i_tr).cue_validity_label = conmat.mats.cue_validity_label{i_tr};

    % event onset times
    conmat.trials(i_tr).event_onset_frames = conmat.mats.event_onset_frames(i_tr);
    conmat.trials(i_tr).event_onset_times = conmat.mats.event_onset_times(i_tr);

    % event position
    conmat.trials(i_tr).event_pos = conmat.mats.event_pos(i_tr);
    conmat.trials(i_tr).event_pos_label = conmat.mats.event_pos_label{i_tr};

    % event type [chroma decrease or increase]
    conmat.trials(i_tr).event_type = conmat.mats.event_type(i_tr);
    conmat.trials(i_tr).event_col = conmat.mats.event_col{i_tr};
        
    % pre-cue times
    conmat.trials(i_tr).pre_cue_frames = conmat.mats.pre_cue_frames(:,i_tr);
    conmat.trials(i_tr).pre_cue_times = conmat.mats.pre_cue_times(:,i_tr);
    
    % post-cue times
    conmat.trials(i_tr).post_cue_times = max(p.postcue_time)/1000;
    
    % post-cue frames
    conmat.trials(i_tr).post_cue_frames = conmat.trials(i_tr).post_cue_times*p.scr_refrate;
end



    

end

