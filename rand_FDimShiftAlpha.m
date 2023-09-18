function [ conmat ] = rand_FDimShiftAlpha(p,RDK,flag_training)
%rand_FDimShiftAlpha randomizes experimental conditions
% move onset only works for constant frequency for all RDKs (i.e. 120)




% set trial number etc
if flag_training~=0
    conmat.totaltrials = numel(p.stim.condition)*numel(p.stim.eventnum)*p.stim.con_repeats_t;
    conmat.totalblocks = 1;
else
    conmat.totaltrials = numel(p.stim.condition)*numel(p.stim.eventnum)*p.stim.con_repeats_e;
    conmat.totalblocks = p.stim.blocknum;
end
conmat.trialsperblock = conmat.totaltrials/conmat.totalblocks;

% matrix with onset times of on framesfor RDKs
t.onframesonset = nan(numel(RDK.RDK),p.scr_refrate*p.stim.time_postcue);
t.onframesonset_times = t.onframesonset; % onset times in s
for i_rdk = 1:numel(RDK.RDK)
    t.mat = ceil(1:p.scr_refrate/RDK.RDK(i_rdk).freq:size(t.onframesonset,2));
    t.onframesonset(i_rdk,t.mat)=1;
    t.onframesonset_times(i_rdk,t.mat)=t.mat./p.scr_refrate;
end
% move
t.movonset_frames=nan(1,p.scr_refrate*p.stim.time_postcue);
t.movonset_times=nan(1,p.scr_refrate*p.stim.time_postcue);
t.mat = 1:p.scr_refrate/RDK.RDK(1).mov_freq:size(t.movonset_frames,2);
t.movonset_frames(t.mat)=1;
t.movonset_times(t.mat)=t.mat./p.scr_refrate;


%% start randomization
% randomize cue (1 = attend to RDK 1 [color]; 2 = attend to RDK 2 [motion])
conmat.mats.cue = [repmat(p.stim.condition(1),1,conmat.totaltrials/numel(p.stim.condition)) ...
    repmat(p.stim.condition(2),1,conmat.totaltrials/numel(p.stim.condition))];

% randomize event numbers per trial
% check if event numbers add up
if mod(conmat.totaltrials,numel(p.stim.eventnum))~=0 || mod(conmat.totaltrials/numel(p.stim.eventnum),2)~=0
    error('rando:eventdistribute', 'Can not distribute event numbers (ratio: [%s]) equally across %1.0f trials',...
        num2str(p.stim.eventnum), conmat.totaltrials);
end
conmat.mats.eventnum = repmat(p.stim.eventnum,1,conmat.totaltrials/numel(p.stim.eventnum));

% randomize eventtype (1 = target; 2 = distractor)
if mod(sum(conmat.mats.eventnum==2),4)~=0 || (mod(sum(conmat.mats.eventnum==2)/4,2)~=0 && sum(conmat.mats.eventnum==2)~=4)
    error('rando:evtypedistribute', 'Can not distribute event types equally across %1.0f trials',...
        conmat.totaltrials);
end
conmat.mats.eventtype = nan(max(p.stim.eventnum),conmat.totaltrials);
conmat.mats.eventtype(:,conmat.mats.eventnum==1) = repmat([1 2; nan nan], 1,sum(conmat.mats.eventnum==1)/2);
conmat.mats.eventtype(:,conmat.mats.eventnum==2) = repmat([1 1 2 2; 1 2 1 2], 1,sum(conmat.mats.eventnum==2)/4);

% determine event RDK and also event class: RDK 1 = color event; RDK 2 = movement event
t.mat = [1 2; 2 1];
conmat.mats.eventRDK = nan(max(p.stim.eventnum),conmat.totaltrials);
for i_r = 1:size(conmat.mats.eventRDK,1)
    for i_c = 1:size(conmat.mats.eventRDK,2)
        if ~isnan(conmat.mats.eventtype(i_r,i_c))
            conmat.mats.eventRDK(i_r,i_c)=t.mat(conmat.mats.eventtype(i_r,i_c),conmat.mats.cue(i_c));
        end
    end
end

% write class of event
conmat.mats.eventclass = repmat({''},max(p.stim.eventnum),conmat.totaltrials);
[conmat.mats.eventclass{conmat.mats.eventRDK==1}] = deal('colorchange');
[conmat.mats.eventclass{conmat.mats.eventRDK==2}] = deal('globalmotion');

% randomize start color of colored RDK
conmat.mats.colRDKcol = nan(1,conmat.totaltrials);
for i_c = 1:2 % cue
    for i_ev = 1:numel(p.stim.eventnum)
        % index cue and eventnum
        t.idx = conmat.mats.cue == i_c & conmat.mats.eventnum == p.stim.eventnum(i_ev);
        % distribute equally and fill up the rest randomly
        t.mat = [repmat(1:size(p.col_colors,1),1,floor(sum(t.idx)/size(p.col_colors,1))), ...
            randsample(1:size(p.col_colors,1), mod(sum(t.idx),size(p.col_colors,1)))];
        conmat.mats.colRDKcol(t.idx) = t.mat(randperm(numel(t.mat)));
    end
end


% randomize start movement direction of non-colored moving RDK
conmat.mats.movRDKmovdir = nan(1,conmat.totaltrials);
for i_c = 1:2 % cue
    for i_ev = 1:numel(p.stim.eventnum)
        t.idx = conmat.mats.cue == i_c & conmat.mats.eventnum == p.stim.eventnum(i_ev);
        t.mat = [repmat(1:size(p.mov_directions,1),1,floor(sum(t.idx)/size(p.mov_directions,1))), ...
            randsample(1:size(p.mov_directions,1), mod(sum(t.idx),size(p.mov_directions,1)))];
        conmat.mats.movRDKmovdir(t.idx) = t.mat(randperm(numel(t.mat)));
    end
end


% randomize event colors
conmat.mats.eventcol = nan(max(p.stim.eventnum),conmat.totaltrials);
% index color events
t.idx1 = (repmat(conmat.mats.cue,2,1)==1 & conmat.mats.eventtype == 1) | ...
    (repmat(conmat.mats.cue,2,1)==2 & conmat.mats.eventtype == 2);
t.cols = 1:size(p.col_colors,1); % what colors are possible?
for i_col = 1:numel(t.cols) % loop across colors
    t.idx2 = repmat(conmat.mats.colRDKcol,2,1)==i_col; % find all trials with specific stimulus colors
    t.idx = t.idx1 & t.idx2; % index color events for trials with specific color
    % create vector for event colors with an (almost) uniform distribution
    t.mat = [repmat(t.cols(~ismember(t.cols, i_col)),1,floor(sum(t.idx,'all')/(numel(t.cols)-1))) ...
        randsample(t.cols(~ismember(t.cols, i_col)), mod(sum(t.idx,'all'),numel(t.cols)-1))];
    conmat.mats.eventcol(t.idx) = t.mat(randperm(numel(t.mat)));
end
% arrayfun(@(x) sum(conmat.mats.eventcol==x, 'all'), unique(conmat.mats.eventcol(~isnan(conmat.mats.eventcol))))


% randomize event directions
conmat.mats.eventdirection = nan(max(p.stim.eventnum),conmat.totaltrials);
% index color events
t.idx1 = (repmat(conmat.mats.cue,2,1)==2 & conmat.mats.eventtype == 1) | ...
    (repmat(conmat.mats.cue,2,1)==1 & conmat.mats.eventtype == 2);
t.mov = 1:size(p.mov_directions,1); % what directions are possible?
for i_mov = 1:numel(t.mov) % loop across movement directions
    t.idx2 = repmat(conmat.mats.movRDKmovdir,2,1)==i_mov; % find all trials with specific movement direction
    t.idx = t.idx1 & t.idx2; % index movement events for trials with specific movement direction
    % create vector for event colors with an (almost) uniform distribution
    t.mat = [repmat(t.mov(~ismember(t.mov, i_mov)),1,floor(sum(t.idx,'all')/(numel(t.mov)-1))) ...
        randsample(t.mov(~ismember(t.mov, i_mov)), mod(sum(t.idx,'all'),numel(t.mov)-1))];
    conmat.mats.eventdirection(t.idx) = t.mat(randperm(numel(t.mat)));
end
% arrayfun(@(x) sum(conmat.mats.eventdirection==x, 'all'), unique(conmat.mats.eventdirection(~isnan(conmat.mats.eventdirection))))


 % pre-allocate possible presentation times [works for two events so far]
conmat.mats.event_onset_frames = nan(max(p.stim.eventnum),conmat.totaltrials);
% position for color events
t.poss_frames_c = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_col-p.stim.event.min_offset));
t.poss_frames_c_1 = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_col-p.stim.event.min_offset-p.stim.event.min_dist-0.01));
t.poss_frames_c_2 = find(p.stim.event.min_onset+p.stim.event.min_dist<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_col-p.stim.event.min_offset));
% position for movement events
t.poss_frames_m = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_mov-p.stim.event.min_offset));
t.poss_frames_m_1 = find(p.stim.event.min_onset<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_mov-p.stim.event.min_offset-p.stim.event.min_dist-0.01));
t.poss_frames_m_2 = find(p.stim.event.min_onset+p.stim.event.min_dist<t.movonset_times & ...
    t.movonset_times<(p.stim.time_postcue-p.stim.event.length_mov-p.stim.event.min_offset));


% for single events first (color)
t.idx = repmat(conmat.mats.eventnum,2,1) == 1 & conmat.mats.eventRDK == 1;
if sum(t.idx(:))<numel(t.poss_frames_c) % if more possible positons than actual events
    t.poss_frames_mat = t.poss_frames_c(randsample(numel(t.poss_frames_c),sum(t.idx(:))));
else
    t.poss_frames_mat = [repmat(t.poss_frames_c,1,floor(sum(t.idx(:))/numel(t.poss_frames_c)))...
        t.poss_frames_m (randsample(numel(t.poss_frames_m),mod(sum(t.idx(:)),numel(t.poss_frames_c))))];
end
% write to frame mat
conmat.mats.event_onset_frames(t.idx)=t.poss_frames_mat;

% for single events first (movement)
t.idx = repmat(conmat.mats.eventnum,2,1) == 1 & conmat.mats.eventRDK == 2;
if sum(t.idx(:))<numel(t.poss_frames_m) % if more possible positons than actual events
    t.poss_frames_mat = t.poss_frames_m(randsample(numel(t.poss_frames_m),sum(t.idx(:))));
else
    t.poss_frames_mat = [repmat(t.poss_frames_m,1,floor(sum(t.idx(:))/numel(t.poss_frames_m)))...
        t.poss_frames_m(randsample(numel(t.poss_frames_m),mod(sum(t.idx(:)),numel(t.poss_frames_m))))];
end
% write to frame mat
conmat.mats.event_onset_frames(t.idx)=t.poss_frames_mat;


% for two events
t.idx = repmat(conmat.mats.eventnum,2,1) == 2;
t.idx2 = find(t.idx(1,:));
t.poss_frames_mat = [];
% loop across all events
for i_ev = 1:numel(t.idx2) % index first events
    % first event
    switch conmat.mats.eventRDK(1,t.idx2(i_ev))
        case 1 % color
            t.poss_frames_mat(1,i_ev) = t.poss_frames_c_1(randsample(numel(t.poss_frames_c_1),1));
        case 2 % movement
            t.poss_frames_mat(1,i_ev) = t.poss_frames_m_1(randsample(numel(t.poss_frames_m_1),1));
    end
    % second event
    switch conmat.mats.eventRDK(2,t.idx2(i_ev))
        case 1 % color
            t.idx3 = find(t.poss_frames_c_2>(t.poss_frames_mat(1,i_ev)+p.stim.event.min_dist*p.scr_refrate));
            t.poss_frames_mat(2,i_ev) = t.poss_frames_c_2(t.idx3(randsample(numel(t.idx3),1)));
        case 2 % movement
            t.idx3 = find(t.poss_frames_m_2>(t.poss_frames_mat(1,i_ev)+p.stim.event.min_dist*p.scr_refrate));
            t.poss_frames_mat(2,i_ev) = t.poss_frames_m_2(t.idx3(randsample(numel(t.idx3),1)));
            % randomly select from pssoble frames (beta distribution random number) --> compensate for righward distribution?
            % t.idx5 = round(betarnd(1,5,1)*(t.idx4(end)-t.idx4(1))+t.idx4(1));
            % t.idx5 = round(betarnd(1.2,3,1)*(t.idx4(end)-t.idx4(1))+t.idx4(1));
    end
end
conmat.mats.event_onset_frames(t.idx)=t.poss_frames_mat;

% extract onset times
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
t.allframes = p.stim.time_precue(1)*p.scr_refrate:p.stim.time_precue(2)*p.scr_refrate;
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

%% randomize all information across experiment
t.tidx = randperm(conmat.totaltrials);
conmat.mats.cue = conmat.mats.cue(:,t.tidx);
conmat.mats.eventnum = conmat.mats.eventnum(:,t.tidx);
conmat.mats.eventtype = conmat.mats.eventtype(:,t.tidx);
conmat.mats.eventclass = conmat.mats.eventclass(:,t.tidx);
conmat.mats.eventRDK = conmat.mats.eventRDK(:,t.tidx);
conmat.mats.colRDKcol = conmat.mats.colRDKcol(:,t.tidx);
conmat.mats.movRDKmovdir = conmat.mats. movRDKmovdir(:,t.tidx);
conmat.mats.eventcol = conmat.mats.eventcol(:,t.tidx);
conmat.mats.eventdirection = conmat.mats.eventdirection(:,t.tidx);
conmat.mats.event_onset_frames = conmat.mats.event_onset_frames(:,t.tidx);
conmat.mats.event_onset_times = conmat.mats.event_onset_times(:,t.tidx);
conmat.mats.pre_cue_frames = conmat.mats.pre_cue_frames(:,t.tidx);
conmat.mats.pre_cue_times = conmat.mats.pre_cue_times(:,t.tidx);

conmat.mats.block = repmat(1:conmat.totalblocks,conmat.trialsperblock,1);
conmat.mats.block = conmat.mats.block(:);

%% write all information into trial structure
% create frame mat, onset time for events

for i_tr = 1:conmat.totaltrials
    % trialnumber
    conmat.trials(i_tr).trialnum = i_tr;
    
    % block number
    conmat.trials(i_tr).blocknum = conmat.mats.block(i_tr);
    
    % cue ((RDK1, RDK2) [1,2])
    conmat.trials(i_tr).cue = conmat.mats.cue(i_tr);
    
    % number of events
    conmat.trials(i_tr).eventnum = conmat.mats.eventnum(i_tr);
    
    % type of events ((target, distractor) [1, 2])
    conmat.trials(i_tr).eventtype = conmat.mats.eventtype(:,i_tr);
    
    % which RDK shows event?
    conmat.trials(i_tr).eventRDK = conmat.mats.eventRDK(:,i_tr);

    % which event class? global motion or color change
    conmat.trials(i_tr).eventclass = conmat.mats.eventclass(:,i_tr);

    % what is the color of the colored RDK in each trial (RDK 1 and RDK 3)
    conmat.trials(i_tr).colRDKcol = conmat.mats.colRDKcol(:,i_tr);

    % what is the movement direction of the grey moving RDK in each trial (RDK 2 and RDK 4)
    conmat.trials(i_tr).movRDKmovdir = conmat.mats.movRDKmovdir(:,i_tr);
    
     % event color 
    conmat.trials(i_tr).eventcol = conmat.mats.eventcol(:,i_tr);

    % eventdirection ((according to RDK.event.direction) [1 2 3 4])
    conmat.trials(i_tr).eventdirection = conmat.mats.eventdirection(:,i_tr);
    
    % event onset frames
    conmat.trials(i_tr).event_onset_frames = conmat.mats.event_onset_frames(:,i_tr);
    
    % event onset times
    conmat.trials(i_tr).event_onset_times = conmat.mats.event_onset_times(:,i_tr);
    
    % pre-cue frames
    conmat.trials(i_tr).pre_cue_frames = conmat.mats.pre_cue_frames(:,i_tr);
    
    % pre-cue times
    conmat.trials(i_tr).pre_cue_times = conmat.mats.pre_cue_times(:,i_tr);
    
    % post-cue times
    conmat.trials(i_tr).post_cue_times = p.stim.time_postcue;
    
    % post-cue frames
    conmat.trials(i_tr).post_cue_frames = p.stim.time_postcue*p.scr_refrate;
end



    

end

