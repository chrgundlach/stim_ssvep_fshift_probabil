function [] = pres_feedback(responses,p,ps, key,RDK)
%PRES_FEEDBACK calculate and display feedback for previous block
%   Returns Percentage hit, false alarms, and reaction time

WaitSecs(0.5);
%% calculate feedback
% get number of all events
t.num_presses=sum(cellfun(@(x) sum(sum(~isnan(x))),{responses.button_presses_t}));    % number of total button presses
summ.targnum = numel(responses);

summ.hits = sum(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)));
summ.errors = sum(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.event_response_type},'UniformOutput',false)));
summ.misses = sum(cell2mat(cellfun(@(x) strcmpi(x,'miss'),{responses.event_response_type},'UniformOutput',false)));
summ.FA = sum(cell2mat(cellfun(@(x) strcmpi(x,'FA'),{responses.event_response_type},'UniformOutput',false)));

t.mat = {responses.event_response_RT};
summ.RT_hit_mean = mean(cell2mat(t.mat(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)))));
summ.RT_hit_std = std(cell2mat(t.mat(cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)))));
summ.RT_error_mean = mean(cell2mat(t.mat(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.event_response_type},'UniformOutput',false)))));
summ.RT_error_std = std(cell2mat(t.mat(cell2mat(cellfun(@(x) strcmpi(x,'error'),{responses.event_response_type},'UniformOutput',false)))));


% behavioral effects separately for experimental conditions
for i_con = 1:numel(p.stim.condition)
    summ.targnum(1+i_con) = sum([responses.condition]==i_con);
    summ.hits(1+i_con) = sum(cell2mat(cellfun(@(x,y) strcmpi(x,'hit') & (y==i_con),...
        {responses.event_response_type},{responses.condition},'UniformOutput',false)));
    summ.errors(1+i_con) = sum(cell2mat(cellfun(@(x,y) strcmpi(x,'hit') & (y==i_con),...
        {responses.event_response_type},{responses.condition},'UniformOutput',false)));
    summ.misses(1+i_con) = sum(cell2mat(cellfun(@(x,y) strcmpi(x,'miss') & (y==i_con),...
        {responses.event_response_type},{responses.condition},'UniformOutput',false)));
    summ.FA(1+i_con) = sum(cell2mat(cellfun(@(x,y) strcmpi(x,'FA') & (y==i_con),...
        {responses.event_response_type},{responses.condition},'UniformOutput',false)));
    t.mat = {responses.event_response_RT};

    summ.RT_hit_mean(1+i_con) = mean(cell2mat(t.mat( ...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)) & ...
        cell2mat({responses.condition}) == i_con...
        )));
    summ.RT_hit_std(1+i_con) = std(cell2mat(t.mat( ...
        cell2mat(cellfun(@(x) strcmpi(x,'hit'),{responses.event_response_type},'UniformOutput',false)) & ...
        cell2mat({responses.condition}) == i_con...
        )));
end

% pixels for shift into 4 quadrants
quadshift = [p.scr_res(1)*(1/4) p.scr_res(2)*(1/4); p.scr_res(1)*(3/4) p.scr_res(2)*(1/4); ...
    p.scr_res(1)*(1/4) p.scr_res(2)*(3/4); p.scr_res(1)*(3/4) p.scr_res(2)*(3/4)];

%% presentation of results
% KbQueueCreate(ps.RespDev, keysOfInterest) %
% KbQueueStart(ps.RespDev);

% output to screen
fprintf('\nResults for events at RDK1 [%1.02f %1.02f %1.02f] and RDK2 [%1.02f %1.02f %1.02f]:\n', ...
    RDK.RDK(1).col(1,1:3)', RDK.RDK(2).col(1,1:3)')
t.textin = {'All:                  '};
t.cue_validity_label = {'valid  ';'valid  ';'invalid';'invalid';'neutral';'neutral'};
for i_con = 1:numel(p.stim.condition)
    t.textin{i_con+1} = sprintf('%s RDK%1.0f:', t.cue_validity_label{i_con}, p.stim.targetRDK(i_con));
end

for i_con = 1:numel(t.textin)
    fprintf(1,...
        '\n%s Hitrate: %06.2f; Hits: %02.0f; Misses: %02.0f; errors: %02.0f; FA: %02.0f; RT: M: %3.0f, Std: %3.0f ms',...
        t.textin{i_con}, ...
        summ.hits(i_con)/summ.targnum(i_con)*100, ...
        summ.hits(i_con), ...
        summ.misses(i_con), ...
        summ.errors(i_con), ...
        summ.FA(i_con), ...
        summ.RT_hit_mean(i_con), ...
        summ.RT_hit_std(i_con))
end
fprintf('\nMit "q" geht es weiter.\n')


% draw text and stimuli (before shifting to the quadrants)
% offscreen window
[ps.offwin,ps.offrect]=Screen('OpenOffscreenWindow',p.scr_num, [0 0 0 0], [0 0 p.scr_res(1)/2 p.scr_res(2)/2], [], [], []);
% get center of offscreen window
[ps.xCenter_off, ps.yCenter_off] = RectCenter(ps.offrect);

text2present=                   [...                % text for feedback
    'P A U S E'...
    sprintf('\n\nHits = %1.0f; Hitrate = %1.2f%%',summ.hits(1), summ.hits(1)/summ.targnum(1)*100)...
    sprintf('\n\nReaktionszeit: M = %1.0fms, Std = %1.0fms',summ.RT_hit_mean(1), summ.RT_hit_std(1))...
    sprintf('\n\n\nMisses = %1.0f; errors = %1.0f',summ.misses(1), summ.errors(1))...
    sprintf('\n\nFA = %1.0f',summ.FA(1)	)...
    ];

% draw text
Screen('TextSize', ps.offwin, 18);
% DrawFormattedText(tx.instruct, INST.text{1}, tx.xCenter_off, p.scr_res(1)/2 * 0.1, p.stim_color);
DrawFormattedText(ps.offwin, text2present, 'center', p.scr_res(1)/2 * 0.2, p.crs.color);

for i_quad = 1:4 % shifst to quadrants
    newpos_stim(:,i_quad) = ...
        CenterRectOnPointd(ps.offrect,quadshift(i_quad,1),quadshift(i_quad,2))';
end
 
% [key.pressed, key.firstPress]=KbQueueCheck;
key.rkey=key.SECRET;
[key.keyisdown,key.secs,key.keycode] = KbCheck; 
while ~(key.keycode(key.rkey)==1)                       % continuously present feedback (wait for q)
    [key.keyisdown,key.secs,key.keycode] = KbCheck;
    Screen('DrawTextures', ps.window, repmat(ps.offwin,1,4),[], newpos_stim, [], [], [], []);
    Screen('Flip', ps.window, 0);                   % flip screen
end
try Screen('Close', ps.offwin); end


end

