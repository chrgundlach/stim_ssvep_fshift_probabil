function [smoothed] = random_walk_smoothing(data_in, varargin)
% random walk smoothing to create near-uniform ditribution of values
%   Detailed explanation goes here


%% preset parameters
if nargin < 1
    help random_walk_smoothing;
    return
end

allpos = min(data_in):max(data_in);             % possible values for data
maxnum_iter = 10000;                            % max number of iterations

%% check for input
if ~(round(numel(varargin)/2) == numel(varargin)/2)
    error('Odd number of input arguments??')
end

for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'allpos' % define possible values for data
            if numel(Value) == 2 % case allpos represents highest and lowest value
                allpos = min(Value):max(Value);
            elseif numel(Value) > 2
                allpos = Value;
            else
                error('too few arguments for parameter allpos')
            end
        case 'iteration' % define maximum number of iterations
            maxnum_iter = Value;                       
    end
end

%% do random walk with replacement
% extract deviation from ideal state
data_edges = [allpos(1)-mean(diff(allpos))/2 allpos+mean(diff(allpos))/2];
% figure; subplot(3,1,1);histogram(data_in, data_edges); subplot(3,1,2); plot(data_in); subplot(3,1,3); plot(diff(data_in))
[bincount,t.edges] = histcounts(data_in,data_edges); % count values for each bin
bincount_diff = bincount-mean(bincount); 
bincount_diff_err_new = mean(abs(bincount_diff)); % error as maximum deviation from ideal near-uniform state
bincount_diff_err_old = inf; % initialize new error value


% start random replacement
smoothed = data_in;
iter_num = 0;
while iter_num <= maxnum_iter & bincount_diff_err_new > 0 & bincount_diff_err_old > bincount_diff_err_new
    bincount_diff_err_old = bincount_diff_err_new;
    
    % do resampling while maxnum_iter wasn't reached and error is larger than 0, error did get smaller previously
    [bincount_diff_sort_v bincount_diff_sort_idx] = sort(bincount_diff,'descend');
    % replace values
    idx_large_count = find(smoothed==allpos(bincount_diff_sort_idx(1))); % find index for largest value
    val_small_count = allpos(bincount_diff_sort_idx(find(bincount_diff_sort_v == bincount_diff_sort_v(end)))); % find value with lowest count
    smoothed(datasample(idx_large_count,1))=datasample(val_small_count,1); % replace with new value
    % update values
    [bincount,t.edges] = histcounts(smoothed,data_edges); % count values for each bin
    bincount_diff = bincount-mean(bincount);
    bincount_diff_err_new = mean(abs(bincount_diff)); % error as maximum deviation from ideal near-uniform state
    iter_num = iter_num+1;
%     figure; subplot(3,1,1);histogram(smoothed, data_edges); subplot(3,1,2); plot(smoothed); subplot(3,1,3); plot(diff(smoothed))
end



end

