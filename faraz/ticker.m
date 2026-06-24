function t = ticker(t,first,last)
% t         either ticker text to initialize new ticker or ticker structure to continue existing
% first     starting iteration number (default = 0)
% last      final iteration number (optional)

% check inputs
if nargin < 3, last = []; end
if nargin < 2, first = 0; end

% check for existing ticker or start new one
if ~isstruct(t)
    if isempty(last)
        t = struct('text',t,'curr',first,'last',last,'tstart',0,'tcurr',0);
        fprintf(1,'%s:  %d',t.text,t.curr);
    else
        t = struct('text',t,'curr',first,'last',last,'tstart',0,'tcurr',0);
        fprintf(1,'%s:  %d / %d',t.text,t.curr,t.last);
    end
    t.tstart = tic;
elseif isstruct(t)
    if ~isempty(t.last)
        if t.tcurr == 0
            fprintf(1,repmat('\b',1, ...
                3 + numel([num2str(t.last) num2str(t.curr)])));
        else
            timestr = sprintf('(%.6g s)',t.tcurr);
            fprintf(1,repmat('\b',1, ...
                4 + numel([num2str(t.last) num2str(t.curr)]) + numel(timestr)));
        end
        t.curr = t.curr + 1;
        t.tcurr = toc(t.tstart);
        fprintf(1,'%d / %d (%.6g s)',t.curr,t.last,t.tcurr);
        if t.curr >= t.last
            fprintf('\n');
        end
    else
        if t.tcurr == 0
            fprintf(1,repmat('\b',1,numel(num2str(t.curr))));
        else
            timestr = sprintf('(%.6g s)',t.tcurr);
            fprintf(1,repmat('\b',1,1 + numel(num2str(t.curr)) + numel(timestr)));
        end
        t.curr = t.curr + 1;
        t.tcurr = toc(t.tstart);
        fprintf(1,'%d (%.6g s)',t.curr,t.tcurr);
    end
end
end