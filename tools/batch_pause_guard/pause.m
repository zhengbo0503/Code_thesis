function varargout = pause(varargin)
%PAUSE - Reject interactive waits only while the batch runner adds this folder.
%   Timed pauses and pause('query'/'on'/'off') retain the builtin behavior.
%   Do not add this folder to a permanent MATLAB search path.
if nargin == 0
    error('run_all_experiments:InteractivePause', ...
        ['An experiment reached an interactive pause() check. This step ' ...
         'was stopped; inspect the preceding diagnostics. Run the experiment ' ...
         'individually, or use InteractivePauses=true, to inspect it interactively.']);
end
if nargout > 0
    [varargout{1:nargout}] = builtin('pause', varargin{:});
else
    builtin('pause', varargin{:});
end
end
