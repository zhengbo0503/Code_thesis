function write_figure_csv(destPath, data, header)
%WRITE_FIGURE_CSV - write a thesis figure CSV natively and publish it
%
%   Usage:
%       write_figure_csv(destPath, values, header)
%       write_figure_csv(destPath, tbl)
%
%   Purpose:
%       Single writer for every CSV consumed by a thesis figure. It writes
%       two identical files, both formatted from the same in-memory array so
%       the published copy is never a re-read of a rounded intermediate:
%           (1) destPath itself, keeping the chapter's native output/ or
%               data/ copy where the experiment has always put it;
%           (2) the same base name inside the thesis figure directory, which
%               is what the .tex files actually read.
%       Values are printed with %.17g, the shortest format guaranteed to
%       round-trip an IEEE double, so a rerun is not silently rounded to the
%       precision of the previously published numbers.
%       Each file is written under a temporary name in its own destination
%       directory and then renamed over the target, so a concurrent LaTeX
%       build can never read a half-written CSV.
%
%   Arguments:
%       (1) destPath - char/string native destination, e.g.
%           "../output/lau.csv". Its base name is reused for the published
%           thesis copy.
%       (2) data - numeric matrix whose columns are the fields, or a table.
%       (3) header - cell array of column names. Required for numeric data;
%           for a table it defaults to the table's VariableNames.
%
%   Outputs:
%       None. Writes the two CSV files described above.
%
%   Notes:
%       The thesis directory may be overridden with the THESIS_FIGS_DIR
%       environment variable, which exists so the export path can be tested
%       against a temporary directory without touching the real thesis.
%
%   Author: Zhengbo Zhou

DEFAULT_THESIS_FIGS = '/Users/cyae/Dropbox/tex/phdthesis/figs';

thesisFigs = getenv('THESIS_FIGS_DIR');
if isempty(thesisFigs)
    thesisFigs = DEFAULT_THESIS_FIGS;
end

if istable(data)
    if nargin < 3 || isempty(header)
        header = data.Properties.VariableNames;
    end
    values = table2array(data);
else
    values = data;
end
values = double(values);

if nargin < 3 && ~exist('header', 'var')
    error('write_figure_csv:MissingHeader', ...
        'A header is required when the data is not a table.');
end
header = cellstr(header);
header = header(:).';

if size(values, 2) ~= numel(header)
    error('write_figure_csv:WidthMismatch', ...
        '%s: header has %d names but the data has %d columns.', ...
        destPath, numel(header), size(values, 2));
end

% Report rather than repair: a non-finite entry is a real experiment result
% and must not be silently replaced, but it must not pass unnoticed either.
nBad = sum(~isfinite(values(:)));
if nBad > 0
    warning('write_figure_csv:NonFinite', ...
        '%s: %d non-finite value(s) written as-is; check the experiment.', ...
        destPath, nBad);
end

local_write_atomic(char(destPath), header, values);

[~, base, ext] = fileparts(char(destPath));
if isfolder(thesisFigs)
    local_write_atomic(fullfile(thesisFigs, [base ext]), header, values);
else
    warning('write_figure_csv:NoThesisDir', ...
        ['Thesis figure directory %s does not exist; wrote only the ' ...
         'native copy %s.'], thesisFigs, destPath);
end

end

% -------------------------------------------------------------------------
function local_write_atomic(filename, header, values)
%LOCAL_WRITE_ATOMIC Write one CSV via a temporary file, then rename it.

outDir = fileparts(filename);
if isempty(outDir)
    outDir = pwd;
end
if ~isfolder(outDir)
    mkdir(outDir);
end

tmp = [tempname(outDir) '.csv'];
fid = fopen(tmp, 'w');
if fid < 0
    error('write_figure_csv:Open', ...
        'Cannot open a temporary file in %s.', outDir);
end

try
    fprintf(fid, '%s\n', strjoin(header, ','));
    if ~isempty(values)
        rowFormat = [repmat('%.17g,', 1, size(values, 2) - 1), '%.17g\n'];
        fprintf(fid, rowFormat, values.');
    end
    fclose(fid);
catch err
    fclose(fid);
    if isfile(tmp)
        delete(tmp);
    end
    rethrow(err);
end

[ok, msg] = movefile(tmp, filename, 'f');
if ~ok
    if isfile(tmp)
        delete(tmp);
    end
    error('write_figure_csv:Replace', ...
        'Could not replace %s: %s', filename, msg);
end

end
