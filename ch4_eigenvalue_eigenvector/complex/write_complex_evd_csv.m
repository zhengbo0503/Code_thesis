function filename = write_complex_evd_csv(results,kind,dataDir)
%WRITE_COMPLEX_EVD_CSV - Export two-sided accuracy results at 17 digits
%   Usage: filename = write_complex_evd_csv(results,kind,dataDir)
%   Purpose: Write one local CSV, in mode-major/case-minor order. Keep new
%       complex outputs separate from the current thesis figure contract.
%   Arguments: results is a completed sweep; kind is varying_kappa or
%       varying_n; dataDir is the output directory.
%   Outputs: filename is the atomically replaced CSV; schema is in README.md.
kind = validatestring(kind,{'varying_kappa','varying_n'});
assert(all(results.completed(:)),'Cannot export an incomplete sweep.');
header = {'mode','seed','n','kappa_target','kappa_prescribed','kappa_actual'};
for metric = {'forward','residual','unitarity'}
    for method = {'mp3','jacobi','mp2','eig'}
        header{end+1} = [metric{1} '_' method{1}]; %#ok<AGROW>
    end
end
for metric = {'info','sweeps','rotations'}
    for method = {'mp3','jacobi','mp2'}
        header{end+1} = [metric{1} '_' method{1}]; %#ok<AGROW>
    end
end
header = [header,{'scond_at','reference_at'}];
[ncase,nmode] = size(results.completed);
values = nan(ncase*nmode,numel(header));
for jm = 1:nmode
    for j = 1:ncase
        if strcmp(kind,'varying_kappa')
            n = results.n; kappa = results.kappas(j);
        else
            n = results.nvalues(j); kappa = results.kappa;
        end
        values((jm-1)*ncase+j,:) = [results.modes(jm),results.seeds(j,jm), ...
            n,kappa,results.prescribedKappas(j,jm),results.actualKappas(j,jm), ...
            reshape(results.forward(j,jm,:),1,[]), ...
            reshape(results.residual(j,jm,:),1,[]), ...
            reshape(results.unitarityV(j,jm,:),1,[]), ...
            reshape(results.info(j,jm,:),1,[]), ...
            reshape(results.sweeps(j,jm,:),1,[]), ...
            reshape(results.rotations(j,jm,:),1,[]), ...
            results.scalecondAt(j,jm),results.boundAt(j,jm)];
    end
end
assert(isreal(values) && all(isfinite(values(:))), ...
    'Accuracy output contains nonreal or nonfinite values.');
filename = fullfile(dataDir,['complex_evals_' kind '_' results.profile '.csv']);
local_write_atomic(filename,header,values);
end

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
