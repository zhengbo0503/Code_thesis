% probe_env - report toolchain availability; must never crash
fprintf('MATLAB   : %s (%s)\n', version, version('-release'));
fprintf('arch     : %s\n', computer('arch'));
have = @(n) ~isempty(which(n));
fprintf('Advanpix : %d\n', have('mp'));
if have('mp')
    try
        mp.Digits(34); x = mp('1')/mp('3');
        fprintf('  mp ok, digits=%d, 1/3=%s\n', mp.Digits(), char(x));
    catch e
        fprintf('  mp FAILED: %s\n', e.message);
    end
end
fprintf('Anymatrix: %d\n', have('anymatrix'));
fprintf('m2tikz   : %d\n', have('matlab2tikz'));
fprintf('PROBE_DONE\n');
