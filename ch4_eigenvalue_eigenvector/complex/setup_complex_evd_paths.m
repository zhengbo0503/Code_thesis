function root = setup_complex_evd_paths()
%SETUP_COMPLEX_EVD_PATHS - Add the complex two-sided and real shared routines
%   Usage: root = setup_complex_evd_paths()
%   Purpose: Resolve paths independently of the working directory.
%   Arguments: None.
%   Outputs: root is this complex experiment directory.
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(fileparts(root),'shared'));
end
