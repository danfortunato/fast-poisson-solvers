fprintf('Initializing path directories...')
root = fileparts([pwd, filesep]);
addpath(root)
addpath([root, '/adi']);
addpath([root, '/cube']);
addpath([root, '/cylinder']);
addpath([root, '/rectangle']);
addpath([root, '/solidsphere']);
addpath([root, '/tests']);
addpath([root, '/transforms']);
addpath([root, '/vis']);
figdir = strcat(root(1:end-5), '/paper/figures/');
fprintf(' done.\n');

mexcodes = {'adi', 'adi_tri'};
for i = 1:length(mexcodes)
    if ~exist(sprintf('%s.%s', mexcodes{i}, mexext), 'file')
        fprintf('Compiling %s.c...', mexcodes{i});
        file = which(sprintf('%s.c', mexcodes{i}));
        outdir = fileparts(file);
        mex('-silent', '-outdir', outdir, file);
        fprintf(' done.\n');
    end
end

clear root mexcodes file outdir
