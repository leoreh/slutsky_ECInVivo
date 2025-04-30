function grph_save(varargin)

% saves the figure handle in specified format and path
% supports jpg, fig, and ai (Adobe Illustrator) formats

% 10 jan 25 LH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fh', gcf, @(x) isa(x, 'matlab.ui.Figure'));
addOptional(p, 'fpath', 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs', @ischar);
addOptional(p, 'fname', 'figure', @ischar);
addOptional(p, 'frmt', 'ai');

parse(p, varargin{:});
fh = p.Results.fh;
fpath = p.Results.fpath;
fname = p.Results.fname;
frmt = p.Results.frmt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure path exists
if ~exist(fpath, 'dir')
    mkdir(fpath);
end

% simplify fname if it is an lme formula
fname = frml2char(fname);

% construct full file path
fullPath = fullfile(fpath, fname);

% assert figure format
if ~iscell(frmt) & ~isstring(frmt)
    frmt = {frmt};
end
nfrmts = length(frmt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifrmt = 1 : nfrmts
    switch frmt{ifrmt}
        case 'jpg'
            print(fh, fullPath, '-djpeg', '-r300');
        case 'fig'
            savefig(fh, fullPath);
        case 'ai'
            % save as eps with additional parameters for AI compatibility
            print(fh, fullPath, '-depsc', '-vector');
    end
end

end

% EOF