function [psd_data, psd_cfg] = mcu_psdOrg(expDsgn)

% organizes psd data from multiple sessions into format compatible with
% fieldtrip's ft_freqstatistics. loads data according to mcu_sessions and
% arranges in cells according to specified dimorder (e.g., subj x freq x time)
%
% INPUT
%   expDsgn     numeric. experimental design preset {1}
%
% OUTPUT
%   dout        cell of psd mats organized in a way convinient for later
%               processing
%
% CALLS
%   mcu_sessions
%   basepaths2vars
%   as_loadConfig
%
% 06 Jun 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd_cfg.expDsgn = expDsgn;
switch psd_cfg.expDsgn

    case 1
        grps = mcu_sessions('wt');
        psd_cfg.fmt = '{state}[mouse x freq x time]';
        flg_emg = true;
        vars = {'psdEmg'};
        nstates = 2;

    case 2
        grps = {'wt_bsl'; 'mcu_bsl'};
        psd_cfg.fmt = '{gen}[state x freq x mouse]';
        flg_emg = false;
        vars = {'psd'};

    otherwise
        error('Unsupported experimental design: %d', expDsgn);
end

ngrps = numel(grps);
flg_norm = true;            % normalize to broadband power

if flg_norm
    psd_cfg.txt_yax = 'Norm. Power';
else
    psd_cfg.txt_yax = 'Power [dB]';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
psdCell = cell(1, numel(grps));

% organize data based on output format
switch expDsgn

    case 1                      % '{state}[mouse x freq x time]'
        
        % load data for each group, in this case mouse 
        for igrp = 1 : ngrps
            basepaths = mcu_sessions(grps{igrp});
            [psdCell{igrp}, psdInfo] = load_psd(basepaths, vars, flg_norm);
        end

        % organize dout
        psd_data = cell(1, nstates);
        for istate = 1 : nstates
            % get data for current state across all groups
            stateData = cellfun(@(x) squeeze(x(istate, :, :)), psdCell, 'UniformOutput', false);

            % concatenate along mouse dimension (3rd dim becomes mouse)
            psd_data{istate} = cat(3, stateData{:});

            % permute to get [mouse x freq x time] format
            psd_data{istate} = permute(psd_data{istate}, [3 1 2]);
        end

    case 2                      % '{gen}[state x freq x mouse]'
        
        % load data for each genotype group
        psd_data = cell(1, ngrps);
        for igrp = 1 : ngrps
            basepaths = mcu_sessions(grps{igrp});
            [psd_data{igrp}, psdInfo] = load_psd(basepaths, vars, flg_norm);
        end

    otherwise
        error('Unsupported output format');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd_cfg.faxis = psdInfo.faxis;
psd_cfg.sstates = psdInfo.sstates;
psd_cfg.snames = psdInfo.snames;
psd_cfg.clr = psdInfo.clr;
psd_cfg.grps = grps;
psd_cfg.fmt = psd_cfg.fmt;
psd_cfg.flg_norm = flg_norm;
psd_cfg.flg_emg = flg_emg;

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psdMat, psdInfo] = load_psd(basepaths, vars, flg_norm)

% load psd
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
psd = catfields([v(:).psd], 'addim', true);
psdMat = psd.psd;

% normalize if requested
if flg_norm
    broadPow = sum(psdMat, 2);
    psdMat = psdMat ./ broadPow;
end

psdInfo = v(1).psd.info;

end