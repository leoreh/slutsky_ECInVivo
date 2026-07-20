function geno = mcu_geno(sbjID)
% MCU_GENO Map subject ids to a genotype categorical.
%
%   geno = MCU_GENO(sbjID)
%
%   SUMMARY:
%       The single subject -> genotype rule for the MCU project. Membership
%       lives in mcu_cfg (miceWT / miceMCU / miceCAG) and the labels in
%       cfg.lbl.grp, so a new cohort is added in one place.
%
%       Only genotypes actually present are kept as categories. This is not
%       cosmetic: fitlme / fitglme reject a fixed-effect predictor holding a
%       declared-but-absent level with "design matrix X must be of full
%       column rank", so a two-genotype table built while three genotypes are
%       defined would fail to fit. Category order still follows cfg.lbl.grp,
%       so Control stays the reference level and coefficient order is
%       unchanged for the two-genotype tables the manuscript already uses.
%
%       A subject in none of the lists is labelled Control with a warning.
%       Previously this was silent, which mislabelled the whole viral cohort
%       as controls.
%
%   INPUTS:
%       sbjID - <cell/categorical/string> [N x 1] subject ids, e.g. the
%               sbjID column produced by v2tbl via get_mname.
%
%   OUTPUTS:
%       geno  - <categorical> [N x 1] genotype, categories ordered as in
%               cfg.lbl.grp but restricted to those present.
%
%   DEPENDENCIES:
%       mcu_cfg.
%
%   HISTORY:
%       Jul 2026 (LH) split out of mcu_tblVivo for the CAG (viral KO) cohort.

cfg = mcu_cfg;

% Work on cellstr so cell, string and categorical inputs behave alike
sbjID = cellstr(string(sbjID(:)));

% Cohort lookup. Control is the fallback, hence the ones
genoIdx = ones(numel(sbjID), 1);
genoIdx(ismember(sbjID, cfg.miceMCU)) = 2;
genoIdx(ismember(sbjID, cfg.miceCAG)) = 3;

% Flag unregistered subjects instead of passing them off as controls
idxUnk = ~ismember(sbjID, [cfg.miceWT, cfg.miceMCU, cfg.miceCAG]);
if any(idxUnk)
    warning('mcu_geno:unknownSbj', ...
        ['Subject(s) not in any mcu_cfg cohort, labelled %s: %s. ', ...
        'Add them to mcu_cfg.'], cfg.lbl.grp{1}, ...
        strjoin(unique(sbjID(idxUnk))', ', '));
end

geno = categorical(genoIdx, 1 : numel(cfg.lbl.grp), cfg.lbl.grp);
geno = removecats(geno);

end     % EOF
