function nDone = fr_migrate(basepaths)
% FR_MIGRATE Rewrites <basename>.fr.mat to the spk_rate schema.
%
%   nDone = FR_MIGRATE(basepaths)
%
%   SUMMARY:
%       Renames the fields of an existing fr struct instead of recomputing
%       it. mfr is carried across verbatim, which is what keeps every result
%       fitted on tblVivo.fr identical across the rename. The fields spk_rate
%       no longer produces (states, gain, ratio, norm, fanoFactor, gini_pop,
%       gini_unit, mfr_std) are dropped. Per-state values now come from
%       spk_states; the two readers left on the dropped fields are
%       spikes/utypes/legacy/selectUnits.m and the Ruggiero2024 scripts, and
%       both stop working on a migrated file.
%
%       Migration is per basepath, so recordings that belong to another
%       project keep the legacy schema until they are passed in.
%
%   INPUTS:
%       basepaths - (Cell) Recording folders to migrate.
%
%   OUTPUTS:
%       nDone     - (Num)  Number of files rewritten.
%
%   DEPENDENCIES:
%       backup_file.
%
%   HISTORY:
%       Created: 260719
%
%   See also: SPK_RATE

nDone = 0;

for iPath = 1 : numel(basepaths)

    basepath = basepaths{iPath};
    [~, basename] = fileparts(basepath);
    frFile = fullfile(basepath, [basename, '.fr.mat']);

    if ~isfile(frFile)
        fprintf('[FR_MIGRATE] no fr file in %s\n', basename);
        continue
    end

    s = load(frFile, 'fr');
    if isfield(s.fr, 'rate')
        continue                    % already on the new schema
    end

    % the mea pipeline writes a different fr struct under the same filename
    % (mea_frPrep: .fr/.frOrig/.t/.uGood). it has no calc_fr schema to
    % rename and must not be touched
    if ~isfield(s.fr, 'strd')
        fprintf('[FR_MIGRATE] %s is not a calc_fr file, skipped\n', basename);
        continue
    end

    frNew.rate = s.fr.strd;
    frNew.t    = s.fr.tstamps(:)';
    frNew.mfr  = s.fr.mfr;

    % info also carries the rename, so a migrated file and a fresh one
    % report their parameters under the same names
    frNew.info.binsize = s.fr.info.binsize;
    frNew.info.winCalc = s.fr.info.winCalc;
    frNew.info.winBL   = s.fr.info.winBL;
    frNew.info.smet    = s.fr.info.smoothMethod;
    frNew.info.migrated = datetime("now");

    backup_file(frFile);
    fr = frNew;                                                  %#ok<NASGU>
    save(frFile, 'fr')
    nDone = nDone + 1;
end

fprintf('[FR_MIGRATE] rewrote %d of %d files\n', nDone, numel(basepaths));

end     % EOF
