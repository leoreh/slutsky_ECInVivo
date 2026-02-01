
import re
import os
from collections import defaultdict

# List of FMA functions (from my previous analysis)
fma_funcs = [
    'ActivityTemplates','AngularVelocity','BrainStates','CCGParameters','CCG_old','CofiringCoefficient','CoherenceBands','CompareDistributions',
    'CountSpikesPerCycle','DefineZone','FieldPSP','FieldShift','FilterLFP','FindDeltaWaves','FindRipples','FindSpikeTrains','FindSpindles',
    'FindSubcycles','FiringCurve','FiringMap','FitCCG','FunctionalClustering','GenerateSpikeTrain','IsInZone','IsolationDistance','LinearVelocity',
    'MTCoherence','MTCoherogram','MTPointSpectrogram','MTPointSpectrum','MTSpectrogram','MTSpectrum','MapStats','NormalizeFields','PETHTransition',
    'PhaseCurve','PhaseMango','PhaseMap','PhasePrecession','QuietPeriods','RadialMaze','RadialMazeTurns','ReactivationStrength','ReconstructPosition',
    'RemoveArtefacts','RippleStats','SelectSpikes','ShortTimeCCG','SpectrogramBands','SurveyCCG','SurveyFiringMaps','SurveyPhasePrecession',
    'SurveyTuningCurves','SyncHist','SyncMap','TestRemapping','TestSkewness','ThresholdSpikes','TimeToPhase','TuneArtefactTimes',
    'CheckChronux','Contiguous','FindField','mtspecgrampt_optimized','AddSpikeTimes','CustomDefaults','GetAngles','GetChannels','GetCurrentSession',
    'GetCustomDefaults','GetEvents','GetEventTypes','GetLFP','GetPositions','GetSpikeAmplitudes','GetSpikeFeatures','GetSpikes','GetSpikeTimes',
    'GetSpikeWaveforms','GetUnits','GetWidebandData','GlobalSettings','SetCurrentSession','BatchInfo','CancelBatch','CleanBatches','DBAddFigure',
    'DBAddVariable','DBConnect','DBCreate','DBCreateTables','DBDuplicate','DBExportGallery','DBExternalStoragePath','DBGetFigures','DBGetValues',
    'DBGetVariables','DBList','DBListFigures','DBListVariables','DBMatchValues','DBMerge','DBRemove','DBRemoveFigures','DBRemoveVariables','DBUse',
    'DebugBatch','GetBatch','ShowBatch','StartBatch','CheckMyM','DBDisplay','GetNextField','GetNextItem','InsertDate','ParseBatch','RunBatch',
    'Accumulate','Array2Matrix','Array2PagedMatrix','BartlettTest','CircularANOVA','CircularConfidenceIntervals','CircularDistribution',
    'CircularMean','CircularRegression','CircularShift','CircularVariance','CompareSlopes','Concentration','ConcentrationTest','ConsolidateIntervals',
    'CountInIntervals','DistanceTransform','ExcludeIntervals','ExtendArray','FindInInterval','FisherTest','InIntervals','Intervals','IsExtremum',
    'IsFirstAfter','IsLastBefore','MatchPairs','mean2str','MultinomialConfidenceIntervals','nansem','npcdf','Restrict','RunningAverage','sem',
    'semedian','SineWavePeaks','SubtractIntervals','ToIntervals','WatsonU2Test','XCorr1','ZeroCrossings','ZeroToOne','clinspace','glinspace',
    'int2zstr','isastring','isdmatrix','isdscalar','isdvector','isimatrix','isiscalar','isivector','islmatrix','islscalar','islvector','isradians',
    'issamples','minmax','mean2str',
    'ChangeBinaryGain','LoadBinary','LoadEvents','LoadParameters','LoadPositions','LoadSpikeAmplitudes','LoadSpikeFeatures','LoadSpikeTimes',
    'LoadSpikeWaveforms','NewEvents','ProcessBinary','ResampleBinary','SaveBinary','SaveEvents','SaveRippleEvents','AdjustAxes','AdjustColorMap',
    'Bright','hsl2hsv','hsv2hsl','Monochrome','MultiPlotXY','own_clim','PiecewiseLinearAxis','PlotCCG','PlotCircularDistribution','PlotColorCurves',
    'PlotColorMap','PlotCSD','PlotDistribution2','PlotHVLines','PlotIntervals','PlotLinkage','PlotMean','PlotPhasePrecession','PlotRepeat',
    'PlotRippleStats','PlotSamples','PlotShortTimeCCG','PlotSlope','PlotSpikeWaveforms','PlotSync','PlotTicks','PlotXY','SideAxes','SplitTitle',
    'SquareSubplot','Subpanel','TableFigure','UIAddLine','UIInPolygon','UISelect'
]

# Case insensitive map
fma_map = {f.lower(): f for f in fma_funcs}
fma_regex = re.compile(r'\b(' + '|'.join(re.escape(f) for f in fma_funcs) + r')\b', re.IGNORECASE)

usage_map = defaultdict(set)

with open('fma_usages.txt', 'r', encoding='utf-8') as f:
    # Skip header
    lines = f.readlines()

current_file = None
for line in lines:
    line = line.strip()
    if not line or line.startswith('----') or line.startswith('Path ') or line.startswith('Searching Batch'):
        continue
    
    # Check if this line describes a file path or a match
    # Select-String output format: Path LineNumber Line
    # But I used `Select-Object Path, LineNumber, Line`, which formats it as table if not piped to string properly?
    # PowerShell `Out-File` with objects creates a table.
    # The output I saw in view_file was:
    # D:\Code\...   202   line content
    
    # I need to parse the fixed width or split strictly.
    # Path is usually long, then LineNumber, then Line.
    # The spacing is dynamic.
    
    # Heuristic: split by space. First token is path (if absolute).
    # But Windows paths contain spaces? No, user path `d:\Code\slutsky_ECInVivo` seems safe.
    # `d:\OneDrive...` has spaces!
    
    # Regex for line:  ^(?P<path>.:\\.*?) \s+ (?P<ln>\d+) \s+ (?P<content>.*)$
    match = re.match(r'^(?P<path>[A-Za-z]:\\.*?)\s+(?P<ln>\d+)\s+(?P<content>.*)$', line)
    if match:
        path = match.group('path')
        content = match.group('content')
        
        # Check for FMA functions in content
        found = fma_regex.findall(content)
        for fn in found:
            real_name = fma_map.get(fn.lower(), fn)
            # Exclude self-references if possible (function definition)
            # But here we are just searching usage.
            usage_map[real_name].add(path)

# Generate Report
with open('FMAToolbox_Usage_Report.md', 'w', encoding='utf-8') as f:
    f.write('# FMAToolbox Usage Report\n\n')
    f.write('This report details all scripts in the codebase capable of calling FMAToolbox functions.\n\n')
    
    # Sort functions by name
    sorted_funcs = sorted(usage_map.keys())
    
    for func in sorted_funcs:
        files = sorted(list(usage_map[func]))
        if not files: continue
        
        f.write(f'## {func}\n')
        f.write(f'**Used in {len(files)} files:**\n')
        for file in files:
            # Make relative path
            rel_path = file.replace(r'D:\Code\slutsky_ECInVivo\\', '').replace(r'D:\Code\slutsky_ECInVivo\d:\Code\slutsky_ECInVivo\\', '') # Fix duplicate prefix if any
            # Fix simple replace
            if 'd:\\code\\slutsky_ecinvivo\\' in file.lower():
                # case insensitive replace
                start = file.lower().find('d:\\code\\slutsky_ecinvivo\\')
                rel_path = file[start+len('d:\\code\\slutsky_ecinvivo\\'):]
            
            f.write(f'- `{rel_path}`\n')
        f.write('\n')
        
    f.write('## Summary\n')
    f.write(f'Total unique functions found: {len(sorted_funcs)}\n')

