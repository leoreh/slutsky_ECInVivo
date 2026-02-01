$exclusion = "packages\FMAToolbox"

$batch1 = "\b(ActivityTemplates|AngularVelocity|BrainStates|CCGParameters|CCG_old|CofiringCoefficient|CoherenceBands|CompareDistributions|CountSpikesPerCycle|DefineZone|FieldPSP|FieldShift|FilterLFP|FindDeltaWaves|FindRipples|FindSpikeTrains|FindSpindles|FindSubcycles|FiringCurve|FiringMap|FitCCG|FunctionalClustering|GenerateSpikeTrain|IsInZone|IsolationDistance|LinearVelocity|MTCoherence|MTCoherogram|MTPointSpectrogram|MTPointSpectrum|MTSpectrogram|MTSpectrum|MapStats|NormalizeFields|PETHTransition|PhaseCurve|PhaseMango|PhaseMap|PhasePrecession|QuietPeriods|RadialMaze|RadialMazeTurns|ReactivationStrength|ReconstructPosition|RemoveArtefacts|RippleStats|SelectSpikes|ShortTimeCCG|SpectrogramBands|SurveyCCG|SurveyFiringMaps|SurveyPhasePrecession|SurveyTuningCurves|SyncHist|SyncMap|TestRemapping|TestSkewness|ThresholdSpikes|TimeToPhase|TuneArtefactTimes)\b"

$batch2 = "\b(CheckChronux|Contiguous|FindField|mtspecgrampt_optimized|AddSpikeTimes|CustomDefaults|GetAngles|GetChannels|GetCurrentSession|GetCustomDefaults|GetEvents|GetEventTypes|GetLFP|GetPositions|GetSpikeAmplitudes|GetSpikeFeatures|GetSpikes|GetSpikeTimes|GetSpikeWaveforms|GetUnits|GetWidebandData|GlobalSettings|SetCurrentSession|BatchInfo|CancelBatch|CleanBatches|DBAddFigure|DBAddVariable|DBConnect|DBCreate|DBCreateTables|DBDuplicate|DBExportGallery|DBExternalStoragePath|DBGetFigures|DBGetValues|DBGetVariables|DBList|DBListFigures|DBListVariables|DBMatchValues|DBMerge|DBRemove|DBRemoveFigures|DBRemoveVariables|DBUse|DebugBatch|GetBatch|ShowBatch|StartBatch|CheckMyM|DBDisplay|GetNextField|GetNextItem|InsertDate|ParseBatch|RunBatch)\b"

$batch3 = "\b(Accumulate|Array2Matrix|Array2PagedMatrix|BartlettTest|CircularANOVA|CircularConfidenceIntervals|CircularDistribution|CircularMean|CircularRegression|CircularShift|CircularVariance|CompareSlopes|Concentration|ConcentrationTest|ConsolidateIntervals|CountInIntervals|DistanceTransform|ExcludeIntervals|ExtendArray|FindInInterval|FisherTest|InIntervals|Intervals|IsExtremum|IsFirstAfter|IsLastBefore|MatchPairs|mean2str|MultinomialConfidenceIntervals|nansem|npcdf|Restrict|RunningAverage|sem|semedian|SineWavePeaks|SubtractIntervals|ToIntervals|WatsonU2Test|XCorr1|ZeroCrossings|ZeroToOne|clinspace|glinspace|int2zstr|isastring|isdmatrix|isdscalar|isdvector|isimatrix|isiscalar|isivector|islmatrix|islscalar|islvector|isradians|issamples|minmax|mean2str)\b"

$batch4 = "\b(ChangeBinaryGain|LoadBinary|LoadEvents|LoadParameters|LoadPositions|LoadSpikeAmplitudes|LoadSpikeFeatures|LoadSpikeTimes|LoadSpikeWaveforms|NewEvents|ProcessBinary|ResampleBinary|SaveBinary|SaveEvents|SaveRippleEvents|AdjustAxes|AdjustColorMap|Bright|hsl2hsv|hsv2hsl|Monochrome|MultiPlotXY|own_clim|PiecewiseLinearAxis|PlotCCG|PlotCircularDistribution|PlotColorCurves|PlotColorMap|PlotCSD|PlotDistribution2|PlotHVLines|PlotIntervals|PlotLinkage|PlotMean|PlotPhasePrecession|PlotRepeat|PlotRippleStats|PlotSamples|PlotShortTimeCCG|PlotSlope|PlotSpikeWaveforms|PlotSync|PlotTicks|PlotXY|SideAxes|SplitTitle|SquareSubplot|Subpanel|TableFigure|UIAddLine|UIInPolygon|UISelect)\b"

$batch5 = "\b(Map|Distance|Frequency|Phase|Sync|Recall|Store|Diff|Filter|Smooth|Bin|Clip|Insert|Interpolate|Match|Threshold|Samples|mM|sz|wrap)\b"

$files = Get-ChildItem -Path . -Recurse -Filter *.m | Where-Object { $_.FullName -notmatch [regex]::Escape($exclusion) }

Write-Host "Searching Batch 1..."
$files | Select-String -Pattern $batch1 | Select-Object Path, LineNumber, Line | Out-File fma_usages.txt -Encoding utf8

Write-Host "Searching Batch 2..."
$files | Select-String -Pattern $batch2 | Select-Object Path, LineNumber, Line | Out-File fma_usages.txt -Append -Encoding utf8

Write-Host "Searching Batch 3..."
$files | Select-String -Pattern $batch3 | Select-Object Path, LineNumber, Line | Out-File fma_usages.txt -Append -Encoding utf8

Write-Host "Searching Batch 4..."
$files | Select-String -Pattern $batch4 | Select-Object Path, LineNumber, Line | Out-File fma_usages.txt -Append -Encoding utf8

Write-Host "Searching Batch 5 (Ambiguous)..."
$files | Select-String -Pattern $batch5 | Select-Object Path, LineNumber, Line | Out-File fma_ambiguous.txt -Encoding utf8
