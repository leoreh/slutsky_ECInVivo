function [basepaths, v] = mcu_sessions(queryStr, vars)

if nargin < 2
    vars = [];
end

if strcmp(queryStr, 'mcu_bsl')

    basepaths = {...
        'F:\Data\lh132\lh132_230413_094013',...
        'F:\Data\lh133\lh133_230413_094013',...
        'F:\Data\lh134\lh134_230504_091744',...
        'F:\Data\lh136\lh136_230519_090043',...
        'F:\Data\lh140\lh140_230619_090023',...
        'F:\Data\lh137\lh137_230516_091852',...
        };

elseif strcmp(queryStr, 'wt_bsl_ripp')

    basepaths = {...
        'F:\Data\lh107\lh107_220518_091200',...
        'F:\Data\lh122\lh122_221223_092656',...
        'F:\Data\lh123\lh123_221219_094508',...
        'F:\Data\lh126\lh126_230111_091208',...
        'F:\Data\lh119\lh119_221114_081305',...
        %         'F:\Data\lh99\lh99_220118_085402',...
        %         'F:\Data\lh96\lh96_220120_090157',...
        %         'F:\Data\lh142\lh142_231005_091832',...
        %         'F:\Data\lh100\lh100_220413_111004',...
        %         'F:\Data\lh130\lh130_230322_084541',...
        };

elseif strcmp(queryStr, 'mcu_wsh')

    basepaths = {...
        'F:\Data\lh132\lh132_230419_090121',...
        'F:\Data\lh133\lh133_230419_090121',...
        'F:\Data\lh134\lh134_230510_090806',...
        'F:\Data\lh136\lh136_230526_090813',...
        'F:\Data\lh140\lh140_230625_090049'};

elseif strcmp(queryStr, 'wt_bsl')

    basepaths = {...
        'F:\Data\lh96\lh96_220120_090157',...
        'F:\Data\lh100\lh100_220413_111004',...
        'F:\Data\lh107\lh107_220518_091200',...
        'F:\Data\lh122\lh122_221223_092656',...
        'F:\Data\lh142\lh142_231005_091832',...
        %         'F:\Data\lh99\lh99_220118_085402',...
        %         'F:\Data\lh130\lh130_230322_084541',...
        %         'F:\Data\lh126\lh126_230111_091208',...
        %         'F:\Data\lh123\lh123_221219_094508',...
        };

elseif strcmp(queryStr, 'wt_wsh')

    basepaths = {...
        'F:\Data\lh96\lh96_220126_085016',...
        'F:\Data\lh100\lh100_220420_095858',...
        'F:\Data\lh107\lh107_220524_091100',...
        'F:\Data\lh122\lh122_221230_090102',...
        'F:\Data\lh142\lh142_231011_092620',...
%         'F:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'wt_bac1')

    basepaths = {...
        'F:\Data\lh96\lh96_220122_090154',...
        'F:\Data\lh100\lh100_220415_100221',...
        'F:\Data\lh107\lh107_220520_093000',...
        'F:\Data\lh122\lh122_221225_091518',...
        'F:\Data\lh142\lh142_231007_095203',...
        %         'F:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'wt_bac2')

    basepaths = {...
        'F:\Data\lh96\lh96_220123_090009',...
        'F:\Data\lh100\lh100_220416_110253',...
        'F:\Data\lh107\lh107_220521_091000',...
        'F:\Data\lh122\lh122_221226_100133',...
        'F:\Data\lh142\lh142_231008_094841',...
        %         'F:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'wt_bac3')

    basepaths = {...
        'F:\Data\lh96\lh96_220124_090127',...
        'F:\Data\lh100\lh100_220418_101641',...
        'F:\Data\lh107\lh107_220522_093000',...
        'F:\Data\lh122\lh122_221228_102653',...
        'F:\Data\lh142\lh142_231009_102856',...
        %         'F:\Data\lh99\lh99_220124_090128',...
        };

elseif contains(queryStr, 'lh')
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [~, basepaths] = getSessionVars('mname', queryStr, 'varsFile', ["session"],...
        'varsName', ["session"], 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);

elseif strcmp(queryStr, 'prePost')

    basepaths{1} = {...
        'F:\Data\lh96\lh96_220119_090035',...
        'F:\Data\lh96\lh96_220120_090157',...
        'F:\Data\lh96\lh96_220126_085016',...
        };

    basepaths{2} = {...
        'F:\Data\lh99\lh99_220117_090735',...
        'F:\Data\lh99\lh99_220118_085402',...
        'F:\Data\lh99\lh99_220124_090128',...
        'F:\Data\lh99\lh99_220125_090041',...
        };

    basepaths{3} = {...
        'F:\Data\lh100\lh100_220413_111004',...
        'F:\Data\lh100\lh100_220420_095858',...
        };

    basepaths{4} = {...
        'F:\Data\lh107\lh107_220517_094900',...
        'F:\Data\lh107\lh107_220518_091200',...
        'F:\Data\lh107\lh107_220524_091100',...
        };

    basepaths{5} = {...
        'F:\Data\lh122\lh122_221221_092941',...
        'F:\Data\lh122\lh122_221223_092656',...
        'F:\Data\lh122\lh122_221230_090102',...
        'F:\Data\lh122\lh122_221231_090112',...
        'F:\Data\lh122\lh122_230101_090144',...
        };

    basepaths{6} = {...
        'F:\Data\lh142\lh142_231005_091832',...
        'F:\Data\lh142\lh142_231011_092620',...
        'F:\Data\lh142\lh142_231012_093127',...
        };

    orgArray{1} = [1, 2, 3];
    orgArray{2} = [1, 2, 3, 4];
    orgArray{3} = [1, 3];
    orgArray{4} = [1, 2, 3];
    orgArray{5} = [1 : 5];
    orgArray{6} = [1, 3, 4];

end

if ~isempty(vars)

    varsName = vars;
    varIdx = find(strcmp(vars, 'sleep_states'));
    if ~isempty(varIdx)
        varsName(varIdx) = "ss";
    end
    varIdx = find(strcmp(vars, 'psdEmg'));
    if ~isempty(varIdx)
        varsName(varIdx) = "psd";
    end
    varIdx = find(strcmp(vars, 'frEmg'));
    if ~isempty(varIdx)
        varsName(varIdx) = "fr";
    end
    varIdx = find(strcmp(vars, 'swv_metrics'));
    if ~isempty(varIdx)
        varsName(varIdx) = "swv";
    end


    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', vars,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
end

end