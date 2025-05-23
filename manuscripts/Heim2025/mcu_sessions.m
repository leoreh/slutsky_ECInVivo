function [basepaths, v] = mcu_sessions(queryStr, vars)

if nargin < 2
    vars = [];
end

if strcmp(queryStr, 'mcu')

    basepaths = {...
        'lh132',...
        'lh133',...
        'lh134',...
        'lh136',...
        'lh140',...
        };

elseif strcmp(queryStr, 'wt')

    basepaths = {...
        'lh96',...
        'lh100',...
        'lh107',...
        'lh122',...
        'lh142',...
        };

elseif strcmp(queryStr, 'eeg')

    basepaths = {...
        'lh96',...
        'lh105',...
        'lh106',...
        };
            % 'lh100',...

elseif strcmp(queryStr, 'mcu_bsl')

    basepaths = {...
        'E:\Data\lh132\lh132_230413_094013',...
        'E:\Data\lh133\lh133_230413_094013',...
        'E:\Data\lh134\lh134_230504_091744',...
        'E:\Data\lh136\lh136_230519_090043',...
        'E:\Data\lh140\lh140_230619_090023',...
        'E:\Data\lh137\lh137_230516_091852',...
        };

elseif strcmp(queryStr, 'wt_bsl_ripp')

    basepaths = {...
        'E:\Data\lh107\lh107_220518_091200',...
        'E:\Data\lh122\lh122_221223_092656',...
        'E:\Data\lh123\lh123_221219_094508',...
        'E:\Data\lh126\lh126_230111_091208',...
        'E:\Data\lh119\lh119_221114_081305',...
        %         'E:\Data\lh99\lh99_220118_085402',...
        %         'E:\Data\lh96\lh96_220120_090157',...
        %         'E:\Data\lh142\lh142_231005_091832',...
        %         'E:\Data\lh100\lh100_220413_111004',...
        %         'E:\Data\lh130\lh130_230322_084541',...
        };

elseif strcmp(queryStr, 'mcu_wsh')

    basepaths = {...
        'E:\Data\lh132\lh132_230419_090121',...
        'E:\Data\lh133\lh133_230419_090121',...
        'E:\Data\lh134\lh134_230510_090806',...
        'E:\Data\lh136\lh136_230526_090813',...
        'E:\Data\lh140\lh140_230625_090049'};

elseif strcmp(queryStr, 'wt_bsl')

    basepaths = {...
        'E:\Data\lh96\lh96_220120_090157',...
        'E:\Data\lh100\lh100_220413_111004',...
        'E:\Data\lh107\lh107_220518_091200',...
        'E:\Data\lh122\lh122_221223_092656',...
        'E:\Data\lh142\lh142_231005_091832',...
        %         'E:\Data\lh99\lh99_220118_085402',...
        %         'E:\Data\lh130\lh130_230322_084541',...
        %         'E:\Data\lh126\lh126_230111_091208',...
        %         'E:\Data\lh123\lh123_221219_094508',...
        };

elseif strcmp(queryStr, 'wt_wsh')

    basepaths = {...
        'E:\Data\lh96\lh96_220126_085016',...
        'E:\Data\lh100\lh100_220420_095858',...
        'E:\Data\lh107\lh107_220524_091100',...
        'E:\Data\lh122\lh122_221230_090102',...
        'E:\Data\lh142\lh142_231011_092620',...
%         'E:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'wt_bac1')

    basepaths = {...
        'E:\Data\lh96\lh96_220122_090154',...
        'E:\Data\lh100\lh100_220415_100221',...
        'E:\Data\lh107\lh107_220520_093000',...
        'E:\Data\lh122\lh122_221225_091518',...
        'E:\Data\lh142\lh142_231007_095203',...
        %         'E:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'mcu_bac1')

    basepaths = {...
        'E:\Data\lh132\lh132_230415_112936',...
        'E:\Data\lh133\lh133_230415_112936',...
        'E:\Data\lh134\lh134_230506_092535',...
        'E:\Data\lh136\lh136_230521_085735',...
        'E:\Data\lh140\lh140_230621_090711',...
        };

elseif strcmp(queryStr, 'wt_bac2')

    basepaths = {...
        'E:\Data\lh96\lh96_220123_090009',...
        'E:\Data\lh100\lh100_220416_110253',...
        'E:\Data\lh107\lh107_220521_091000',...
        'E:\Data\lh122\lh122_221226_100133',...
        'E:\Data\lh142\lh142_231008_094841',...
        %         'E:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'wt_bac3')

    basepaths = {...
        'E:\Data\lh96\lh96_220124_090127',...
        'E:\Data\lh100\lh100_220418_101641',...
        'E:\Data\lh107\lh107_220522_093000',...
        'E:\Data\lh122\lh122_221228_102653',...
        'E:\Data\lh142\lh142_231009_102856',...
        %         'E:\Data\lh99\lh99_220124_090128',...
        };

elseif strcmp(queryStr, 'ra')

    basepaths = {...
        'E:\Data\Colleagues\RA\MCU1',...
        'E:\Data\Colleagues\RA\MCU2',...
        'E:\Data\Colleagues\RA\MCU3_20211203_084720',...
        'E:\Data\Colleagues\RA\MCU4',...
        'E:\Data\Colleagues\RA\MCU5',...
        };

elseif contains(queryStr, 'lh')
    xlsname = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Data summaries\sessionList.xlsx';
    [~, basepaths] = getSessionVars('mname', queryStr, 'varsFile', ["session"],...
        'varsName', ["session"], 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);

elseif strcmp(queryStr, 'prePost')

    basepaths{1} = {...
        'E:\Data\lh96\lh96_220119_090035',...
        'E:\Data\lh96\lh96_220120_090157',...
        'E:\Data\lh96\lh96_220126_085016',...
        };

    basepaths{2} = {...
        'E:\Data\lh99\lh99_220117_090735',...
        'E:\Data\lh99\lh99_220118_085402',...
        'E:\Data\lh99\lh99_220124_090128',...
        'E:\Data\lh99\lh99_220125_090041',...
        };

    basepaths{3} = {...
        'E:\Data\lh100\lh100_220413_111004',...
        'E:\Data\lh100\lh100_220420_095858',...
        };

    basepaths{4} = {...
        'E:\Data\lh107\lh107_220517_094900',...
        'E:\Data\lh107\lh107_220518_091200',...
        'E:\Data\lh107\lh107_220524_091100',...
        };

    basepaths{5} = {...
        'E:\Data\lh122\lh122_221221_092941',...
        'E:\Data\lh122\lh122_221223_092656',...
        'E:\Data\lh122\lh122_221230_090102',...
        'E:\Data\lh122\lh122_221231_090112',...
        'E:\Data\lh122\lh122_230101_090144',...
        };

    basepaths{6} = {...
        'E:\Data\lh142\lh142_231005_091832',...
        'E:\Data\lh142\lh142_231011_092620',...
        'E:\Data\lh142\lh142_231012_093127',...
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
    varIdx = find(strcmp(vars, 'sleep_statesEmg'));
    if ~isempty(varIdx)
        varsName(varIdx) = "ssEmg";
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


    [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', vars,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""]);
end

end