classdef reduct_displayer < handle
    % REDUCT_DISPLAYER Interactive GUI for dimensionality reduction and clustering of waveforms
    %
    % REDUCT_DISPLAYER creates an interactive GUI that allows users to:
    %   1. Perform dimensionality reduction on waveforms
    %   2. Cluster the data
    %   3. Visualize and compare waveforms across different clusters and labels
    %   4. Export results for further analysis
    %
    % Syntax:
    %   app = reduct_displayer(ied_obj)
    %   app = reduct_displayer(ied_obj, wave_labels)
    %   app = reduct_displayer(ied_obj, wave_labels, Name, Value)
    %
    % Required Input:
    %   ied_obj - IED.data object containing the waveforms to analyze.
    %             only the accepted waveforms will be used.
    %
    % Optional Input:
    %   wave_labels - label data for each waveform. A cell array with the following structure:
    %       1. Label type name (text scalar)
    %       2. Label meanings for each unique label value
    %          (vector of unique, at the same order as values in 3 - unique(...,'stable'))
    %       3. Label assignments for each waveform 
    %          (vector of labels for each waveform, at the same order as ied_obj.pos(ied_obj.accepted))
    %
    % Name-Value Pairs:
    %   Any property of the GUI elements can be set, including:
    %   - 'reduction_function' - {'PCA','t-SNE','UMAP','cMDS'}
    %   - 'reduction_source' - {'Base Waveforms','Distance Mat'} 
    %   - 'cluster_method' - {'Manual','by Label','kmeans','dbscan','hierarchical','None'}
    %   - 'cluster_source' - {'Reduced Data 2d','Distance Mat','Base Waveforms'}
    %   See properties (public) section for complete list of settable parameters
    %
    % Examples:
    %   % Basic usage with just IED object
    %   app = reduct_displayer(ied_obj);
    %
    %   % With labels and custom parameters
    %   labels = {
    %     'Type',  {'Spike','Sharp'}, categorical(repmat("Spike",100,1));
    %     'Side',  {'Left','Right'},  categorical(repmat("Left",100,1))
    %   };
    %   app = reduct_displayer(ied_obj, labels, ...
    %                         'reduction_function', 't-SNE', ...
    %                         'cluster_method', 'kmeans', ...
    %                         'kmean_nClusters', 5);
    %
    % See also: IED.data, tsne, umap, kmeans, dbscan, linkage
    %
    % To Do:
    %   * Implement create code from app.
    %   * Implement function hints & tab-completion.
    %   * Finish documentations.
    %   * let user edit ied from event_marker (remove & add events).
    %   * implement 3d option.
    %   * simplify wave_lb checks - remove unnecessary conditions.
    %   * implement change label / cluster via context menu on line /
    %     main_scatter.
    %   * make context menu when displaying waveform on new axe show only
    %     relevant options (listen to ev_app prop).

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%  properties %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % User Settable app components
    properties (Access = public)

        %%%%% Reduction Tab %%%%
        % general elements
        reduction_function              matlab.ui.control.DropDown        % Dimensionality reduction method: 'PCA','t-SNE','UMAP','cMDS'
        reduction_source                matlab.ui.control.DropDown        % Input data type: 'Base Waveforms','Distance Mat'
            reduction_dist_fun          matlab.ui.control.EditField       % Distance function for distance matrix calculation, e.g. 'wv2dist', 
            % or any function handle that takes in waveforms and returns a distance matrix, as a string accepted by str2func
        reduction_wv_win               matlab.ui.control.NumericEditField % Window size around waveform in ms, >0
        reduction_all_norm             matlab.ui.container.TreeNode      % If to apply Z-score normalization and alignment to waveforms before reduction, True/False as a string
        reduction_z_norm               matlab.ui.container.TreeNode      % Z-score normalize waveforms before reduction: True/False as a string
        reduction_align_norm           matlab.ui.container.TreeNode      % Try to Align waveforms to biggest sink before reduction: True/False as a string

        % tsne components
        tsne_dist_fun                  matlab.ui.control.DropDown        % t-SNE distance metric: 'euclidean','cosine',etc, or custom function as a string accepted by str2func
        tsne_exaggeration              matlab.ui.control.NumericEditField % t-SNE exaggeration factor, >1
        tsne_perplexity               matlab.ui.control.NumericEditField % t-SNE perplexity parameter, >0
        tsne_Algorithm                 matlab.ui.control.DropDown        % t-SNE algorithm: 'exact','barneshut'
        tsne_LearnRate                matlab.ui.control.NumericEditField % t-SNE learning rate, >0
        tsne_max_iter                 matlab.ui.control.NumericEditField % Maximum iterations, integer >0

        % umap components
        umap_dist_fun                  matlab.ui.control.DropDown        % UMAP distance metric: 'euclidean','cosine',etc, or custom function as a string accepted by str2func
        umap_n_neighbors              matlab.ui.control.NumericEditField % Number of neighbors, integer >0
        umap_spread                    matlab.ui.control.NumericEditField % Spread of points, >0
        umap_min_dist                 matlab.ui.control.NumericEditField % how tightly UMAP is allowed to pack points together, >0
        umap_init                      matlab.ui.control.DropDown        % Initialization: 'spectral','random'
        umap_fast_approx               matlab.ui.control.DropDown        % Use fast approximation: 'True','False'

        %%%%% Cluster Tab %%%%
        % general elements
        cluster_source                 matlab.ui.control.DropDown        % Clustering input: 'Reduced Data 2d','Distance Mat','Base Waveforms'
            cluster_dist_fun           matlab.ui.control.EditField       % Distance function for distance matrix calculation, e.g. 'wv2dist', 
            % or any function handle that takes in waveforms and returns a distance matrix, as a string accepted by str2func
        cluster_wv_win                matlab.ui.control.NumericEditField % Window size around waveform in ms, >0
        cluster_norm_tree             matlab.ui.container.CheckBoxTree
            cluster_all_norm          matlab.ui.container.TreeNode      % If to apply Z-score normalization and alignment to waveforms before clustering, True/False as a string
            cluster_z_norm            matlab.ui.container.TreeNode      % Z-score normalize waveforms before clustering: True/False as a string
            cluster_align_norm        matlab.ui.container.TreeNode      % Try to Align waveforms to biggest sink before clustering: True/False as a string
        cluster_method                 matlab.ui.control.DropDown        % Clustering method: 'Manual','by Label','kmeans','dbscan','hierarchical','None'

        % kmean components
        kmean_nClusters               matlab.ui.control.NumericEditField % Number of clusters, integer >0
        kmean_dist_fun                 matlab.ui.control.DropDown        % k-means distance: 'sqeuclidean','cityblock',etc
        kmean_start                    matlab.ui.control.DropDown        % Initialization: 'plus','cluster','sample','uniform'
        kmean_MaxIter                 matlab.ui.control.NumericEditField % Maximum iterations, integer >0
        kmean_online                   matlab.ui.control.DropDown        % Online phase: 'On','Off'
        kmean_replicates              matlab.ui.control.NumericEditField % Number of replicates, integer >0

        % dbscan components
        dbscan_epsilon                matlab.ui.control.NumericEditField % Neighborhood distance, >0
        dbscan_minpts                 matlab.ui.control.NumericEditField % Minimum points in neighborhood, integer >0
        dbscan_dist_fun                matlab.ui.control.DropDown        % DBSCAN distance: 'euclidean','cityblock',etc, or custom function as a string accepted by str2func

        % hierarchical components 
        HC_link_method                 matlab.ui.control.DropDown        % Linkage method: 'weighted','average','complete',etc
        HC_clust_critera               matlab.ui.control.DropDown        % Clustering criterion: 'inconsistent','distance'
        HC_clust_cutoff               matlab.ui.control.NumericEditField % Cutoff value for clustering, >0
        HC_clust_MaxClust             matlab.ui.control.NumericEditField % Maximum number of clusters, integer >0

        % manual clustering components
        % button is in private

        % by_label components
        label2clust_by                 matlab.ui.control.DropDown        % If clustering by label, which label to use - unrecommended to use during construction

        %%%%% Waveform Display Tab %%%%
        % general elements
        disp_wv_win                    matlab.ui.control.NumericEditField % Display window around waveform in ms, >0
        disp_all_norm                  matlab.ui.container.TreeNode      % If to apply Z-score normalization and alignment to waveforms before display, True/False as a string
        disp_z_norm                    matlab.ui.container.TreeNode      % Z-score normalize waveforms before display: True/False as a string
        disp_align_norm                matlab.ui.container.TreeNode      % Try to Align waveforms to biggest sink before display: True/False as a string
        disp_type                      matlab.ui.control.DropDown        % Display type: 'Full Cluster Average','Clusters as Labels','Overlays'
    end

    % Backbone app components
    properties (Hidden)
        %%%%% Figure Base %%%%%
        UIFigure                        matlab.ui.Figure
        main_grid                       matlab.ui.container.GridLayout
        TabGroup                        matlab.ui.container.TabGroup

        % reduction results
        reduction_axe                   matlab.ui.control.UIAxes
        main_scatter                    matlab.graphics.chart.primitive.Scatter % displaying the reduction & clustering results
        label_scatter                   matlab.graphics.chart.primitive.Scatter % under main_scatter, displaying what label each point belong to
        cluster_legend_scatter          matlab.graphics.chart.primitive.Scatter % invisible scatter that is used for creating legend
        label_legend_scatter            matlab.graphics.chart.primitive.Scatter % invisible scatter that is used for creating legend

        % context menues
        single_wv_menu                  matlab.ui.container.ContextMenu % menu for displaying only 1 waveform
            show_on_signal_new_win          matlab.ui.container.Menu
            show_on_signal_jump2loc         matlab.ui.container.Menu
            wv_on_other_disp            matlab.ui.container.Menu
        cluster_disp_menu               matlab.ui.container.ContextMenu % menu for interacting with cluster_disp_axes
            open_new_win                    matlab.ui.container.Menu

        %%%%% Reduction Tab %%%%
        reduction_manager               matlab.ui.container.Tab
        reduction_grid                  matlab.ui.container.GridLayout
        tsne_grid                       matlab.ui.container.GridLayout
        umap_grid                       matlab.ui.container.GridLayout

        % trees & buttons
        reduction_norm_tree             matlab.ui.container.CheckBoxTree
        umap_rerun_button               matlab.ui.control.Button

        %%%%% Cluster Tab %%%%
        cluster_manager                 matlab.ui.container.Tab
        cluster_grid                    matlab.ui.container.GridLayout

        % subgrids
        kmean_grid                      matlab.ui.container.GridLayout
        dbscan_grid                     matlab.ui.container.GridLayout
        HC_grid                         matlab.ui.container.GridLayout
        clust_by_label_grid             matlab.ui.container.GridLayout

        % trees & buttons
        start_marking_button            matlab.ui.control.StateButton

        % hold the ROIs when doing manual clustering
        manual_ROIs                     matlab.graphics.Graphics

        %%%%% Waveform Display Tab %%%%
        WaveformsDisplayTab             matlab.ui.container.Tab
        display_grid                    matlab.ui.container.GridLayout

        % trees & buttons
        disp_norm_tree                  matlab.ui.container.CheckBoxTree

        % hold the axes that display the waveforms
        cluster_disp_axes               matlab.graphics.Graphics

        
        %%%%% Export Tab %%%%
        ExportOptionsTab                matlab.ui.container.Tab
        export_grid                     matlab.ui.container.GridLayout

        % trees & buttons
        GenerateCodeButton              matlab.ui.control.Button
        ExportDataButton                matlab.ui.control.Button
        SaveAppButton                   matlab.ui.control.Button
    end

    % public data holders
    properties (Access = public)
        % the basic data to work on
        ied_obj IED.data

        % the data after reduction, in 2d - each col is a diffrent dim
        reduced_data

        % indecies for each row in base_data, marking what cluster it belong to
        % as a vector, with same number of rows as reduced_data
        cluster_idx

        % data labels - some conditions of intrest to place on each waveform.
        % cell label_types X 2, col 1 is label type name, col 2 is logical
        % vector - true for each waveform that answer the condition.
        wave_labels = {}

        % seed for each reduction & cluster run. Note that it does not
        % effect, which is random.
        rng_seed = 0
    end

    % private data holders
    properties (Hidden)
        % distance matrix & parameters, to avoid unneeded recalculation if
        % reduction_source is "Distance Mat"
        reduce_dist_mat = [];
        reduce_dist_parms = {};
        
        % distance matrix & parameters, to avoid unneeded recalculation if
        % cluster_source is "Distance Mat"
        cluster_dist_mat = [];
        cluster_dist_parms = {};

        % flag if to continue showing instructions for manual clustering.
        manual_instruct_show = true;
    end

    % holder for event marker displayment
    properties (Hidden)
        ev_app event_marker_gui
        app_ev_close_listener event.listener
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%  Methods %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%% Calculation methods %%%%%%%%%%
    methods (Hidden)

        function reduce_data(app)
            % perform reduction according to existing parameters in app.

            % make it reproducable in each run - in UMAP it is decided by
            % default seed, instead.
            rng(app.rng_seed)

            % display process dlg
            progress_dlg = uiprogressdlg(app.UIFigure,"Indeterminate","on","Message","Collecting data","Cancelable","off","Title","Reduce Data");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            %%%% collect data %%%%
            if app.reduction_source.Value == "Distance Mat"
                data2reduce = app.reduce_dist_mat; % note that getting dist mat only calculate if needed, see get
            elseif app.reduction_source.Value == "Base Waveforms"
                data2reduce = app.collect_wv(app.reduction_grid);
            end
            
            % deal with nan - set aside, and reenter later
            nan_rows = any(isnan(data2reduce),2);
            data2reduce(nan_rows,:) = [];

            %%%% switch reduction methods %%%%
            progress_dlg.Message = "Performing reduction";
            switch app.reduction_function.Value
                case 'PCA'
                    [~,scores] = pca(data2reduce,'Centered',false);
                    app.reduced_data = scores(:,1:2);

                case 't-SNE'
                    % convert distance function if needed, so tsne can use it
                    dist_fun = app.tsne_dist_fun.Value;
                    if contains(dist_fun,'@')
                        dist_fun = str2func(dist_fun);
                    end
                    % run tsne
                    app.reduced_data = tsne(data2reduce,...
                        'Distance',dist_fun,...
                        'Exaggeration', app.tsne_exaggeration.Value,...
                        'Perplexity',   app.tsne_perplexity.Value,...
                        'Algorithm',    app.tsne_Algorithm.Value,...
                        'LearnRate',    app.tsne_LearnRate.Value, ...
                        'Options',      struct('MaxIter',app.tsne_max_iter.Value), ...
                        'Verbose',      0);
                case 'UMAP'
                    % validate right toolbox is used. Would be nicer if
                    % UMAP would actully be a package, instead of a folder
                    % group.
                    path_shadows_solver("umap\")
                    
                    % convert distance function if needed, so tsne can use it
                    dist_fun = app.umap_dist_fun.Value;
                    if contains(dist_fun,'@')
                        dist_fun = str2func(dist_fun);
                    end
                    app.reduced_data = run_umap(data2reduce, ...
                        'metric', dist_fun,...
                        'n_neighbors',  app.umap_n_neighbors.Value, ...
                        'spread',       app.umap_spread.Value, ...
                        'min_dist',     app.umap_min_dist.Value, ...
                        'init',         app.umap_init.Value, ...
                        'fast_approximation', app.umap_fast_approx.Value == "True", ...
                        'verbose',      'none', ...
                        'randomize',    true); % change to false to use rng default instead.  
                case 'cMDS'
                    scores = cmdscale(data2reduce);
                    app.reduced_data = scores(:,1:2);
            end

            % restore nan rows - push nan into reduced results
            reduced_data_with_nans = nan(numel(nan_rows),size(app.reduced_data,2));
            reduced_data_with_nans(~nan_rows,:) = app.reduced_data;
            app.reduced_data = reduced_data_with_nans;

            %%%% Change Graphics %%%%
            progress_dlg.Message = "Updating Graphics";

            % make sure limits change to match new plot
            [app.reduction_axe.XLimMode, app.reduction_axe.YLimMode] = deal('auto');

            % note that we assume that for any reduction method, the number
            % of points is constant, and is equal to sum(app.ied_obj.accepted)
            % therefore we can update the x - y location without changing
            % clusters index
            app.main_scatter.XData = app.reduced_data(:,1);
            app.main_scatter.YData = app.reduced_data(:,2);

            close(progress_dlg)
            %%%% Recluster %%%%
            if app.cluster_source.Value == "Reduced Data 2d"
                % source data has changed - we should rerun the clustering
                % algorithem
                app.cluster_data()
            end

            % update label_scatter
            app.label_scatter.XData = app.reduced_data(:,1);
            app.label_scatter.YData = app.reduced_data(:,2);
        end


        function cluster_data(app)
            % preform clustering by existing parameters in app.
            % Change graphics to display new clusters.
            %   Note:
            %       app.cluster_idx usage assume that clusters are positive
            %       intergers without "jumps".

            % make it reproducable in each run
            rng(app.rng_seed)

            % display process dlg
            progress_dlg = uiprogressdlg(app.UIFigure,"Indeterminate","on","Message","Collecting data","Cancelable","off","Title","Cluster Data");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            %%%% collect data %%%%
            if app.cluster_source.Value == "Reduced Data 2d"
                data2cluster = app.reduced_data;
            elseif app.cluster_source.Value == "Distance Mat"
                data2cluster = app.cluster_dist_mat; % note that getting dist mat only calculate if needed, see get
            elseif app.cluster_source.Value == "Base Waveforms"
                data2cluster = app.collect_wv(app.cluster_grid);
            end

            % deal with nan - set aside, and reenter later
            nan_rows = any(isnan(data2cluster),2);
            data2cluster(nan_rows,:) = [];

            %%%% switch clustering methods %%%%
            progress_dlg.Message = "Performing clustering";
            switch app.cluster_method.Value
                case "kmeans"
                    app.cluster_idx = kmeans(data2cluster, app.kmean_nClusters.Value,...
                        'Distance',     app.kmean_dist_fun.Value, ...
                        'start',        app.kmean_start.Value, ...
                        'MaxIter',      app.kmean_MaxIter.Value, ...
                        'OnlinePhase',  app.kmean_online.Value, ...
                        'Replicates',   app.kmean_replicates.Value);

                case "dbscan"
                    % convert distance function if needed, so dbscan can use it
                    dist_fun = app.dbscan_dist_fun.Value;
                    if contains(dist_fun,'@')
                        dist_fun = str2func(dist_fun);
                    end
                    % run dbscan
                    app.cluster_idx = dbscan(data2cluster, app.dbscan_epsilon.Value, app.dbscan_minpts.Value, ...
                        "Distance", dist_fun);
                    app.cluster_idx(app.cluster_idx == -1) = max(max(app.cluster_idx),0)+1; % make outliers a seperate cluster

                case "hierarchical"
                    % can only be done on distance mat.
                    % data need to be transformed 2 vector so "linkage" will know that.
                    data2cluster = squareform(data2cluster);

                    % Since distance mat already calculated, we will
                    % perform clustering by linking & clustering seperetly.
                    link_mat = linkage(data2cluster,app.HC_link_method.Value);

                    % split behavior by if we are using "MaxClust" or not
                    if isempty(app.HC_clust_MaxClust.Value) || isinf(app.HC_clust_MaxClust.Value)
                        app.cluster_idx = cluster(link_mat, ...
                            "criterion",    app.HC_clust_critera.Value,...
                            "cutoff",       app.HC_clust_cutoff.Value);
                    else
                        app.cluster_idx = cluster(link_mat, ...
                            "maxclust",     app.HC_clust_MaxClust.Value);
                    end

                case "Manual"
                    % if a previos manual clustering exist, get it. Else,
                    % do nothing.
                    if any(isvalid(app.manual_ROIs))
                        app.StartMarkingButtonValueChanged() % now start marking is unselected, so using roi to caluculate
                    end
                case "None"
                    app.cluster_idx = ones(sum(app.ied_obj.accepted),1);
                case "by Label"
                    app.cluster_idx = categorical(app.label2clust_by.Value,unique(app.label2clust_by.Value,'stable'));
            end

            % restore nan rows - push nan into clustered results
            cluster_with_nans = nan(numel(nan_rows), 1);
            cluster_with_nans(~nan_rows,:) = app.cluster_idx;
            app.cluster_idx = cluster_with_nans;

            %%%% Change Graphics %%%%
            progress_dlg.Message = "Updating Graphics";

            % repaint dots
            warning("off","MATLAB:handle_graphics:exceptions:SceneNode")
            app.main_scatter.CData = app.cluster_idx;
            warning("on","MATLAB:handle_graphics:exceptions:SceneNode")

            % refresh legend
            app.refresh_colormap_legend()

            close(progress_dlg)
        end

    end

    %%%%%%%%%% Tab Management %%%%%%%%%%
    methods (Hidden)

        function build_reduce_tab(app)
            % remove & add elements from reduce tab according to changes in:
            % reduction_source
            % reduction_function
            %
            % the main idea is that there are base elements (identfied by
            % tag, 'base_elem', that

            %%%% create base status - everything but base are hidden %%%%

            % display ALL elements
            all_child = app.reduction_grid.Children;
            app.reduction_grid.ColumnWidth = repelem({'fit'},numel(all_child));
            % hide existing non base-elements
            child2hide = all_child(~ismember(get(all_child,'Tag'),'base_elem'));
            cols2hide = unique(arrayfun(@(x) x.Layout.Column, child2hide));
            app.reduction_grid.ColumnWidth(cols2hide) = {0};

            % make all properties that can be enabled, enabled
            enable_able = findobj(all_child,'-property','Enable');
            set(enable_able,'Enable','on')

            %%%% change other fields according to existing parameters %%%%
            % This must happen here, so the right elements will be shown

            % change other fields according to existing parameters
            app.reduction_source.Enable = 'on';
            if app.reduction_function.Value == "cMDS"
                app.reduction_source.Value = "Distance Mat";
                app.reduction_source.Enable = 'off';
            end

            %%%% Reshow elements as needed %%%%

            % collect columns of elements to reshow
            col2show = [];
            switch app.reduction_source.Value
                case "Base Waveforms"
                    % No need to do anything
                case "Distance Mat"
                    col2show = [col2show app.reduction_dist_fun.Layout.Column];
            end
            switch app.reduction_function.Value
                case "PCA"
                    % No need to do anything
                case "t-SNE"
                    col2show = [col2show app.tsne_grid.Layout.Column];
                case "UMAP"
                    col2show = [col2show app.umap_grid.Layout.Column];
                case "cMDS"
                    % No need to do anything
            end

            % reshow elements
            app.reduction_grid.ColumnWidth(unique(col2show)) = {'fit'};
        end


        function build_cluster_tab(app,event)
            % remove & add elements from cluster tab according to changes in:
            % cluster_source
            % cluster_method
            %
            % the main idea is that there are base elements (identfied by
            % tag, 'base_elem', that

            %%%% create base status - everything but base are hidden %%%%

            % display ALL elements
            all_child = app.cluster_grid.Children;
            app.cluster_grid.ColumnWidth = repelem({'fit'},numel(all_child));
            % hide existing non base-elements
            child2hide = all_child(~ismember(get(all_child,'Tag'),'base_elem'));
            cols2hide = unique(arrayfun(@(x) x.Layout.Column, child2hide));
            app.cluster_grid.ColumnWidth(cols2hide) = {0};

            % make all properties that can be enabled, enabled
            enable_able = findobj(all_child,'-property','Enable');
            set(enable_able,'Enable','on')

            %%%% change other fields & disable according to existing parameters %%%%
            % This must happen here, so the right elements will be shown

            % with focuse on "cluster_method"
            if app.cluster_method.Value == "Manual"
                % manual only works on reduced data, since this is what displayed
                app.cluster_source.Value = "Reduced Data 2d";
                app.cluster_source.Enable = 'off';
                cols2hide = [app.cluster_wv_win.Layout.Column, app.cluster_norm_tree.Layout.Column];
                app.cluster_grid.ColumnWidth(cols2hide) = {0};

            elseif app.cluster_method.Value == "hierarchical"
                % hierarchical is seperated to 3 stages: creating distance mat,
                % linkage & cluster. So source data is always distance matrix.
                app.cluster_source.Value = "Distance Mat";
                app.cluster_source.Enable = 'off';

            elseif app.cluster_method.Value == "dbscan"
                if exist("event","Var") && ... % for startup calling
                        ismember(event.Source,[app.dbscan_dist_fun, app.cluster_source]) && ...
                        ismember(event.Value,["precomputed","Distance"])
                    app.cluster_source.Value = "Distance Mat";
                    app.dbscan_dist_fun.Value = "precomputed";
                end
                if app.dbscan_dist_fun.Value == "precomputed" && app.cluster_source.Value == "Distance Mat"
                    app.cluster_source.Enable = 'off';
                end
                % in dbscan, if input is distance mat the distance function
                % must be "precomputed", and vise versa.
                % The "enable" condition on "cluster_source" is in order to
                % let user change the dist_fun & exit this state.

            elseif app.cluster_method.Value == "None"
                % disable everything but cluster_method, as they are irrelevant
                set([app.cluster_dist_fun, app.cluster_wv_win, app.cluster_norm_tree, app.cluster_source],'Enable','off')

            elseif app.cluster_method.Value == "by Label"
                % use a dropdown to select what label to work use.
                % Everything else, hide.
                cols2hide = [app.cluster_wv_win.Layout.Column, app.cluster_norm_tree.Layout.Column];
                app.cluster_grid.ColumnWidth(cols2hide) = {0};
            end

            % with focuse on specific methods argumants
            if ~isempty(app.HC_clust_MaxClust.Value) && ~isinf(app.HC_clust_MaxClust.Value)
                % hierarchical clustering can't accept both
                % "MaxClust" & "critera - cutoff".
                % Use empty state of HC_clust_MaxClust to diff between what
                % user wants.
                app.HC_clust_critera.Enable = "off";
                app.HC_clust_cutoff.Enable = "off";
            end

            %%%% Reshow elements as needed %%%%

            % collect columns of elements to reshow
            col2show = [];
            switch app.cluster_source.Value
                case "Reduced Data 2d"
                    % No need to do anything
                case "Base Waveforms"
                    % No need to do anything
                case "Distance Mat"
                    col2show = [col2show app.cluster_dist_fun.Layout.Column];
            end

            switch app.cluster_method.Value
                case "Manual"
                    col2show = [col2show app.start_marking_button.Layout.Column];
                case "kmeans"
                    col2show = [col2show app.kmean_grid.Layout.Column];
                case "dbscan"
                    col2show = [col2show app.dbscan_grid.Layout.Column];
                case "hierarchical"
                    col2show = [col2show app.HC_grid.Layout.Column];
                case "by Label"
                    col2show = [col2show app.clust_by_label_grid.Layout.Column];
            end

            % reshow elements
            app.cluster_grid.ColumnWidth(unique(col2show)) = {'fit'};
        end


        function update_wv_disp(app)
            % update the waveforms displays by requested label

            % delete old graphics
            delete(app.cluster_disp_axes)

            % collect base info
            [base_wv, t_wv] = app.collect_wv(app.display_grid);
            t_wv = t_wv*1000;
            base_wv = double(base_wv);
            if iscategorical(app.cluster_idx)
                % in case of categorical labels, the "value" of each
                % category is defined by its position in app.cluster_idx.
                % This value is what determine the color later
                clust_ids = unique(app.cluster_idx,'stable');
            else
                % non-categorical data, the value is diractly the cluster
                % index, also determining the color
                clust_ids = unique(app.cluster_idx);
            end
            clust_ids(ismissing(clust_ids)) = []; % missing clusters are not really clusters, just placeholders

            nClusts = numel(clust_ids);
            base_elem_num = numel(findobj(app.display_grid,"Tag","base_elem"));

            % if there are to many clusters - offer to skip plotting
            if nClusts > 20
                answer = uiconfirm(app.UIFigure,...
                    "More than 20 diffrent clusters. " + ...
                    "Plotting all of them may take long and and cause memory issues, continue?",...
                    "Too many clusters","Icon","warning",...
                    "Options",["Yes","No"],"CancelOption","No","DefaultOption","No");
                if answer == "No"
                    return
                end
            end
            % collect labels data
            label_data = app.disp_type.Value;
            poss_labels = unique(label_data{2},"stable");
            nLabels = numel(poss_labels);

            % collect colors so they will match the clusters in main scatter
            colors = colormap(app.reduction_axe);

            % add graphics
            app.display_grid.ColumnWidth = [repelem({'fit'},base_elem_num), repelem({'fit'},nClusts)];
            app.cluster_disp_axes = gobjects(nClusts,1);
            for iClust = nClusts:-1:1
                % add axis
                app.cluster_disp_axes(iClust) = uiaxes(app.display_grid);
                app.cluster_disp_axes(iClust).Layout.Column = iClust + base_elem_num;
                app.cluster_disp_axes(iClust).Layout.Row = [1 2];
                hold(app.cluster_disp_axes(iClust),"on")
                app.cluster_disp_axes(iClust).ContextMenu = app.cluster_disp_menu;

                %%%% How to plot %%%%

                % collect which values belong to this cluster & how it is
                % called
                in_cluster = app.cluster_idx == clust_ids(iClust);
                if isnumeric(clust_ids)
                    clust_name = sprintf('Cluster %d',clust_ids(iClust));
                else
                    clust_name = string(clust_ids(iClust));
                end

                % deal with special "Overlays" case
                if isscalar(label_data{1}) && label_data{1} == "Overlays"
                    % plot every waveform seperetly, on top of each other
                    p = plot(app.cluster_disp_axes(iClust),t_wv,base_wv(in_cluster,:)');
                    title(app.cluster_disp_axes(iClust), clust_name,'Interpreter','none')

                    % tag every line so it will have its sample position
                    pos_indx = app.ied_obj.filter_accpeted(in_cluster);
                    wv_samp = num2cell(app.ied_obj.pos(pos_indx));
                    [p.UserData] = deal(wv_samp{:});
                    wv_samp = cellfun(@num2str,wv_samp,'UniformOutput',false);
                    [p.DisplayName] = deal(wv_samp{:});

                    % add context menu to lines
                    [p.ContextMenu] = deal(app.single_wv_menu);
                    [p.ButtonDownFcn] = deal(@(src,evt) app.clust_disp_overlays_lines_ButtonDownFcn(src,evt));

                    legend(app.cluster_disp_axes(iClust),'off')
                    continue
                end

                % create a unique color for each label
                cluster_color = rgb2hsv(colors(iClust,:));
                labels_colors = repmat(cluster_color,[nLabels,1]);
                % labels_colors(:,2) = normalize(labels_colors(:,2) + linspace(-0.15,0.15,nLabels),'range');
                if nLabels ~= 1
                    labels_colors(:,2) = linspace(0.25,1,nLabels);
                end
                labels_colors = hsv2rgb(labels_colors);

                % stdplot, sepereted by labels
                for iLabel = 1:nLabels
                    in_label = ismember(label_data{2},poss_labels(iLabel));
                    wv2plot = base_wv(in_label & in_cluster,:);

                    label_line(iLabel) = IED.utils.stdshade(app.cluster_disp_axes(iClust),wv2plot,0.5,labels_colors(iLabel,:),t_wv);
                end
                legend(app.cluster_disp_axes(iClust), label_line, string(label_data{1}),'Location','best',...
                    'Interpreter','none',...
                    "ItemHitFcn",@(~,evt) app.cluster_disp_axes_legend_hide_on_click(evt),...
                    "ContextMenu",[])

                title(app.cluster_disp_axes(iClust), clust_name,'Interpreter','none')
            end

            % make all axes have same ylimits
            all_ylims = vertcat(app.cluster_disp_axes.YLim);
            [app.cluster_disp_axes.YLim] = deal([ min(all_ylims(:,1)), max(all_ylims(:,2)) ]);
        end

    end

    %%%%%%%%%% Graphics Controls %%%%%%%%%%
    methods (Hidden)

        function refresh_colormap_legend(app)
            % refresh everything that help create a "legend" - sepereting
            % diffrent elements.
            % Include:
            %   * recreate invisible scatter for legend
            %   * update colormap by number of unique colors needed
            %   * update labels datatips
            % recreate invisible scatters to match current datasets

            % remove existing invisible sactters & non scatter elements (i.e. text)
            all_scatters = app.reduction_axe.Children;
            delete(all_scatters(~ismember(all_scatters, [app.main_scatter app.label_scatter])))

            % find how many legend entries
            if iscategorical(app.cluster_idx)
                % in case of categorical labels, the "value" of each
                % category is defined by its position in app.cluster_idx.
                % This value is what determine the color later
                clust_ids = unique(app.cluster_idx,'stable');
            else
                % non-categorical data, the value is diractly the cluster
                % index, also determining the color
                clust_ids = unique(app.cluster_idx);
            end
            clust_ids(ismissing(clust_ids)) = []; % missing clusters are not really clusters, just placeholders
            nClusts = numel(clust_ids);
            nLabels = numel(app.disp_type.Value{1});

            if nLabels == 1
                % if only 1 label - do not display any labels.
                % this will cause the iLabel loop to be skipped.
                nLabels = 0;
            end
            % change colormap so clusters & labels will be distinct
            colors = distinguishable_colors(nClusts + nLabels,[1 1 1;0 0 0]);
            colormap(app.reduction_axe, colors);

            % add clusters invisible scatter
            for iClust = 1:nClusts
                if isnumeric(clust_ids)
                    clust_name = sprintf('Cluster %d',clust_ids(iClust));
                else
                    clust_name = string(clust_ids(iClust));
                end
                scatter(app.reduction_axe, nan, nan, [], colors(iClust,:),'o','filled',...
                    'DisplayName',clust_name);
            end

            % add labels invisible scatter
            for iLabel = 1:nLabels
                scatter(app.reduction_axe, nan, nan, [], colors(iLabel+nClusts,:),'o',...
                    'DisplayName',string(app.disp_type.Value{1}(iLabel)));
            end

            % label_scatter.CData is indexing rows in colormap. Number of
            % clusters determine what is the right row to index - fix to
            % match
            warning("off","MATLAB:handle_graphics:exceptions:SceneNode")
            app.label_scatter.CData = app.label_scatter.CData - min(app.label_scatter.CData) + nClusts + 1;
            warning("on","MATLAB:handle_graphics:exceptions:SceneNode")

            legend(app.reduction_axe,"Location","best",...
                'Interpreter','none',...
                "ItemHitFcn",@(~,evt) app.reduction_axe_legend_on_click(evt), ...
                "ContextMenu",[])
        end


        function add_new_cluster_menu(app, evt)
            % Let user mark a new ROI
            if evt.Button == 3
                ROI_colors = distinguishable_colors(numel(app.manual_ROIs)+1);
                app.manual_ROIs(end+1) = drawfreehand(app.reduction_axe,"Color",ROI_colors(end,:));
            end
        end

    end

    %%%%%%%%%% Save & Load %%%%%%%%
    methods
        function sobj = saveobj(app)
            % Force the object into a structure,
            % that can be constracted on load

            % colllect setable elements
            dashbord_elem = [?matlab.ui.control.DropDown, ?matlab.ui.control.NumericEditField, ?matlab.ui.control.EditField, ?matlab.ui.container.TreeNode];
            class_info = ?reduct_displayer;
            props_in_dashbord = arrayfun(@(x) ~isempty(x.Validation) && ismember(x.Validation.Class,dashbord_elem), class_info.PropertyList);
            setable_props = class_info.PropertyList(props_in_dashbord);

            % convert them into a struct by class definition
            sobj = struct();
            for iProp = 1:numel(setable_props)

                % collect prop info
                relv_class = setable_props(iProp).Validation.Class;
                prop_name = setable_props(iProp).Name;

                % choose what to save based on class
                switch relv_class
                    case ?matlab.ui.control.DropDown
                        % if editable, take current value. Else, take the
                        % relevant item by ValueIndex (to pass ItemsData)

                        if app.(prop_name).Editable
                            sobj.(prop_name) = app.(prop_name).Value;
                        else
                            if isMATLABReleaseOlderThan("R2023b")
                                value_indx= app.value_index_substetute(app.(prop_name).Value, app.(prop_name));
                                sobj.(prop_name) = app.(prop_name).Items{value_indx};
                            else
                                sobj.(prop_name) = app.(prop_name).Items{app.(prop_name).ValueIndex};
                            end
                        end

                    case {?matlab.ui.control.NumericEditField, ?matlab.ui.control.EditField}
                        % just save the value of the edit field
                        sobj.(prop_name) = app.(prop_name).Value;

                    case ?matlab.ui.container.TreeNode
                        % check if it is checked in parent tree.
                        % Save logical flag to indicate checked.

                        relv_tree = ancestor(app.(prop_name),'matlab.ui.container.CheckBoxTree');
                        sobj.(prop_name) = ismember(app.(prop_name),relv_tree.CheckedNodes);
                end
            end
            
            % save ROIs
            sobj.manual_ROIs = app.manual_ROIs;
            
            % save ied_obj, wave_labels & rng_seed
            % sobj.ied_obj = app.ied_obj.file_loc; % memory saving option
            sobj.ied_obj = app.ied_obj;
            sobj.wave_labels = [app.disp_type.Items(1:(end-2))' vertcat(app.disp_type.ItemsData{1:(end-2)})];
            sobj.wave_labels(:,1) = extractBefore(sobj.wave_labels(:,1),' seperated');
            sobj.rng_seed = app.rng_seed;

        end
    end

    methods (Static)
        function app = loadobj(sobj)
            % unpack saved struct to recreate app

            % not name value fields
            ied_obj = sobj.ied_obj;
            wave_labels = sobj.wave_labels;
            rng_seed = sobj.rng_seed;
            if isfield(sobj, 'manual_ROIs')
                manual_ROIs = sobj.manual_ROIs;
            else
                manual_ROIs = gobjects(0,0);
                sobj.manual_ROIs = gobjects(0,0);
            end
            sobj = rmfield(sobj,{'ied_obj', 'wave_labels', 'rng_seed', 'manual_ROIs'});

            % collect name value fields
            name_vals = namedargs2cell(sobj);

            % create app
            try
                app = reduct_displayer(ied_obj, wave_labels, name_vals{:});
            catch err
                if err.identifier == "reduct_displayer:wave_labels"
                    childs = allchild(groot);
                    delete(childs(1)) % close failed attempt

                    warning('Failed to init due to wave_labels. Setting all detection to accepted and try again')
                    ied_obj.accepted = true(size(ied_obj.accepted));
                    app = reduct_displayer(ied_obj, wave_labels, name_vals{:});
                end
            end
            % insert manual_ROIs & reset seed
            app.manual_ROIs = manual_ROIs;
            app.rng_seed = rng_seed;
        end
    end

    %%%%%%%%%% Set & Get for Properties %%%%%%%%

    methods
        function set.rng_seed(app,val)
            % refresh graphics after setting new seed

            % set new seed
            app.rng_seed = val;

            % refresh graphics
            app.reduce_data();
            app.cluster_data();
        end


        function val = get.reduce_dist_mat(app)
            % when taking out distance mat,
            % check if there was any change to parameters, if so, recalculate.
            % Else, use the already calculated matrix.

            checked_nodes = app.reduction_norm_tree.CheckedNodes;
            dist_fun = app.reduction_dist_fun.Value;
            wv_win = app.reduction_wv_win.Value;
            if isempty(app.reduce_dist_parms) || ... % inital case
                    any(~ismember(checked_nodes, app.reduce_dist_parms{1})) || ... % no new checks
                    ~all(ismember(app.reduce_dist_parms{1}, checked_nodes)) || ... % no new unchecks
                    ~strcmp(dist_fun, app.reduce_dist_parms{2}) || ... % same function
                    ~ismembertol(wv_win,app.reduce_dist_parms{3}) % same window size

                % save new parameters
                app.reduce_dist_parms{1} = checked_nodes;
                app.reduce_dist_parms{2} = dist_fun;
                app.reduce_dist_parms{3} = wv_win;

                % calculate new distance matrix
                base_wv = app.collect_wv(app.reduction_grid);
                dist_fun = str2func(dist_fun);
                app.reduce_dist_mat = round(dist_fun(base_wv),6); % deal with double percision error
            end

            % return new value
            val = app.reduce_dist_mat;
        end


        function val = get.cluster_dist_mat(app)
            % when taking out distance mat,
            % check if there was any change to parameters, if so, recalculate.
            % Else, use the already calculated matrix.

            checked_nodes = app.cluster_norm_tree.CheckedNodes;
            dist_fun = app.cluster_dist_fun.Value;
            wv_win = app.cluster_wv_win.Value;
            if isempty(app.cluster_dist_parms) || ... % inital case
                    any(~ismember(checked_nodes, app.cluster_dist_parms{1})) || ... % no new checks
                    ~all(ismember(app.cluster_dist_parms{1}, checked_nodes)) || ... % no new unchecks
                    ~strcmp(dist_fun, app.cluster_dist_parms{2}) || ... % same function
                    ~ismembertol(wv_win,app.cluster_dist_parms{3}) % same window size

                % save new parameters
                app.cluster_dist_parms{1} = checked_nodes;
                app.cluster_dist_parms{2} = dist_fun;
                app.cluster_dist_parms{3} = wv_win;

                % calculate new distance matrix
                base_wv = app.collect_wv(app.cluster_grid);
                dist_fun = str2func(dist_fun);
                app.cluster_dist_mat = round(dist_fun(base_wv),6); % deal with double percision error
            end

            % return new value
            val = app.cluster_dist_mat;
        end
    end


    %%%%%%%%%% Helper function %%%%%%%%
    methods (Access = public)
        function [base_wv,t_stamps] = collect_wv(app, parent_grid)
            % this function collect the waveforms from the ied object and
            % normilaise them as requested.
            % The decision acording to what to normilaze them is done by
            % checking the relevent elements, that are children of
            % parent_grid.

            % collect elements that are relevant to reduction
            switch parent_grid
                case app.reduction_grid
                    z_node = app.reduction_z_norm;
                    align_node = app.reduction_align_norm;
                    wv_win = app.reduction_wv_win.Value;
                    checked_nodes = app.reduction_norm_tree.CheckedNodes;
                case app.cluster_grid
                    z_node = app.cluster_z_norm;
                    align_node = app.cluster_align_norm;
                    wv_win = app.cluster_wv_win.Value;
                    checked_nodes = app.cluster_norm_tree.CheckedNodes;
                case app.display_grid
                    z_node = app.disp_z_norm;
                    align_node = app.disp_align_norm;
                    wv_win = app.disp_wv_win.Value;
                    checked_nodes = app.disp_norm_tree.CheckedNodes;
            end

            [base_wv, t_stamps]= app.ied_obj.extract_discharges(wv_win./1000);
            base_wv = double(base_wv);

            % normilased waveforms
            if ismember(checked_nodes, z_node)
                base_wv = zscore(base_wv,1,2);
            end
            if ismember(checked_nodes, align_node)
                base_wv = align_waveforms(base_wv);
                base_wv(isnan(base_wv)) = 0;
            end
        end
        
        function add_label_types(app, wave_labels)
            % add new label types to the dropdown menu, so they can be displayed on the reduction & waveforms plots.
            % wave_labels is a cell matrix with 3 columns:
            % 1) label type name
            % 2) meaning of each label (by order of appearance)
            % 3) which label each pos (event in the ied object) relate to
            %
            % In case the label type already exist (by name), it will be overwritten. It is not recommended to overwrite
            % Default types 'Full Cluster Average', 'Clusters as Labels', 'Overlays' or 'None', as results may be cause unexpected behavior.
            %
            % See also: validate_labels

            
            % validate new labels
            app.validate_labels(wave_labels)
            
            % make 2 & 3rd cols a col vec
            wave_labels(:,2) = cellfun(@(x) x(:), wave_labels(:,2), 'UniformOutput', false);
            wave_labels(:,3) = cellfun(@(x) x(:), wave_labels(:,3), 'UniformOutput', false);

            % detect overwrites
            old_labels_types = extractBefore(app.disp_type.Items, " seperated");
            old_labels_types(cellfun(@isempty, old_labels_types)) = [];
            labels2overwrite = strcmp(old_labels_types,wave_labels(:,1));
            app.disp_type.Items(labels2overwrite) = [];
            app.disp_type.ItemsData(labels2overwrite) = [];

            % push into dropdown menu
            app.disp_type.Items = [(wave_labels(:,1)+" seperated")', app.disp_type.Items];
            app.disp_type.ItemsData = [num2cell(wave_labels(:,2:3),2)', app.disp_type.ItemsData];

            app.label2clust_by.Items = [(wave_labels(:,1)+" seperated")', app.label2clust_by.Items];
            app.label2clust_by.ItemsData = [wave_labels(:,3)', app.label2clust_by.ItemsData];
        end

        function validate_labels(app, wave_labels)
            % validate new labels.
            % wave_labels is a cell matrix with 3 columns:
            % 1) label type name (cellstr)
            % 2) meaning of each label (by order of appearance) (with same number of elements as there are unique labels in col 3)
            % 3) which label each pos (event in the ied object) relate to (vector with the same number of elements as in ied.pos(ied.accepted))

            % validate labels structure
            if ~iscell(wave_labels) || ~(ismatrix(wave_labels) && (size(wave_labels,2) == 3))

                error('reduct_displayer:wave_labels',...
                    ['wave_labels must be cell matrix with 3 columns:\n' ...
                    '1) label type name\n2) which label each pos relate to\n3) each label meaning'])

            elseif ~iscellstr(wave_labels(:,1))

                error('reduct_displayer:wave_labels','wave_labels 1st col must be cellstr, how to call this label type')

            elseif  (~all( ...
                    cellfun(@(x,y) numel(x)==numel(unique(y)),...
                    wave_labels(:,2), wave_labels(:,3)) ...
                    ))

                error('reduct_displayer:wave_labels', ...
                    ['wave_labels 2st col must have same numebr of elements as there are unique labels in col 3;' ...
                    'what each label means (by order of appearance).\n' ...
                    'Recommended creating via unique(labels,"stable")'])

            elseif  ~all( cellfun(@(x) isvector(x) && numel(x)==sum(app.ied_obj.accepted), wave_labels(:,3)) )

                error('reduct_displayer:wave_labels',...
                    ['wave_labels 3st col must be a vector with the same number of elements as in ied.pos(ied.accepted),' ...
                    'which label each event belongs to'])

            end
        end
    end

    methods (Hidden)
        function [value_indx] = value_index_substetute(~, val2search, dropdown_obj)
            % replase value_indx with manual search for older versions
            
            % make serach constant
            if ~iscell(val2search)
                val2search = {val2search};
            end
            if ~isempty(dropdown_obj.ItemsData) && ~iscell(dropdown_obj.ItemsData{1})
                dropdown_obj.ItemsData = num2cell(dropdown_obj.ItemsData);
            end

            if isempty(dropdown_obj.ItemsData)
                % use simple item comparison instead of item data.
                % this assume that val2search is a valid dropdown_obj.Items
                [~,value_indx] = ismember(val2search,dropdown_obj.Items);
            
            else
                % match cells of dropdown_obj.ItemsData to find the
                % matching index.
                for iItem = 1:numel(dropdown_obj.Items)
                    % match all cells
                    for iCell = numel(val2search):-1:1
                        if isempty(val2search{iCell})
                            match_idx(iCell) = isempty(dropdown_obj.ItemsData{iItem}{iCell});
                        else
                            if iscellstr(val2search{iCell})
                                match_idx(iCell) = all(string(dropdown_obj.ItemsData{iItem}{iCell}) == string(val2search{iCell}));
                            else
                                match_idx(iCell) = all(dropdown_obj.ItemsData{iItem}{iCell} == val2search{iCell});
                            end
                            
                        end
                    end
                    
                    % return if fully matched
                    if all(match_idx)
                        value_indx = iItem;
                        break
                    end
                end
            end
        end
    
    end

    %%%%%%%%%% Callbacks that handle component events %%%%%%%%%%
    methods (Hidden)

        % Code that executes after component creation
        function startupFcn(app, ied_obj, wave_labels, names, vals)
            % insert data into the app - call reduction method,

            arguments
                app
                ied_obj IED.data
                wave_labels = {};
            end

            arguments (Repeating)
                names char {mustBeTextScalar}
                vals
            end


            %%%% write in properites %%%%

            %%% ied obj
            app.ied_obj = ied_obj;

            %%% wave_labels

            % create defult waveforms disp_type ItemData & label2clust_by
            app.disp_type.ItemsData = {...
                { {'Cluster'}, true(sum(ied_obj.accepted),1) }, ... % Full Cluster Average case
                { {''}, []}, ... % Clusters as Labels case
                { "Overlays", [] }... % Overlays case
                };
            app.disp_type.Items = {'Full Cluster Average', 'Clusters as Labels', 'Overlays'};
            app.label2clust_by.ItemsData = {...
                true(sum(ied_obj.accepted),1) % None case
                };
            % validate labels & add them in
            if ~isempty(wave_labels)
                app.add_label_types(wave_labels)
            end

            %%% Override default dashbord values with user requests

            % validate only the right elements are refernced %

            % colllect setable elements
            dashbord_elem = [?matlab.ui.control.DropDown, ?matlab.ui.control.NumericEditField, ?matlab.ui.control.EditField, ?matlab.ui.container.TreeNode];
            class_info = ?reduct_displayer;
            props_in_dashbord = arrayfun(@(x) ~isempty(x.Validation) && ismember(x.Validation.Class,dashbord_elem), class_info.PropertyList);
            setable_props = class_info.PropertyList(props_in_dashbord);

            % find which elements are to be set
            [~,props2set] = ismember(names, {setable_props.Name});
            if any(props2set == 0)
                problem_names = join(names(props2set == 0),', ');
                error("Unrecognized name-value arguments: {%s}",problem_names{:})
            end

            % place value in the right place %
            for iName = 1:numel(props2set)

                % collect prop_data
                relv_class = setable_props(props2set(iName)).Validation.Class;
                val2use = vals{iName};
                if isa(val2use,"function_handle")
                    val2use = func2str(val2use);
                end
                prop_name = setable_props(props2set(iName)).Name;

                % set value by class
                switch relv_class
                    case ?matlab.ui.control.DropDown
                        % only accept scalar texts
                        validateattributes(val2use,{'char','string'},{'scalartext'},'reduct_displayer',prop_name)

                        % find pos of element in drop down items, than set to that pos
                        drop_down_items = app.(prop_name).Items;
                        [~,val_pos] = ismember(val2use, drop_down_items);

                        % if the value exist, set to it, else error
                        if val_pos ~= 0
                            if isMATLABReleaseOlderThan("R2023b")
                                if isempty(app.(prop_name).ItemsData)
                                    app.(prop_name).Value = val2use;
                                else
                                    app.(prop_name).Value = app.(prop_name).ItemsData{val_pos};
                                end
                            else
                                app.(prop_name).ValueIndex = val_pos;
                            end
                        else
                            if app.(prop_name).Editable
                                app.(prop_name).Value = val2use;
                            else
                                error('%s must be one of: %s',prop_name, string(join(drop_down_items,', ')) )
                            end
                        end

                    case ?matlab.ui.control.NumericEditField

                        if isempty(val2use) && (isMATLABReleaseOlderThan('R2023a') || app.(prop_name).AllowEmpty)
                            if isMATLABReleaseOlderThan('R2023a')
                                app.(prop_name).Value = inf;
                            else
                                app.(prop_name).Value = val2use;
                            end
                            continue
                        else
                            % only accept scalar numerics
                            validateattributes(val2use,{'numeric'},{'scalar'},'reduct_displayer',prop_name)
                        end

                        % validate that value is in limits, if not write it in
                        field_lims = app.(prop_name).Limits;
                        if val2use < field_lims(1) || val2use > field_lims(2)
                            warning("In %s, value must be between [%s]. Skipping",prop_name,join(string(field_lims),', '))
                        else
                            app.(prop_name).Value = val2use;
                        end

                    case ?matlab.ui.control.EditField
                        % only accept scalar texts
                        validateattributes(val2use,{'char','string'},{'scalartext'},'reduct_displayer',prop_name)

                        % write in
                        app.(prop_name).Value = val2use;

                    case ?matlab.ui.container.TreeNode
                        % only accept scalar texts
                        validateattributes(val2use,{'numeric','logical'},{'binary'},'reduct_displayer',prop_name)

                        % find relevant tree
                        relv_tree = ancestor(app.(prop_name),'matlab.ui.container.CheckBoxTree');

                        % check / uncheck box
                        if val2use
                            relv_tree.CheckedNodes = [relv_tree.CheckedNodes app.(prop_name)];
                        else
                            relv_tree.CheckedNodes(ismember(relv_tree.CheckedNodes, app.(prop_name))) = [];
                        end

                end
            end

            %%% Preset graphics  %%%
            % create empty main scatter, so other functions can use it
            hold(app.reduction_axe,"on")

            % create main scatter, and its proparties
            app.main_scatter = scatter(app.reduction_axe,nan,nan,36,"filled","o","MarkerEdgeColor",'k','LineWidth',1);
            app.main_scatter.ContextMenu = app.single_wv_menu;
            app.main_scatter.ButtonDownFcn = @(src,evt) app.main_scatter_ButtonDownFcn(src,evt);
            app.main_scatter.Annotation.LegendInformation.IconDisplayStyle = "off"; % prevent legend
            app.cluster_idx = ones(sum(ied_obj.accepted),1);
            

            % create empty label scatter
            app.label_scatter = scatter(app.reduction_axe,...
                nan(sum(ied_obj.accepted),1),nan(sum(ied_obj.accepted),1),...
                108,'filled',"o","PickableParts","none");
            uistack(app.label_scatter,"bottom")
            [~,~, app.label_scatter.CData] = unique(app.disp_type.Value{2},"stable");
            if numel(unique(app.disp_type.Value{2})) == 1
                % if current disp_type include a single label for all
                % (i.e. "Full Cluster Average"), do not show the underlining labels.
                app.label_scatter.Visible = "off";
            end
            app.label_scatter.Annotation.LegendInformation.IconDisplayStyle = "off"; % prevent legend
            
            % create informative datatips for main_scatter
            app.main_scatter.DataTipTemplate.DataTipRows(end-1) = []; % remove "size" datatip
            app.main_scatter.DataTipTemplate.DataTipRows(end).Label = "Cluster Number"; % change color to cluster
            label_datatip = dataTipTextRow('Label', app.disp_type.Value{1}(app.label_scatter.CData));
            app.main_scatter.DataTipTemplate.DataTipRows(end+1) = label_datatip;
            sample_datatip = dataTipTextRow('Detection Position [sample]', app.ied_obj.pos(app.ied_obj.accepted),'%d');
            app.main_scatter.DataTipTemplate.DataTipRows(end+1) = sample_datatip;

            %%% Calculate reduction & clusters %%%
            app.reduce_data();
            if app.cluster_source.Value ~= "Reduced Data 2d"
                app.cluster_data();
            end

            %%% Build tabs %%%
            app.build_reduce_tab()
            app.build_cluster_tab()
        end


        %%%%% Clustering Tab %%%%

        % Callback function: HC_clust_MaxClust, HC_clust_critera,
        % ...and 16 other components
        function refresh_clust(app, event)
            % Rebuilt cluster tab & recluster.
            % Order is important, as "build_cluster_tab" affect clustering
            % parameters.

            app.build_cluster_tab(event);
            app.cluster_data();
        end

        
        % Value changed function: StartMarkingButton
        function StartMarkingButtonValueChanged(app)
            % let user mark clusters manually.
            % Note: this only work in 2d.

            value = app.start_marking_button.Value;

            % display instructions
            if value && app.manual_instruct_show
                answer = uiconfirm(app.UIFigure,...
                    ["Use right click anywhere on the axis background to start drawing a shape. " ...
                    "Hold the mouse till finish marking a cluster. " ...
                    "After relesing, you can edit / delete the existing shape (right click it), or add more clusters using right click on the axis. " ...
                    "When finished, deselect the 'Start Marking' button, and the clusters will be reevaluate. " ...
                    "Cancle by marking 0 clusters."...
                    "Note:" "1) In case of overlaps, the top ROI is used.","2) Points on the edge of the shape count as in it"],"Instructions",...
                    "Options",["Got it","Do not show again"],...
                    "CancelOption","Got it","DefaultOption","Got it","Icon","info");
                if answer == "Do not show again"
                    app.manual_instruct_show = false;
                end
            end

            if value
                % delete existing ROIs saved in memory
                delete(app.manual_ROIs)
                app.manual_ROIs = gobjects(0);

                % add right click create ROI to axes
                % app.reduction_axe.ContextMenu = app.manual_cluster_axe_menu;
                app.reduction_axe.ButtonDownFcn = @(~,evt) app.add_new_cluster_menu(evt);
            else
                % remove right click create ROI from axes
                % app.reduction_axe.ContextMenu = gobjects(0);
                app.reduction_axe.ButtonDownFcn = '';


                % if no ROIs - cancle operation (don't change anything)
                if isempty(app.manual_ROIs)
                    return
                end

                % make ROIs invisible
                [app.manual_ROIs.Visible] = deal("off");
                [app.manual_ROIs.HandleVisibility] = deal("off");

                % remove all existing cluster idx
                app.cluster_idx = nan(size(app.cluster_idx));

                % convert all ROIs to clusters
                % note that order is important - later drawn ROIs take
                % precedence over eariler
                xq = app.reduced_data(:,1);
                yq = app.reduced_data(:,2);
                for iROI = 1:numel(app.manual_ROIs)
                    % find what points are on / in ROI
                    ploy_pos = app.manual_ROIs(iROI).Position;
                    [p_in,p_on] = inpolygon(xq, yq, ploy_pos(:,1), ploy_pos(:,2));

                    % change those points to new cluster index
                    app.cluster_idx(p_in | p_on) = iROI;
                end



                % place all points not in a cluster in a single cluster
                app.cluster_idx(isnan(app.cluster_idx)) = numel(app.manual_ROIs) + 1;

                % an empty ROI will create a "jump" in idx, causing weird
                % results. fix that.
                [~,~,app.cluster_idx] = unique(app.cluster_idx,"stable");

                % repaint dots
                app.main_scatter.CData = app.cluster_idx;

                % refresh legend
                app.refresh_colormap_legend()
            end
        end


        %%%%% Reduction Tab %%%%

        % Callback function: reduction_dist_fun, reduction_function,
        % ...and 8 other components
        function refresh_reduction(app, ~)
            % Rebuilt reduction tab & reduce.
            % Order is important, as "build_reduce_tab" affect reduction
            % parameters.
            % If app.cluster_source.Value is "Reduced Data 2d", reduce_data
            % will also cause reclustering.

            app.build_reduce_tab();
            app.reduce_data();
        end


        %%%%% Waveform Display Tab %%%%

        % Value changed function: disp_wv_win, Selection change function: TabGroup
        function refresh_wv_disp_tab(app)
            % replot waveforms on cluster_disp_axes

            % plot waveforms only when waveform tab is open
            selectedTab = app.TabGroup.SelectedTab;
            if selectedTab == app.WaveformsDisplayTab
                app.update_wv_disp()
            end

        end


        % Value changed function: disp_type
        function disp_typeValueChanged(app)
            % switch display to a new label setup.
            % Note that value{1} is each label meaning, value{2} is what
            % label match each waveform.

            % collect data
            label_data = app.disp_type.Value;
            if isMATLABReleaseOlderThan("R2023b")
                [value_indx] = app.value_index_substetute(label_data, app.disp_type);
            else
                value_indx = app.disp_type.ValueIndex;
            end
            
            value_title = app.disp_type.Items(value_indx);

            % deal with special "Overlays" case
            if value_title == "Overlays"
                app.label_scatter.Visible = "off";
                app.refresh_colormap_legend()
                app.update_wv_disp()
                return
            end

            % deal with special "Clusters as Labels" case
            if value_title == "Clusters as Labels"
                % switch item data to current clusters
                [app.disp_type.ItemsData{value_indx}{2}, label_data{2}] = deal(app.cluster_idx);
                [app.disp_type.ItemsData{value_indx}{1}, label_data{1}] = deal(compose({'Cluster %d'},unique(app.cluster_idx,"stable")'));
            end

            % update label_scatter cdata to match new labels
            [~,~, app.label_scatter.CData] = unique(label_data{2},'stable');

            % update datatips for new label
            tip2edit = strcmp('Label',{app.main_scatter.DataTipTemplate.DataTipRows.Label});
            app.main_scatter.DataTipTemplate.DataTipRows(tip2edit).Value = label_data{1}(app.label_scatter.CData);

            % if only 1 label, hide label_scatter
            if numel(label_data{1}) == 1
                app.label_scatter.Visible = "off";
            else
                app.label_scatter.Visible = "on";
            end

            % refresh legend
            app.refresh_colormap_legend()

            % refresh clusters waveforms display
            app.update_wv_disp()
        end


        %%%%% Export Tab %%%%
        
        % Button pushed function: SaveAppButton
        function SaveAppButtonPushed(app)
            % ask user where to  save, than save results
            [file,path] = uiputfile('.mat','Choose where to save the app');
            if ~isnumeric(file) % file is numeric only if user cancled
                save(fullfile(path,file),"app","-v7.3")
            end

        end


        % Button pushed function: ExportDataButton
        function ExportDataButtonPushed(app)
            % ask user where to  save, than save results
            [file,path] = uiputfile('.mat','Choose where to save the result data');
            if ~isnumeric(file) % file is numeric only if user cancled
                reduced_data = app.reduced_data; %#ok<PROP> % I save it outside of app, this is a usefull name
                cluster_idx = app.cluster_idx; %#ok<PROP> % I save it outside of app, this is a usefull name
                wave_lb = [app.disp_type.Items' vertcat(app.disp_type.ItemsData{:})];
                save(fullfile(path,file),"reduced_data","cluster_idx","wave_lb")
            end

        end


        % Button pushed function: GenerateCodeButton
        function GenerateCodeButtonPushed(app)
            % create in editor code to match what the app does
            uialert(app.UIFigure,'Code is incomplete! Aborting','Incomplete code','Icon','error')
            error('Code is incomplete! Aborting')

            % open file editor
            file_object = matlab.desktop.editor.newDocument();

            % write intro
            generation_time = datetime("today","Format","uuMMdd HH:mm:ss");
            file_object.Text = sprintf("%% Code generation from 'reduct_displayer' on %s\n",generation_time);

            %%% Data about code version

            % move to class defenition
            current_path = pwd;
            cd(fileparts(which("reduct_displayer")))

            % collect HEAD, and if any changes were done
            [cmd_status,git_HEAD] = system("git rev-parse --short HEAD");
            if cmd_status ~= 0
                git_HEAD = 'HEAD colllaction failed';
            else
                git_HEAD(end) = []; % remove final new line
            end
            [cmd_status, git_changes] = system('git status --porcelain');
            if cmd_status ~= 0
                git_changes = ', Status colllaction failed';
            elseif ~isempty(git_changes)
                git_changes = ', However, folder was changed since last commit';
            end

            % return to user loc
            cd(current_path)

            file_object.Text = file_object.Text + sprintf("Git commit head: %s%s\n",git_HEAD,git_changes) + newline;
            file_object.Text = file_object.Text + ...
                "%%%%%%%% Inset Animal data & Comments here %%%%%%%%" + newline;
            %%% Load data process

            % create section header
            file_object.Text = file_object.Text + ...
                "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" + newline +...
                "%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%" + newline + ...
                "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" + newline;

            % colllect ied_object
            file_object.Text = file_object.Text + "% Load IED object" + newline + ...
                sprintf("load('%s')",app.ied_obj.file_loc) + newline;

            %




        end


        %%%%% Scatters, Lines  & their legends %%%%

        % Button Down function of waveforms lines when "disp_type" is "Overlays"
        function clust_disp_overlays_lines_ButtonDownFcn(app,src,evt)
            relevant_fig = ancestor(src,{'figure','uifigure'});
            switch relevant_fig.SelectionType
                case 'alt'
                    app.show_on_signal_new_win.MenuSelectedFcn = @(~,~) app.show_on_signal_new_win_MenuSelected(src.UserData);
                    app.show_on_signal_jump2loc.MenuSelectedFcn = @(~,~) app.show_on_signal_jump2loc_MenuSelected(src.UserData);
                    app.wv_on_other_disp.MenuSelectedFcn = @(~,~) app.wv_on_other_disp_MenuSelected(src.UserData,src);
                    app.wv_on_other_disp.Visible = "on";
                case 'normal'
                    if isempty(app.main_scatter.ZData)
                        % 2d graph
                        datatip(src,evt.IntersectionPoint(1),evt.IntersectionPoint(2));
                    else
                        % 3d graph
                        datatip(src,evt.IntersectionPoint(1),evt.IntersectionPoint(2),evt.IntersectionPoint(3));
                    end
            end
        end


        % Button Down function: main_scatter
        function main_scatter_ButtonDownFcn(app,~,evt)

            switch app.UIFigure.SelectionType
                case 'alt' % right mouse click
                    % insert clicked point into menues functions

                    % collect all data points
                    x_data = app.main_scatter.XData;
                    y_data = app.main_scatter.YData;
                    z_data = app.main_scatter.ZData;
                    if isempty(z_data)
                        % this code always assume 3d - in case of 2d, just use 0.
                        z_data = zeros(size(y_data));
                    end
                    data_points = [x_data; y_data; z_data]';

                    % find the sample relating to this point in the original signal
                    [~,p_idx] = min(pdist2(data_points, evt.IntersectionPoint));
                    all_accepted = find(app.ied_obj.accepted);
                    relevant_sample = app.ied_obj.pos(all_accepted(p_idx));

                    % pass the sample info to uimenu function
                    app.show_on_signal_new_win.MenuSelectedFcn = @(~,~) app.show_on_signal_new_win_MenuSelected(relevant_sample);
                    app.show_on_signal_jump2loc.MenuSelectedFcn = @(~,~) app.show_on_signal_jump2loc_MenuSelected(relevant_sample);
                    label_data = app.disp_type.Value;
                    if isscalar(label_data{1}) && label_data{1} == "Overlays"
                        app.wv_on_other_disp.Visible = "on";
                        app.wv_on_other_disp.MenuSelectedFcn = @(~,~) app.wv_on_other_disp_MenuSelected(relevant_sample,app.main_scatter);
                    else
                        app.wv_on_other_disp.Visible = "off";
                    end

                    % open ui context menu
                    % cp = app.UIFigure.CurrentPoint;
                    % open(app.main_scatter_menu,cp)
                case 'normal' % left click
                    if isempty(app.main_scatter.ZData)
                        % 2d graph
                        datatip(app.main_scatter,evt.IntersectionPoint(1),evt.IntersectionPoint(2));
                    else
                        % 3d graph
                        datatip(app.main_scatter,evt.IntersectionPoint(1),evt.IntersectionPoint(2),evt.IntersectionPoint(3));
                    end
            end
        end

        
        % Item hit function: legend of reduction_axe
        function reduction_axe_legend_on_click(app,evt)
            % Managed item click interactions - hide element on
            % left click, change color on right click

            switch evt.SelectionType
                case 'normal'
                    % toggle visability of clicked element

                    % force limits not to change
                    [app.reduction_axe.XLimMode, app.reduction_axe.YLimMode] = deal('manual');

                    % switch visability of clicked item
                    evt.Peer.Visible = ~evt.Peer.Visible;

                    % make sure only what need to be viable is seen.
                    % every point that match any hidden category will be given value NaN.
                    % Identify points to remove by colormap (see also refresh_colormap_legend)

                    % reset all points to show
                    app.main_scatter.XData = app.reduced_data(:,1);
                    app.label_scatter.XData = app.reduced_data(:,1);

                    % collect all legend-hidden invisable scatters colors
                    hidden_inv_scatter = findobj(app.reduction_axe,'Visible','off','Type','scatter');
                    hidden_inv_scatter(ismember(hidden_inv_scatter,app.label_scatter)) = []; % When only 1 label, this scatter is invisable as well
                    if isempty(hidden_inv_scatter)
                        % no need to hide anything
                    else
                        hidden_colors = vertcat(hidden_inv_scatter.CData);
                        c_map = colormap(app.reduction_axe);
                        [~,c_idx] = ismembertol(hidden_colors,c_map,"ByRows",true);

                        % turn points that point to the same colors to NaN
                        points2remove = any(ismember([double(app.main_scatter.CData), double(app.label_scatter.CData)],c_idx),2);
                        app.main_scatter.XData(points2remove) = nan;
                        app.label_scatter.XData(points2remove) = nan;
                    end
                case 'alt'
                    % let user pick a new color for clicked element
                    current_color = evt.Peer.CData;
                    new_color = uisetcolor(current_color,evt.Peer.DisplayName + " color");
                    c_map = colormap(app.reduction_axe);
                    [~,c_idx] = ismembertol(current_color,c_map,"ByRows",true);
                    c_map(c_idx,:) = new_color;
                    colormap(app.reduction_axe,c_map)
                    evt.Peer.CData = new_color;
            end
        end


        % Item hit function: legend of cluster_disp_axes
        function cluster_disp_axes_legend_hide_on_click(~,evt)
            % toggle visability of clicked element (and its shade)

            switch evt.SelectionType
                case 'normal'

                    % find axe and what elements to work on
                    relv_ax = evt.Peer.Parent;
                    [~,line_loc] = ismember(evt.Peer,relv_ax.Children);

                    % toggle visibility - shade always follow line loc
                    [relv_ax.Children([line_loc, line_loc+1]).Visible] = ...
                        deal(~evt.Peer.Visible);
                case 'alt'
                    % let user pick a new color for clicked element
                    current_color = evt.Peer.Color;
                    new_color = uisetcolor(current_color, evt.Peer.DisplayName + " color");

                    % find the corresponding patch - 1 after the clicked line in axe children
                    line_pos = find(evt.Peer.Parent.Children == evt.Peer);
                    relev_patch = evt.Peer.Parent.Children(line_pos+1);

                    % change colors
                    evt.Peer.Color = new_color;
                    relev_patch.FaceColor = new_color;
            end
        end
    

        %%%%% Context Menues  %%%%
        % Menu Selected function: show_on_signal_new_win
        function show_on_signal_new_win_MenuSelected(app,relevant_sample)
            % open process dlg
            % display process dlg
            progress_dlg = uiprogressdlg(app.UIFigure,"Indeterminate","on","Message","Preparing Display","Cancelable","off","Title","Opening signal view");
            cl = onCleanup(@() close(progress_dlg)); % remove progress_dlg block on any form of exit (i.e. error)

            % convert all accepted detections to events
            detect_pos = app.ied_obj.pos(app.ied_obj.accepted);
            detect_t = detect_pos./app.ied_obj.fs;
            event_win = app.disp_wv_win.Value/1000;
            event_t = (detect_t(:) + [-event_win/2 event_win/2]);%'; % start-end of each event, each event is a col

            % collect labels info
            label_data = app.disp_type.Value;
            poss_labels = unique(label_data{2},"stable");
            if iscategorical(app.cluster_idx)
                poss_clusters = unique(app.cluster_idx,"stable");
            else
                poss_clusters = unique(app.cluster_idx);
            end
            poss_clusters(ismissing(poss_clusters)) = []; % missing clusters are not really clusters, just placeholders
            % all_colors = colormap(app.reduction_axe);

            % seperate events by label & clusters, into event struct for event_marker
            eve_dict = cell_dict(string.empty(), {});
            nClust = numel(poss_clusters);
            nLabel = numel(poss_labels);
            % event_colors = nan(nClust + nLabel,3);
            for iCluster = 1:nClust
                relv_events = app.cluster_idx == poss_clusters(iCluster);
                eve_dict{"Cluster "+iCluster} = event_t(relv_events,:);
                % event_colors(iCluster,:) = all_colors(iCluster,:);
            end
            if nLabel > 1
                for iLabel = 1:nLabel
                    relv_events = ismember(label_data{2},poss_labels(iLabel));
                    if ~isKey(eve_dict, label_data{1}(iLabel))
                        eve_dict{label_data{1}(iLabel)} = event_t(relv_events,:);
                    else
                        new_name = matlab.lang.makeUniqueStrings([label_data{1}{iLabel}; eve_dict.keys], 1); % we only want to change & use the first element
                        eve_dict{new_name(1)} = event_t(relv_events,:);
                    end
                    % event_colors(nClust + iLabel,:) = all_colors(nClust + iLabel,:);
                end
            end

            % open event_marker
            % ev_app_win = multi_window_lock(data_tab.data{:}, ...
            %     "fs", data_tab.fs, "labels", data_tab.Properties.RowNames, ...
            %     "init_time", init_time, "events_time", eve_dict);
            ev_app_win = reduct_displayer_IED_disp(...
                "reduct_app", app, "ied", app.ied_obj, ...
                "pos2start",  find(app.ied_obj.pos == relevant_sample), ...
                "temp_events", eve_dict, ...
                "only_init_accepted", true);

            % hold 1st opened app, for "time jumping" (instead of opening a new window every time).
            % if window was closed, remove the option to "time jump"
            if isempty(app.ev_app) || ~isvalid(app.ev_app)|| ~isvalid(app.ev_app.UIFigure)
                app.ev_app = ev_app_win;
                app.show_on_signal_jump2loc.Visible = "on";
                app.app_ev_close_listener = addlistener(app.ev_app.UIFigure,'ObjectBeingDestroyed',...
                    @(~,~) set(app.show_on_signal_jump2loc,'Visible',"off"));
            end
        end
        
        function show_on_signal_jump2loc_MenuSelected(app,relevant_sample)
            % time2center = relevant_sample./app.ied_obj.fs;
            % app.ev_app.model.move_time(time2center)
            pos2jump = find(app.ied_obj.pos == relevant_sample);
            app.ev_app.current_pos = pos2jump;
        end
        

        % Menu Selected function: wv_on_other_disp 
        function wv_on_other_disp_MenuSelected(app,relevant_sample,context_opener)
            % make the relevant waveform emphasized on the not clicked
            % axis, by 1) Marking it with an arrow (emphasized on reduction
            % axe) or 2) on relevant cluster_disp_axes (unimplamented!)

            % transform relevant_sample to main_scatter datapoint
            [~,point_idx] = ismember(relevant_sample, app.ied_obj.pos(app.ied_obj.accepted));
            
            if context_opener == app.main_scatter
                % open the relevant cluster overlay in a new figure, hide
                % all other signals, toggle plotbrowser

                % collect point info
                cluster_data = app.cluster_idx(point_idx);
                app.UIFigure.CurrentObject = app.cluster_disp_axes(cluster_data);
                [new_fig] = app.open_new_winMenuSelected();
                all_lines = findobj(new_fig,'-depth',2,'Type','Line');
                line2keep = [all_lines.UserData] == relevant_sample;
                [all_lines(~line2keep).Visible] = deal('off');
                warning('off','PlotTools:FunctionDeprecation')
                plotbrowser(new_fig,"on")
                warning('on','PlotTools:FunctionDeprecation')
                drawnow
                plotedit(new_fig,'off')
            else%if isa(context_opener,'matlab.graphics.chart.primitive.Line')
                % add arrow to reduction axe marking the relevant point in
                % main_scatter

                % Delete existing arrow, if any
                existing_arrow = findobj(app.reduction_axe,'-depth',1,'Type','text');
                delete(existing_arrow)

                % collect point info
                x_data = app.main_scatter.XData(point_idx);
                y_data = app.main_scatter.YData(point_idx);
                if ~isempty(app.main_scatter.ZData)
                    z_data = app.main_scatter.ZData(point_idx);
                    text(app.reduction_axe,x_data,y_data,z_data,'\leftarrow','Color','k','FontSize',50,'Interpreter','tex')
                else
                    text(app.reduction_axe,x_data,y_data,'\leftarrow','Color','k','FontSize',50,'PickableParts','none','Interpreter','tex')
                end
            end

        end
    

        % Menu selected function: open_new_win
        function [new_fig] = open_new_winMenuSelected(app)
            % open clicked axe in full window, unrealated to the app

            % collect what axe was clicked on
            current_axe = app.UIFigure.CurrentObject;

            % open a new figure
            new_fig = figure("Color",'w');

            % copy & reposition axe
            new_ax = copyobj(current_axe,new_fig);
            new_ax.Units = 'normalized';
            new_ax.Position = [0.0650    0.0550    0.8525    0.8965];

            % create legend for every line
            all_new_lines = findobj(new_ax,"Type","line");
            legend(new_ax,all_new_lines,"Location","best",...
                'Interpreter','none',...
                "ItemHitFcn",@(~,evt) app.cluster_disp_axes_legend_hide_on_click(evt),...
                "ContextMenu",[])

            % import interactivity if "Overlay" - note that copyobj also copy UserData
            label_data = app.disp_type.Value;
            if isscalar(label_data{1}) && label_data{1} == "Overlays"
                legend(new_ax,'off')
                % create a context menu - it will only refer to the
                % original context menu
                contex_menu_new = copyobj(app.single_wv_menu,new_fig);
                for iChild = 1:numel(contex_menu_new.Children)
                    contex_menu_new.Children(iChild).MenuSelectedFcn = @(src,evt) app.single_wv_menu.Children(iChild).MenuSelectedFcn(src,evt);
                end
                [contex_menu_new.Children.Visible] = deal('on');

                [all_new_lines.ContextMenu] = deal(contex_menu_new);
                [all_new_lines.ButtonDownFcn] = deal(@(src,evt) app.clust_disp_overlays_lines_ButtonDownFcn(src,evt));
            end
        end


    end

    %%%%%%%% Component initialization & destruction %%%%%%%%%
    methods (Hidden)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [1 1 1];
            app.UIFigure.Position = [100 100 714 612];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = @(~,~) app.app_close_fun();

            % Create main_grid
            app.main_grid = uigridlayout(app.UIFigure);
            app.main_grid.ColumnWidth = {'1x'};
            app.main_grid.RowHeight = {'0.4x', '1x'};
            app.main_grid.BackgroundColor = [1 1 1];

            % Create reduction_axe
            app.reduction_axe = uiaxes(app.main_grid);
            title(app.reduction_axe, 'Title')
            xlabel(app.reduction_axe, 'X')
            ylabel(app.reduction_axe, 'Y')
            zlabel(app.reduction_axe, 'Z')
            app.reduction_axe.Layout.Row = 2;
            app.reduction_axe.Layout.Column = 1;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.main_grid);
            app.TabGroup.SelectionChangedFcn = @(~,~) app.refresh_wv_disp_tab();
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%% Reduction Tab %%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%% Base Elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Create reduction_manager
            app.reduction_manager = uitab(app.TabGroup);
            app.reduction_manager.Title = 'Reduction Manager';
            app.reduction_manager.BackgroundColor = [1 1 1];

            % Create reduction_grid
            app.reduction_grid = uigridlayout(app.reduction_manager);
            app.reduction_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.reduction_grid.RowHeight = {'fit', '1x'};
            app.reduction_grid.Scrollable = 'on';
            app.reduction_grid.BackgroundColor = [1 1 1];

            % Create ReductionWaveformwindowmsLabel
            ReductionWaveformwindowmsLabel = uilabel(app.reduction_grid);
            ReductionWaveformwindowmsLabel.Tag = 'base_elem';
            ReductionWaveformwindowmsLabel.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel.Layout.Row = 1;
            ReductionWaveformwindowmsLabel.Layout.Column = 1;
            ReductionWaveformwindowmsLabel.Text = {'Reduction'; 'Waveform window'; '[ms]'};

            % Create reduction_norm_tree
            app.reduction_norm_tree = uitree(app.reduction_grid, 'checkbox');
            app.reduction_norm_tree.Tag = 'base_elem';
            app.reduction_norm_tree.Layout.Row = [1 2];
            app.reduction_norm_tree.Layout.Column = 2;

            % Create reduction_all_norm
            app.reduction_all_norm = uitreenode(app.reduction_norm_tree);
            app.reduction_all_norm.Text = 'All Normalisations';

            % Create reduction_z_norm
            app.reduction_z_norm = uitreenode(app.reduction_all_norm);
            app.reduction_z_norm.Text = 'Z Score Waveforms';

            % Create reduction_align_norm
            app.reduction_align_norm = uitreenode(app.reduction_all_norm);
            app.reduction_align_norm.Text = 'Align Waveforms';

            % Assign Checked Nodes
            app.reduction_norm_tree.CheckedNodesChangedFcn = @(~,evt) app.refresh_reduction(evt);

            % Create ReductiononDropDownLabel
            ReductiononDropDownLabel = uilabel(app.reduction_grid);
            ReductiononDropDownLabel.Tag = 'base_elem';
            ReductiononDropDownLabel.HorizontalAlignment = 'center';
            ReductiononDropDownLabel.Layout.Row = 1;
            ReductiononDropDownLabel.Layout.Column = 3;
            ReductiononDropDownLabel.Text = 'Reduction on:';

            % Create DistMatDistancefunctionEditFieldLabel
            DistMatDistancefunctionEditFieldLabel = uilabel(app.reduction_grid);
            DistMatDistancefunctionEditFieldLabel.HorizontalAlignment = 'center';
            DistMatDistancefunctionEditFieldLabel.Layout.Row = 1;
            DistMatDistancefunctionEditFieldLabel.Layout.Column = 4;
            DistMatDistancefunctionEditFieldLabel.Text = {'Dist Mat'; 'Distance function'};

            % Create RecductionfunctionLabel
            reduction_function_Label = uilabel(app.reduction_grid);
            reduction_function_Label.Tag = 'base_elem';
            reduction_function_Label.HorizontalAlignment = 'center';
            reduction_function_Label.Layout.Row = 1;
            reduction_function_Label.Layout.Column = 5;
            reduction_function_Label.Text = 'Reduction function:';
            
            % Create reduction_wv_win
            app.reduction_wv_win = uieditfield(app.reduction_grid, 'numeric');
            app.reduction_wv_win.Limits = [0 Inf];
            app.reduction_wv_win.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.reduction_wv_win.Tag = 'base_elem';
            app.reduction_wv_win.HorizontalAlignment = 'center';
            app.reduction_wv_win.Layout.Row = 2;
            app.reduction_wv_win.Layout.Column = 1;
            app.reduction_wv_win.Value = 50;

            % Create reduction_source
            app.reduction_source = uidropdown(app.reduction_grid);
            app.reduction_source.Items = {'Base Waveforms', 'Distance Mat'};
            app.reduction_source.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.reduction_source.Tag = 'base_elem';
            app.reduction_source.Layout.Row = 2;
            app.reduction_source.Layout.Column = 3;
            app.reduction_source.Value = 'Base Waveforms';

            % Create reduction_dist_fun
            app.reduction_dist_fun = uieditfield(app.reduction_grid, 'text');
            app.reduction_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.reduction_dist_fun.HorizontalAlignment = 'center';
            app.reduction_dist_fun.Layout.Row = 2;
            app.reduction_dist_fun.Layout.Column = 4;
            app.reduction_dist_fun.Value = 'wv2dist';

            % Create reduction_function
            app.reduction_function = uidropdown(app.reduction_grid);
            app.reduction_function.Items = {'PCA', 't-SNE', 'UMAP', 'cMDS'};
            app.reduction_function.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.reduction_function.Tag = 'base_elem';
            app.reduction_function.Layout.Row = 2;
            app.reduction_function.Layout.Column = 5;
            app.reduction_function.Value = 'PCA';


            %%%%%%%%%%%%%%%%%%%%%% t-SNE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Create tsne_grid
            app.tsne_grid = uigridlayout(app.reduction_grid);
            app.tsne_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.tsne_grid.RowHeight = {'fit', '1x'};
            app.tsne_grid.Layout.Row = [1 2];
            app.tsne_grid.Layout.Column = 6;
            app.tsne_grid.BackgroundColor = [1 1 1];

            % Create DistMatDistancefunctionLabel
            DistMatDistancefunctionLabel = uilabel(app.tsne_grid);
            DistMatDistancefunctionLabel.HorizontalAlignment = 'center';
            DistMatDistancefunctionLabel.Layout.Row = 1;
            DistMatDistancefunctionLabel.Layout.Column = 1;
            DistMatDistancefunctionLabel.Text = {'t-SNE'; 'Distance function'};

            % Create tsne_dist_fun
            app.tsne_dist_fun = uidropdown(app.tsne_grid);
            app.tsne_dist_fun.Items = {'euclidean', 'seuclidean', 'fasteuclidean', 'fastseuclidean', 'cityblock', 'chebychev', 'minkowski', 'mahalanobis', 'cosine', 'correlation', 'spearman', 'hamming', 'jaccard'};
            app.tsne_dist_fun.Editable = 'on';
            app.tsne_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_dist_fun.BackgroundColor = [1 1 1];
            app.tsne_dist_fun.Layout.Row = 2;
            app.tsne_dist_fun.Layout.Column = 1;
            app.tsne_dist_fun.Value = 'euclidean';

            % Create ReductionWaveformwindowmsLabel_2
            ReductionWaveformwindowmsLabel_2 = uilabel(app.tsne_grid);
            ReductionWaveformwindowmsLabel_2.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_2.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_2.Layout.Column = 2;
            ReductionWaveformwindowmsLabel_2.Text = {'t-SNE'; 'Exaggeration'};

            % Create tsne_exaggeration
            app.tsne_exaggeration = uieditfield(app.tsne_grid, 'numeric');
            app.tsne_exaggeration.Limits = [1 Inf];
            app.tsne_exaggeration.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_exaggeration.HorizontalAlignment = 'center';
            app.tsne_exaggeration.Layout.Row = 2;
            app.tsne_exaggeration.Layout.Column = 2;
            app.tsne_exaggeration.Value = 4;

            % Create ReductionWaveformwindowmsLabel_3
            ReductionWaveformwindowmsLabel_3 = uilabel(app.tsne_grid);
            ReductionWaveformwindowmsLabel_3.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_3.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_3.Layout.Column = 3;
            ReductionWaveformwindowmsLabel_3.Text = {'t-SNE'; 'Perplexity'};

            % Create tsne_perplexity
            app.tsne_perplexity = uieditfield(app.tsne_grid, 'numeric');
            app.tsne_perplexity.Limits = [0 Inf];
            app.tsne_perplexity.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_perplexity.HorizontalAlignment = 'center';
            app.tsne_perplexity.Layout.Row = 2;
            app.tsne_perplexity.Layout.Column = 3;
            app.tsne_perplexity.Value = 30;
            
            % Create tSNEAlgorithmDropDownLabel
            tSNEAlgorithmDropDownLabel = uilabel(app.tsne_grid);
            tSNEAlgorithmDropDownLabel.HorizontalAlignment = 'center';
            tSNEAlgorithmDropDownLabel.Layout.Row = 1;
            tSNEAlgorithmDropDownLabel.Layout.Column = 4;
            tSNEAlgorithmDropDownLabel.Text = {'t-SNE'; 'Algorithm'};

            % Create tsne_Algorithm
            app.tsne_Algorithm = uidropdown(app.tsne_grid);
            app.tsne_Algorithm.Items = {'exact', 'barneshut'};
            app.tsne_Algorithm.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_Algorithm.Layout.Row = 2;
            app.tsne_Algorithm.Layout.Column = 4;
            app.tsne_Algorithm.Value = 'exact';

            % Create ReductionWaveformwindowmsLabel_4
            ReductionWaveformwindowmsLabel_4 = uilabel(app.tsne_grid);
            ReductionWaveformwindowmsLabel_4.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_4.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_4.Layout.Column = 5;
            ReductionWaveformwindowmsLabel_4.Text = {'t-SNE'; 'LearnRate'};

            % Create tsne_LearnRate
            app.tsne_LearnRate = uieditfield(app.tsne_grid, 'numeric');
            app.tsne_LearnRate.Limits = [0 Inf];
            app.tsne_LearnRate.RoundFractionalValues = 'on';
            app.tsne_LearnRate.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_LearnRate.HorizontalAlignment = 'center';
            app.tsne_LearnRate.Layout.Row = 2;
            app.tsne_LearnRate.Layout.Column = 5;
            app.tsne_LearnRate.Value = 500;

            

            % Create tsne_max_iter
            tsne_max_iter_label = uilabel(app.tsne_grid);
            tsne_max_iter_label.HorizontalAlignment = 'center';
            tsne_max_iter_label.Layout.Row = 1;
            tsne_max_iter_label.Layout.Column = 6;
            tsne_max_iter_label.Text = {'t-SNE'; 'Max Iter'};

            app.tsne_max_iter = uieditfield(app.tsne_grid, 'numeric');
            app.tsne_max_iter.Limits = [1 Inf];
            app.tsne_max_iter.RoundFractionalValues = 'on';
            app.tsne_max_iter.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.tsne_max_iter.HorizontalAlignment = 'center';
            app.tsne_max_iter.Layout.Row = 2;
            app.tsne_max_iter.Layout.Column = 6;
            app.tsne_max_iter.Value = 1000;
            

            %%%%%%%%%%%%%%%%%%%%%% UMAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % create umap_grid
            app.umap_grid = uigridlayout(app.reduction_grid);
            app.umap_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.umap_grid.RowHeight = {'fit', '1x'};
            app.umap_grid.Layout.Row = [1 2];
            app.umap_grid.Layout.Column = 7;
            app.umap_grid.BackgroundColor = [1 1 1];
            
            % Create umap_dist_fun
            umap_distfun_label = uilabel(app.umap_grid);
            umap_distfun_label.HorizontalAlignment = 'center';
            umap_distfun_label.Layout.Row = 1;
            umap_distfun_label.Layout.Column = 1;
            umap_distfun_label.Text = {'UMAP'; 'Distance function'};

            app.umap_dist_fun = uidropdown(app.umap_grid);
            app.umap_dist_fun.Items = {'euclidean', 'seuclidean', 'fasteuclidean', 'fastseuclidean', 'cityblock', 'chebychev', 'minkowski', 'mahalanobis', 'cosine', 'correlation', 'spearman', 'hamming', 'jaccard'};
            app.umap_dist_fun.Editable = 'on';
            app.umap_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_dist_fun.BackgroundColor = [1 1 1];
            app.umap_dist_fun.Layout.Row = 2;
            app.umap_dist_fun.Layout.Column = 1;
            app.umap_dist_fun.Value = 'euclidean';
            
            % create umap_n_neighbors
            umap_n_neighbors_label = uilabel(app.umap_grid);
            umap_n_neighbors_label.HorizontalAlignment = 'center';
            umap_n_neighbors_label.Layout.Row = 1;
            umap_n_neighbors_label.Layout.Column = 2;
            umap_n_neighbors_label.Text = {'UMAP'; 'n Neighbors'};

            app.umap_n_neighbors = uieditfield(app.umap_grid, 'numeric');
            app.umap_n_neighbors.Limits = [0 Inf];
            app.umap_n_neighbors.RoundFractionalValues = 'on';
            app.umap_n_neighbors.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_n_neighbors.HorizontalAlignment = 'center';
            app.umap_n_neighbors.Layout.Row = 2;
            app.umap_n_neighbors.Layout.Column = 2;
            app.umap_n_neighbors.Value = 15;
            
            % create umap_spread
            umap_spread_label = uilabel(app.umap_grid);
            umap_spread_label.HorizontalAlignment = 'center';
            umap_spread_label.Layout.Row = 1;
            umap_spread_label.Layout.Column = 3;
            umap_spread_label.Text = {'UMAP'; 'Spread'};

            app.umap_spread = uieditfield(app.umap_grid, 'numeric');
            app.umap_spread.Limits = [0 Inf];
            app.umap_spread.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_spread.HorizontalAlignment = 'center';
            app.umap_spread.Layout.Row = 2;
            app.umap_spread.Layout.Column = 3;
            app.umap_spread.Value = 1;
            
            % create umap_min_dist
            umap_min_dist_label = uilabel(app.umap_grid);
            umap_min_dist_label.HorizontalAlignment = 'center';
            umap_min_dist_label.Layout.Row = 1;
            umap_min_dist_label.Layout.Column = 4;
            umap_min_dist_label.Text = {'UMAP'; 'Min Distance'};

            app.umap_min_dist = uieditfield(app.umap_grid, 'numeric');
            app.umap_min_dist.Limits = [0 10];
            app.umap_min_dist.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_min_dist.HorizontalAlignment = 'center';
            app.umap_min_dist.Layout.Row = 2;
            app.umap_min_dist.Layout.Column = 4;
            app.umap_min_dist.Value = 0.3;

            % create umap_init
            umap_init_label = uilabel(app.umap_grid);
            umap_init_label.HorizontalAlignment = 'center';
            umap_init_label.Layout.Row = 1;
            umap_init_label.Layout.Column = 5;
            umap_init_label.Text = {'UMAP'; 'Start'};

            app.umap_init = uidropdown(app.umap_grid);
            app.umap_init.Items = {'spectral', 'random'};
            app.umap_init.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_init.BackgroundColor = [1 1 1];
            app.umap_init.Layout.Row = 2;
            app.umap_init.Layout.Column = 5;
            app.umap_init.Value = 'spectral';
            
            % create umap_fast_approx
            umap_fast_approx_label = uilabel(app.umap_grid);
            umap_fast_approx_label.HorizontalAlignment = 'center';
            umap_fast_approx_label.Layout.Row = 1;
            umap_fast_approx_label.Layout.Column = 6;
            umap_fast_approx_label.Text = {'UMAP'; 'Fast Approximation'};

            app.umap_fast_approx = uidropdown(app.umap_grid);
            app.umap_fast_approx.Items = {'True', 'False'};
            app.umap_fast_approx.ValueChangedFcn = @(~,evt) app.refresh_reduction(evt);
            app.umap_fast_approx.BackgroundColor = [1 1 1];
            app.umap_fast_approx.Layout.Row = 2;
            app.umap_fast_approx.Layout.Column = 6;
            app.umap_fast_approx.Value = 'False';

            % create umap_rerun_button
            % Create SaveAppButton
            app.umap_rerun_button = uibutton(app.umap_grid, 'push');
            app.umap_rerun_button.Layout.Row = [1 2];
            app.umap_rerun_button.Layout.Column = 7;
            app.umap_rerun_button.Text = 'Rerun';
            app.umap_rerun_button.ButtonPushedFcn = @(~,evt) app.refresh_reduction(evt);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%% Cluster Tab %%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%% Base Elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Create cluster_manager
            app.cluster_manager = uitab(app.TabGroup);
            app.cluster_manager.Title = 'Cluster Manager';
            app.cluster_manager.BackgroundColor = [1 1 1];

            % Create cluster_grid
            app.cluster_grid = uigridlayout(app.cluster_manager);
            app.cluster_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.cluster_grid.RowHeight = {'fit', '1x'};
            app.cluster_grid.Scrollable = 'on';
            app.cluster_grid.BackgroundColor = [1 1 1];

            % Create ClusterMethodLabel
            ClusterMethodLabel = uilabel(app.cluster_grid);
            ClusterMethodLabel.Tag = 'base_elem';
            ClusterMethodLabel.HorizontalAlignment = 'center';
            ClusterMethodLabel.Layout.Row = 1;
            ClusterMethodLabel.Layout.Column = 5;
            ClusterMethodLabel.Text = 'Cluster Method:';

            % Create StartMarkingButton
            app.start_marking_button = uibutton(app.cluster_grid, 'state');
            app.start_marking_button.ValueChangedFcn = @(~,~) app.StartMarkingButtonValueChanged();
            app.start_marking_button.Text = 'Start Marking';
            app.start_marking_button.Layout.Row = [1 2];
            app.start_marking_button.Layout.Column = 9;

            % Create cluster_method
            app.cluster_method = uidropdown(app.cluster_grid);
            app.cluster_method.Items = {'Manual', 'by Label', 'kmeans', 'dbscan', 'hierarchical', 'None'};
            app.cluster_method.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.cluster_method.Tag = 'base_elem';
            app.cluster_method.Layout.Row = 2;
            app.cluster_method.Layout.Column = 5;
            app.cluster_method.Value = 'kmeans';

            % Create ClusterOnLabel
            ClusterOnLabel = uilabel(app.cluster_grid);
            ClusterOnLabel.Tag = 'base_elem';
            ClusterOnLabel.HorizontalAlignment = 'center';
            ClusterOnLabel.Layout.Row = 1;
            ClusterOnLabel.Layout.Column = 1;
            ClusterOnLabel.Text = 'Cluster On:';

            % Create cluster_source
            app.cluster_source = uidropdown(app.cluster_grid);
            app.cluster_source.Items = {'Reduced Data 2d', 'Distance Mat', 'Base Waveforms'};
            app.cluster_source.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.cluster_source.Tag = 'base_elem';
            app.cluster_source.Layout.Row = 2;
            app.cluster_source.Layout.Column = 1;
            app.cluster_source.Value = 'Reduced Data 2d';

            % Create cluster_norm_tree
            app.cluster_norm_tree = uitree(app.cluster_grid, 'checkbox');
            app.cluster_norm_tree.Tag = 'base_elem';
            app.cluster_norm_tree.Layout.Row = [1 2];
            app.cluster_norm_tree.Layout.Column = 4;

            % Create cluster_all_norm
            app.cluster_all_norm = uitreenode(app.cluster_norm_tree);
            app.cluster_all_norm.Text = 'All Normalisations';

            % Create cluster_z_norm
            app.cluster_z_norm = uitreenode(app.cluster_all_norm);
            app.cluster_z_norm.Text = 'Z Score Waveforms';

            % Create cluster_align_norm
            app.cluster_align_norm = uitreenode(app.cluster_all_norm);
            app.cluster_align_norm.Text = 'Align Waveforms';

            % Assign Checked Nodes
            app.cluster_norm_tree.CheckedNodesChangedFcn = @(~,evt) app.refresh_clust(evt);

            % Create cluster_wv_win
            app.cluster_wv_win = uieditfield(app.cluster_grid, 'numeric');
            app.cluster_wv_win.Limits = [0 Inf];
            app.cluster_wv_win.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.cluster_wv_win.Tag = 'base_elem';
            app.cluster_wv_win.HorizontalAlignment = 'center';
            app.cluster_wv_win.Layout.Row = 2;
            app.cluster_wv_win.Layout.Column = 3;
            app.cluster_wv_win.Value = 50;

            % Create ClusterWaveformwindowmsEditFieldLabel
            ClusterWaveformwindowmsEditFieldLabel = uilabel(app.cluster_grid);
            ClusterWaveformwindowmsEditFieldLabel.Tag = 'base_elem';
            ClusterWaveformwindowmsEditFieldLabel.HorizontalAlignment = 'center';
            ClusterWaveformwindowmsEditFieldLabel.Layout.Row = 1;
            ClusterWaveformwindowmsEditFieldLabel.Layout.Column = 3;
            ClusterWaveformwindowmsEditFieldLabel.Text = {'Cluster'; 'Waveform window'; '[ms]'};

            % Create cluster_dist_fun
            app.cluster_dist_fun = uieditfield(app.cluster_grid, 'text');
            app.cluster_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.cluster_dist_fun.HorizontalAlignment = 'center';
            app.cluster_dist_fun.Layout.Row = 2;
            app.cluster_dist_fun.Layout.Column = 2;
            app.cluster_dist_fun.Value = 'wv2dist';

            % Create DistMatDistancefunctionEditField_2Label
            DistMatDistancefunctionEditField_2Label = uilabel(app.cluster_grid);
            DistMatDistancefunctionEditField_2Label.HorizontalAlignment = 'center';
            DistMatDistancefunctionEditField_2Label.Layout.Row = 1;
            DistMatDistancefunctionEditField_2Label.Layout.Column = 2;
            DistMatDistancefunctionEditField_2Label.Text = {'Dist Mat'; 'Distance function'};

            % create clust_by_label_grid
            app.clust_by_label_grid = uigridlayout(app.cluster_grid);
            app.clust_by_label_grid.ColumnWidth = {'fit'};
            app.clust_by_label_grid.RowHeight = {'fit', '1x'};
            app.clust_by_label_grid.Layout.Row = [1 2];
            app.clust_by_label_grid.Layout.Column = 10;
            app.clust_by_label_grid.BackgroundColor = [1 1 1];

            % Create label2clust_by_Label
            label2clust_by_Label = uilabel(app.clust_by_label_grid);
            label2clust_by_Label.HorizontalAlignment = 'center';
            label2clust_by_Label.Layout.Row = 1;
            label2clust_by_Label.Layout.Column = 2;
            label2clust_by_Label.Text = 'Label 2 Cluster by:';

            % create label2clust_by
            app.label2clust_by = uidropdown(app.clust_by_label_grid);
            app.label2clust_by.Items = {'None'};
            app.label2clust_by.ValueChangedFcn =  @(~,evt) app.refresh_clust(evt);
            app.label2clust_by.Layout.Row = 2;
            app.label2clust_by.Layout.Column = 2;
            app.label2clust_by.Value = 'None';

            % Create dbscan_grid
            app.dbscan_grid = uigridlayout(app.cluster_grid);
            app.dbscan_grid.ColumnWidth = {'fit', 'fit', 'fit'};
            app.dbscan_grid.RowHeight = {'fit', '1x'};
            app.dbscan_grid.Layout.Row = [1 2];
            app.dbscan_grid.Layout.Column = 7;
            app.dbscan_grid.BackgroundColor = [1 1 1];

            % Create ReductionWaveformwindowmsLabel_9
            ReductionWaveformwindowmsLabel_9 = uilabel(app.dbscan_grid);
            ReductionWaveformwindowmsLabel_9.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_9.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_9.Layout.Column = 1;
            ReductionWaveformwindowmsLabel_9.Text = {'dbscan'; 'epsilon'};

            % Create dbscan_epsilon
            app.dbscan_epsilon = uieditfield(app.dbscan_grid, 'numeric');
            app.dbscan_epsilon.Limits = [0 Inf];
            app.dbscan_epsilon.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.dbscan_epsilon.HorizontalAlignment = 'center';
            app.dbscan_epsilon.Layout.Row = 2;
            app.dbscan_epsilon.Layout.Column = 1;
            app.dbscan_epsilon.Value = 3;

            % Create ReductionWaveformwindowmsLabel_10
            ReductionWaveformwindowmsLabel_10 = uilabel(app.dbscan_grid);
            ReductionWaveformwindowmsLabel_10.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_10.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_10.Layout.Column = 2;
            ReductionWaveformwindowmsLabel_10.Text = {'dbscan'; 'minpts'};

            % Create dbscan_minpts
            app.dbscan_minpts = uieditfield(app.dbscan_grid, 'numeric');
            app.dbscan_minpts.Limits = [1 Inf];
            app.dbscan_minpts.RoundFractionalValues = 'on';
            app.dbscan_minpts.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.dbscan_minpts.HorizontalAlignment = 'center';
            app.dbscan_minpts.Layout.Row = 2;
            app.dbscan_minpts.Layout.Column = 2;
            app.dbscan_minpts.Value = 4;

            % Create DistMatDistancefunctionLabel_2
            DistMatDistancefunctionLabel_2 = uilabel(app.dbscan_grid);
            DistMatDistancefunctionLabel_2.Tag = 'base_elem';
            DistMatDistancefunctionLabel_2.HorizontalAlignment = 'center';
            DistMatDistancefunctionLabel_2.Layout.Row = 1;
            DistMatDistancefunctionLabel_2.Layout.Column = 3;
            DistMatDistancefunctionLabel_2.Text = {'dbscan'; 'Distance function'};

            % Create dbscan_dist_fun
            app.dbscan_dist_fun = uidropdown(app.dbscan_grid);
            app.dbscan_dist_fun.Items = {'precomputed', 'euclidean', 'squaredeuclidean', 'seuclidean', 'mahalanobis', 'cityblock', 'minkowski', 'chebychev', 'cosine', 'correlation', 'hamming', 'jaccard', 'spearman'};
            app.dbscan_dist_fun.Editable = 'on';
            app.dbscan_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.dbscan_dist_fun.BackgroundColor = [1 1 1];
            app.dbscan_dist_fun.Layout.Row = 2;
            app.dbscan_dist_fun.Layout.Column = 3;
            app.dbscan_dist_fun.Value = 'euclidean';

            % Create HC_grid
            app.HC_grid = uigridlayout(app.cluster_grid);
            app.HC_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit'};
            app.HC_grid.RowHeight = {'fit', '1x'};
            app.HC_grid.Layout.Row = [1 2];
            app.HC_grid.Layout.Column = 8;
            app.HC_grid.BackgroundColor = [1 1 1];

            % Create ReductionWaveformwindowmsLabel_12
            ReductionWaveformwindowmsLabel_12 = uilabel(app.HC_grid);
            ReductionWaveformwindowmsLabel_12.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_12.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_12.Layout.Column = 3;
            ReductionWaveformwindowmsLabel_12.Text = {'HC Cluster'; 'cutoff'};

            % Create HC_clust_cutoff
            app.HC_clust_cutoff = uieditfield(app.HC_grid, 'numeric');
            app.HC_clust_cutoff.Limits = [0 Inf];
            app.HC_clust_cutoff.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.HC_clust_cutoff.HorizontalAlignment = 'center';
            app.HC_clust_cutoff.Layout.Row = 2;
            app.HC_clust_cutoff.Layout.Column = 3;
            app.HC_clust_cutoff.Value = 1;

            % Create HC_clust_critera
            app.HC_clust_critera = uidropdown(app.HC_grid);
            app.HC_clust_critera.Items = {'inconsistent', 'distance'};
            app.HC_clust_critera.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.HC_clust_critera.Layout.Row = 2;
            app.HC_clust_critera.Layout.Column = 2;
            app.HC_clust_critera.Value = 'inconsistent';

            % Create kmeanDistancefunctionLabel_3
            kmeanDistancefunctionLabel_3 = uilabel(app.HC_grid);
            kmeanDistancefunctionLabel_3.HorizontalAlignment = 'center';
            kmeanDistancefunctionLabel_3.Layout.Row = 1;
            kmeanDistancefunctionLabel_3.Layout.Column = 2;
            kmeanDistancefunctionLabel_3.Text = {'HC Cluster'; 'Criterion'};

            % Create kmeanDistancefunctionLabel_2
            kmeanDistancefunctionLabel_2 = uilabel(app.HC_grid);
            kmeanDistancefunctionLabel_2.HorizontalAlignment = 'center';
            kmeanDistancefunctionLabel_2.Layout.Row = 1;
            kmeanDistancefunctionLabel_2.Layout.Column = 1;
            kmeanDistancefunctionLabel_2.Text = {'HC Linkage'; 'Method'};

            % Create HC_link_method
            app.HC_link_method = uidropdown(app.HC_grid);
            app.HC_link_method.Items = {'weighted', 'average', 'centroid', 'complete', 'median', 'single', 'ward'};
            app.HC_link_method.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.HC_link_method.Layout.Row = 2;
            app.HC_link_method.Layout.Column = 1;
            app.HC_link_method.Value = 'weighted';

            % Create HC_clust_MaxClust
            app.HC_clust_MaxClust = uieditfield(app.HC_grid, 'numeric');
            app.HC_clust_MaxClust.Limits = [1 Inf];
            app.HC_clust_MaxClust.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.HC_clust_MaxClust.Layout.Row = 2;
            app.HC_clust_MaxClust.Layout.Column = 4;
            if isMATLABReleaseOlderThan('R2023a')
                % no "AllowEmpty" before R2023, use "inf" instead.
                app.HC_clust_MaxClust.Value = inf;
            else
                % allow empty instead of inf - more intuitive
                app.HC_clust_MaxClust.AllowEmpty = 'on';
                app.HC_clust_MaxClust.Value = [];
            end


            % Create ReductionWaveformwindowmsLabel_11
            ReductionWaveformwindowmsLabel_11 = uilabel(app.HC_grid);
            ReductionWaveformwindowmsLabel_11.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_11.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_11.Layout.Column = 4;
            ReductionWaveformwindowmsLabel_11.Text = {'HC Cluster'; 'MaxClust'};

            % Create kmean_grid
            app.kmean_grid = uigridlayout(app.cluster_grid);
            app.kmean_grid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            app.kmean_grid.RowHeight = {'fit', '1x'};
            app.kmean_grid.Layout.Row = [1 2];
            app.kmean_grid.Layout.Column = 6;
            app.kmean_grid.BackgroundColor = [1 1 1];

            % Create ReductionWaveformwindowmsLabel_5
            ReductionWaveformwindowmsLabel_5 = uilabel(app.kmean_grid);
            ReductionWaveformwindowmsLabel_5.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_5.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_5.Layout.Column = 4;
            ReductionWaveformwindowmsLabel_5.Text = {'kmean'; 'MaxIter'};

            % Create kmean_MaxIter
            app.kmean_MaxIter = uieditfield(app.kmean_grid, 'numeric');
            app.kmean_MaxIter.Limits = [0 Inf];
            app.kmean_MaxIter.RoundFractionalValues = 'on';
            app.kmean_MaxIter.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_MaxIter.HorizontalAlignment = 'center';
            app.kmean_MaxIter.Layout.Row = 2;
            app.kmean_MaxIter.Layout.Column = 4;
            app.kmean_MaxIter.Value = 1000;

            % Create ReductionWaveformwindowmsLabel_6
            ReductionWaveformwindowmsLabel_6 = uilabel(app.kmean_grid);
            ReductionWaveformwindowmsLabel_6.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_6.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_6.Layout.Column = 1;
            ReductionWaveformwindowmsLabel_6.Text = {'kmean'; 'nClusters'};

            % Create kmean_nClusters
            app.kmean_nClusters = uieditfield(app.kmean_grid, 'numeric');
            app.kmean_nClusters.Limits = [1 Inf];
            app.kmean_nClusters.RoundFractionalValues = 'on';
            app.kmean_nClusters.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_nClusters.HorizontalAlignment = 'center';
            app.kmean_nClusters.Layout.Row = 2;
            app.kmean_nClusters.Layout.Column = 1;
            app.kmean_nClusters.Value = 3;

            % Create ReductionWaveformwindowmsLabel_7
            ReductionWaveformwindowmsLabel_7 = uilabel(app.kmean_grid);
            ReductionWaveformwindowmsLabel_7.HorizontalAlignment = 'center';
            ReductionWaveformwindowmsLabel_7.Layout.Row = 1;
            ReductionWaveformwindowmsLabel_7.Layout.Column = 6;
            ReductionWaveformwindowmsLabel_7.Text = {'kmean'; 'Replicates'};

            % Create kmean_replicates
            app.kmean_replicates = uieditfield(app.kmean_grid, 'numeric');
            app.kmean_replicates.Limits = [1 Inf];
            app.kmean_replicates.RoundFractionalValues = 'on';
            app.kmean_replicates.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_replicates.HorizontalAlignment = 'center';
            app.kmean_replicates.Layout.Row = 2;
            app.kmean_replicates.Layout.Column = 6;
            app.kmean_replicates.Value = 1;

            % Create kmeanStartDropDownLabel
            kmeanStartDropDownLabel = uilabel(app.kmean_grid);
            kmeanStartDropDownLabel.HorizontalAlignment = 'center';
            kmeanStartDropDownLabel.Layout.Row = 1;
            kmeanStartDropDownLabel.Layout.Column = 3;
            kmeanStartDropDownLabel.Text = {'kmean'; 'Start'};

            % Create kmean_start
            app.kmean_start = uidropdown(app.kmean_grid);
            app.kmean_start.Items = {'plus', 'cluster', 'sample', 'uniform'};
            app.kmean_start.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_start.Layout.Row = 2;
            app.kmean_start.Layout.Column = 3;
            app.kmean_start.Value = 'plus';

            % Create kmeanDistancefunctionLabel
            kmeanDistancefunctionLabel = uilabel(app.kmean_grid);
            kmeanDistancefunctionLabel.HorizontalAlignment = 'center';
            kmeanDistancefunctionLabel.Layout.Row = 1;
            kmeanDistancefunctionLabel.Layout.Column = 2;
            kmeanDistancefunctionLabel.Text = {'kmean'; 'Distance function'};

            % Create kmean_dist_fun
            app.kmean_dist_fun = uidropdown(app.kmean_grid);
            app.kmean_dist_fun.Items = {'sqeuclidean', 'cityblock', 'cosine', 'correlation', 'hamming'};
            app.kmean_dist_fun.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_dist_fun.Layout.Row = 2;
            app.kmean_dist_fun.Layout.Column = 2;
            app.kmean_dist_fun.Value = 'sqeuclidean';

            % Create kmeanOnlinePhaseDropDownLabel
            kmeanOnlinePhaseDropDownLabel = uilabel(app.kmean_grid);
            kmeanOnlinePhaseDropDownLabel.HorizontalAlignment = 'center';
            kmeanOnlinePhaseDropDownLabel.Layout.Row = 1;
            kmeanOnlinePhaseDropDownLabel.Layout.Column = 5;
            kmeanOnlinePhaseDropDownLabel.Text = {'kmean'; 'Online Phase'};

            % Create kmean_online
            app.kmean_online = uidropdown(app.kmean_grid);
            app.kmean_online.Items = {'On', 'Off'};
            app.kmean_online.ValueChangedFcn = @(~,evt) app.refresh_clust(evt);
            app.kmean_online.Layout.Row = 2;
            app.kmean_online.Layout.Column = 5;
            app.kmean_online.Value = 'On';

            % Create WaveformsDisplayTab
            app.WaveformsDisplayTab = uitab(app.TabGroup);
            app.WaveformsDisplayTab.Title = 'Waveforms Display';
            app.WaveformsDisplayTab.BackgroundColor = [1 1 1];
            app.WaveformsDisplayTab.Scrollable = 'on';

            % Create display_grid
            app.display_grid = uigridlayout(app.WaveformsDisplayTab);
            app.display_grid.ColumnWidth = {'fit', 'fit'};
            app.display_grid.RowHeight = {'fit', '1x'};
            app.display_grid.Scrollable = 'on';
            app.display_grid.BackgroundColor = [1 1 1];

            % Create DisplayTypeLabel
            DisplayTypeLabel = uilabel(app.display_grid);
            DisplayTypeLabel.HorizontalAlignment = 'center';
            DisplayTypeLabel.Layout.Row = 1;
            DisplayTypeLabel.Layout.Column = 2;
            DisplayTypeLabel.Text = 'Display Type:';

            % Create disp_type
            app.disp_type = uidropdown(app.display_grid);
            app.disp_type.Items = {'Full Cluster Average', 'Clusters as Labels', 'Overlays'};
            app.disp_type.ValueChangedFcn = @(~,~) app.disp_typeValueChanged();
            app.disp_type.Tag = 'base_elem';
            app.disp_type.Layout.Row = 2;
            app.disp_type.Layout.Column = 3;
            app.disp_type.Value = 'Full Cluster Average';

            % Create cluster_norm_tree
            app.disp_norm_tree = uitree(app.display_grid, 'checkbox');
            app.disp_norm_tree.Tag = 'base_elem';
            app.disp_norm_tree.Layout.Row = [1 2];
            app.disp_norm_tree.Layout.Column = 2;

            % Create cluster_all_norm
            app.disp_all_norm = uitreenode(app.disp_norm_tree);
            app.disp_all_norm.Text = 'All Normalisations';

            % Create cluster_z_norm
            app.disp_z_norm = uitreenode(app.disp_all_norm);
            app.disp_z_norm.Text = 'Z Score Waveforms';

            % Create cluster_align_norm
            app.disp_align_norm = uitreenode(app.disp_all_norm);
            app.disp_align_norm.Text = 'Align Waveforms';

            % Assign Checked Nodes
            app.disp_norm_tree.CheckedNodesChangedFcn =  @(~,~) app.refresh_wv_disp_tab();

            % Create DisplayWaveformwindowmsEditFieldLabel
            DisplayWaveformwindowmsEditFieldLabel = uilabel(app.display_grid);
            DisplayWaveformwindowmsEditFieldLabel.HorizontalAlignment = 'center';
            DisplayWaveformwindowmsEditFieldLabel.Layout.Row = 1;
            DisplayWaveformwindowmsEditFieldLabel.Layout.Column = 1;
            DisplayWaveformwindowmsEditFieldLabel.Text = {'Display'; 'Waveform window'; '[ms]'};

            % Create disp_wv_win
            app.disp_wv_win = uieditfield(app.display_grid, 'numeric');
            app.disp_wv_win.Limits = [0 Inf];
            app.disp_wv_win.ValueChangedFcn = @(~,~) app.refresh_wv_disp_tab();
            app.disp_wv_win.Tag = 'base_elem';
            app.disp_wv_win.HorizontalAlignment = 'center';
            app.disp_wv_win.Layout.Row = 2;
            app.disp_wv_win.Layout.Column = 1;
            app.disp_wv_win.Value = 50;

            % Create ExportOptionsTab
            app.ExportOptionsTab = uitab(app.TabGroup);
            app.ExportOptionsTab.Title = 'Export Options';
            app.ExportOptionsTab.BackgroundColor = [1 1 1];

            % Create export_grid
            app.export_grid = uigridlayout(app.ExportOptionsTab);
            app.export_grid.ColumnWidth = {'1x', '1x', '1x'};
            app.export_grid.RowHeight = {'1x'};
            app.export_grid.BackgroundColor = [1 1 1];

            % Create SaveAppButton
            app.SaveAppButton = uibutton(app.export_grid, 'push');
            app.SaveAppButton.Layout.Row = 1;
            app.SaveAppButton.Layout.Column = 1;
            app.SaveAppButton.Text = 'Save App';
            app.SaveAppButton.ButtonPushedFcn = @(~,~) app.SaveAppButtonPushed();

            % Create ExportDataButton
            app.ExportDataButton = uibutton(app.export_grid, 'push');
            app.ExportDataButton.Layout.Row = 1;
            app.ExportDataButton.Layout.Column = 2;
            app.ExportDataButton.Text = 'Export Data';
            app.ExportDataButton.ButtonPushedFcn = @(~,~) app.ExportDataButtonPushed();

            % Create GenerateCodeButton
            app.GenerateCodeButton = uibutton(app.export_grid, 'push');
            app.GenerateCodeButton.Layout.Row = 1;
            app.GenerateCodeButton.Layout.Column = 3;
            app.GenerateCodeButton.Text = 'Generate Code';
            app.GenerateCodeButton.ButtonPushedFcn = @(~,~) app.GenerateCodeButtonPushed();

            % Create cluster_disp_menu
            app.cluster_disp_menu = uicontextmenu(app.UIFigure);

            % Create open_new_win
            app.open_new_win = uimenu(app.cluster_disp_menu);
            app.open_new_win.MenuSelectedFcn = @(~,~) app.open_new_winMenuSelected();
            app.open_new_win.Text = 'Open in New Window';

            % Create main_scatter_menu
            app.single_wv_menu = uicontextmenu(app.UIFigure);
            % app.main_scatter_menu.ContextMenuOpeningFcn = @(s,evt) app.main_scatter_menu_opened(s,evt);

            % Create show_on_signal_new_win
            app.show_on_signal_new_win = uimenu(app.single_wv_menu);
            app.show_on_signal_new_win.Text = 'Show On Signal (new window)';

            app.show_on_signal_jump2loc = uimenu(app.single_wv_menu);
            app.show_on_signal_jump2loc.Text = 'Show On Signal (jump to location)';
            app.show_on_signal_jump2loc.Visible = "off";

            % Create show_on_signal_new_win
            app.wv_on_other_disp = uimenu(app.single_wv_menu);
            app.wv_on_other_disp.Text = 'Mark waveform on other axe';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    
        % close request fun
        function app_close_fun(app)
            % delete noisy leftovers of the function before closing
            
            % delete listener for events subwindow
            if ~isempty(app.app_ev_close_listener)
                delete(app.app_ev_close_listener)
            end

            delete(app.UIFigure)
        end
    end

    %%%%%%% App creation and deletion %%%%%%%%
    methods (Access = public)

        % Construct app
        function app = reduct_displayer(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            % registerApp(app, app.UIFigure)

            % Execute the startup function
            % runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            startupFcn(app, varargin{:})

            if nargout == 0
                clear app
            end
        end


        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end



end