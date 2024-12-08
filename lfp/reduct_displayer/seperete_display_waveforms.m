function [map_fig,waveforms_figs] = seperete_display_waveforms(reduced_data, marks_vec, label_idx, ied)
    % seperete_display_waveforms - Create graphics displaying waveforms for each cluster / label.
    % This function creates visual representations of waveform data for different
    % clusters and marks within the data.
    %
    % USAGE:
    %   [waveforms_figs] = seperete_display_waveforms(reduced_data, marks_vec, label_idx, ied)
    %
    % INPUT:
    %   reduced_data - Numeric matrix of shape (N, 2), where N is the number of data
    %                  points (events). It contains the reduced data for visualization (e.g., PCA or t-SNE).
    %   marks_vec -    Logical vector of shape (N, 1), indicating marked (true) or
    %                  unmarked (false) events for each data point.
    %   label_idx -    Vector (numeric, string, or categorical) of shape (N, 1), containing
    %                  labels for each data point. These labels identify the cluster to which
    %                  each data point belongs.
    %   ied -          IED.data object, containing each waveform data.
    %
    % OUTPUT:
    %   waveforms_figs - An array of figure handles, each displaying waveforms for a cluster.

    unique_labels = unique(label_idx);
    label_idx_cat = categorical(label_idx);

    colors = distinguishable_colors(numel(unique_labels));
    map_fig = figure('Name','Full Map','Color','w','WindowState','maximized');
    tl = tiledlayout(1,3,'Parent',map_fig);

    map_marked_true = nexttile(tl);
    gscatter(map_marked_true,reduced_data(marks_vec,1),reduced_data(marks_vec,2),label_idx_cat(marks_vec),colors,'O',5,'filled');
    title(map_marked_true,'Marked True IED')
    legend(map_marked_true,'off')

    map_marked_false = nexttile(tl);
    gscatter(map_marked_false,reduced_data(~marks_vec,1),reduced_data(~marks_vec,2),label_idx_cat(~marks_vec),colors,'O',5,'filled');
    title('Marked Not IED (False)')
    legend(map_marked_false,'off')

    map_all_marked = nexttile(tl);
    gscatter(map_all_marked,reduced_data(:,1),reduced_data(:,2),label_idx_cat,colors,'O',5,'filled');
    title('All Detections')
    legend(map_all_marked,"show",'Interpreter','none')
    linkaxes([map_marked_true map_marked_false map_all_marked])

    for iLabel = numel(unique_labels):-1:1
        current_label = unique_labels(iLabel);
        clust_idx = label_idx == current_label;

        if isnumeric(current_label)
            fig_name = sprintf('Cluster %d',current_label);
        else
            fig_name = string(current_label);
        end

        %%%%%%%% Marked True / Marked False %%%%%%%%
        % collect true - false data
        only_marked_true_clust = clust_idx &  marks_vec;
        only_marked_false_clust = clust_idx &  ~marks_vec;

        % plot true - false data
        [marked_true_events, t_true] = ied.extract_discharges(150/1000, only_marked_true_clust);
        [marked_false_events,t_false] = ied.extract_discharges(150/1000, only_marked_false_clust);

        %%% create graphics
        waveforms_figs(iLabel) = figure('Name',fig_name,'Color','w','WindowState','maximized');
        tl = tiledlayout(2,3,"Parent",waveforms_figs(iLabel));

        % cluster map true
        ax_cluster_true = nexttile(tl);
        gscatter(reduced_data(marks_vec,1), reduced_data(marks_vec,2), clust_idx(marks_vec),'kg','o','filled')
        title("Reduced Data Marked True")
        legend(ax_cluster_true,'off')

        % all waveforms true
        ax_wv_true = nexttile(tl);
        p = plot(ax_wv_true,t_true.*1000,marked_true_events');
        xlabel('Time [ms]')
        ylabel('Voltage')
        title("Waveforms Marked True")

        % avg waveforms true
        ax_wv_avg_true = nexttile(tl);
        IED.utils.stdshade(ax_wv_avg_true,marked_true_events,0.5,'k',t_true.*1000);
        xlabel('Time [ms]')
        ylabel('Voltage')
        title("Mean + STD Waveforms Marked True")

        % cluster map false
        ax_cluster_false = nexttile(tl);
        gscatter(reduced_data(~marks_vec,1), reduced_data(~marks_vec,2), clust_idx(~marks_vec),'kr','o','filled')
        title("Reduced Data Marked False")
        legend(ax_cluster_false,'off')

        % waveforms true
        ax_wv_false = nexttile(tl);
        plot(ax_wv_false,t_false.*1000,marked_false_events')
        xlabel('Time [ms]')
        ylabel('Voltage')
        title("All Waveforms Marked False")

        % avg waveforms true
        ax_wv_avg_false = nexttile(tl);
        IED.utils.stdshade(ax_wv_avg_false,marked_false_events,0.5,'k',t_false.*1000);
        xlabel('Time [ms]')
        ylabel('Voltage')
        title("Mean + STD Waveforms Marked False")

        % link axes and set limits
        linkaxes([ax_cluster_true, ax_cluster_false])
        linkaxes([ax_wv_true, ax_wv_avg_true, ax_wv_false, ax_wv_avg_false])
        if ~isempty(t_true)
            xlim(ax_wv_true,[min(t_true), max(t_true)].*1000)
        else
            xlim(ax_wv_false,[min(t_false), max(t_false)].*1000)
        end

        % title figure
        sgtitle(waveforms_figs(iLabel),fig_name,'Interpreter','none');
    end

    %%%%%% match the Ylimits of all axes displaying waveforms
    
    % find biggest ylimits
    biggesy_ylims = [nan nan];
    for iFig = waveforms_figs
        current_ylims = iFig.Children.Children(2).YLim; % Should be "All Waveforms Marked False"
        biggesy_ylims(1) = min(biggesy_ylims(1), current_ylims(1));
        biggesy_ylims(2) = max(biggesy_ylims(2), current_ylims(2));
    end

    % fix ylims for all figures, increase font size
    for iFig = waveforms_figs
        fontsize(iFig,"increase")
        fontsize(iFig,"increase")
        fontsize(iFig,"increase")
        iFig.Children.Children(2).YLim = biggesy_ylims; % because of linkaxes we only need to change 1 axes
    end

end