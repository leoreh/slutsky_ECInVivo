classdef event_map_time_win_fixed < event_map


methods
    function create_global_map(app, map_ax_parent, map_ax_row, map_ax_col)
        create_global_map@event_map(app, map_ax_parent, map_ax_row, map_ax_col)

        app.ROI_time.InteractionsAllowed = "translate";
    end

    function map_lims_changed(app)
        map_lims_changed@event_map(app)
        
        % only allow translation, overriding super behavior
        app.ROI_time.InteractionsAllowed = "translate";
    end
end

end