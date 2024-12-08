classdef event_marker_views < handle
% abstract class, defining what methods any view must have
properties (Abstract)
    app_base
    model
end
methods (Abstract)
    set_view
    % add constant view elements to app
    
    data_loaded_responce
    % to be called outside of view object, in a listener responce to model
    % "data_loaded" event.
    % define how the view should be updated when data_sources enter /
    % removed from model.
    % note that the model expect "all change at once" and not "1 data
    % source at a time", i.e. this function should treat all data_sources
    % in the model everytime it is called.
    % if view have no responce, create an empty function.

    events_changed_responve
    % to be called outside of view object, in a listener responce to model
    % "events_changed" event.
    % define how the view should be updated when base events are changed in
    % the model - event times are added / removed, event types are added /
    % removed
end

end