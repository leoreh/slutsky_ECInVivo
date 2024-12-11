# Targets
This repository contains code to create MATLAB's app for creating a GUI for marking events on time series data.
The guiding principle is modularity and extensibility - therefore everything is trying to be organized under the MVC pattern, 
with a central app manager that place all the pieces together. Composite is key.

# How to use
Either inherit an existing app manager, the basic one is "event_marker_gui", or create your own.
You can also copy an existing app manager and just modify the GUI layout and/or views and/or controllers.
The interface is the important part - 
Each view & controller should save the app manager & model as properties during construction, 
and have a later "set_controller / set_view" method to construct the controller/view parts.

### NOTE:
This is a work in progress, documentation is sparse. Some architectural changes are needed, but the core functionality is there.
See "event_marker_gui_v2" for a suggested improvement.

### TO DO
* Make a better unified and sensible architecture.
* Solve resizing bug - current process is recommanding all axes to have the same YTick.