function [labels_pos,body_parts] = read_label_file(labels_file)
% Takes deeplabcut output CSV file, and read it.
% INPUTS:
%       labels_file - The path to the CSV file.
% OUTPUTS:
%       labels_pos  - double 4d matrix. Position of each bodypart in each
%                     frame. Ordered: frame X mouse X body part X (x,y,likelihood)
%       body_parts  - cellstr, body parts names, order matching labels_pos
%                     3rd dimension.


labels_tab_cell = readcell(labels_file,"Range",'1:4');
% get individuals
[individuals,~,ic] = unique(labels_tab_cell(2,2:end),'stable');
nDataTypes = numel(unique(labels_tab_cell(4,2:end),'stable'));

coords = readmatrix(labels_file,"Range",'B5');
coords = [coords,nan(size(coords,1),size(labels_tab_cell,2)-size(coords,2))];

% orginaze cordinates in frame X mouse X label X (x,y,likelihood)
for iIndividual = numel(individuals):-1:1
    mouse_cords = coords(:,ic == iIndividual);
    
    for iType = nDataTypes:-1:1
        celled_cords{iIndividual}(:,:,iType) = mouse_cords(:,iType:nDataTypes:end);
    end

end
labels_pos = cat(4,celled_cords{:});
labels_pos = permute(labels_pos,[1 4 2 3]);

% get all bodyparts names by order
body_parts = unique(labels_tab_cell(3,2:end),'stable');


end