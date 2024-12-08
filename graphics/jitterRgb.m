function rgb_jittered = jitterRgb(rgbClr)

% recieves rgb triplet/s and changes it slightly via its hsv values
%
% INPUT:
%   rgb       numeric triplet of rgb. if mat will treat each row a
%             different triplet
%
% 31 mar 24 LH      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
jitFctr = [0.05, -0.05, 0.05];      % [hue, saturation, brightness]

nclrs = size(rgbClr, 1);

% convert to hsv
hsvClr = rgb2hsv(rgbClr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsv_jittered = nan(nclrs, 3);
for iclr = 1 : nclrs
    for ival = 1 : 3
        hsv_jittered(iclr, ival) = hsvClr(iclr, ival) + jitFctr(ival) * (iclr - 1);
    end
end

% ensure values are in range 0-1
hsv_jittered(:, 1) = mod(hsv_jittered(:, 1), 1);
hsv_jittered(:, 2) = min(max(hsv_jittered(:, 2), 0), 1);
hsv_jittered(:, 3) = min(max(hsv_jittered(:, 3), 0), 1);

% convert back to rgb
rgb_jittered = hsv2rgb(hsv_jittered);

end