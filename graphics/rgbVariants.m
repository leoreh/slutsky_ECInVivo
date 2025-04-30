function rgb_variants = rgbVariants(baseColor, k)
% generates k rgb triplets of same hue but different saturation/brightness
%
% INPUT:
% baseColor  - 1x3 rgb triplet of base color
% k         - number of variants to generate
%
% OUTPUT:
% rgb_variants - kx3 matrix of rgb triplets
%
% 01 feb 25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert to hsv
hsvColor = rgb2hsv(baseColor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create even spacing for saturation and value
sat = linspace(0.05, 1, k);           % maintain strong colors
val = linspace(0.05, 1, k);           % bright colors for visibility

% initialize output
hsv_variants = zeros(k, 3);

% create variants by varying saturation and value but keeping same hue
for i = 1:k
    hsv_variants(i, :) = [hsvColor(1), sat(i), val(i)];
end

% convert back to rgb
rgb_variants = hsv2rgb(hsv_variants);

end