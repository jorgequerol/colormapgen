function varargout = colormapgen(varargin)
%  colormapgen Generates a colormap with linear luminance, linear hue and 
%              constant saturation.
%     RGB = colormapgen() generates "mauricio" colormap.
%  
%     RGB = colormapgen(M,G1,G2,H1,H2,S) generates a colormap with M
%     colors, with luminance or grayscale [0,1] from G1 to G2, with hue
%     [0,360] from H1 to H2, and with saturation S [0,1].
%
%     colormapgen() plots a sample of the generated colormap.
%
%     Example:
%     mauricio = colormapgen(64,0.2,0.9,300,-60,0.9);
%


    if (nargin <= 0)
        M = 64;
        G1 = 0.05;
        G2 = 0.85;
        H1 = 360;
        H2 = -60;
        s = 0.95;
    else
        M = varargin{1};
        G1 = varargin{2};
        G2 = varargin{3};
        H1 = varargin{4};
        H2 = varargin{5};
        s = varargin{6};
    end
    
    m = 4*M;
    h = linspace(0,360,m)';
    l = linspace(0,1,m)';

    % Get F(G,H) -> L obtaining true G values
    [HH,LL] = meshgrid(h, l); clear h l;
    SS = s*ones(m);
    HSL = cat(3, HH, SS, LL);
    RGB = hsl2rgb(HSL); clear HSL;
    GG = rgb2gray(RGB); clear RGB;
    try
        F = scatteredInterpolant(GG(:),HH(:),LL(:));
    catch
        % Older MATLAB versions do not support "scatteredInterpolant"
        F = TriScatteredInterp(GG(:),HH(:),LL(:));
    end
    clear HH SS LL GG;

    % Generate colormap
    G = linspace(G1,G2,M)';
    H = mod(linspace(H1,H2,M)',360);
    S = s*ones(M,1); clear s;
    L = F(G,H); clear G F;
    RGB = squeeze(abs(hsl2rgb(permute([H S L],[1 3 2]))));  clear H S L;
    
    if(nargout > 0)
        varargout{1} = RGB;
    else
        figure(gcf());
        imagesc(linspace(0,1,M)');
        colormap(gcf, RGB);
        set(gca, 'YDir', 'Normal', 'XTick', []);
    end
end

function RGB = hsl2rgb(HSL)
% hsl2rgb Converts HSL to RGB colors

% This function was derived from Pascal Getreuer's "colorspace".

    H = HSL(:,:,1);
    S = HSL(:,:,2);
    L = HSL(:,:,3);
    Delta = S.*min(L,1-L);
    m0 = L-Delta;
    m2 = L+Delta;
    N = size(H);
    H = min(max(H(:),0),360)/60;
    m0 = m0(:);
    m2 = m2(:);
    F = H - round(H/2)*2;
    M = [m0, m0 + (m2-m0).*abs(F), m2];
    Num = length(m0);
    j = [2 1 0;1 2 0;0 2 1;0 1 2;1 0 2;2 0 1;2 1 0]*Num;
    k = floor(H) + 1;
    RGB = reshape([M(j(k,1)+(1:Num).'),M(j(k,2)+(1:Num).'),M(j(k,3)+(1:Num).')],[N,3]);
end
