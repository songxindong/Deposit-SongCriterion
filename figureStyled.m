function h_fig = figureStyled(style, varargin)
% return a new figure handle and configure figure format properties to certain preset style
% h_fig:        figure handle
% 'Style':      style name, choose from: 
%               'Science��,  'PNAS',     'Elsevier',
%               'Nature',   'JNeurosci',   
% 'Scaling':    zoom in scaling factor, should be a non-zero small integer
% 'Position':   position of the figure, should be a 4-number vector in cm
% 'Column':     column # in the figure, typically, this number is 1, 1.5, 2
%               this property can overwrite position's width
% 
global p

%% Input Parser
p = inputParser;
validStyles =    	{   'Nature',...
                        'Science',...
                        'PNAS',...
                        'Elsevier', ...
                        'Thesis'};
checkStyle =        @(x) any(validatestring(x, validStyles));
defaultScaling =    1;
defaultPosition =   [   1       1       10      10]; 
                    %   left    bottom  width   height   
addRequired(    p,  'Style',    checkStyle);
addParameter(   p,  'Scaling',  defaultScaling,     @isnumeric);
addParameter(   p,  'Position', defaultPosition,    @isnumeric);
addParameter(   p,  'Column',   NaN,                @isnumeric);

parse(p, style, varargin{:});

%% Scaling
%   Check more in 
%   http://undocumentedmatlab.com/blog/graphic-sizing-in-matlab-r2015b
%   For most LCD screens, get(0, 'ScreenPixelsPerInch') =       96
%   For Dell U2713HM,   the physical number is 2560/23.49 =	    109
%   For BENQ BL3200PT,  the physical number is 2560/27.9 =      91.76
%   For Thinkpad Helix, the physical number is                  190
%   For Thinkpad X1yoga,the physical number is 2560/(14*cos(atan(9/16))) = 209.8
set(0,  'ScreenPixelsPerInch',  96* p.Results.Scaling);
set(0,  'DefaultAxesFontName',  'Arial');

%% Consistent Settings
h_fig = figure;
    set(h_fig,  'Units',                'Centimeters'); % always use Centimeters
    set(h_fig,  'defaultAxesUnits',     'Centimeters'); % for measurement consistency
    set(h_fig,  'defaultAxesFontSize',  7);         % Font Size is 7pt
	set(h_fig,  'defaultTextFontSize',  7);         % Font Size is 7pt
    set(h_fig,  'Color',                [1 1 1]);   % Background color: pure white
    set(h_fig,  'defaultAxesBox',       'on');      % box on
    set(h_fig,  'defaultAxesLayer',     'top');    	% put axes on the top
    set(h_fig,  'defaultAxesNextPlot',  'add');     % hold on for more components
    set(h_fig,  'defaultAxesColorOrder',[0 0 0]);   % default color: pure black            
    % set(h_fig,  'defaultAxesClipping',  'off');

%% Style related Settings
switch style
    case 'Nature'           % as of 2019
        % check https://www.nature.com/documents/nature-final-artwork.pdf
        % and https://www.nature.com/documents/nature-extended-data.pdf
        set(h_fig,  'defaultTextFontName',  'Helvetica');        
        set(h_fig,  'defaultAxesFontName',  'Helvetica');
            % "All text should be in a sans-serif typeface, preferably Helvetica or Arial."      
            % Separate panels in multi-part gures should be labelled with 8 pt bold,
            %   upright (not italic) a, b, c��s...
            % Maximum text size for all other text: 7pt.
            % Minimal text size: 5pt.
        % Figure Size & Unit & Resolution
        switch p.Results.Column
            case 1;     width = 8.9;
            case 1.5;   width = 12.0;    % 12.0 or 13.6 cm
            case 2;     width = 18.3;
            otherwise;  width = p.Results.Position(3);
        end
            % Provide fles at about the size they are to be printed. 
            % Nature��s standard figure sizes are 89 mm (single column) and
            % 183 mm (double column). The full depth of the page is 247 mm.
            % Figures can also be 1.5 columns where necessary �C 120 mm or 136 mm.
       % Line Size 
        set(h_fig,  'defaultAxesLineWidth', 0.75);
        set(h_fig,  'defaultLineLineWidth', 0.75);
            % Lines and strokes should be set between 0.25 and 1 pt.

    case 'Science'          % as of 2014
        % check http://www.sciencemag.org/site/feature/contribinfo/prep/prep_revfigs.xhtml
        % Font
        set(h_fig,  'defaultTextFontName',  'Helvetica');        
        set(h_fig,  'defaultAxesFontName',  'Helvetica');
            % "Use a sans-serif font whenever possible (we prefer Helvetica)."      
            % body text and figure legend in Science are 8pt
            % Science Signaling requires this to be 7pt
            % figure part letter should be bold 9pt A B C 
        % Figure Size & Unit & Resolution
        switch p.Results.Column
            case 1;     width = 5.5;
            case 1.5;   width = 9.1;    % estimated
            case 2;     width = 12.0;
            otherwise;  width = p.Results.Position(3);
        end
            % The width of figures, when printed, will usually be 
            % 5.5 cm (2.25 inches or 1 column) or 
            % 12.0 cm (4.75 inches or 2 columns)
            % initial submission requires 300dpi
            % revised submission requires 1200dpi
       % Line Size 
        set(h_fig,  'defaultAxesLineWidth', 1);
        set(h_fig,  'defaultLineLineWidth', 1);
            % prepare to  have 50% or even 33% reduction
            % the minimum linewidth after reduction should be 0.5 pt   

    case 'PNAS'             % as of 2015
        % check http://www.pnas.org/site/misc/digitalart.pdf
        set(h_fig,  'defaultTextFontName',  'Arial');        
        set(h_fig,  'defaultAxesFontName',  'Arial');
            % Use standard fonts such as Arial, Helvetica, Times, Symbol...    
            % Size any text in your figure to at least 6-8 points.             
            % figure part letter should be bold 9pt A B C  
            % PNAS main test                font size is 9 pt
            % PNAS materials and methods    font size is 7 pt
        % Figure Size & Unit & Resolution
        switch p.Results.Column
            case 1;     width = 8.7;
            case 1.5;   width = 11.4;
            case 2;     width = 17.8;
            otherwise;  width = p.Results.Position(3);
        end;
            % on of the following widths: 1, 1.5, or 2 columns 
            % 1     column,     20.5    picas   3.42"   8.7     cm
            % 1.5   columns,    27      picas   4.5"    11.4    cm
            % 2     columns     42.125  picas   7"      17.8    cm                        
            % initial submission requires 300dpi
            % revised submission requires 1200dpi 
    	% Line Size
        set(h_fig,  'defaultAxesLineWidth', 1);
        set(h_fig,  'defaultLineLineWidth', 1);
            % It is preferable to have graph lines at least 0.25 pt wide    
            
	case 'Elsevier'
        % as of 2016
        % check https://www.elsevier.com/authors/author-schemas/artwork-and-media-instructions
        % Font                                  (as in artwork-sizing) 
        set(h_fig,  'defaultTextFontName',  'Arial');   % in overview    
        set(h_fig,  'defaultAxesFontName',  'Arial');   % in overview
            % Arial (or Helvetica), Courier, Symbol, Times(or Times New Roman) 
            % As a general rule, the lettering on the artwork should have a 
            % finished, printed size of 7 pt for normal text and no smaller
            % than 6 pt for subscript and superscript characters. Smaller 
            % lettering will yield text that is hardly legible. This is a 
            % rule-of-thumb rather than a strict rule. There are instances 
            % where other factors in the artwork (e.g., tints and shadings)
            % dictate a finished size of perhaps 10 pt.
        % Figure Size & Unit & Resolution       (as in artwork-sizing)
    	switch p.Results.Column
            case 1;     width = 9.0;
            case 1.5;   width = 14.0;
            case 2;     width = 19.0;
            otherwise;  width = p.Results.Position(3); 
    	end;
            % 1     column,                             9.0     cm
            % 1.5   columns,                            14.0    cm
            % 2     columns,                            19.0    cm   
        % Line Size                 (as in artwork-sizing & artwork-types)
        set(h_fig,  'defaultAxesLineWidth', 1);
        set(h_fig,  'defaultLineLineWidth', 1);
            % When Elsevier decides on the size of a line art graphic, in 
            % addition to the lettering, there are several other factors to 
            % assess. 
            % Line weights range from 0.10 pt to 1.5 pt
    case 'Thesis'
        set(h_fig,  'defaultTextFontName',  'Arial');   % in overview    
        set(h_fig,  'defaultAxesFontName',  'Arial');   % in overview    	
        switch p.Results.Column
            case 1;     width = (8.5-1.0-1.5)*2.54 /2;
            case 2;     width = (8.5-1.0-1.5)*2.54;
            otherwise;  width = p.Results.Position(3); 
    	end;
        set(h_fig,  'defaultAxesLineWidth', 1);
        set(h_fig,  'defaultLineLineWidth', 1);
        
            
    case 'JNeurosci'
        
    otherwise
                
end

%% Position
posi =      p.Results.Position;
posi(3) =   width;
set(h_fig,  'Position', posi);

%% Axes things         
%         set(h_axes, 'TickDir',              'Out');
%         set(h_axes, 'TickLength',           [0.015 0.02]);
%         h = xlabel('');
%         set(h,      'VerticalAlignment',    'Middle');
%               Tightness:  'Bottom' > 'Baseline' > 'Middle' > 'Cap' > 'Top'
%         h = ylabel('');
%         set(h,      'VerticalAlignment',    'Baseline');
%               Tightness:  'Bottom' < 'Baseline' < 'Middle' < 'Cap' < 'Top'

