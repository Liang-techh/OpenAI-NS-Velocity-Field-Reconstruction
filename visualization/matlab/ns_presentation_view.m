function app=ns_presentation_view(varargin)
%NS_PRESENTATION_VIEW Reproduce the selected full-field viewing preset.
% The slice is moved onto the axis. Frozen coefficients remain unchanged.
app=ns_explorer([],varargin{:});app.SetControl('radius',.621);
app.SetControl('seedz',-.221);app.SetControl('extent',.938);
app.SetControl('colorScale',.932);app.SetControl('slice',0);app.SetTime(.5);
end
