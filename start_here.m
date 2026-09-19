function app=start_here(mode)
%START_HERE Launch the frozen-field visualization from the repository root.
% start_here('full'), start_here('core'), or start_here('presentation').
if nargin<1,mode='presentation';end
prepare_matlab;
addpath(fullfile(fileparts(mfilename('fullpath')),'visualization','matlab'));
switch mode
    case 'core',app=ns_core_explorer;
    case 'full',app=ns_explorer;
    case 'presentation',app=ns_presentation_view;
    otherwise,error('ns:Mode','Use full, core or presentation.');
end
end
