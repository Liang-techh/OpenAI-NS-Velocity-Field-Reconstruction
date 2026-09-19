function T=ns_core_timeseries(dataFile,candidateId,varargin)
%NS_CORE_TIMESERIES Geometry of a fixed physical observation window over time.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
if nargin<2||isempty(candidateId),candidateId='ST054-Q2';end
p=inputParser;addParameter(p,'Times',linspace(.25,.75,13));addParameter(p,'RadiusMax',.35);addParameter(p,'ZHalfSpan',.60);addParameter(p,'Plot',true);parse(p,varargin{:});o=p.Results;
D=load(dataFile);if ~iscell(D.models),D.models=num2cell(D.models);end
ids=cellfun(@(m)m.id,D.models,'UniformOutput',false);ix=find(strcmp(ids,candidateId),1);
if isempty(ix),error('ns:Candidate','Unknown candidate.');end
names={'weighted_r_rms','weighted_z_rms','axial_to_radial_extent','inward_fraction','positive_swirl_fraction','bipolar_fraction','axial_pressure_toward_fraction','median_swirl_over_inflow','median_axial_over_inflow','sampled_residual_max'};
a=nan(numel(o.Times),numel(names));
for i=1:numel(o.Times)
 M=ns_core_metrics(D.models{ix},o.Times(i),'RadiusMax',o.RadiusMax,'ZHalfSpan',o.ZHalfSpan);
 for j=1:numel(names),a(i,j)=M.(names{j});end
end
T=array2table([o.Times(:) a],'VariableNames',[{'time'},names]);
if o.Plot
 groups={1:2,3,4:7,8:9,10};labels={'Enstrophy-weighted extents','Window-dependent extent ratio','Direction fractions','Velocity-component ratios','Sampled residual maximum'};
 for k=1:numel(groups)
  figure('Name',[candidateId ' | ' labels{k}]);plot(o.Times(:),a(:,groups{k}),'-o');grid on;xlabel('Physical time');ylabel(labels{k});
  legend(strrep(names(groups{k}),'_',' '),'Location','best');title([candidateId ' | fixed physical window'],'Interpreter','none');
 end
end
end
