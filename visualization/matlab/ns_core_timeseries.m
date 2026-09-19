function T=ns_core_timeseries(dataFile,candidateId,varargin)
%NS_CORE_TIMESERIES Fixed-observation-window diagnostics, not target errors.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
if nargin<2||isempty(candidateId),candidateId='ST054-Q2';end
p=inputParser;addParameter(p,'Times',linspace(.25,.75,13));addParameter(p,'RadiusMax',.35);
addParameter(p,'ZHalfSpan',.60);addParameter(p,'Plot',true);addParameter(p,'Nxy',49);addParameter(p,'Nz',65);parse(p,varargin{:});o=p.Results;
D=load(dataFile);if ~iscell(D.models),D.models=num2cell(D.models);end
ids=cellfun(@(q)q.id,D.models,'UniformOutput',false);ix=find(strcmp(ids,candidateId),1);
if isempty(ix),error('ns:Candidate','Unknown candidate');end
m=D.models{ix};times=o.Times(:);fields={'weighted_r_rms','weighted_z_rms','axial_to_radial_extent','inward_fraction','positive_swirl_fraction','bipolar_outflow_fraction','axial_pressure_toward_fraction','median_swirl_over_inflow','median_axial_over_inflow'};
values=NaN(numel(times),numel(fields));
for i=1:numel(times)
    M=ns_core_metrics(m,times(i),'RadiusMax',o.RadiusMax,'ZHalfSpan',o.ZHalfSpan,'Nxy',o.Nxy,'Nz',o.Nz);
    for j=1:numel(fields),values(i,j)=M.(fields{j});end
end
T=array2table([times values],'VariableNames',[{'time'},fields]);
if o.Plot
    groups={1:2,3,4:7,8:9};titles={'Local vorticity-weighted RMS extents','Local axial/radial RMS ratio','Finite-cylinder direction fractions','Median component ratios (inward samples)'};
    for k=1:4
        figure('Name',[candidateId ' | ' titles{k}],'Color','w');plot(times,values(:,groups{k}),'-o');grid on;
        xlabel('Physical time');legend(strrep(fields(groups{k}),'_',' '),'Location','best');title(titles{k});
    end
end
end
