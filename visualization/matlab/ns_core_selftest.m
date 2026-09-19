function receipt=ns_core_selftest(dataFile,outDir)
%NS_CORE_SELFTEST Native tests for local sampling, metrics, actual GUI controls.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
if nargin<2,outDir=fullfile('outputs','matlab_core');end
if ~exist(outDir,'dir'),mkdir(outDir);end
D=load(dataFile);if ~iscell(D.models),D.models=num2cell(D.models);end
errors=[];
for i=1:numel(D.models)
    m=D.models{i};V=ns_core_volume(m,.537,17,25,.25,.4);ix=[12 258 612 1742 3900];
    x=[V.X(ix)',V.Y(ix)',V.Z(ix)'];[u,p,w,r]=ns_evaluate(m,x,.537);
    e=max([max(abs(u(:,1)-V.U(ix)')),max(abs(u(:,2)-V.V(ix)')),max(abs(u(:,3)-V.W(ix)')),max(abs(p-V.P(ix)')),max(abs(vecnorm(w,2,2)-V.vorticity(ix)')),max(abs(vecnorm(r,2,2)-V.residual(ix)'))]);
    errors(end+1)=e;assert(e<1e-9);assert(max(abs(V.divergence(:)))<1e-9); %#ok<AGROW>
    M=ns_core_metrics(m,.537,'RadiusMax',.25,'ZHalfSpan',.4,'Volume',V);
    assert(isfinite(M.axial_to_radial_extent)&&M.sample_count< numel(V.U));assert(~M.pde_validated);
    [lines,stats]=ns_core_lines(V,.12,.1,24,1.2);assert(stats.valid_count>0);
    for k=1:numel(lines),q=lines{k};if size(q,1)>1,assert(sum(vecnorm(diff(q),2,2))<=1.2+1e-8);end;end
    bad=false;try,ns_core_lines(V,.5,.1,24,1);catch,bad=true;end;assert(bad);
end
T=ns_core_timeseries(dataFile,'ST054-Q2','Times',[.25,.537,.75],'Plot',false,'Nxy',17,'Nz',25);assert(height(T)==3);assert(all(isfinite(T.axial_to_radial_extent)));
a=ns_core_explorer(dataFile,'Visible','off');guard=onCleanup(@()a.Close());assert(isempty(getappdata(a.Figure,'LastRenderError')));
s0=a.Snapshot();a.SetControl('rmax',.2);a.SetControl('zmax',.35);s1=a.Snapshot();assert(s1.grid_spacing(1)<s0.grid_spacing(1));
a.SetTime(.41937);a.SetScalar('Axial pressure force');a.SetMode('Streamlines + surface');assert(isempty(getappdata(a.Figure,'LastRenderError')));
cb=a.TimeSlider.ValueChangingFcn;cb(a.TimeSlider,struct('Value',.527));a.TimeSlider.Value=.527;cb=a.TimeSlider.ValueChangedFcn;cb(a.TimeSlider,struct('Value',.527));st=a.Snapshot();assert(abs(st.t-.527)<1e-12);
cb=a.PlayButton.ButtonPushedFcn;cb(a.PlayButton,[]);pause(1);cb(a.PlayButton,[]);st=a.Snapshot();assert(st.t>.527);
a.SetCandidate('ST054-M3');a.SetCandidate('ST054-Q2');a.SetControl('rmax',.35);a.SetControl('zmax',.6);a.SetScalar('Speed');a.SetMode('Streamlines');a.SetTime(.5);
assert(isempty(getappdata(a.Figure,'LastRenderError')));a.Figure.Visible='on';drawnow;pause(1);exportapp(a.Figure,fullfile(outDir,'core_explorer.png'));
receipt=struct('matlab_version',version,'local_sampling_pass',true,'metrics_pass',true,'arc_caps_pass',true,'time_series_pass',true,'gui_pass',true,'time_callbacks_pass',true,'playback_pass',true,'maximum_reference_error',max(errors),'pde_validated',false,'scope','Software/finite diagnostic tests only, not source geometry or PDE certification');
fid=fopen(fullfile(outDir,'core_receipt.json'),'w');closer=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(receipt),'char');fprintf('%s\n',jsonencode(receipt));clear guard
end
