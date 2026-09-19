function receipt=ns_release_test(outDir)
%NS_RELEASE_TEST Actual packaged MATLAB numerical/UI checks, not NS acceptance.
if nargin<1,outDir='outputs/matlab';end;if ~isfolder(outDir),mkdir(outDir);end
receipt=ns_selftest([],true,outDir);assert(receipt.numerical_pass&&receipt.gui_pass);
ns_callback_test([],outDir);
D=load(fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat'));if ~iscell(D.models),D.models=num2cell(D.models);end
for k=1:numel(D.models)
 m=D.models{k};V=ns_core_volume(m,.517,19,25,.18,.27);ix=[3 247 789 5000];
 [u,p,om,r]=ns_evaluate(m,[V.X(ix)' V.Y(ix)' V.Z(ix)'],.517); %#ok<ASGLU>
 assert(max(abs(u(:,1)-V.U(ix)'))<1e-10);assert(max(abs(p-V.P(ix)'))<1e-10);
 assert(max(abs(sqrt(sum(r.^2,2))-V.residual(ix)'))<1e-10);
 M=ns_core_metrics(m,.5,'Nr',33,'Nz',49);assert(M.weighted_r_rms>0&&M.weighted_z_rms>0&&~M.pde_validated);
 T=ns_core_timeseries([],m.id,'Times',[.25,.5,.75],'Plot',false);assert(height(T)==3);
end
a=ns_core_explorer([],'Visible','off');cleanup=onCleanup(@()a.Close());
assert(isempty(getappdata(a.Figure,'LastRenderError')));a.SetTime(.631);a.SetControl('rmax',.18);a.SetControl('zmax',.28);a.SetControl('seedZ',.10);a.SetControl('slice',0);
assert(isempty(getappdata(a.Figure,'LastRenderError')));st=a.Snapshot();assert(abs(st.t-.631)<1e-12);
a.Figure.Visible='on';drawnow;pause(1);exportapp(a.Figure,fullfile(outDir,'core_explorer.png'));clear cleanup
receipt.core_numerical_pass=true;receipt.core_gui_pass=true;receipt.release_pde_validated=false;
fid=fopen(fullfile(outDir,'release_matlab_receipt.json'),'w');cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(receipt),'char');disp(jsonencode(receipt));
end
