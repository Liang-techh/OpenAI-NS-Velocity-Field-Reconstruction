function app=ns_core_explorer(dataFile,varargin)
%NS_CORE_EXPLORER Dense local sampling, continuous time and measured geometry.
% Frozen spectral models only. Display sliders never change coefficients.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
p=inputParser;addParameter(p,'Visible','on');parse(p,varargin{:});D=load(dataFile);
if ~isfield(D,'models'),error('ns:Data','Missing frozen models');end
if ~iscell(D.models),D.models=num2cell(D.models);end
models=D.models;ids=cellfun(@(q)q.id,models,'UniformOutput',false);ix=find(strcmp(ids,'ST054-Q2'),1);if isempty(ix),ix=1;end
m=models{ix};if ~strcmp(m.mode,'spectral')||m.pde_validated,error('ns:Model','Expected an unvalidated frozen spectral field');end
s=struct('t',.5,'rmax',.35,'zmax',.60,'slice',0,'seedR',.14,'seedZ',.12,'count',200,'arc',1.5,'iso',.35,'alpha',.20,'colorScale',.25);
V=[];cache='';busy=false;closed=false;last=tic;sl=struct();labs=struct();
fig=uifigure('Name','NS central core | exact local re-sampling','Position',[30 30 1480 960],'Visible',p.Results.Visible);
fig.CloseRequestFcn=@closeApp;
g=uigridlayout(fig,[5 3]);g.RowHeight={42,'1x',190,72,30};
g.ColumnWidth={290,'1x','1x'};g.Padding=[10 10 10 10];
h=uilabel(g,'Text','CENTRAL CORE | frozen ST054 field; geometry diagnostics are not NS acceptance','FontSize',17,'FontWeight','bold');h.Layout.Row=1;h.Layout.Column=[1 3];
panel=uipanel(g,'Title','Observation and seeding','Scrollable','on');panel.Layout.Row=[2 3];panel.Layout.Column=1;
cg=uigridlayout(panel,[15 2]);cg.RowHeight=repmat({32},1,15);cg.ColumnWidth={128,'1x'};cg.Padding=[6 6 6 6];
label('Candidate',1);candidate=uidropdown(cg,'Items',ids,'Value',m.id,'ValueChangedFcn',@changeCandidate);place(candidate,1,2);
label('Scalar / line color',2);scalar=uidropdown(cg,'Items',{'Speed','Vorticity','Swirl velocity','Axial velocity','Pressure','Axial pressure force','Full NS residual'},'Value','Speed','ValueChangedFcn',@refresh);place(scalar,2,2);
label('3-D mode',3);mode=uidropdown(cg,'Items',{'Streamlines','Streamlines + surface','Vorticity surface'},'Value','Streamlines','ValueChangedFcn',@refresh);place(mode,3,2);
sl.rmax=control('Core radius',4,[.12 .8],'rmax');sl.zmax=control('Core z half-span',5,[.15 1.2],'zmax');
sl.slice=control('Slice offset y',6,[-.8 .8],'slice');sl.seedR=control('Seed radius',7,[.01 .7],'seedR');sl.seedZ=control('Seed z half-span',8,[.01 1.1],'seedZ');
sl.count=control('Seed count',9,[20 320],'count');sl.arc=control('Total arc cap',10,[.2 4],'arc');
sl.iso=control('|omega| level',11,[.01 3],'iso');sl.alpha=control('Surface opacity',12,[.03 .6],'alpha');
sl.colorScale=control('Fixed color gain',13,[.05 2],'colorScale');
quality=uidropdown(cg,'Items',{'Fast','Balanced','Fine'},'Value','Balanced','ValueChangedFcn',@refresh);label('Local quality',14);place(quality,14,2);
btn=uibutton(cg,'Text','Full-field window','ButtonPushedFcn',@(src,e)ns_explorer(dataFile));place(btn,15,1);
btn=uibutton(cg,'Text','Time-series plots','ButtonPushedFcn',@series);place(btn,15,2);
ax3=uiaxes(g);ax3.Layout.Row=2;ax3.Layout.Column=2;view(ax3,-37,22);grid(ax3,'on');colormap(ax3,parula(256));b3=colorbar(ax3);
ax2=uiaxes(g);ax2.Layout.Row=2;ax2.Layout.Column=3;colormap(ax2,parula(256));b2=colorbar(ax2);
tab=uitable(g,'ColumnName',{'Local diagnostic','Value'},'ColumnEditable',[false false]);tab.Layout.Row=3;tab.Layout.Column=[2 3];
fg=uigridlayout(g,[2 4]);fg.Layout.Row=4;fg.Layout.Column=[1 3];fg.RowHeight={25,35};fg.ColumnWidth={100,'1x',160,110};fg.Padding=[0 0 0 0];
play=uibutton(fg,'Text','Play','ButtonPushedFcn',@toggle);place(play,1,1);
l=uilabel(fg,'Text','Drag time: low-detail preview. Release: locally re-sample. Streamlines are instantaneous, not particle paths.');place(l,1,2);
clock=uilabel(fg,'Text','t = 0.50000','FontWeight','bold');place(clock,1,3);
btn=uibutton(fg,'Text','Export PNG','ButtonPushedFcn',@png);place(btn,1,4);
time=uislider(fg,'Limits',[m.tmin m.tmax],'Value',s.t,'MajorTicks',linspace(m.tmin,m.tmax,5),...
    'ValueChangingFcn',@(src,e)setTimeInternal(e.Value,true),'ValueChangedFcn',@(src,e)setTimeInternal(e.Value,false));time.Layout.Row=2;time.Layout.Column=[1 4];
status=uilabel(g,'Text','Preparing local data...','FontSize',11);status.Layout.Row=5;status.Layout.Column=[1 3];
timerObject=timer('ExecutionMode','fixedSpacing','Period',.2,'BusyMode','drop','TimerFcn',@tick);
app=struct('Figure',fig,'Axes3D',ax3,'AxesSlice',ax2,'Table',tab,'TimeSlider',time,'PlayButton',play,'Status',status,...
 'SetTime',@setTime,'SetControl',@setControl,'SetCandidate',@setCandidate,'SetMode',@setMode,'SetScalar',@setScalar,...
 'Snapshot',@snapshot,'Refresh',@()render(false),'Close',@()closeApp());
render(false);
    function place(h,r,c),h.Layout.Row=r;h.Layout.Column=c;end
    function label(text,row),h=uilabel(cg,'Text',text,'FontSize',11);place(h,row,1);end
    function h=control(text,row,limits,key)
        l=uilabel(cg,'Text',sprintf('%s: %.3g',text,s.(key)),'FontSize',11);place(l,row,1);labs.(key)={l,text};
        h=uislider(cg,'Limits',limits,'Value',s.(key),'MajorTicks',[],'ValueChangedFcn',@(src,e)changeControl(key,e.Value));place(h,row,2);
    end
    function clampSeeds()
        s.seedR=min(s.seedR,.95*s.rmax);s.seedZ=min(s.seedZ,.95*s.zmax);s.slice=max(-s.rmax,min(s.rmax,s.slice));
        keys={'seedR','seedZ','slice'};
        for j=1:numel(keys),k=keys{j};sl.(k).Value=s.(k);q=labs.(k);q{1}.Text=sprintf('%s: %.3g',q{2},s.(k));end
    end
    function changeControl(key,value)
        if strcmp(key,'count'),value=round(value);end
        s.(key)=value;q=labs.(key);q{1}.Text=sprintf('%s: %.3g',q{2},value);clampSeeds();render(false);
    end
    function setControl(key,value)
        if ~isfield(sl,key)||value<sl.(key).Limits(1)||value>sl.(key).Limits(2),error('ns:Control','Invalid control');end
        sl.(key).Value=value;changeControl(key,value);
    end
    function setTime(value),stop(timerObject);play.Text='Play';time.Value=value;setTimeInternal(value,false);end
    function setTimeInternal(value,preview)
        s.t=max(m.tmin,min(m.tmax,value));clock.Text=sprintf('t = %.5f',s.t);
        if preview&&toc(last)<.12,return,end
        render(preview);
    end
    function setCandidate(value),candidate.Value=value;changeCandidate();end
    function changeCandidate(varargin),stop(timerObject);play.Text='Play';m=models{find(strcmp(ids,candidate.Value),1)};cache='';render(false);end
    function setMode(value),mode.Value=value;render(false);end
    function setScalar(value),scalar.Value=value;render(false);end
    function refresh(varargin),render(false);end
    function q=snapshot(),q=s;q.candidate=m.id;q.mode=mode.Value;q.scalar=scalar.Value;q.pde_validated=false;if ~isempty(V),q.grid_spacing=[V.x(2)-V.x(1),V.z(2)-V.z(1)];end;end
    function toggle(varargin)
        if strcmp(timerObject.Running,'on'),stop(timerObject);play.Text='Play';render(false);
        else,if s.t>=m.tmax,s.t=m.tmin;end;play.Text='Pause';start(timerObject);end
    end
    function tick(varargin)
        if closed,return,end
        s.t=min(m.tmax,s.t+.005);time.Value=s.t;render(true);
        if s.t>=m.tmax,stop(timerObject);play.Text='Play';render(false);end
    end
    function series(varargin),stop(timerObject);play.Text='Play';ns_core_timeseries(dataFile,m.id,'RadiusMax',s.rmax,'ZHalfSpan',s.zmax);end
    function png(varargin)
        stop(timerObject);play.Text='Play';render(false);[name,path]=uiputfile('*.png','Export core view','core_view.png');if isequal(name,0),return,end;exportapp(fig,fullfile(path,name));
    end
    function render(preview)
        if closed||busy,return,end
        busy=true;guard=onCleanup(@unlock);started=tic;
        try
            if preview,nxy=25;nz=33;else
                switch quality.Value,case 'Fast',nxy=33;nz=45;case 'Fine',nxy=65;nz=89;otherwise,nxy=49;nz=65;end
            end
            key=sprintf('%s %.12g %.12g %.12g %d %d',m.id,s.t,s.rmax,s.zmax,nxy,nz);
            if ~strcmp(key,cache),V=ns_core_volume(m,s.t,nxy,nz,s.rmax,s.zmax);cache=key;end
            count=round(s.count);if preview,count=min(count,48);end
            [lines,stats]=ns_core_lines(V,s.seedR,s.seedZ,count,s.arc);
            [az,el]=view(ax3);cla(ax3);hold(ax3,'on');
            if contains(mode.Value,'surface')&&s.iso>min(V.vorticity(:))&&s.iso<max(V.vorticity(:))
                surf=isosurface(V.X,V.Y,V.Z,V.vorticity,s.iso);
                if ~isempty(surf.faces),h=patch(ax3,surf,'FaceColor',[.3 .62 .72],'EdgeColor','none','FaceAlpha',s.alpha);isonormals(V.X,V.Y,V.Z,V.vorticity,h);camlight(ax3,'headlight');lighting(ax3,'gouraud');end
            end
            switch scalar.Value
                case 'Speed',values=V.speed;limits=[0 .8];case 'Vorticity',values=V.vorticity;limits=[0 5];
                case 'Swirl velocity',values=V.swirl;limits=[-.7 .7];case 'Axial velocity',values=V.W;limits=[-.7 .7];
                case 'Pressure',values=V.P;limits=[-.2 .2];case 'Axial pressure force',values=V.PforceZ;limits=[-.5 .5];otherwise,values=V.residual;limits=[0 .06];
            end
            if contains(mode.Value,'Streamlines')
                I=griddedInterpolant({V.y,V.x,V.z},values,'linear','none');q=[];cc=[];
                for j=1:numel(lines),v=lines{j};if size(v,1)<2,continue,end;q=[q;v;NaN(1,3)];cc=[cc;I(v(:,2),v(:,1),v(:,3));NaN];end %#ok<AGROW>
                if ~isempty(q),surface(ax3,[q(:,1)';q(:,1)'],[q(:,2)';q(:,2)'],[q(:,3)';q(:,3)'],[cc';cc'],'FaceColor','none','EdgeColor','interp','LineWidth',.65);end
            end
            xlabel(ax3,'x');ylabel(ax3,'y');zlabel(ax3,'z');axis(ax3,'equal');grid(ax3,'on');xlim(ax3,[-s.rmax s.rmax]);ylim(ax3,[-s.rmax s.rmax]);zlim(ax3,[-s.zmax s.zmax]);view(ax3,az,el);hold(ax3,'off');caxis(ax3,limits*s.colorScale);b3.Label.String=scalar.Value;
            title(ax3,sprintf('%s | t=%.5f\nLocally re-sampled instantaneous flow',m.id,s.t),'Interpreter','none');
            nn=151;if preview,nn=61;end
            x=linspace(-s.rmax,s.rmax,nn);z=linspace(-s.zmax,s.zmax,nn);[X,Z]=meshgrid(x,z);points=[X(:),s.slice+zeros(numel(X),1),Z(:)];
            [u,pres,om,res,pf]=ns_evaluate(m,points,s.t);rad=hypot(points(:,1),points(:,2));
            switch scalar.Value
                case 'Speed',v=vecnorm(u,2,2);case 'Vorticity',v=vecnorm(om,2,2);case 'Swirl velocity',v=(-points(:,2).*u(:,1)+points(:,1).*u(:,2))./max(rad,eps);
                case 'Axial velocity',v=u(:,3);case 'Pressure',v=pres;case 'Axial pressure force',v=pf(:,3);otherwise,v=vecnorm(res,2,2);
            end
            cla(ax2);imagesc(ax2,x,z,reshape(v,size(X)));ax2.YDir='normal';hold(ax2,'on');k=1:max(1,round(nn/16)):nn;ux=reshape(u(:,1),size(X));uz=reshape(u(:,3),size(X));quiver(ax2,X(k,k),Z(k,k),.15*ux(k,k),.15*uz(k,k),0);hold(ax2,'off');axis(ax2,'equal');xlim(ax2,[-s.rmax s.rmax]);ylim(ax2,[-s.zmax s.zmax]);xlabel(ax2,'x');ylabel(ax2,'z');title(ax2,sprintf('%s | exact slice y=%.3g',scalar.Value,s.slice));caxis(ax2,limits*s.colorScale);b2.Label.String=scalar.Value;
            M=ns_core_metrics(m,s.t,'RadiusMax',s.rmax,'ZHalfSpan',s.zmax,'Volume',V);
            tab.Data={'Local enstrophy-weighted r RMS',M.weighted_r_rms;'Local enstrophy-weighted z RMS',M.weighted_z_rms;'Local axial / radial RMS ratio',M.axial_to_radial_extent;'Inward / positive swirl fractions',sprintf('%.4f / %.4f',M.inward_fraction,M.positive_swirl_fraction);'Bipolar / axial-pressure fractions (|z| > .03)',sprintf('%.4f / %.4f',M.bipolar_outflow_fraction,M.axial_pressure_toward_fraction);'Median |u_theta| / |u_r|',M.median_swirl_over_inflow;'Median |u_z| / |u_r|',M.median_axial_over_inflow;'Forward streamline median turn (radians)',stats.median_turn_angle;'Forward streamline mean end/start radius',stats.mean_end_radius_fraction;'Local sampled max |R| (NOT acceptance)',M.residual_max};
            clock.Text=sprintf('t = %.5f',s.t);status.Text=sprintf('%dx%dx%d LOCAL samples | %.2f s | Window-dependent geometry, not source correspondence. Original NS target NOT met.',nxy,nxy,nz,toc(started));
            setappdata(fig,'LastRenderError','');drawnow limitrate nocallbacks;last=tic;
        catch err,setappdata(fig,'LastRenderError',err.message);status.Text=['Render error: ' err.message];warning('ns:CoreRender','%s',getReport(err,'basic'));end
        clear guard
    end
    function unlock(),busy=false;end
    function closeApp(varargin),if closed,return,end;closed=true;if isvalid(timerObject),stop(timerObject);delete(timerObject);end;if isvalid(fig),delete(fig);end;end
end
