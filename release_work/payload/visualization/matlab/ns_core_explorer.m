function app=ns_core_explorer(dataFile,varargin)
%NS_CORE_EXPLORER True local resampling with explicit physical window controls.
% All geometry is frozen; sliders change time or visualization only.
if nargin<1||isempty(dataFile),dataFile=fullfile(fileparts(mfilename('fullpath')),'data','st054_models.mat');end
p=inputParser;addParameter(p,'Visible','on');parse(p,varargin{:});D=load(dataFile);
if ~iscell(D.models),D.models=num2cell(D.models);end
models=D.models;ids=cellfun(@(m)m.id,models,'UniformOutput',false);m=models{find(strcmp(ids,'ST054-Q2'),1)};
s=struct('t',.5,'rmax',.35,'zmax',.60,'seedR',.16,'seedZ',.18,'count',120,'arc',1.2,'slice',0);
fig=uifigure('Name','NS core explorer | frozen research field','Position',[40 40 1440 920],'Visible',p.Results.Visible);
g=uigridlayout(fig,[4 3]);g.ColumnWidth={265,'1x','1x'};g.RowHeight={38,'1x',155,75};
h=uilabel(g,'Text','CENTRAL CORE | exact local resampling; physical aspect ratio; not PDE acceptance','FontSize',16,'FontWeight','bold');h.Layout.Row=1;h.Layout.Column=[1 3];
cg=uigridlayout(g,[11 2]);cg.Layout.Row=2;cg.Layout.Column=1;cg.RowHeight=repmat({38},1,11);cg.ColumnWidth={125,'1x'};
lab('Candidate',1);dd=uidropdown(cg,'Items',ids,'Value',m.id,'ValueChangedFcn',@switchModel);place(dd,1,2);
lab('Slice quantity',2);sc=uidropdown(cg,'Items',{'Speed','Vorticity','Axial pressure force','Full NS residual'},'Value','Vorticity','ValueChangedFcn',@refresh);place(sc,2,2);
sl.rmax=ctrl('Core radius',3,[.08 1.2],'rmax');sl.zmax=ctrl('Core half-height',4,[.08 1.5],'zmax');
sl.seedR=ctrl('Seed radius max',5,[.01 .6],'seedR');sl.seedZ=ctrl('Seed height span',6,[.01 1.0],'seedZ');
sl.count=ctrl('Line count',7,[20 240],'count');sl.arc=ctrl('Arc length cap',8,[.15 4],'arc');sl.slice=ctrl('Slice offset y',9,[-1 1],'slice');
b=uibutton(cg,'Text','Time-series','ButtonPushedFcn',@series);place(b,10,1);
b=uibutton(cg,'Text','Export PNG','ButtonPushedFcn',@png);place(b,10,2);
h=uilabel(cg,'Text','Seeds are clipped to the sampled box.','WordWrap','on');h.Layout.Row=11;h.Layout.Column=[1 2];
ax=uiaxes(g);ax.Layout.Row=2;ax.Layout.Column=2;view(ax,-35,22);colorbar(ax);
az=uiaxes(g);az.Layout.Row=2;az.Layout.Column=3;colorbar(az);
tab=uitable(g,'ColumnName',{'Metric','Value'},'ColumnEditable',[false false]);tab.Layout.Row=3;tab.Layout.Column=[1 3];
f=uigridlayout(g,[2 3]);f.Layout.Row=4;f.Layout.Column=[1 3];f.RowHeight={24,30};f.ColumnWidth={80,'1x',200};
bp=uibutton(f,'Text','Play','ButtonPushedFcn',@play);place(bp,1,1);
status=uilabel(f,'Text','Loading...');place(status,1,2);timeLabel=uilabel(f,'Text','');place(timeLabel,1,3);
ts=uislider(f,'Limits',[m.tmin m.tmax],'Value',s.t,'MajorTicks',linspace(m.tmin,m.tmax,5),'ValueChangingFcn',@(a,e)setValue('t',e.Value,true),'ValueChangedFcn',@(a,e)setValue('t',e.Value,false));ts.Layout.Row=2;ts.Layout.Column=[1 3];
V=[];key=[];busy=false;closed=false;last=tic;
tm=timer('ExecutionMode','fixedSpacing','Period',.20,'BusyMode','drop','TimerFcn',@tick);fig.CloseRequestFcn=@closeApp;
app=struct('Figure',fig,'TimeSlider',ts,'SetTime',@(t)setValue('t',t,false),'SetControl',@(n,v)setValue(n,v,false),'Close',@closeApp,'Snapshot',@()s,'Refresh',@()render(false));render(false);
 function place(a,r,c),a.Layout.Row=r;a.Layout.Column=c;end
 function lab(t,r),h=uilabel(cg,'Text',t);place(h,r,1);end
 function h=ctrl(t,r,limits,n)
  lab(t,r);h=uislider(cg,'Limits',limits,'Value',s.(n),'MajorTicks',[],'ValueChangingFcn',@(a,e)setValue(n,e.Value,true),'ValueChangedFcn',@(a,e)setValue(n,e.Value,false));place(h,r,2);
 end
 function setValue(n,v,preview)
  if strcmp(n,'t'),validateattributes(v,{'numeric'},{'scalar','>=',m.tmin,'<=',m.tmax});ts.Value=v;
  elseif ~isfield(sl,n),error('ns:Control','Unknown core control.');
  else,validateattributes(v,{'numeric'},{'scalar','>=',sl.(n).Limits(1),'<=',sl.(n).Limits(2)});sl.(n).Value=v;end
  if strcmp(n,'count'),v=round(v);end;s.(n)=v;
  if preview&&toc(last)<.18,return,end;render(preview);
 end
 function switchModel(varargin),m=models{find(strcmp(ids,dd.Value),1)};key=[];render(false);end
 function refresh(varargin),render(false);end
 function play(varargin)
  if strcmp(tm.Running,'on'),stop(tm);bp.Text='Play';render(false);else,if s.t>=m.tmax,s.t=m.tmin;end;bp.Text='Pause';start(tm);end
 end
 function tick(varargin),s.t=min(m.tmax,s.t+.01);ts.Value=s.t;render(true);if s.t>=m.tmax,stop(tm);bp.Text='Play';render(false);end;end
 function series(varargin),stop(tm);bp.Text='Play';ns_core_timeseries(dataFile,m.id,'RadiusMax',s.rmax,'ZHalfSpan',s.zmax);end
 function png(varargin)
  stop(tm);bp.Text='Play';render(false);[n,pth]=uiputfile('*.png','Export actual core view');if ~isequal(n,0),exportapp(fig,fullfile(pth,n));end
 end
 function render(preview)
  if busy||closed,return,end;busy=true;guard=onCleanup(@unlock);start=tic;
  try
   if preview,n=25;nz=33;else,n=65;nz=97;end
   next=[s.t s.rmax s.zmax n nz];if isempty(V)||~isequal(key,next),V=ns_core_volume(m,s.t,n,nz,s.rmax,s.zmax);key=next;end
   [aa,ee]=view(ax);cla(ax);hold(ax,'on');count=s.count;if preview,count=min(36,count);end
   angles=(0:ceil(count/20)-1)*2*pi/ceil(count/20);
   [rr,zz,th]=ndgrid(min(s.seedR,.95*s.rmax)*[.25 .5 .75 1],linspace(-min(s.seedZ,.95*s.zmax),min(s.seedZ,.95*s.zmax),5),angles);
   ix=round(linspace(1,numel(rr),count));x=rr(ix).*cos(th(ix));y=rr(ix).*sin(th(ix));z=zz(ix);
   maxv=min(5000,ceil(s.arc/(.2*min([diff(V.x),diff(V.z)])))+12);
   F=stream3(V.X,V.Y,V.Z,V.U,V.V,V.W,x(:),y(:),z(:),[.2,maxv]);B=stream3(V.X,V.Y,V.Z,-V.U,-V.V,-V.W,x(:),y(:),z(:),[.2,maxv]);
   I=griddedInterpolant({V.y,V.x,V.z},V.speed,'linear','none');pts=[];col=[];
   for k=1:numel(F)
    q=[flipud(trim(B{k},s.arc/2));trim(F{k},s.arc/2)];if size(q,1)<2,continue,end
    pts=[pts;q;NaN(1,3)];col=[col;I(q(:,2),q(:,1),q(:,3));NaN]; %#ok<AGROW>
   end
   if ~isempty(pts),surface(ax,[pts(:,1)';pts(:,1)'],[pts(:,2)';pts(:,2)'],[pts(:,3)';pts(:,3)'],[col';col'],'FaceColor','none','EdgeColor','interp','LineWidth',.65);end
   axis(ax,'equal');grid(ax,'on');xlim(ax,[-s.rmax s.rmax]);ylim(ax,[-s.rmax s.rmax]);zlim(ax,[-s.zmax s.zmax]);view(ax,aa,ee);caxis(ax,[0 .3]);xlabel(ax,'x');ylabel(ax,'y');zlabel(ax,'z');title(ax,[m.id ' | instantaneous streamlines'],'Interpreter','none');hold(ax,'off');
   if preview,nn=51;else,nn=121;end
   xv=linspace(-s.rmax,s.rmax,nn);zv=linspace(-s.zmax,s.zmax,nn);[X,Z]=meshgrid(xv,zv);Y=s.slice+zeros(size(X));
   [u,~,om,res,pf]=ns_evaluate(m,[X(:),Y(:),Z(:)],s.t);
   switch sc.Value,case 'Speed',c=sqrt(sum(u.^2,2));lim=[0 .3];case 'Vorticity',c=sqrt(sum(om.^2,2));lim=[0 2];case 'Axial pressure force',c=pf(:,3);lim=[-.1 .1];otherwise,c=sqrt(sum(res.^2,2));lim=[0 .06];end
   cla(az);imagesc(az,xv,zv,reshape(c,size(X)));az.YDir='normal';axis(az,'equal');xlim(az,[-s.rmax s.rmax]);ylim(az,[-s.zmax s.zmax]);caxis(az,lim);xlabel(az,'x');ylabel(az,'z');title(az,sprintf('%s | y=%.3g',sc.Value,s.slice));
   if ~preview
    M=ns_core_metrics(m,s.t,'RadiusMax',s.rmax,'ZHalfSpan',s.zmax);
    names={'weighted_r_rms','weighted_z_rms','axial_to_radial_extent','inward_fraction','positive_swirl_fraction','bipolar_fraction','axial_pressure_toward_fraction','median_swirl_over_inflow','sampled_residual_max'};
    rows=cell(numel(names),2);for j=1:numel(names),rows{j,1}=strrep(names{j},'_',' ');rows{j,2}=M.(names{j});end;tab.Data=rows;
   end
   timeLabel.Text=sprintf('Physical t = %.5f',s.t);status.Text=sprintf('%dx%dx%d local nodes | %.2f s | diagnostics only; PDE gates not met',n,n,nz,toc(start));setappdata(fig,'LastRenderError','');drawnow limitrate nocallbacks;last=tic;
  catch err,setappdata(fig,'LastRenderError',err.message);status.Text=['Render error: ' err.message];warning('ns:Core','%s',err.message);end
  clear guard
 end
 function q=trim(q,L),if size(q,1)<2,return,end;lastq=find([0;cumsum(sqrt(sum(diff(q).^2,2)))]<=L,1,'last');q=q(1:lastq,:);end
 function unlock(),busy=false;end
 function closeApp(varargin),if closed,return,end;closed=true;if isvalid(tm),stop(tm);delete(tm);end;if isvalid(fig),delete(fig);end;end
end
