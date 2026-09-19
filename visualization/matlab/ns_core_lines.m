function [lines,stats]=ns_core_lines(V,seedR,seedZ,count,cap)
%NS_CORE_LINES Instantaneous streamlines, not material trajectories.
% Positive/negative integration constructs the displayed lines. Metrics use
% forward portions only. Finite arc/domain stops are reported, not hidden.
validateattributes(count,{'numeric'},{'scalar','integer','>=',1,'<=',400});
validateattributes(cap,{'numeric'},{'scalar','finite','positive'});
if seedR<=0||seedR>=V.x(end)||seedZ<0||seedZ>=V.z(end),error('ns:Seeds','Seeds must be inside the local box.');end
angles=(0:max(7,ceil(count/20)-1))*2*pi/max(8,ceil(count/20));
[r,z,a]=ndgrid(seedR*[.25 .5 .75 1],linspace(-seedZ,seedZ,5),angles);
take=unique(round(linspace(1,numel(r),count)));
sx=r(take).*cos(a(take));sy=r(take).*sin(a(take));sz=z(take);
ds=min([diff(V.x),diff(V.y),diff(V.z)]);step=.15;
maxv=min(5000,ceil(cap/(step*ds))+20);
F=stream3(V.X,V.Y,V.Z,V.U,V.V,V.W,sx(:),sy(:),sz(:),[step maxv]);
B=stream3(V.X,V.Y,V.Z,-V.U,-V.V,-V.W,sx(:),sy(:),sz(:),[step maxv]);
lines=cell(numel(F),1);turn=[];span=[];contraction=[];lengths=[];
for k=1:numel(F)
    fw=trimArc(F{k},cap/2);bw=trimArc(B{k},cap/2);
    if isempty(fw),lines{k}=bw;continue,end
    if isempty(bw),lines{k}=fw;else,lines{k}=[flipud(bw);fw(2:end,:)];end
    if size(fw,1)<3,continue,end
    rad=hypot(fw(:,1),fw(:,2));th=unwrap(atan2(fw(:,2),fw(:,1)));
    if min(rad)<1e-7,continue,end
    turn(end+1)=th(end)-th(1);span(end+1)=max(fw(:,3))-min(fw(:,3)); %#ok<AGROW>
    contraction(end+1)=rad(end)/rad(1);lengths(end+1)=sum(vecnorm(diff(fw),2,2)); %#ok<AGROW>
end
stats=struct('valid_count',numel(turn),'seed_count',numel(F),'median_turn_angle',median(turn),...
 'mean_turn_angle',mean(turn),'mean_axial_span',mean(span),'mean_end_radius_fraction',mean(contraction),...
 'mean_forward_arc',mean(lengths),'scope','Forward instantaneous streamline segments only; finite seed, arc and domain dependent; not particle trajectories');
end
function q=trimArc(q,cap)
if isempty(q)||size(q,1)<2,return,end
d=[0;cumsum(vecnorm(diff(q),2,2))];last=find(d<=cap,1,'last');
if last<size(q,1)&&d(last+1)>d(last)
    a=(cap-d(last))/(d(last+1)-d(last));q=[q(1:last,:);q(last,:)+a*(q(last+1,:)-q(last,:))];
else,q=q(1:last,:);end
end
