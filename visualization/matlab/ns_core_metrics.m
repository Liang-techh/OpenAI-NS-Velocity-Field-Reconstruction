function M=ns_core_metrics(m,t,varargin)
%NS_CORE_METRICS Cylinder-restricted, trapezoid-weighted local diagnostics.
% Local extents depend on the chosen observation cylinder, not a detected
% core boundary or a similarity score against an unavailable source field.
p=inputParser;addParameter(p,'RadiusMax',.35);addParameter(p,'ZHalfSpan',.60);
addParameter(p,'Nxy',49);addParameter(p,'Nz',65);addParameter(p,'Volume',[]);
addParameter(p,'BipolarExclusion',.03);parse(p,varargin{:});o=p.Results;
V=o.Volume;if isempty(V),V=ns_core_volume(m,t,o.Nxy,o.Nz,o.RadiusMax,o.ZHalfSpan);end
if abs(V.time-t)>1e-12,error('ns:MetricTime','Volume time differs from requested time.');end
R=hypot(V.X,V.Y);Z=V.Z;inside=R<=o.RadiusMax & abs(Z)<=o.ZHalfSpan;
Ur=(V.X.*V.U+V.Y.*V.V)./max(R,eps);Ut=(-V.Y.*V.U+V.X.*V.V)./max(R,eps);
wx=weights(V.x);wy=weights(V.y);wz=weights(V.z);
W=reshape(wy,[],1,1).*reshape(wx,1,[],1).*reshape(wz,1,1,[]);
valid=inside & isfinite(V.vorticity) & isfinite(V.speed);
w=W.*V.vorticity.^2;w(~valid)=0;mass=sum(w(:));
if mass<=0,error('ns:NoEnstrophy','No finite nonzero vorticity in observation cylinder.');end
M=struct('time',t,'radius_max',o.RadiusMax,'z_half_span',o.ZHalfSpan,'weight','vorticity squared times Cartesian trapezoid volume weights',...
    'sample_count',nnz(valid),'bipolar_exclusion',o.BipolarExclusion,'pde_validated',false,'source_correspondence_verified',false);
M.weighted_r_rms=sqrt(sum(w(:).*R(:).^2)/mass);M.weighted_z_rms=sqrt(sum(w(:).*Z(:).^2)/mass);
M.axial_to_radial_extent=M.weighted_z_rms/M.weighted_r_rms;
q=valid&R>1e-8;M.inward_fraction=mean(Ur(q)<0);M.positive_swirl_fraction=mean(Ut(q)>0);
q=valid&abs(Z)>o.BipolarExclusion;M.bipolar_outflow_fraction=mean(Z(q).*V.W(q)>0);
M.axial_pressure_toward_fraction=mean(Z(q).*V.PforceZ(q)<0);
q=valid&R>1e-8&Ur<0&abs(Ur)>1e-5;
M.ratio_sample_count=nnz(q);M.mean_swirl_over_inflow=mean(abs(Ut(q)./Ur(q)));
M.median_swirl_over_inflow=median(abs(Ut(q)./Ur(q)));M.mean_axial_over_inflow=mean(abs(V.W(q)./Ur(q)));
M.median_axial_over_inflow=median(abs(V.W(q)./Ur(q)));
M.speed_max=max(V.speed(valid));M.residual_max=max(V.residual(valid));
M.scope='Finite local-cylinder geometry; no target correspondence, continuous bound or PDE acceptance. Ratios exclude near-zero radial inflow.';
end
function w=weights(x)
x=x(:);d=diff(x);w=[d(1)/2;(d(1:end-1)+d(2:end))/2;d(end)/2];
end
