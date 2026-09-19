function M=ns_core_metrics(m,t,varargin)
%NS_CORE_METRICS Finite-cylinder diagnostics, not fitted target-image scores.
% Moments use actual enstrophy and cylindrical volume weights 2*pi*r*dr*dz.
% Window choices change these diagnostics; they are not intrinsic core boundaries.
p=inputParser;
addParameter(p,'RadiusMax',.35);addParameter(p,'ZHalfSpan',.60);
addParameter(p,'Nr',65);addParameter(p,'Nz',97);addParameter(p,'ExcludedMidplane',.03);
parse(p,varargin{:});o=p.Results;
validateattributes(o.RadiusMax,{'numeric'},{'scalar','positive','<=',2});
validateattributes(o.ZHalfSpan,{'numeric'},{'scalar','positive','<=',2});
validateattributes(o.Nr,{'numeric'},{'scalar','integer','>=',9});
validateattributes(o.Nz,{'numeric'},{'scalar','integer','>=',9});
r=linspace(0,o.RadiusMax,o.Nr);z=linspace(-o.ZHalfSpan,o.ZHalfSpan,o.Nz);
[R,Z]=meshgrid(r,z);pts=[R(:),zeros(numel(R),1),Z(:)];
[u,~,om,res,pf]=ns_evaluate(m,pts,t);en=sum(om.^2,2);
wr=ones(1,o.Nr);wr([1 end])=.5;wz=ones(o.Nz,1);wz([1 end])=.5;
w=2*pi*R.*(wz*wr)*(r(2)-r(1))*(z(2)-z(1));weight=w(:).*en;mass=sum(weight);
if ~isfinite(mass)||mass<=1e-30,error('ns:Core','No resolved nonzero enstrophy in this window.');end
M.time=t;M.radius_max=o.RadiusMax;M.z_half_span=o.ZHalfSpan;M.sample_shape=[o.Nr o.Nz];
M.weighted_r_rms=sqrt(sum(weight.*R(:).^2)/mass);
M.weighted_z_rms=sqrt(sum(weight.*Z(:).^2)/mass);
M.axial_to_radial_extent=M.weighted_z_rms/M.weighted_r_rms;
axisMask=R(:)>1e-8;away=axisMask & abs(Z(:))>o.ExcludedMidplane;
M.inward_fraction=mean(u(axisMask,1)<0);M.positive_swirl_fraction=mean(u(axisMask,2)>0);
M.bipolar_fraction=mean(sign(Z(away)).*u(away,3)>0);
M.axial_pressure_toward_fraction=mean(sign(Z(away)).*pf(away,3)<0);
ratioMask=axisMask & u(:,1)<-1e-5;
M.ratio_valid_count=sum(ratioMask);M.median_swirl_over_inflow=median(abs(u(ratioMask,2)./u(ratioMask,1)));
M.median_axial_over_inflow=median(abs(u(ratioMask,3)./u(ratioMask,1)));
M.speed_max=max(sqrt(sum(u.^2,2)));M.sampled_residual_max=max(sqrt(sum(res.^2,2)));
M.pde_validated=false;M.scope='Finite fixed-cylinder enstrophy moments and sample fractions; no target match or continuum certificate.';
end
