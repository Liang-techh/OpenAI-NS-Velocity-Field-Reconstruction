function V=ns_core_volume(m,t,nxy,nz,rmax,zmax)
%NS_CORE_VOLUME Re-sample a local physical box; never zoom a coarse full grid.
% The box is square in x/y. Cylindrical masks are applied by ns_core_metrics.
% Only the frozen spectral models are accepted; no missing p/R is fabricated.
if ~strcmp(m.mode,'spectral'),error('ns:CoreMode','Use ns_explorer for imported sampled grids.');end
validateattributes(nxy,{'numeric'},{'scalar','integer','>=',9,'<=',129});
validateattributes(nz,{'numeric'},{'scalar','integer','>=',9,'<=',193});
validateattributes(rmax,{'numeric'},{'scalar','finite','positive','<=',2});
validateattributes(zmax,{'numeric'},{'scalar','finite','positive','<=',2});
x=linspace(-rmax,rmax,nxy);z=linspace(-zmax,zmax,nz);
[X,Y]=meshgrid(x,x);s=X(:).^2+Y(:).^2;xx=X(:);yy=Y(:);
D=ns_cylindrical(m,s,z,t,true);sz=[nxy nxy nz];
V.x=x;V.y=x;V.z=z;[V.X,V.Y,V.Z]=meshgrid(x,x,z);
V.U=reshape(xx.*D.A-yy.*D.B,sz);V.V=reshape(yy.*D.A+xx.*D.B,sz);V.W=reshape(D.C,sz);
V.P=reshape(D.P,sz);V.PforceZ=reshape(-D.Pz,sz);V.PforceR=reshape(-2*sqrt(s).*D.Ps,sz);
V.speed=sqrt(V.U.^2+V.V.^2+V.W.^2);V.radius=hypot(V.X,V.Y);
V.swirl=reshape(sqrt(s).*D.B,sz);V.axial=V.W;
V.vorticity=reshape(sqrt(s.*D.Bz.^2+s.*(D.Az-2*D.Cs).^2+4*(D.B+s.*D.Bs).^2),sz);
V.residual=reshape(sqrt(s.*(D.Ra.^2+D.Rb.^2)+D.Rc.^2),sz);
V.divergence=reshape(D.divergence,sz);V.time=t;V.approximate=false;
V.domain='Local physical box; no change to frozen coefficients';
end
