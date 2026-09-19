function V=ns_core_volume(m,t,nxy,nz,rmax,zmax)
%NS_CORE_VOLUME Re-evaluate a physical local box, not a zoomed coarse volume.
% Tensor evaluation uses the unchanged finite basis, with no time interpolation.
if ~strcmp(m.mode,'spectral'),error('ns:Core','Core diagnostics require a spectral model.');end
validateattributes(t,{'numeric'},{'scalar','real','finite','>=',m.tmin,'<=',m.tmax});
validateattributes(nxy,{'numeric'},{'scalar','integer','>=',9,'<=',129});
validateattributes(nz,{'numeric'},{'scalar','integer','>=',9,'<=',193});
validateattributes(rmax,{'numeric'},{'scalar','real','finite','>',0,'<=',2});
validateattributes(zmax,{'numeric'},{'scalar','real','finite','>',0,'<=',2});
x=linspace(-rmax,rmax,nxy);z=linspace(-zmax,zmax,nz);
[xx,yy]=meshgrid(x,x);s=xx(:).^2+yy(:).^2;
D=ns_cylindrical(m,s,z,t,true);xx=xx(:);yy=yy(:);shape=[nxy nxy nz];
V.x=x;V.y=x;V.z=z;[V.X,V.Y,V.Z]=meshgrid(x,x,z);
V.U=reshape(xx.*D.A-yy.*D.B,shape);V.V=reshape(yy.*D.A+xx.*D.B,shape);V.W=reshape(D.C,shape);
V.P=reshape(D.P,shape);V.PforceZ=reshape(-D.Pz,shape);
V.speed=sqrt(V.U.^2+V.V.^2+V.W.^2);V.swirl=reshape(sqrt(s).*D.B,shape);
V.vorticity=reshape(sqrt(s.*D.Bz.^2+s.*(D.Az-2*D.Cs).^2+4*(D.B+s.*D.Bs).^2),shape);
V.residual=reshape(sqrt(s.*(D.Ra.^2+D.Rb.^2)+D.Rc.^2),shape);
V.divergence=reshape(D.divergence,shape);V.time=t;V.approximate=false;
end
