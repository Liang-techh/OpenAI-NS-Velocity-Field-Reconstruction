"""Build and exercise an installed wheel away from the source checkout."""
from pathlib import Path
import subprocess, sys, tempfile
ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as d:
    tmp=Path(d);wheels=tmp/'wheel';site=tmp/'site'
    subprocess.run([sys.executable,'-m','pip','wheel','--no-deps','--no-build-isolation',str(ROOT),'-w',str(wheels)],check=True)
    wheel=next(wheels.glob('*.whl'))
    subprocess.run([sys.executable,'-m','pip','install','--no-deps','--target',str(site),str(wheel)],check=True)
    code="import sys;sys.path.insert(0,sys.argv[1]);from ns_reconstruction import load_candidate,verify_integrity;" \
         "m=verify_integrity();f=load_candidate();u,p=f.fields([[.1,0,.1]],.5);assert u.shape==(1,3) and p.shape==(1,);" \
         "assert not m['pde_validated'];print('Installed-wheel API and frozen data: PASS')"
    subprocess.run([sys.executable,'-c',code,str(site)],cwd=tmp,check=True)
