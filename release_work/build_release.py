"""Build a clean publication; never edit the source research repository."""
from pathlib import Path
import argparse, hashlib, json, re, shutil, sys
import numpy as np
from scipy.io import loadmat
SCIENCE='c77492a48e9c0f13d4d51244987c57c28519409b'
VISUAL='5d49ccf07d097a4bffeb35f6398936493b062d00'
MAT_SHA='59d43fdc40cd4df1675c2191f85810b9351802b2f4b1e0ac25a3bccfd9914c62'
ORIGINAL={'ST054-Q2':'d772949621e0f7703a9d6b36f28828532ba3e5ec678dec14db5ce14eb8b851d5','ST054-M3':'2b2b571986297b52bc884a2ee4808a7ade1ade86f38f938604f2a64118882df3'}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(root,name,obj):
 p=root/name;p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(obj,indent=2)+'\n' if not isinstance(obj,str) else obj,encoding='utf-8')

def build(science,visual,out,staging):
 if out.exists() and any(out.iterdir()):raise ValueError('Refuse nonempty destination')
 out.mkdir(parents=True,exist_ok=True);copied={}
 def copy(src,name):
  dst=out/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst);copied[name]=sha(src)
 for p in (staging/'payload').rglob('*'):
  if p.is_file():copy(p,p.relative_to(staging/'payload').as_posix())
 for name in ('spacetime.py','hybrid_basis.py','asymmetric_basis.py','validate.py'):
  copy(science/'experiments/root_st030'/name,'runtime/'+name)
 copy(science/'experiments/root_st054/results.json','evidence/source_results.json')
 for p in visual.glob('*.m'):copy(p,'visualization/matlab/'+p.name)
 for name in ('st054_models.mat','st054_models.json'):copy(visual/'data'/name,'visualization/matlab/data/'+name)
 assert sha(out/'visualization/matlab/data/st054_models.mat')==MAT_SHA
 for name in ('matlab_test_receipt.json','callback_receipt.json','matlab_explorer.png'):
  p=visual/'tests/output'/name
  if p.exists():copy(p,'evidence/source_matlab/'+name)
 # New core code in the publication payload takes precedence over imported names.
 for p in (staging/'payload/visualization/matlab').glob('*.m'):copy(p,'visualization/matlab/'+p.name)
 p=out/'visualization/matlab/ns_explorer.m';text=p.read_text()
 before="cg=uigridlayout(panel,[18 2]);cg.ColumnWidth={124,'1x'};cg.RowHeight=repmat({32},1,18);"
 after="cg=uigridlayout(panel,[19 2]);cg.ColumnWidth={124,'1x'};cg.RowHeight=repmat({32},1,19);"
 assert before in text;text=text.replace(before,after)
 anchor="mat=uibutton(cg,'Text','Export view MAT','ButtonPushedFcn',@exportMat);place(mat,18,2);"
 assert anchor in text
 text=text.replace(anchor,anchor+"\ncore=uibutton(cg,'Text','Core explorer','ButtonPushedFcn',@(~,~)ns_core_explorer(dataFile));place(core,19,1);\nseries=uibutton(cg,'Text','Core time series','ButtonPushedFcn',@(~,~)ns_core_timeseries(dataFile,m.id));place(series,19,2);")
 p.write_text(text)
 sys.path.insert(0,str(out/'runtime'));from spacetime import Family
 from validate import validate
 data=loadmat(out/'visualization/matlab/data/st054_models.mat',simplify_cells=True)
 source=json.loads((out/'evidence/source_results.json').read_text())
 registry={'default_candidate':'ST054-Q2','selection_scope':'Lower L2 among the two complete ST054 candidates; M3 is the lower-peak alternative. Not a global-best claim.','pde_validated':False,'candidates':{}}
 parity={}
 for m in data['models']:
  ident=m['id'];assert m['original_raw_sha256']==ORIGINAL[ident]
  raw=np.asarray(m['raw_coefficients'],float);assert raw.shape==(2594,)
  obj={'schema':'root_st001_compact_spacetime_v1','nr':9,'nz':12,'nt':8,'basis_kind':'hybrid9_axial_opposite3_v1','coefficients':raw.tolist(),'metadata':{'id':ident,'science_commit':SCIENCE,'source_mat_sha256':MAT_SHA,'original_research_json_sha256':ORIGINAL[ident],'scope':'Exact raw vector from native-tested MAT; new JSON metadata and file hash, not original research JSON bytes.'},'pde_validated':False,'paper_exact':False,'blowup_proved':False,'field_identity_claim':False}
  name=f'candidates/{ident}/candidate.json';write(out,name,obj);f,r=Family.load(out/name);error=np.zeros(3)
  for i,t in enumerate(m['ref_times']):
   u,p=f.fields(r,m['ref_points'],t);res=f.analytic_residual(r,m['ref_points'],t)
   error=np.maximum(error,[np.max(abs(u-m['ref_velocity'][i])),np.max(abs(p-m['ref_pressure'][i])),np.max(abs(res-m['ref_residual'][i]))])
  assert error.max()<1e-8;parity[ident]=error.tolist()
  registry['candidates'][ident]={'file':name,'sha256':sha(out/name),'original_research_json_sha256':ORIGINAL[ident],'coefficient_array_sha256':hashlib.sha256(raw.astype('<f8').tobytes()).hexdigest(),'pde_validated':False,'paper_exact':False,'source_correspondence_verified':False,'blowup_proved':False}
  for seed in (9175491,9175492):
   dst=out/f'evidence/replay/{ident}_{seed}.json';dst.parent.mkdir(parents=True,exist_ok=True)
   report=validate(out/name,dst,seed)
   ref=next(x for x in source['paired_results'] if x['id']==ident and x['seed']==seed)
   fine=[x for x in report['spatial_refinement'] if x['space_step']==.005]
   assert abs(max(x['momentum_max'] for x in fine)-ref['max'])<1e-7
   assert abs(max(x['momentum_L2'] for x in fine)-ref['volume_L2'])<1e-7
   assert report['gates']==source['common_gates'];assert not report['all_numeric_gates_pass']
 write(out,'candidates/registry.json',registry)
 write(out,'evidence/release_reference_parity.json',{'max_absolute_u_p_R_errors':parity,'tolerance':1e-8,'pde_validated':False})
 write(out,'provenance.json',{'source_repository':'Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1','science_commit':SCIENCE,'visualization_commit':VISUAL,'previous_target_main':'065b67e9a3e8e22d697b49ccc67163db60a79783','backup_branch':'archive/pre-st054-replacement-2026-09-19','copied_input_sha256':copied,'source_mat_sha256':MAT_SHA,'scope':'User-requested publication replacement. No refit or change to the source repository. Candidate JSON is newly serialized; original and release hashes differ. Two launcher buttons added; focused core tools are new.'})
 write(out,'repository-metadata.json',{'description':'Independent Navier-Stokes research with reproducible ST054 candidates, full-residual validation, and interactive MATLAB visualization.','github_about_updated':False,'reason':'The available connector does not expose a repository-settings write action.'})
 for p in out.rglob('*'):
  if p.is_file() and p.suffix in ('.md','.py','.m','.json','.yml','.txt'):
   assert not re.search('[\u3400-\u9fff]',p.read_text(encoding='utf-8')),str(p)
 for p in out.rglob('__pycache__'):shutil.rmtree(p)
 inventory={p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file()}
 write(out,'manifest.json',{'schema':'publication_sha256_v1','files':dict(sorted(inventory.items())),'pde_validated':False})
 print(json.dumps({'published_files':len(inventory)+1,'full_replays':4,'parity':parity},indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--science',type=Path,required=True);p.add_argument('--visual',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--staging',type=Path,default=Path(__file__).parent);a=p.parse_args();build(a.science,a.visual,a.out,a.staging)
