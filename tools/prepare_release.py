"""One-time pinned setup for a source checkout; initialized bundles stay offline."""
from pathlib import Path
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
marker=ROOT/'src/ns_reconstruction/data/manifest.json'
if marker.is_file():
    sys.path.insert(0,str(ROOT/'src'))
    from ns_reconstruction import verify_integrity
    verify_integrity()
    records=json.loads((ROOT/'evidence/import_receipt.json').read_text())['imports']
    for row in records:
        path=ROOT/row['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row['release_sha256']:
            raise ValueError(f"Existing imported asset changed or missing: {row['path']}")
    print('Initialized bundle verified; no download needed. Scientific target still unmet.')
else:
    subprocess.run([sys.executable,str(ROOT/'tools/build_release.py')],check=True)
