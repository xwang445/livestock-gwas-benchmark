from pathlib import Path
import csv,hashlib,sys,tempfile,unittest
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import workflow
class Workflow(unittest.TestCase):
 def test_protects_source_and_repository(self):
  with tempfile.TemporaryDirectory() as td:
   source=Path(td)/'historical';source.mkdir()
   for target in [source/'new',ROOT/'figures'/'new',ROOT/'derived_data'/'new']:
    with self.assertRaises(ValueError):workflow.new_output(target,source)
    self.assertFalse(target.exists())
 def test_existing_outputs_require_explicit_resume(self):
  with tempfile.TemporaryDirectory() as td:
   target=Path(td)/'run';target.mkdir();(target/'sentinel').write_text('preserve')
   with self.assertRaises(ValueError):workflow.new_output(target)
   self.assertEqual(workflow.new_output(target,resume=True),target.resolve())
   self.assertEqual((target/'sentinel').read_text(),'preserve')
 def test_archived_sources_are_intact_and_inactive(self):
  for row in csv.DictReader((ROOT/'audit/legacy_relocation.tsv').read_text().splitlines(),delimiter='\t'):
   self.assertFalse((ROOT/row['old_path']).exists())
   self.assertEqual(hashlib.sha256((ROOT/row['archived_path']).read_bytes()).hexdigest(),row['sha256'])
 def test_summary_dispatch_uses_verified_profile(self):
  with tempfile.TemporaryDirectory() as td:
   source=Path(td)/'source';source.mkdir();out=Path(td)/'out'
   with patch.object(sys,'argv',['workflow.py','summarize','--profile','mac','--out',str(out)]),patch.object(workflow,'resolve_profile',return_value={'benchmark_root':str(source),'inputs':{}}) as resolver,patch.object(workflow,'run') as run:
    workflow.main()
   resolver.assert_called_once_with(ROOT,'mac',None,check_inputs=True)
   call=run.call_args.args
   self.assertEqual(call[:2],('scripts/benchmark_core.py','summarize'))
   self.assertIn(str(source),call)
   self.assertFalse((source/'derived_data').exists())
