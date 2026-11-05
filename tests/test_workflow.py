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
   with patch.object(sys,'argv',['workflow.py','summarize','--profile','mac','--simulation',str(source/'sim'),'--genotype-prefix',str(source/'panel'),'--out',str(out)]),patch.object(workflow,'resolve_profile',return_value={'benchmark_root':str(source),'inputs':{}}) as resolver,patch.object(workflow,'run') as run:
    workflow.main()
   resolver.assert_called_once_with(ROOT,'mac',None,check_inputs=True,stage='summarize',historical=False)
   call=run.call_args_list[0].args
   self.assertEqual(call[0],'scripts/revision_2026-11-05/association/derive_association_tables.py')
   self.assertIn(str(source),call)
   self.assertFalse((source/'derived_data').exists())
 def test_historical_strict_route_requires_explicit_flag(self):
  with tempfile.TemporaryDirectory() as td:
   source=Path(td)/'source';source.mkdir();out=Path(td)/'out'
   with patch.object(sys,'argv',['workflow.py','summarize','--historical','--profile','hpcc','--out',str(out)]),patch.object(workflow,'resolve_profile',return_value={'benchmark_root':str(source),'inputs':{}}),patch.object(workflow,'run') as run:
    workflow.main()
   self.assertEqual(run.call_args.args[:3],('scripts/benchmark_core.py','summarize','--historical'))
 def test_current_manhattan_route_uses_new_inputs_and_does_not_precreate_output(self):
  with tempfile.TemporaryDirectory() as td:
   source=Path(td)/'source';source.mkdir();out=Path(td)/'out'
   with patch.object(sys,'argv',['workflow.py','manhattan','--profile','hpcc','--out',str(out)]),patch.object(workflow,'resolve_profile',return_value={'benchmark_root':str(source),'inputs':{}}),patch.object(workflow,'run') as run:
    workflow.main()
   self.assertIn('--require-complete',run.call_args.args)
   self.assertIn('--gcta-dir',run.call_args.args)
   self.assertFalse(out.exists())
