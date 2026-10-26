"""Generation must fail when an already admitted operator is refused."""
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("emit", Path(__file__).with_name("emit.py"))
emit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(emit)


class ExpectedProducts(unittest.TestCase):
    def test_refused_ported_pair_fails_without_shrinking_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = Path.cwd()
            try:
                os.chdir(tmp)
                root = Path("tools/fortad")
                root.mkdir(parents=True)
                entry = {"stem": "triangle", "active": ["x"], "outputs": ["y"]}
                (root / "catalogue.json").write_text(json.dumps([entry]))
                manifest = root / "ported.json"
                manifest.write_text('["triangle"]\n')
                with patch.object(emit, "run", return_value=False), \
                     patch.object(emit, "write_test") as write_test, \
                     contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(emit.main(), 1)
                    write_test.assert_not_called()
                self.assertEqual(json.loads(manifest.read_text()), ["triangle"])
            finally:
                os.chdir(old_cwd)


if __name__ == "__main__":
    unittest.main()
