import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python" / "reseq" / "plot_data_stats.py"


class PlotDataStatsTest(unittest.TestCase):
    def test_help_exits_successfully(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            Path(tmp_dir, "DataStats.py").write_text("# test stub\n", encoding="utf-8")
            env = dict(os.environ)
            env["RESEQ_PYMODS"] = tmp_dir
            env["MPLCONFIGDIR"] = tmp_dir

            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--help"],
                capture_output=True,
                text=True,
                check=False,
                env=env,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("usage:", result.stdout)


if __name__ == "__main__":
    unittest.main()
