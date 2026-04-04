import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "python" / "reseq" / "prepare_names.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def write_fastq(path: Path, header: str, sequence: str = "ACGT", quality: str = "!!!!") -> None:
    path.write_text(f"{header}\n{sequence}\n+\n{quality}\n", encoding="utf-8")


def write_fastq_gz(path: Path, header: str, sequence: str = "TGCA", quality: str = "####") -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"{header}\n{sequence}\n+\n{quality}\n")


class ReSeqPrepareNamesTest(unittest.TestCase):
    def test_rewrites_first_file_headers_for_plain_inputs(self):
        tmp_path = Path(tempfile.mkdtemp(prefix="reseq-prepare-names-"))
        self.addCleanup(shutil.rmtree, tmp_path)

        file1 = tmp_path / "reads_1.fastq"
        file2 = tmp_path / "reads_2.fastq"
        write_fastq(file1, "@INST:1:FCID:2:2104:15343:197393 1:N:0:ATCACG")
        write_fastq(file2, "@INST:1:FCID:2:2104:15343:197393 2:N:0:ATCACG")

        result = run_script(str(file1), str(file2))

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines()[0], "@INST:1:FCID:2:2104:15343:197393")

    def test_supports_gzip_inputs(self):
        tmp_path = Path(tempfile.mkdtemp(prefix="reseq-prepare-names-"))
        self.addCleanup(shutil.rmtree, tmp_path)

        file1 = tmp_path / "reads_1.fastq.gz"
        file2 = tmp_path / "reads_2.fastq.gz"
        write_fastq_gz(file1, "@INST:1:FCID:2:2104:15343:197393 1:N:0:ATCACG")
        write_fastq_gz(file2, "@INST:1:FCID:2:2104:15343:197393 2:N:0:ATCACG")

        result = run_script(str(file1), str(file2))

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines()[0], "@INST:1:FCID:2:2104:15343:197393")

    def test_requires_exactly_two_inputs(self):
        result = run_script()

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("usage", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
