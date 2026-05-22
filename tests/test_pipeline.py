from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
import gzip
import os
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ontyper import demultiplex, report
from ontyper.genotypes import extract, identify, visualize


DATA_DIR = Path(__file__).parent / "data"
SUBSET_FASTA = DATA_DIR / "pipeline_subset_large.fasta"
DEMULTIPLEX_FILE = DATA_DIR / "pipeline_demultiplex_large.csv"

EXPECTED_DEMULTIPLEX = {
    "STMock03_Ssa408Uos.fas": 1030,
    "STMock12_Ssa408Uos.fas": 1022,
    "STMock03_Ssa410Uos.fas": 1063,
    "STMock08_Ssa410Uos.fas": 1001,
}

EXPECTED_GENOTYPES = {
    "STMock03_Ssa408Uos": {"13x": 265, "15x": 254, "21x": 100, "nsequences": 1030},
    "STMock12_Ssa408Uos": {"11x": 237, "20x": 113, "21x": 129, "nsequences": 1022},
    "STMock03_Ssa410Uos": {"22x": 156, "25x": 96, "30x": 129, "nsequences": 1063},
    "STMock08_Ssa410Uos": {"22x": 149, "32x": 103, "34x": 74, "nsequences": 1001},
}

EXPECTED_DETAILED = {
    ("STMock03_Ssa408Uos", "15x"): {
        "count": 254,
        "first_ids": [
            "e0589df2-9065-4922-9ee7-71bc36a97d13",
            "02ec0c61-48ff-4700-815c-b4b33e062d05",
            "4615b1b5-77a4-4694-97f6-25ed6ec6c3b0",
        ],
    },
    ("STMock12_Ssa408Uos", "21x"): {
        "count": 129,
        "first_ids": [
            "57ff1fc5-a611-4b78-8e57-423d51c3c745",
            "de120f76-63f7-4604-8817-9b3f1c4c1ae4",
            "b35edc7f-24b6-46a5-b684-44dfa47c3e5c",
        ],
    },
    ("STMock03_Ssa410Uos", "22x"): {
        "count": 156,
        "first_ids": [
            "6132ffd0-4553-4732-a33b-15020364a801",
            "0f0f67fd-4002-4c8a-a18d-d889340a336c",
            "a8cb31cd-5993-4ec2-86d5-697ca10c6c66",
        ],
    },
    ("STMock08_Ssa410Uos", "32x"): {
        "count": 103,
        "first_ids": [
            "cc5d4b88-ed58-4ff2-833d-a2ae6a6e1b15",
            "880351c4-472a-43ce-8326-757cc23d996f",
            "386ecd78-d4cf-4b95-a20b-19acf6f4034c",
        ],
    },
}


def count_fasta_records(path: Path) -> int:
    count = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count



def fasta_ids(path: Path) -> list[str]:
    ids = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids



def parse_identify_output(text: str) -> dict[str, dict[str, int]]:
    lines = [line for line in text.splitlines() if line.strip()]
    header = [part.strip() for part in lines[0].split("\t")]
    repeat_bins = header[1:-1]

    table = {}
    for line in lines[1:]:
        parts = [part.strip() for part in line.split("\t")]
        counts = {repeat_bins[i]: int(parts[i + 1]) for i in range(len(repeat_bins))}
        counts["nsequences"] = int(parts[-1])
        table[parts[0]] = counts
    return table



def parse_detailed_output(path: Path) -> dict[tuple[str, str], list[str]]:
    detailed = {}
    with open(path) as handle:
        next(handle)
        for line in handle:
            locus, repeat_bin, names = line.rstrip("\n").split("\t")
            detailed[(locus, repeat_bin)] = [name.split()[0] for name in names.split(",") if name]
    return detailed



def write_fastq_from_fasta(source_fasta: Path, target_fastq: Path) -> Path:
    opener = gzip.open if str(target_fastq).endswith(".gz") else open

    header = None
    sequence_lines = []
    with open(source_fasta) as infile, opener(target_fastq, "wt") as out:
        for line in infile:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_lines)
                    out.write(f"@{header}\n{sequence}\n+\n{'I' * len(sequence)}\n")
                header = line[1:]
                sequence_lines = []
                continue
            sequence_lines.append(line)

        if header is not None:
            sequence = "".join(sequence_lines)
            out.write(f"@{header}\n{sequence}\n+\n{'I' * len(sequence)}\n")

    return target_fastq



def run_demultiplex(outdir: Path, infile: Path = SUBSET_FASTA) -> Path:
    args = SimpleNamespace(
        infile=str(infile),
        demfile=str(DEMULTIPLEX_FILE),
        threads=1,
        mode="single",
        outdir=str(outdir),
        cutoff="",
        delimiter=",",
        window=150,
        chunksize=1000,
        primer_mismatches=5,
        tag_mismatches=0,
        single_mismatches=5,
    )

    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        assert demultiplex.main(args) == 0
    return outdir



def run_identify(indir: Path, workdir: Path) -> tuple[dict[str, dict[str, int]], Path, Path]:
    output_file = workdir / "genotypes.tsv"
    args = SimpleNamespace(
        ntimes="50",
        demfile=str(DEMULTIPLEX_FILE),
        indir=str(indir),
        method="strict",
        delimiter=",",
        relative=False,
        detailed=True,
        outfile=str(output_file),
    )

    old_cwd = Path.cwd()
    try:
        workdir.mkdir(parents=True, exist_ok=True)
        os.chdir(workdir)
        with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
            assert identify.main(args) == 0
    finally:
        os.chdir(old_cwd)

    detailed_file = workdir / "genotypes-detailed.tsv"
    assert output_file.exists()
    assert detailed_file.exists()
    return parse_identify_output(output_file.read_text()), detailed_file, output_file



def test_demultiplex_recovers_expected_sample_sizes(tmp_path):
    demux_dir = run_demultiplex(tmp_path / "demultiplexed")

    produced = {path.name: count_fasta_records(path) for path in demux_dir.glob("*.fas")}
    assert produced == EXPECTED_DEMULTIPLEX



def test_demultiplex_recovers_expected_sample_sizes_from_fastq_gz(tmp_path):
    fastq_file = write_fastq_from_fasta(SUBSET_FASTA, tmp_path / "pipeline_subset_large.fastq.gz")
    demux_dir = run_demultiplex(tmp_path / "demultiplexed", fastq_file)

    produced = {path.name: count_fasta_records(path) for path in demux_dir.glob("*.fas")}
    assert produced == EXPECTED_DEMULTIPLEX



def test_identify_recovers_reference_genotypes(tmp_path):
    demux_dir = run_demultiplex(tmp_path / "demultiplexed")
    genotype_table, detailed_file, _ = run_identify(demux_dir, tmp_path / "identify")
    detailed = parse_detailed_output(detailed_file)

    for sample, expected in EXPECTED_GENOTYPES.items():
        observed = genotype_table[sample]
        assert observed["nsequences"] == expected["nsequences"]
        for repeat_bin, count in expected.items():
            if repeat_bin == "nsequences":
                continue
            assert observed[repeat_bin] == count

    for key, expected in EXPECTED_DETAILED.items():
        ids = detailed[key]
        assert len(ids) == expected["count"]
        assert ids[:3] == expected["first_ids"]



def test_extract_returns_the_requested_genotype_bin(tmp_path):
    demux_dir = run_demultiplex(tmp_path / "demultiplexed")
    _, detailed_file, _ = run_identify(demux_dir, tmp_path / "identify")

    output_file = tmp_path / "stmock03_ssa410uos_22x.fas"
    args = SimpleNamespace(
        infile=str(SUBSET_FASTA),
        genotypes=str(detailed_file),
        locus="STMock03_Ssa410Uos",
        number="22x",
        outfile=str(output_file),
    )

    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        assert extract.main(args) == 0

    ids = fasta_ids(output_file)
    assert len(ids) == 156
    assert ids[:3] == EXPECTED_DETAILED[("STMock03_Ssa410Uos", "22x")]["first_ids"]



def test_extract_accepts_plain_fastq_input(tmp_path):
    fastq_file = write_fastq_from_fasta(SUBSET_FASTA, tmp_path / "pipeline_subset_large.fastq")
    demux_dir = run_demultiplex(tmp_path / "demultiplexed", fastq_file)
    _, detailed_file, _ = run_identify(demux_dir, tmp_path / "identify")

    output_file = tmp_path / "stmock03_ssa410uos_22x_from_fastq.fas"
    args = SimpleNamespace(
        infile=str(fastq_file),
        genotypes=str(detailed_file),
        locus="STMock03_Ssa410Uos",
        number="22x",
        outfile=str(output_file),
    )

    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        assert extract.main(args) == 0

    ids = fasta_ids(output_file)
    assert len(ids) == 156
    assert ids[:3] == EXPECTED_DETAILED[("STMock03_Ssa410Uos", "22x")]["first_ids"]



def test_visualize_writes_a_pdf_figure(tmp_path):
    demux_dir = run_demultiplex(tmp_path / "demultiplexed")
    _, _, genotype_file = run_identify(demux_dir, tmp_path / "identify")

    output_file = tmp_path / "stmock03_ssa410uos_distribution.pdf"
    args = SimpleNamespace(
        genotypes=str(genotype_file),
        locus="STMock03_Ssa410Uos",
        outfile=str(output_file),
    )

    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        assert visualize.main(args) == 0

    assert output_file.exists()
    assert output_file.stat().st_size > 1000
    with open(output_file, "rb") as handle:
        assert handle.read(4) == b"%PDF"



def test_report_writes_a_single_html_summary(tmp_path):
    demux_dir = run_demultiplex(tmp_path / "demultiplexed")
    _, _, genotype_file = run_identify(demux_dir, tmp_path / "identify")

    output_file = tmp_path / "report.html"
    args = SimpleNamespace(
        demultiplexed=str(demux_dir),
        genotypes=str(genotype_file),
        outfile=str(output_file),
    )

    with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
        assert report.main(args) == 0

    html = output_file.read_text()
    assert "ontyper report" in html
    assert "Interactive summary of demultiplexing and genotype identification results." in html
    assert "STMock03_Ssa410Uos" in html
    assert 'const REPORT_DATA =' in html
    assert 'Repeat count histogram' in html
