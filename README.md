# ontyper

Microsatellite genotyping toolkit for Oxford Nanopore sequencing data.

## Installation

First you need to clone the repository:

```bash
git clone git@github.com:reslp/ontyper.git
```

### As a global CLI tool (recommended)

```bash
# Install globally using uv tool (from local directory)
cd /path/to/ontyper
uv tool install .

# Or install directly from a git repository (once published)
# uv tool install git+https://github.com/username/ontyper.git

# The 'ontyper' command is now available globally
ontyper --help

# To upgrade later
uv tool upgrade ontyper

# To uninstall
uv tool uninstall ontyper
```

### For development

```bash
# Using pip (editable install)
pip install -e .

# Using uv in a virtual environment
uv venv
source .venv/bin/activate  # or .venv-mba/bin/activate on some systems
uv pip install -e .
```

## Usage

### Demultiplex reads

```bash
ontyper demultiplex -i reads.fasta -d demultiplex.txt -o output_dir/
ontyper demultiplex -i reads.fastq.gz -d demultiplex.txt -o output_dir/
```

### Identify genotypes

```bash
ontyper genotypes identify -d demultiplex.txt -i demultiplexed_dir/ > genotypes.tsv

# Or write the table directly to a file
ontyper genotypes identify -d demultiplex.txt -i demultiplexed_dir/ -o genotypes.tsv

# With detailed output (creates genotypes-detailed.tsv)
ontyper genotypes identify -d demultiplex.txt -i demultiplexed_dir/ --detailed -o genotypes.tsv
```

### Extract sequences for a specific genotype

```bash
ontyper genotypes extract -i reads.fasta -g genotypes-detailed.tsv -l TTMock01_Tth445 -n 20x
ontyper genotypes extract -i reads.fastq.gz -g genotypes-detailed.tsv -l TTMock01_Tth445 -n 20x
```

### Visualize genotype distributions

```bash
ontyper genotypes visualize -g genotypes.tsv -l TTMock01_Tth445
```

### Build an interactive HTML report

```bash
ontyper report -d demultiplexed_dir/ -g genotypes.tsv -o report.html
```

## Command Structure

```
ontyper
├── demultiplex              # Demultiplex ONT reads by primer/tag
├── report                   # Build an interactive single-file HTML report
└── genotypes                # Genotype-related commands
    ├── identify             # Identify microsatellite repeat genotypes
    ├── extract              # Extract sequences for a specific genotype bin
    └── visualize            # Plot repeat-count distributions for one sample/locus
```

## Development

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest
```

## License

MIT
