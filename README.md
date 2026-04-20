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
```

### Identify genotypes

```bash
ontyper genotypes identify -d demultiplex.txt -i demultiplexed_dir/ > genotypes.tsv

# With detailed output (creates genotypes-detailed.tsv)
ontyper genotypes identify -d demultiplex.txt -i demultiplexed_dir/ --detailed > genotypes.tsv
```

### Extract sequences for a specific genotype

```bash
ontyper genotypes extract -i reads.fasta -g genotypes-detailed.tsv -l TTMock01_Tth445 -n 20x
```

## Command Structure

```
ontyper
├── demultiplex              # Demultiplex ONT reads by primer/tag
└── genotypes                # Genotype-related commands
    ├── identify             # Identify microsatellite repeat genotypes
    └── extract              # Extract sequences for a specific genotype bin
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
