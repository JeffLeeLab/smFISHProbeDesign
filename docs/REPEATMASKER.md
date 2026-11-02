# RepeatMasker Setup for ProbeDesign

This guide explains how to install and configure RepeatMasker for automatic repeat masking in ProbeDesign.

## Overview

RepeatMasker screens DNA sequences for interspersed repeats and low complexity DNA sequences. When integrated with ProbeDesign, it automatically masks repeat regions to avoid designing probes in those areas.

## Quick Start (Human/Mouse)

```bash
# 0. Activate the environment
mamba activate probedesign

# 1. Install RepeatMasker
mamba install -n probedesign -c conda-forge -c bioconda repeatmasker

# 2. Download the required database partition (~8.9GB compressed, ~56GB extracted)
cd $CONDA_PREFIX/share/RepeatMasker/Libraries/famdb
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz
gunzip dfam39_full.7.h5.gz

# 3. Test it works
probedesign design input.fa --repeatmask --probes 32
```

All commands assume the `probedesign` environment is active (`mamba activate probedesign`). `$CONDA_PREFIX` then points at it, e.g. `~/miniforge3/envs/probedesign`. Install RepeatMasker into that environment, never into `base`.

## Installation Details

### Step 1: Install RepeatMasker via Conda

```bash
mamba install -n probedesign -c conda-forge -c bioconda repeatmasker
```

This installs:
- RepeatMasker 4.2.2
- RMBlast (search engine)
- TRF (Tandem Repeat Finder)
- HMMER
- Root database partition (dfam39_full.0.h5, ~77MB)

### Step 2: Download Species-Specific Database

RepeatMasker requires the [Dfam](https://www.dfam.org/) database. The database is split into partitions by taxonomic group. You only need to download the partition(s) for your species of interest.

#### Database Partitions

| Partition | Taxonomic Group | Size (compressed) | Species Examples |
|-----------|-----------------|-------------------|------------------|
| 0 (root) | Root families | 17 MB | **Required for all** |
| 7 | **Mammalia** | **8.9 GB** | **Human, Mouse, Rat** |
| 1-6, 8-16 | Other groups | 4-11 GB each | Plants, insects, fish, etc. |

**For human and mouse sequences, you need partition 0 (included) + partition 7.**

#### Download Commands

```bash
# Navigate to RepeatMasker's famdb directory
cd $CONDA_PREFIX/share/RepeatMasker/Libraries/famdb

# Download Mammalia partition (for human, mouse, rat, etc.)
curl -O https://www.dfam.org/releases/current/families/FamDB/dfam39_full.7.h5.gz

# Extract (creates ~56GB file - ensure sufficient disk space!)
gunzip dfam39_full.7.h5.gz
```

**Disk Space Requirements:**
- Compressed download: ~8.9 GB
- Extracted database: ~56 GB
- Total needed during extraction: ~65 GB free space

### Step 3: Verify Installation

```bash
# Check RepeatMasker version
$CONDA_PREFIX/share/RepeatMasker/RepeatMasker -v

# Test with human species (requires partition 7)
PATH="$CONDA_PREFIX/bin:$PATH" \
  $CONDA_PREFIX/share/RepeatMasker/RepeatMasker \
  -species human -help 2>&1 | head -10
```

## Usage with ProbeDesign

### Automatic Repeat Masking

Use the `--repeatmask` flag to automatically run RepeatMasker:

```bash
# Human sequences (default)
probedesign design input.fa --repeatmask --probes 32

# Mouse sequences
probedesign design input.fa --repeatmask --species mouse --probes 32
```

This will:
1. Run RepeatMasker on your input FASTA
2. Mask repeat regions with N's
3. Design probes avoiding those regions

### Manual Repeat Masking

If you prefer to run RepeatMasker separately or use pre-masked files:

```bash
# Option 1: Use a pre-masked file
probedesign design input.fa --repeatmask-file input_repeatmasked.fa --probes 32

# Option 2: Run RepeatMasker manually first
RepeatMasker -species human input.fa
probedesign design input.fa --repeatmask-file input.fa.masked --probes 32
```

**Note:** You cannot use both `--repeatmask` and `--repeatmask-file` together.

### Combined with Other Masking

RepeatMasker can be combined with pseudogene and genome masking:

```bash
probedesign design input.fa \
  --repeatmask \
  --pseudogene-mask \
  --genome-mask \
  --index-dir bowtie_indexes \
  --probes 32
```

## Supported Species

| Species | Flag Value | Database Partition |
|---------|------------|-------------------|
| Human | `--species human` | Partition 7 (Mammalia) |
| Mouse | `--species mouse` | Partition 7 (Mammalia) |
| Rat | `--species rat` | Partition 7 (Mammalia) |

For other species, check the [Dfam taxonomy](https://www.dfam.org/browse) to determine which partition is needed.

## Troubleshooting

### "partition is absent" Error

```
Taxon "human" is in partition 7 of the current FamDB however,
this partition is absent.
```

**Solution:** Download partition 7 as described above.

### "No module named 'h5py'" Error

RepeatMasker needs Python with h5py. ProbeDesign handles this automatically by setting the PATH, but if running RepeatMasker manually:

```bash
# Run with conda Python in PATH
PATH="$CONDA_PREFIX/bin:$PATH" RepeatMasker -species human input.fa
```

### "No space left on device" During Extraction

The extracted database is ~56GB. Ensure you have at least 65GB free space before extraction.

```bash
# Check available space
df -h $CONDA_PREFIX/

# If needed, clean up and retry
rm -f dfam39_full.7.h5  # Remove partial extraction
gunzip dfam39_full.7.h5.gz
```

### RepeatMasker Not Found

```bash
# Install via conda
mamba install -n probedesign -c conda-forge -c bioconda repeatmasker

# Verify installation
$CONDA_PREFIX/share/RepeatMasker/RepeatMasker -v
```

`find_repeatmasker()` in `masking.py` first checks a hard-coded Homebrew-cask path (`/opt/homebrew/Caskroom/miniforge/base/share/RepeatMasker/RepeatMasker`), then `RepeatMasker` on `PATH`. With the environment active, the `PATH` lookup finds the conda install.

## Database Updates

Dfam databases are updated periodically (~annually). To update:

1. Check current version: https://www.dfam.org/releases/
2. Download new partition files (naming convention: `dfam{version}_full.{partition}.h5.gz`)
3. Replace existing `.h5` files in the famdb directory

## Performance Notes

- **First run**: RepeatMasker builds species-specific libraries on first use (~1-2 minutes)
- **Typical runtime**: 30 seconds to 2 minutes for typical FISH target sequences
- **Memory usage**: ~2-4 GB RAM during analysis
- **For batch processing**: Consider pre-masking all sequences and using `--repeatmask-file`

## Alternative: Web-Based RepeatMasker

If you cannot install RepeatMasker locally (e.g., disk space constraints), use the web service:

1. Go to [RepeatMasker Web Server](http://www.repeatmasker.org/cgi-bin/WEBRepeatMasker)
2. Upload your FASTA file
3. Select species and options
4. Download the `.masked` output file
5. Use with ProbeDesign: `probedesign design input.fa --repeatmask-file input.fa.masked`

## References

- [Dfam Database](https://www.dfam.org/)
- [RepeatMasker Documentation](https://www.repeatmasker.org/webrepeatmaskerhelp.html)
- [FamDB GitHub](https://github.com/Dfam-consortium/FamDB)
