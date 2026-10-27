# Varaser

Varaser erases germline variants from a BAM file by restoring the reference base at user-specified variant positions. Reads overlapping those positions are edited so that the sequence, CIGAR, and tags are consistent with the reference, which reduces the individual-identifying information in the BAM.

## Files

- `Varaser.py`: current release version (2026-09-06)

## Usage

```bash
python3 Varaser.py <bam> <vcf> <prefix> <outdir> -N <thread> -R <ref.fasta> [-T fastq] [--force-sort]
```

| Argument | Description |
|---|---|
| `bam` | Input BAM (coordinate-sorted, indexed) |
| `vcf` | VCF of germline variants to erase |
| `prefix` | Prefix for output file names |
| `outdir` | Output directory (created if missing) |
| `-N` | Number of threads (default: 1) |
| `-R`, `--reference` | Reference FASTA, required. Used to recalculate NM/MD of modified reads; use the same reference the BAM was aligned to |
| `-T` | Output type; `fastq` also writes FASTQ files (default: `bam`) |
| `--force-sort` | Sort the final BAM even if it appears to be sorted already |

## Outputs

- `<outdir>/<prefix>.Varaser.sorted.bam` (+ `.bai`): BAM with the variants erased
- `<outdir>/<prefix>.runtime.txt`: total run time (`Run`) and erasing time (`Erase`), in minutes
- With `-T fastq`: `<prefix>_1.fastq.gz` and `<prefix>_2.fastq.gz`

## Processing rules

- Only reads overlapping a target variant are modified; all other reads are kept unchanged.
- SNV/MNP: the base is replaced with the reference base only when the read carries exactly the alternate allele at that position.
- Insertion/deletion: the CIGAR is updated and the sequence and base qualities are adjusted. Deleted reference bases are restored with the mean base quality of the read.
- Spliced reads (N operations): intron lengths are kept in their original order. If the number of N operations does not match, the read is left unmodified and counted in the summary.
- Tags: original tags are preserved with their types. Only `NM` and `MD` of modified reads are recalculated against the reference.
- Tags not recalculated: `AS`, `XS`, `XA`, `SA`, `NH`, `HI`, and other tags keep their original values.
- Known limitation: `MC` (mate CIGAR) and `TLEN` of the mate of a modified read are not updated. To update them, run `samtools fixmate -m` on the output, then sort by coordinate again.

## Dependencies

- Python 3: `pysam`, `pandas`, `sortedcontainers`
- External tools: `samtools` (validated with 1.21), `bgzip` and `tabix` (htslib)

## Notes

- For normally oriented pairs, TLEN was identical before and after `fixmate -m`. For same-strand pairs (both mates on the same strand), `fixmate` can differ by 1 bp.
