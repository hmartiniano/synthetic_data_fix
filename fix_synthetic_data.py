"""Add standard QC FORMAT fields (GQ, DP, AD) to synthetic VCF files.

Fields added:
  - GQ: Genotype Quality (default: 99)
  - DP: Read Depth (default: 30)
  - AD: Allelic Depth (e.g. 30,0 for hom-ref, 15,15 for het, 0,30 for hom-alt)

Optionally, a fraction (0% to 5% per variant, random) of
genotypes can be simulated as having insufficient quality (low GQ, low DP,
or skewed heterozygous allele balance).
"""

from __future__ import annotations

import argparse
import gzip
import io
import random
import sys
from typing import Dict, List, Optional, TextIO

GQ_HEADER = '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n'
DP_HEADER = '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
AD_HEADER = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments.

    Parameters
    ----------
    args : Optional[List[str]]
        List of command-line arguments. If None, uses sys.argv[1:].

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.

    """
    parser = argparse.ArgumentParser(
        description="Add QC FORMAT fields (GQ, DP, AD) to synthetic VCF files."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_vcf",
        default="-",
        help="Input VCF file path (.vcf or .vcf.gz) or '-' for stdin (default: -).",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_vcf",
        default="-",
        help="Output VCF path (.vcf or .vcf.gz) or '-' for stdout (default: -).",
    )
    parser.add_argument(
        "--simulate-low-quality",
        action="store_true",
        default=False,
        help=(
            "Simulate a random fraction (0 to max) of genotypes per variant "
            "with insufficient quality (failing GQ, DP, or AB thresholds)."
        ),
    )
    parser.add_argument(
        "--max-low-quality-fraction",
        type=float,
        default=0.05,
        help=(
            "Max fraction of low-quality genotypes per variant when "
            "simulation is enabled (default: 0.05)."
        ),
    )
    parser.add_argument(
        "--default-gq",
        type=int,
        default=99,
        help="Default Genotype Quality (GQ) for high-quality calls (default: 99).",
    )
    parser.add_argument(
        "--default-dp",
        type=int,
        default=30,
        help="Default read depth (DP) for high-quality calls (default: 30).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible runs.",
    )
    return parser.parse_args(args)


def build_baseline_mapping(gq: int = 99, dp: int = 30) -> Dict[str, str]:
    """Pre-build constant-time lookup table for common genotypes.

    Parameters
    ----------
    gq : int
        Default genotype quality.
    dp : int
        Default total depth.

    Returns
    -------
    Dict[str, str]
        Dictionary mapping GT tokens to GT:GQ:DP:AD formatted strings.

    """
    half_dp = dp // 2
    rem_dp = dp - half_dp

    return {
        # Diploid phased
        "0|0": f"0|0:{gq}:{dp}:{dp},0",
        "0|1": f"0|1:{gq}:{dp}:{half_dp},{rem_dp}",
        "1|0": f"1|0:{gq}:{dp}:{rem_dp},{half_dp}",
        "1|1": f"1|1:{gq}:{dp}:0,{dp}",
        # Diploid unphased
        "0/0": f"0/0:{gq}:{dp}:{dp},0",
        "0/1": f"0/1:{gq}:{dp}:{half_dp},{rem_dp}",
        "1/0": f"1/0:{gq}:{dp}:{rem_dp},{half_dp}",
        "1/1": f"1/1:{gq}:{dp}:0,{dp}",
        # Haploid / hemizygous
        "0": f"0:{gq}:{dp}:{dp}",
        "1": f"1:{gq}:{dp}:0,{dp}",
        # Missing
        "./.": "./.:.:0:0,0",
        ".|.": ".|.:.:0:0,0",
        ".": ".:.:0:0",
    }


def format_call_fallback(gt: str, gq: int, dp: int) -> str:
    """Format an uncommon or multi-allelic genotype string.

    Parameters
    ----------
    gt : str
        Genotype string (e.g. '0/2', '2/2', './1').
    gq : int
        Genotype quality.
    dp : int
        Total read depth.

    Returns
    -------
    str
        Formatted GT:GQ:DP:AD string.

    """
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(sep)

    # Missing allele handling
    if all(a == "." for a in alleles):
        return f"{gt}:.:0:{','.join(['0'] * len(alleles))}"

    # Extract integer allele indices if possible
    try:
        int_alleles = [int(a) for a in alleles if a != "."]
    except ValueError:
        return f"{gt}:{gq}:{dp}:{dp},0"

    max_allele = max(int_alleles) if int_alleles else 0
    ad_counts = [0] * (max_allele + 1)

    if len(int_alleles) == 1:
        # Haploid
        ad_counts[int_alleles[0]] = dp
    elif len(int_alleles) == 2:
        # Diploid
        a1, a2 = int_alleles
        if a1 == a2:
            ad_counts[a1] = dp
        else:
            half = dp // 2
            ad_counts[a1] = half
            ad_counts[a2] = dp - half
    else:
        # Polyploid
        share = dp // len(int_alleles)
        for a in int_alleles:
            ad_counts[a] += share

    ad_str = ",".join(str(c) for c in ad_counts)
    return f"{gt}:{gq}:{dp}:{ad_str}"


def make_low_quality_call(gt: str, rng: random.Random) -> str:
    """Generate an insufficient quality call (failing GQ, DP, or AB thresholds).

    Parameters
    ----------
    gt : str
        Genotype string (e.g. '0/1', '0/0', '1|1').
    rng : random.Random
        Random generator instance.

    Returns
    -------
    str
        GT:GQ:DP:AD call with failing quality parameters.

    """
    if gt in ("./.", ".|.", "."):
        return f"{gt}:.:0:0,0"

    # Modes of failure:
    # 0: Low GQ (< 20, e.g. 5-18)
    # 1: Low DP (< 10, e.g. 2-8)
    # 2: Skewed allele balance for hets (AB < 0.2 or > 0.8)
    alleles = gt.replace("|", "/").split("/")
    is_het = len(alleles) > 1 and len(set(alleles)) > 1

    modes = [0, 1]
    if is_het:
        modes.append(2)

    mode = rng.choice(modes)

    if mode == 0:
        # Low GQ (threshold min_gq is 20)
        low_gq = rng.randint(5, 18)
        dp = rng.randint(20, 35)
        half = dp // 2
        if is_het:
            return f"{gt}:{low_gq}:{dp}:{half},{dp - half}"
        if "1" in gt:
            return f"{gt}:{low_gq}:{dp}:0,{dp}"
        return f"{gt}:{low_gq}:{dp}:{dp},0"

    if mode == 1:
        # Low DP (threshold min_dp is 10)
        low_dp = rng.randint(2, 8)
        gq = rng.randint(10, 30)
        half = low_dp // 2
        if is_het:
            return f"{gt}:{gq}:{low_dp}:{half},{low_dp - half}"
        if "1" in gt:
            return f"{gt}:{gq}:{low_dp}:0,{low_dp}"
        return f"{gt}:{gq}:{low_dp}:{low_dp},0"

    # Skewed AB for het (normal is ~0.5; failing is < 0.2 or > 0.8)
    gq = rng.randint(30, 99)
    dp = rng.randint(25, 40)
    if rng.random() < 0.5:
        # ALT ratio < 0.2
        alt_reads = rng.randint(1, max(1, int(dp * 0.15)))
        ref_reads = dp - alt_reads
    else:
        # ALT ratio > 0.8
        alt_reads = rng.randint(int(dp * 0.85), dp - 1)
        ref_reads = dp - alt_reads

    return f"{gt}:{gq}:{dp}:{ref_reads},{alt_reads}"


def fix_synthetic_vcf(
    input_file: TextIO,
    output_file: TextIO,
    simulate_low_quality: bool = False,
    max_low_quality_fraction: float = 0.05,
    default_gq: int = 99,
    default_dp: int = 30,
    seed: Optional[int] = None,
) -> None:
    """Stream and fix synthetic VCF data with GQ, DP, and AD FORMAT tags.

    Parameters
    ----------
    input_file : TextIO
        Input text stream reading VCF data.
    output_file : TextIO
        Output text stream for writing enriched VCF.
    simulate_low_quality : bool
        Whether to introduce a random fraction of low-quality calls.
    max_low_quality_fraction : float
        Maximum fraction of low-quality calls per variant.
    default_gq : int
        Default genotype quality.
    default_dp : int
        Default read depth.
    seed : Optional[int]
        Random seed for reproducible low-quality simulation.

    """
    rng = random.Random(seed)
    baseline_mapping = build_baseline_mapping(default_gq, default_dp)

    header_injected = False

    for line in input_file:
        if line.startswith("##"):
            output_file.write(line)
            if line.startswith("##FORMAT=<ID=GT,") and not header_injected:
                output_file.write(GQ_HEADER)
                output_file.write(DP_HEADER)
                output_file.write(AD_HEADER)
                header_injected = True
            continue

        if line.startswith("#CHROM"):
            if not header_injected:
                output_file.write(GQ_HEADER)
                output_file.write(DP_HEADER)
                output_file.write(AD_HEADER)
                header_injected = True
            output_file.write(line)
            continue

        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 10:
            output_file.write(line)
            continue

        format_field = parts[8]

        # Only transform if FORMAT is strictly GT
        if format_field == "GT":
            parts[8] = "GT:GQ:DP:AD"
            samples = parts[9:]
            num_samples = len(samples)

            if not simulate_low_quality or max_low_quality_fraction <= 0:
                # Fast O(1) path
                enriched_samples = [
                    baseline_mapping.get(s)
                    or format_call_fallback(s, default_gq, default_dp)
                    for s in samples
                ]
            else:
                # Random fraction of low-quality calls in [0, max_fraction]
                lq_rate = rng.uniform(0.0, max_low_quality_fraction)
                num_lq = int(num_samples * lq_rate)

                lq_indices = (
                    set(rng.sample(range(num_samples), num_lq)) if num_lq > 0 else set()
                )

                enriched_samples = []
                for idx, s in enumerate(samples):
                    if idx in lq_indices:
                        enriched_samples.append(make_low_quality_call(s, rng))
                    else:
                        enriched = baseline_mapping.get(s) or format_call_fallback(
                            s, default_gq, default_dp
                        )
                        enriched_samples.append(enriched)

            output_file.write("\t".join(parts[:9] + enriched_samples) + "\n")
        else:
            # If FORMAT already contains other fields, leave line intact
            output_file.write(line)


def main() -> None:
    """CLI entrypoint for fixing synthetic VCF files."""
    args = parse_args()

    # Open input stream
    if args.input_vcf == "-":
        input_stream: TextIO = sys.stdin
    elif args.input_vcf.endswith(".gz"):
        input_stream = io.TextIOWrapper(gzip.open(args.input_vcf, "rb"))
    else:
        input_stream = open(args.input_vcf, "r", encoding="utf-8")

    # Open output stream
    if args.output_vcf == "-":
        output_stream: TextIO = sys.stdout
    elif args.output_vcf.endswith(".gz"):
        output_stream = io.TextIOWrapper(gzip.open(args.output_vcf, "wb"))
    else:
        output_stream = open(args.output_vcf, "w", encoding="utf-8")

    try:
        fix_synthetic_vcf(
            input_stream,
            output_stream,
            simulate_low_quality=args.simulate_low_quality,
            max_low_quality_fraction=args.max_low_quality_fraction,
            default_gq=args.default_gq,
            default_dp=args.default_dp,
            seed=args.seed,
        )
    finally:
        if input_stream is not sys.stdin:
            input_stream.close()
        if output_stream is not sys.stdout:
            output_stream.close()


if __name__ == "__main__":
    main()
