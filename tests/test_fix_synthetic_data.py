"""Unit and integration tests for fix_synthetic_data.py."""

from __future__ import annotations

import io
import random

from fix_synthetic_data import (
    build_baseline_mapping,
    fix_synthetic_vcf,
    format_call_fallback,
    make_low_quality_call,
    parse_args,
)


def test_build_baseline_mapping() -> None:
    """Test standard baseline mapping dictionary."""
    mapping = build_baseline_mapping(gq=99, dp=30)

    assert mapping["0/0"] == "0/0:99:30:30,0"
    assert mapping["0|0"] == "0|0:99:30:30,0"
    assert mapping["0/1"] == "0/1:99:30:15,15"
    assert mapping["1|0"] == "1|0:99:30:15,15"
    assert mapping["1/1"] == "1/1:99:30:0,30"
    assert mapping["0"] == "0:99:30:30"
    assert mapping["1"] == "1:99:30:0,30"
    assert mapping["./."] == "./.:.:0:0,0"
    assert mapping[".|."] == ".|.:.:0:0,0"
    assert mapping["."] == ".:.:0:0"


def test_format_call_fallback_multi_allelic() -> None:
    """Test formatting fallback for multi-allelic and hemizygous variants."""
    res_0_2 = format_call_fallback("0/2", gq=99, dp=30)
    assert res_0_2 == "0/2:99:30:15,0,15"

    res_2_2 = format_call_fallback("2|2", gq=99, dp=30)
    assert res_2_2 == "2|2:99:30:0,0,30"

    res_hemi = format_call_fallback("2", gq=99, dp=30)
    assert res_hemi == "2:99:30:0,0,30"

    res_missing = format_call_fallback("./.", gq=99, dp=30)
    assert res_missing == "./.:.:0:0,0"


def test_make_low_quality_call() -> None:
    """Test generating failing QC calls (low GQ, low DP, or skewed AB)."""
    rng = random.Random(42)

    # Test het calls over multiple iterations
    for _ in range(30):
        call = make_low_quality_call("0/1", rng)
        parts = call.split(":")
        assert len(parts) == 4
        gt, gq_str, dp_str, ad_str = parts
        assert gt == "0/1"
        gq = int(gq_str)
        dp = int(dp_str)
        ad_vals = [int(x) for x in ad_str.split(",")]

        # Must fail at least one threshold: GQ < 20, DP < 10, or AB out of [0.2, 0.8]
        ab = ad_vals[1] / sum(ad_vals) if sum(ad_vals) > 0 else 0.0
        fails_qc = (gq < 20) or (dp < 10) or (ab < 0.2 or ab > 0.8)
        assert fails_qc, f"Call {call} did not fail any QC threshold"


def test_fix_synthetic_vcf_baseline() -> None:
    """Test streaming enrichment without low-quality simulation."""
    input_vcf = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
        "chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0|0\t0|1\t1|1\n"
        "chrX\t200\t.\tC\tT\t100\tPASS\t.\tGT\t0\t1\t./.\n"
    )

    in_stream = io.StringIO(input_vcf)
    out_stream = io.StringIO()

    fix_synthetic_vcf(
        in_stream,
        out_stream,
        simulate_low_quality=False,
        default_gq=99,
        default_dp=30,
    )

    output = out_stream.getvalue()
    lines = output.splitlines()

    # Check header
    assert (
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">' in lines
    )
    assert '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">' in lines
    assert (
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">' in lines
    )

    # Check variant 1
    v1 = [line for line in lines if line.startswith("chr1\t100")][0].split("\t")
    assert v1[8] == "GT:GQ:DP:AD"
    assert v1[9] == "0|0:99:30:30,0"
    assert v1[10] == "0|1:99:30:15,15"
    assert v1[11] == "1|1:99:30:0,30"

    # Check variant 2 (haploid + missing)
    v2 = [line for line in lines if line.startswith("chrX\t200")][0].split("\t")
    assert v2[8] == "GT:GQ:DP:AD"
    assert v2[9] == "0:99:30:30"
    assert v2[10] == "1:99:30:0,30"
    assert v2[11] == "./.:.:0:0,0"


def test_fix_synthetic_vcf_with_simulation() -> None:
    """Test low quality call simulation (0-5% per variant)."""
    num_samples = 200
    sample_ids = "\t".join(f"S{i}" for i in range(num_samples))
    sample_gts = "\t".join("0/1" for _ in range(num_samples))

    input_vcf = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_ids}\n"
        f"chr1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t{sample_gts}\n"
    )

    in_stream = io.StringIO(input_vcf)
    out_stream = io.StringIO()

    fix_synthetic_vcf(
        in_stream,
        out_stream,
        simulate_low_quality=True,
        max_low_quality_fraction=0.05,
        default_gq=99,
        default_dp=30,
        seed=123,
    )

    output = out_stream.getvalue()
    v1 = [line for line in output.splitlines() if line.startswith("chr1\t100")][
        0
    ].split("\t")
    calls = v1[9:]
    assert len(calls) == num_samples

    # Count calls with failing QC
    failing_calls = 0
    for call in calls:
        gt, gq_str, dp_str, ad_str = call.split(":")
        gq = int(gq_str)
        dp = int(dp_str)
        ad = [int(x) for x in ad_str.split(",")]
        ab = ad[1] / sum(ad) if sum(ad) > 0 else 0.0
        if (gq < 20) or (dp < 10) or (ab < 0.2 or ab > 0.8):
            failing_calls += 1

    # In 200 samples with max 5%, failing count must be between 0 and 10 (<= 5%)
    assert 0 <= failing_calls <= int(num_samples * 0.05)


def test_parse_args() -> None:
    """Test CLI argument parsing."""
    args = parse_args(
        ["-i", "input.vcf", "-o", "output.vcf.gz", "--simulate-low-quality"]
    )
    assert args.input_vcf == "input.vcf"
    assert args.output_vcf == "output.vcf.gz"
    assert args.simulate_low_quality is True
    assert args.max_low_quality_fraction == 0.05
    assert args.default_gq == 99
    assert args.default_dp == 30
