from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NamingCase:
    label: str
    filenames: tuple[str, ...]
    paired: bool
    expected: frozenset[str]


CASES: tuple[NamingCase, ...] = (
    NamingCase(
        "pe_R1R2",
        ("sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_R1R2_001",
        ("sampleA_R1_001.fastq.gz", "sampleA_R2_001.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_lower_r1r2",
        ("sampleA_r1.fastq.gz", "sampleA_r2.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_read1read2",
        ("sampleA_read1.fastq.gz", "sampleA_read2.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_1_2",
        ("sampleA_1.fastq.gz", "sampleA_2.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_forward_reverse",
        ("sampleA_forward.fastq.gz", "sampleA_reverse.fastq.gz"),
        True,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "pe_multi",
        ("s1_R1.fastq.gz", "s1_R2.fastq.gz", "s2_R1.fastq.gz", "s2_R2.fastq.gz"),
        True,
        frozenset({"s1", "s2"}),
    ),
    NamingCase(
        "se_plain",
        ("sampleA.fastq.gz",),
        False,
        frozenset({"sampleA"}),
    ),
    NamingCase(
        "se_keeps_suffix",
        ("sampleA_R1.fastq.gz",),
        False,
        frozenset({"sampleA_R1"}),
    ),
)
