"""
Assemble and sort some COVID reads...
"""

import subprocess
from typing import List, Optional

from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchFile, file_glob

from .docs import metadata


@small_task
def bbduk(
    read1: LatchFile,
    read2: LatchFile,
    sample_name: str,
    contaminants: Optional[LatchFile],
) -> List[LatchFile]:

    _bbduk_cmd = [
        "bbduk.sh",
        "in1=",
        read1.local_path,
        "in2=",
        read2.local_path,
        "out1=",
        f"{sample_name}_1_filtered.fastq",
        "out2=",
        f"{sample_name}_2_filtered.fastq",
        "threads=",
        "31",
    ]

    if contaminants is not None:
        _bbduk_cmd.extend(
            [
                "ref=",
                contaminants.local_path,
            ]
        )

    subprocess.run(_bbduk_cmd)

    return file_glob("*filtered.fastq", "latch:///bbduk_outputs")


@workflow(metadata)
def assemble_and_sort(
    read1: LatchFile,
    read2: LatchFile,
    sample_name: str = "BBDuk_Sample",
    contaminants: Optional[LatchFile] = None,
) -> List[LatchFile]:
    """Low complexity read filtering

    HighComplexity
    ---

    Low complexity filtering for short read datasets with
    [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) [^1].

    ---
    [^1]: Bushnell B, Rood J, Singer E (2017) BBMerge –
    Accurate paired shotgun read merging via overlap.
    PLOS ONE 12(10): e0185056. https://doi.org/10.1371/journal.pone.0185056

    """

    return bbduk(
        read1=read1, read2=read2, sample_name=sample_name, contaminants=contaminants
    )


LaunchPlan(
    assemble_and_sort,
    "Crohn's disease gut microbiome (SRR579292)",
    {
        "read1": LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
        "read2": LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
        "sample_name": "SRR579292",
        "contaminants": LatchFile("s3://latch-public/test-data/4318/adapters.fasta"),
    },
)
