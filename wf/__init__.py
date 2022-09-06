"""
Low complexity filtering for reads
"""

import re
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from latch import message, small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import LatchDir, LatchFile

from .docs import metadata


# From: https://github.com/latch-verified/bulk-rnaseq/blob/64a25531e1ddc43be0afffbde91af03754fb7c8c/wf/__init__.py
def _capture_output(command: List[str]) -> Tuple[int, str]:
    captured_stdout = []

    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=1,
        universal_newlines=True,
    ) as process:
        assert process.stdout is not None
        for line in process.stdout:
            print(line)
            captured_stdout.append(line)
        process.wait()
        returncode = process.returncode

    return returncode, "\n".join(captured_stdout)


@small_task
def bbduk(
    read1: LatchFile,
    read2: LatchFile,
    sample_name: str,
    contaminants: Optional[LatchFile],
) -> LatchDir:

    output_dir_name = f"{sample_name}_bbduk_outputs"
    output_dir = Path(output_dir_name).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _bbduk_cmd = [
        "bbduk.sh",
        f"in1={read1.local_path}",
        f"in2={read2.local_path}",
        f"out1={str(output_dir)}/{sample_name}_1_filtered.fastq",
        f"out2={str(output_dir)}/{sample_name}_2_filtered.fastq",
        "threads=31",
    ]

    if contaminants is not None:
        _bbduk_cmd.extend(
            [
                f"ref={contaminants.local_path}",
            ]
        )

    return_code, stdout = _capture_output(_bbduk_cmd)

    version = re.findall("Version.*", stdout)[0]
    running_cmd = re.findall("Executing.*", stdout)[0]

    message(
        "info",
        {
            "title": f"Executing bbduk {version}",
            "body": running_cmd,
        },
    )

    if return_code != 0:
        errors = re.findall("Exception.*", stdout[1])
        for error in errors:
            message(
                "error",
                {
                    "title": f"An error was raised while running bbduk for {sample_name}:",
                    "body": error,
                },
            )
        raise RuntimeError

    return LatchDir(str(output_dir), f"latch:///{output_dir_name}/")


@workflow(metadata)
def highcomplexity(
    read1: LatchFile,
    read2: LatchFile,
    sample_name: str = "BBDuk_Sample",
    contaminants: Optional[LatchFile] = None,
) -> LatchDir:
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
    highcomplexity,
    "Crohn's disease gut microbiome (SRR579292)",
    {
        "read1": LatchFile("s3://latch-public/test-data/4318/SRR579292_1.fastq"),
        "read2": LatchFile("s3://latch-public/test-data/4318/SRR579292_2.fastq"),
        "sample_name": "SRR579292",
        "contaminants": LatchFile("s3://latch-public/test-data/4318/adapters.fasta"),
    },
)
