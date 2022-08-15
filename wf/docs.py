from latch.types import LatchAuthor, LatchMetadata, LatchParameter

metadata = LatchMetadata(
    display_name="HighComplexity",
    documentation="https://github.com/jvfe/highcomplex_latch/blob/main/README.md",
    author=LatchAuthor(
        name="jvfe",
        github="https://github.com/jvfe",
    ),
    repository="https://github.com/jvfe/highcomplex_latch",
    license="MIT",
)

metadata.parameters = {
    "read1": LatchParameter(
        display_name="Read 1",
        description="Paired-end read 1 file to be assembled.",
        batch_table_column=True,  # Show this parameter in batched mode.
        section_title="Data",
    ),
    "read2": LatchParameter(
        display_name="Read 2",
        description="Paired-end read 2 file to be assembled.",
        batch_table_column=True,  # Show this parameter in batched mode.
    ),
    "sample_name": LatchParameter(
        display_name="Sample name",
        description="Sample name (will define output file names)",
    ),
    "contaminants": LatchParameter(
        display_name="Contaminants", description="FASTA file with contaminant sequences"
    ),
}
