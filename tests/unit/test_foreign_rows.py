"""The guard that keeps re-selection stages from dropping appended genomes."""

from __future__ import annotations

import logging

import pytest

from repgenr.core.context import WorkdirContext
from repgenr.core.errors import UserInputError
from repgenr.core.integrity import refuse_foreign_rows
from repgenr.core.manifest import GenomeRecord

_LOG = logging.getLogger("test")


def test_guard_names_the_rows_and_the_override(workdir) -> None:
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert_many(
        [GenomeRecord("GCF_1", "a.fasta", "gtdb"), GenomeRecord("SRR1", "b.fasta", "sra")]
    )
    with pytest.raises(UserInputError, match="SRR1.*--drop-foreign|--drop-foreign.*SRR1"):
        refuse_foreign_rows(ctx, "metadata", drop_foreign=False, logger=_LOG)
    refuse_foreign_rows(ctx, "metadata", drop_foreign=True, logger=_LOG)  # a warning, no error


def test_guard_is_quiet_without_foreign_rows(workdir) -> None:
    ctx = WorkdirContext(workdir, create=True)
    ctx.manifest.upsert_many([GenomeRecord("GCF_1", "a.fasta", "gtdb")])
    refuse_foreign_rows(ctx, "ingest", drop_foreign=False, logger=_LOG)
