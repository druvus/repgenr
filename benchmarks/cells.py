"""Benchmark cell matrix for the scaling/bias audit.

A cell describes one measured run. Cells use the stateless data-channel steps
wherever possible (``dereplicate-chunk`` via fofn, ``phylo-build`` via
``--genomes-dir``) so the genome sets are read in place with zero staging;
chunked-dereplication cells copy ``genomes/`` into a real workdir, and dense
sourmash cells run a small in-process driver (there is no CLI switch to
disable the branchwater sparse path).
"""

from __future__ import annotations

from dataclasses import dataclass, field

DEREP_TOOLS = ("skder", "galah", "sourmash")  # sourmash = sparse (branchwater installed)
DEREP_SIZES = (100, 1000, 5000)
SCENARIOS = ("balanced", "clonal")


@dataclass(frozen=True)
class Cell:
    id: str
    kind: str            # derep_step | derep_dense | derep_stage | tree_step
    tool: str
    set_name: str        # directory name under sets/
    subset: int | None = None   # run on the first N genomes of the set
    extra_args: tuple[str, ...] = field(default=())
    timeout_s: int = 14400

    @property
    def n(self) -> int:
        return self.subset if self.subset is not None else int(self.set_name.split("_")[1])


def derep_cells() -> list[Cell]:
    cells = []
    for scenario in SCENARIOS:
        for n in DEREP_SIZES:
            for tool in DEREP_TOOLS:
                cells.append(Cell(
                    id=f"derep-{tool}-{scenario}-{n}",
                    kind="derep_step", tool=tool,
                    set_name=f"{scenario}_{n}_clustered",
                ))
            cells.append(Cell(
                id=f"derep-sourmashdense-{scenario}-{n}",
                kind="derep_dense", tool="sourmash",
                set_name=f"{scenario}_{n}_clustered",
            ))
    return cells


def chunked_cells() -> list[Cell]:
    """Chunked vs single-pass at n=5000 plus the order-sensitivity pair (B2)."""
    return [
        Cell(
            id=f"chunked-{tool}-clonal-5000-{order}",
            kind="derep_stage", tool=tool,
            set_name=f"clonal_5000_{order}",
            extra_args=("--process-size", "1000"),
        )
        for tool in ("skder", "sourmash")
        for order in ("clustered", "random")
    ]


def tree_cells() -> list[Cell]:
    cells = [
        Cell(id=f"tree-mashtree-balanced-{n}", kind="tree_step",
             tool="mashtree", set_name=f"balanced_{n}_clustered")
        for n in (100, 1000, 5000)
    ]
    # sourmash tree builder: pure-Python O(n^3) NJ -- NEVER schedule 5000;
    # 100/500/1000 give the n^3 fit that tests its declared 10000 limit.
    cells += [
        Cell(id=f"tree-sourmashtb-balanced-{n}", kind="tree_step",
             tool="sourmash", set_name="balanced_1000_clustered",
             subset=(None if n == 1000 else n))
        for n in (100, 500, 1000)
    ]
    # ML builders need an MSA: subset via the sibeliaz aligner (capped small).
    cells += [
        Cell(id="tree-fasttree-balanced-100", kind="tree_step", tool="fasttree",
             set_name="balanced_1000_clustered", subset=100,
             extra_args=("--aligner", "sibeliaz"), timeout_s=7200),
        Cell(id="tree-fasttree-balanced-200", kind="tree_step", tool="fasttree",
             set_name="balanced_1000_clustered", subset=200,
             extra_args=("--aligner", "sibeliaz"), timeout_s=7200),
        Cell(id="tree-iqtree-balanced-100", kind="tree_step", tool="iqtree",
             set_name="balanced_1000_clustered", subset=100,
             extra_args=("--aligner", "sibeliaz"), timeout_s=7200),
    ]
    return cells


def all_cells() -> list[Cell]:
    return derep_cells() + chunked_cells() + tree_cells()


def tiers(cells: list[Cell]) -> dict[str, list[Cell]]:
    """smoke (n<=100) -> mid (n<=1000, unchunked) -> heavy (the rest)."""
    out: dict[str, list[Cell]] = {"smoke": [], "mid": [], "heavy": []}
    for cell in cells:
        if cell.n <= 100:
            out["smoke"].append(cell)
        elif cell.n <= 1000 and cell.kind != "derep_stage":
            out["mid"].append(cell)
        else:
            out["heavy"].append(cell)
    return out
