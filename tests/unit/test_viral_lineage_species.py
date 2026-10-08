"""NCBI Virus species from the lineage binomial, and normalised segment labels."""

from __future__ import annotations

import json
import logging
import zlib
from pathlib import Path

import pytest

from repgenr.core.errors import WorkdirError
from repgenr.viral.ncbi_virus import (
    UNLABELLED_SEGMENTS,
    normalise_segment,
    parse_report,
    read_records,
    species_from_lineage,
    write_records,
)

_ARENA = ["Viruses", "Riboviria", "Arenaviridae", "Mammarenavirus"]
_HANTA = ["Viruses", "Hantaviridae", "Orthohantavirus"]


def _line(acc: str, organism: str, lineage: list[str] | None, segment=None, length=3400):
    # A taxid per organism name, as in NCBI (records of one taxid share a species).
    virus: dict = {"organism_name": organism, "tax_id": zlib.crc32(organism.encode())}
    if lineage is not None:
        virus["lineage"] = [{"name": n, "tax_id": i} for i, n in enumerate(lineage)]
    row = {"accession": acc, "length": length, "completeness": "COMPLETE", "virus": virus}
    if segment is not None:
        row["segment"] = segment
    return json.dumps(row)


# --- species --------------------------------------------------------------


def test_species_is_the_binomial_in_the_lineage() -> None:
    lineage = [*_ARENA, "Mammarenavirus brazilense", "Sabia virus"]
    assert species_from_lineage(lineage, ["Mammarenavirus"]) == "Mammarenavirus brazilense"
    hyphen = [*_ARENA, "Mammarenavirus dhati-welelense", "Dhati Welel virus"]
    assert species_from_lineage(hyphen, ["Mammarenavirus"]) == "Mammarenavirus dhati-welelense"


@pytest.mark.parametrize(
    ("below", "genus"),
    [
        # Live (review of #227): older names below the binomial that start
        # with the genus but are not binomials.
        (["Hepatovirus ahepa", "Hepatovirus A"], "Hepatovirus"),
        (["Rotavirus alphagastroenteritidis", "Rotavirus A"], "Rotavirus"),
        (["Cypovirus altineae", "Cypovirus 1"], "Cypovirus"),
        (["Enterovirus alphacoxsackie", "Enterovirus A"], "Enterovirus"),
        (["Orthohantavirus hantanense", "Orthohantavirus sp. 'x'"], "Orthohantavirus"),
        # An earlier binomial below the current one: the shallowest is the species.
        (["Orthobunyavirus cacheense", "Orthobunyavirus maguariense"], "Orthobunyavirus"),
    ],
)
def test_names_below_the_binomial_are_not_the_species(below, genus) -> None:
    lineage = ["Viruses", genus, *below, "Some virus strain 1"]
    assert species_from_lineage(lineage, [genus]) == below[0]


def test_no_binomial_gives_none() -> None:
    lineage = [*_HANTA, "unclassified Orthohantavirus", "unclassified Hantavirus", "Academ virus"]
    assert species_from_lineage(lineage, ["Orthohantavirus"]) is None
    assert species_from_lineage([], []) is None
    assert species_from_lineage([*_HANTA, "Orthohantavirus sp."], ["Orthohantavirus"]) is None


def test_subgenus_does_not_hide_the_binomial() -> None:
    """Live (Embecovirus, 889 records): the last single-word '-virus' name is
    the subgenus, so 'Betacoronavirus gravedinis' never matched it and every
    record kept its organism name. The genus is the binomial's first word."""
    lineage = [
        "Viruses",
        "Coronaviridae",
        "Orthocoronavirinae",
        "Betacoronavirus",
        "Embecovirus",
        "Betacoronavirus gravedinis",
        "Betacoronavirus 1",
        "Bovine coronavirus",
    ]
    (rec,) = parse_report([_line("U00735.2", "Bovine coronavirus", lineage)])
    assert rec.species == "Betacoronavirus-gravedinis"
    assert rec.genus == "Betacoronavirus" and rec.family == "Coronaviridae"


def test_one_taxid_takes_one_species(caplog) -> None:
    """Live: Mudanjiang phlebovirus (taxid 2983973) under 'Phlebovirus
    baishanense' and under 'unclassified Phlebovirus'; Hepatovirus A (12092)
    under 'Hepatovirus ahepa' and under 'Hepatovirus fejalco' > 'ahepa'."""

    def line(acc, taxid, lineage, organism):
        row = json.loads(_line(acc, organism, lineage))
        row["virus"]["tax_id"] = taxid
        return json.dumps(row)

    phlebo = ["Viruses", "Phenuiviridae", "Phlebovirus"]
    hav = ["Viruses", "Picornaviridae", "Hepatovirus"]
    lines = [
        line(
            "P1",
            2983973,
            [*phlebo, "Phlebovirus baishanense", "Mudanjiang phlebovirus"],
            "Mudanjiang phlebovirus",
        ),
        line(
            "P2",
            2983973,
            [*phlebo, "unclassified Phlebovirus", "Mudanjiang phlebovirus"],
            "Mudanjiang phlebovirus",
        ),
        line("H1", 12092, [*hav, "Hepatovirus ahepa", "Hepatovirus A"], "Hepatovirus A"),
        line("H2", 12092, [*hav, "Hepatovirus ahepa", "Hepatovirus A"], "Hepatovirus A"),
        line(
            "H3",
            12092,
            [*hav, "Hepatovirus fejalco", "Hepatovirus ahepa", "Hepatovirus A"],
            "Hepatovirus A",
        ),
    ]
    with caplog.at_level(logging.INFO, logger="t"):
        recs = parse_report(lines, logging.getLogger("t"))
    by_acc = {r.accession: r.species for r in recs}
    assert by_acc["P1"] == by_acc["P2"] == "Phlebovirus-baishanense"
    assert by_acc["H1"] == by_acc["H2"] == by_acc["H3"] == "Hepatovirus-ahepa"
    assert "12092 -> Hepatovirus ahepa" in caplog.text


def test_records_without_a_binomial_are_logged(caplog) -> None:
    lines = [
        _line("D.1", "Academ virus", [*_HANTA, "unclassified Orthohantavirus", "Academ virus"])
    ]
    with caplog.at_level(logging.INFO, logger="t"):
        (rec,) = parse_report(lines, logging.getLogger("t"))
    assert rec.species == "Academ-virus" and rec.genus == "Orthohantavirus"
    assert "1 of 1 record(s)" in caplog.text and "organism name" in caplog.text


def test_parse_report_takes_the_species_from_the_lineage() -> None:
    """Live (Orthohantavirus, Mammarenavirus): strain-level and earlier organism
    names gave one species token per name."""
    recs = parse_report(
        [
            _line(
                "A.1",
                "Hantaan virus CGAa1011",
                [*_HANTA, "Orthohantavirus hantanense", "Hantaan virus CGAa1011"],
            ),
            _line(
                "B.1",
                "Argentinian mammarenavirus",
                [*_ARENA, "Mammarenavirus juninense", "Argentinian mammarenavirus"],
            ),
            _line("C.1", "Mammarenavirus juninense", [*_ARENA, "Mammarenavirus juninense"]),
            _line("D.1", "Academ virus", [*_HANTA, "unclassified Orthohantavirus", "Academ virus"]),
            _line("E.1", "Lone virus", None),
        ]
    )
    by_acc = {r.accession: r for r in recs}
    assert by_acc["A.1"].species == "Orthohantavirus-hantanense"
    assert by_acc["A.1"].organism == "Hantaan virus CGAa1011"
    assert by_acc["B.1"].species == by_acc["C.1"].species == "Mammarenavirus-juninense"
    assert by_acc["B.1"].organism == "Argentinian mammarenavirus"
    assert by_acc["D.1"].species == "Academ-virus"
    assert by_acc["E.1"].species == "Lone-virus" and by_acc["E.1"].lineage == []
    assert by_acc["C.1"].lineage[-1] == "Mammarenavirus juninense"


def test_records_round_trip_with_the_lineage(tmp_path: Path) -> None:
    recs = parse_report(
        [_line("A.1", "Sabia virus", [*_ARENA, "Mammarenavirus brazilense", "Sabia virus"])]
    )
    write_records(tmp_path / "r.json", recs)
    assert read_records(tmp_path / "r.json") == recs


def test_records_without_a_lineage_key_are_refused(tmp_path: Path) -> None:
    """Records written before this change carry the organism name as species."""
    legacy = {
        "accession": "A.1",
        "taxid": "1",
        "organism": "Sabia virus",
        "family": "Arenaviridae",
        "genus": "Mammarenavirus",
        "species": "Sabia-virus",
        "length": 7000,
        "completeness": "COMPLETE",
        "segment": "L",
        "isolate": "",
    }
    path = tmp_path / "virus_records.json"
    path.write_text(json.dumps([legacy]))
    with pytest.raises(WorkdirError, match="--force vmetadata"):
        read_records(path)
    # Records with a lineage but no species source (an unreleased
    # intermediate format) are refused the same way.
    path.write_text(json.dumps([{**legacy, "lineage": ["Viruses"]}]))
    with pytest.raises(WorkdirError, match="--force vmetadata"):
        read_records(path)


# --- segments -------------------------------------------------------------


@pytest.mark.parametrize(
    ("label", "expected"),
    [
        ("M", "M"),
        ("M; medium", "M"),
        ("middle", "M"),
        ("medium", "M"),
        ("Small", "S"),
        ("S; small", "S"),
        ("L; large", "L"),
        ("LARGE", "L"),
        ("S RNA", "S"),
        ("RNA L", "L"),
        ("RNA 1", "1"),
        ("RNA1", "1"),
        ("dsRNA 5", "5"),
        ("Genome segment 8", "8"),
        ("DNA-A", "A"),
        ("DNA A", "A"),
        ("DNA_B", "B"),
        ("b", "B"),
        ("component B", "B"),
        ("circular DNA 1", "1"),
        ("circular DNA 2", "2"),
        ("S1", "S1"),
        ("NSP4", "NSP4"),
        ("4", "4"),
        ("RNA", "RNA"),
        ("", ""),
        (None, ""),
    ],
)
def test_normalise_segment(label, expected) -> None:
    assert normalise_segment(label) == expected


def test_distinct_numbered_segments_stay_distinct() -> None:
    """Taking the first word alone would make 'RNA 1' and 'RNA 2' one segment."""
    assert len({normalise_segment(x) for x in ("RNA 1", "RNA 2", "RNA 3")}) == 3
    assert normalise_segment("DNA-A") != normalise_segment("DNA-B")


def test_unlabelled_segments() -> None:
    for label in ("ANONYMOUS", "Unknown", ""):
        assert normalise_segment(label) in UNLABELLED_SEGMENTS
    assert normalise_segment("M") not in UNLABELLED_SEGMENTS


def test_parse_report_keeps_the_raw_segment_label() -> None:
    (rec,) = parse_report([_line("A.1", "x", _ARENA, segment="M; medium")])
    assert rec.segment == "M; medium"


def test_species_resolution_logs_an_organism_name(caplog) -> None:
    from repgenr.viral.selection import _resolve_species_targets

    recs = parse_report(
        [
            _line(
                "B.1",
                "Argentinian mammarenavirus",
                [*_ARENA, "Mammarenavirus juninense", "Argentinian mammarenavirus"],
            ),
            _line("C.1", "Mammarenavirus juninense", [*_ARENA, "Mammarenavirus juninense"]),
        ]
    )
    log = logging.getLogger("t")
    with caplog.at_level(logging.INFO, logger="t"):
        out = _resolve_species_targets(recs, {"species": ["Argentinian mammarenavirus"]}, log)
    assert out == {"species": ["Argentinian mammarenavirus", "Mammarenavirus-juninense"]}
    assert "organism name" in caplog.text
    # A binomial passes unchanged; an unknown value is left for the no-match error.
    same = {"species": ["Mammarenavirus juninense", "nothing"]}
    assert _resolve_species_targets(recs, same, log) == same
    assert _resolve_species_targets(recs, {"genus": ["x"]}, log) == {"genus": ["x"]}


def test_species_resolution_takes_the_union() -> None:
    """A value that is a species token (a record without a binomial) and the
    organism name of records under a binomial selects both."""
    from repgenr.viral.selection import _record_matches, _resolve_species_targets

    hav = ["Viruses", "Picornaviridae", "Hepatovirus"]
    recs = parse_report(
        [
            _line("H1", "Hepatovirus A", [*hav, "Hepatovirus ahepa", "Hepatovirus A"]),
            _line("H2", "Hepatovirus A", [*hav, "unclassified Hepatovirus", "Hepatovirus A"]),
        ]
    )
    # One organism name, so one taxid and one species; a second taxid
    # without a binomial keeps the organism name as species.
    assert {r.species for r in recs} == {"Hepatovirus-ahepa"}
    recs[1].species = "Hepatovirus-A"
    out = _resolve_species_targets(recs, {"species": ["Hepatovirus A"]}, logging.getLogger("t"))
    assert sorted(out["species"]) == ["Hepatovirus A", "Hepatovirus-ahepa"]
    assert all(_record_matches(r, out) for r in recs)


# --- NCBI Taxonomy lookup -------------------------------------------------------


def _taxonomy_line(taxid: int, family: str, genus: str, species: str | None) -> str:
    classification = {
        "family": {"id": 1, "name": family},
        "genus": {"id": 2, "name": genus},
    }
    if species:
        classification["species"] = {"id": 3, "name": species}
    return json.dumps(
        {
            "query": [str(taxid)],
            "taxonomy": {"tax_id": taxid, "classification": classification},
        }
    )


_BUNYA = ["Viruses", "Peribunyaviridae", "Orthobunyavirus"]


def _maguari_lines() -> list[str]:
    """Live (Orthobunyavirus): the report lineage nests Maguari and Playas
    under 'Orthobunyavirus cacheense' > 'Orthobunyavirus maguariense', while
    NCBI Taxonomy has 'Orthobunyavirus maguariense' as their species."""

    def line(acc, taxid, organism, lineage):
        row = json.loads(_line(acc, organism, lineage))
        row["virus"]["tax_id"] = taxid
        return json.dumps(row)

    nested = ["Orthobunyavirus cacheense", "Orthobunyavirus maguariense"]
    return [
        line("MAG.1", 11575, "Maguari virus", [*_BUNYA, *nested, "Maguari virus"]),
        line("PLA.1", 273344, "Playas virus", [*_BUNYA, *nested, "Playas virus"]),
        line("CVV.1", 35304, "Cache Valley virus", [*_BUNYA, nested[0], "Cache Valley virus"]),
    ]


def _fake_taxonomy_runner(lines: list[str], calls: list[list[str]] | None = None):
    def runner(caps, cmd, *, stdout_path, **kw):
        if calls is not None:
            calls.append([str(c) for c in cmd])
        Path(stdout_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
        return 0

    return runner


def test_lookup_taxonomy_reads_the_classification(tmp_path: Path) -> None:
    from repgenr.viral.ncbi_virus import TaxonClass, lookup_taxonomy

    calls: list[list[str]] = []
    out = lookup_taxonomy(
        ["11575", "273344", "11575", "999"],
        tmp_path,
        logger=logging.getLogger("t"),
        runner=_fake_taxonomy_runner(
            [
                _taxonomy_line(
                    11575, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
                ),
                _taxonomy_line(
                    273344, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
                ),
                _taxonomy_line(999, "Peribunyaviridae", "Orthobunyavirus", None),
            ],
            calls,
        ),
    )
    assert out["11575"] == TaxonClass(
        "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
    )
    assert out["999"].species == ""
    # One call for every distinct taxid, read from an input file; no files left.
    assert len(calls) == 1 and "--inputfile" in calls[0]
    assert not list(tmp_path.iterdir())


def test_lookup_taxonomy_failure_returns_nothing(tmp_path: Path, caplog) -> None:
    from repgenr.core.errors import ToolExecutionError
    from repgenr.viral.ncbi_virus import lookup_taxonomy

    def failing(caps, cmd, **kw):
        raise ToolExecutionError(["datasets"], 1)

    with caplog.at_level(logging.WARNING, logger="t"):
        out = lookup_taxonomy(["1"], tmp_path, logger=logging.getLogger("t"), runner=failing)
    assert out == {} and "lookup" in caplog.text and "lineage" in caplog.text


def test_taxonomy_species_is_preferred_over_the_lineage(tmp_path: Path) -> None:
    from repgenr.viral.ncbi_virus import lookup_taxonomy

    lines = _maguari_lines()
    taxonomy = lookup_taxonomy(
        ["11575", "273344"],
        tmp_path,
        logger=logging.getLogger("t"),
        runner=_fake_taxonomy_runner(
            [
                _taxonomy_line(
                    11575, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
                ),
                _taxonomy_line(
                    273344, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
                ),
            ]
        ),
    )
    recs = {r.accession: r for r in parse_report(lines, None, taxonomy)}
    assert recs["MAG.1"].species == recs["PLA.1"].species == "Orthobunyavirus-maguariense"
    assert recs["MAG.1"].species_source == "taxonomy"
    # Not looked up: the shallowest lineage binomial is the offline fallback.
    assert recs["CVV.1"].species == "Orthobunyavirus-cacheense"
    assert recs["CVV.1"].species_source == "lineage"
    # Without any lookup the lineage rule applies to all.
    offline = {r.accession: r for r in parse_report(lines)}
    assert offline["MAG.1"].species == "Orthobunyavirus-cacheense"


def test_a_species_absent_from_the_lineage_comes_from_the_taxonomy() -> None:
    """Live: Murutucu virus's current species, 'Orthobunyavirus maritubaense',
    is absent from its report lineage."""
    from repgenr.viral.ncbi_virus import TaxonClass

    (line,) = [
        _line(
            "MUR.1", "Murutucu virus", [*_BUNYA, "Orthobunyavirus caraparuense", "Murutucu virus"]
        )
    ]
    taxid = str(json.loads(line)["virus"]["tax_id"])
    taxonomy = {
        taxid: TaxonClass("Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maritubaense")
    }
    (rec,) = parse_report([line], None, taxonomy)
    assert rec.species == "Orthobunyavirus-maritubaense"
    assert rec.organism == "Murutucu virus"


def test_organism_at_genus_level_is_not_given_the_family_as_genus() -> None:
    """Live (Phlebovirus, 6 records): an organism that is the genus itself
    had the family as its genus."""
    (rec,) = parse_report(
        [_line("G.1", "Phlebovirus", ["Viruses", "Phenuiviridae", "Phlebovirus"])]
    )
    assert rec.genus == "Phlebovirus" and rec.family == "Phenuiviridae"
    (rec,) = parse_report([_line("G.2", "Odd agent", ["Viruses", "Phenuiviridae", "Odd agent"])])
    assert rec.genus == "NA" and rec.family == "Phenuiviridae"


def test_fetch_resolves_the_species_through_the_taxonomy(tmp_path: Path) -> None:
    """vmetadata looks the taxids up once after the download and stores the
    resolved species, so vgenome needs no network."""
    import zipfile

    from repgenr.viral import ncbi_virus

    lines = _maguari_lines()

    def download(caps, cmd, **kw):
        zip_path = Path(cmd[cmd.index("--filename") + 1])
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/data_report.jsonl", "\n".join(lines) + "\n")
            zf.writestr(
                "ncbi_dataset/data/genomic.fna",
                "".join(f">{a} x\nACGT\n" for a in ("MAG.1", "PLA.1", "CVV.1")),
            )
        return 0

    calls: list[list[str]] = []
    taxonomy = _fake_taxonomy_runner(
        [
            _taxonomy_line(
                11575, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
            ),
            _taxonomy_line(
                273344, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus maguariense"
            ),
            _taxonomy_line(
                35304, "Peribunyaviridae", "Orthobunyavirus", "Orthobunyavirus cacheense"
            ),
        ],
        calls,
    )
    recs = ncbi_virus.fetch(
        "Orthobunyavirus",
        tmp_path,
        logger=logging.getLogger("t"),
        runner=download,
        taxonomy_runner=taxonomy,
    )
    assert {r.accession: r.species for r in recs} == {
        "MAG.1": "Orthobunyavirus-maguariense",
        "PLA.1": "Orthobunyavirus-maguariense",
        "CVV.1": "Orthobunyavirus-cacheense",
    }
    assert len(calls) == 1
    write_records(tmp_path / "virus_records.json", recs)
    assert {r.species_source for r in read_records(tmp_path / "virus_records.json")} == {"taxonomy"}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["download.fa", "virus_records.json"]
