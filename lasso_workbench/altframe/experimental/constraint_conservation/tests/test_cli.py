import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from lasso_workbench.altframe.experimental.constraint_conservation import cli


def test_cli_smoke(tmp_path, monkeypatch):
    seq = Seq("ATG" + "GCT" * 15 + "TAA")
    record = SeqRecord(seq, id="TEST", name="TEST", description="")
    record.annotations["molecule_type"] = "DNA"
    feature = SeqFeature(
        FeatureLocation(0, len(seq), strand=1),
        type="CDS",
        qualifiers={"gene": ["testGene"]},
    )
    record.features = [feature]

    gbk_path = tmp_path / "test.gbk"
    SeqIO.write(record, gbk_path, "genbank")

    hits_path = tmp_path / "hits.tsv"
    hits_path.write_text(
        "\t".join(
            [
                "gene_name",
                "match_type",
                "orf_strand",
                "orf_frame",
                "gene_start",
                "gene_end",
                "orf_start",
                "orf_end",
                "gbk_file",
                "record_id",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "testGene",
                "out_of_frame",
                "+",
                "1",
                "0",
                str(len(seq)),
                "0",
                str(len(seq)),
                gbk_path.name,
                "TEST",
            ]
        )
        + "\n"
    )

    out_dir = tmp_path / "out"
    argv = [
        "prog",
        "--gbk-dir",
        str(tmp_path),
        "--hits-file",
        str(hits_path),
        "--output-dir",
        str(out_dir),
        "--min-genomes",
        "1",
        "--max-candidates",
        "5",
        "--null-iterations",
        "3",
        "--geom-bins",
        "10",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    result = cli.main()
    assert result == 0
    out_path = out_dir / "altframe_constraint_conservation.tsv"
    assert out_path.exists()
    contents = out_path.read_text().strip().splitlines()
    assert len(contents) >= 2
    header = contents[0].split("\t")
    assert "obs_survival" in header
    assert "null_survival_mean" in header
