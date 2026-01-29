from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from lasso_workbench.altframe.gene_extraction import extract_gene_instances


def test_extract_gene_instances_from_gbk(tmp_path):
    seq = Seq("ATG" + "GCT" * 10 + "TAA")
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

    instances = extract_gene_instances(tmp_path, {"testGene"}, "gene")
    assert "testGene" in instances
    assert len(instances["testGene"]) == 1
    inst = instances["testGene"][0]
    assert inst.gene_start == 0
    assert inst.gene_end == len(seq)
    assert inst.gene_strand == "+"
    assert inst.gene_genomic == str(seq)


def test_extract_gene_instances_with_locus_tag(tmp_path):
    seq = Seq("ATG" + "GCT" * 4 + "TAA")
    record = SeqRecord(seq, id="TEST_TAG", name="TEST_TAG", description="")
    record.annotations["molecule_type"] = "DNA"
    feature = SeqFeature(
        FeatureLocation(0, len(seq), strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["LT0001"]},
    )
    record.features = [feature]

    gbk_path = tmp_path / "tag.gbk"
    SeqIO.write(record, gbk_path, "genbank")

    instances = extract_gene_instances(tmp_path, {"LT0001"}, "locus_tag")
    assert "LT0001" in instances
    assert len(instances["LT0001"]) == 1


def test_extract_gene_instances_respects_codon_start(tmp_path):
    seq = Seq("G" + "ATG" + "GCT" + "TAA")
    record = SeqRecord(seq, id="TEST_CODON", name="TEST_CODON", description="")
    record.annotations["molecule_type"] = "DNA"
    feature = SeqFeature(
        FeatureLocation(0, len(seq), strand=1),
        type="CDS",
        qualifiers={"gene": ["codonGene"], "codon_start": ["2"]},
    )
    record.features = [feature]

    gbk_path = tmp_path / "codon_start.gbk"
    SeqIO.write(record, gbk_path, "genbank")

    instances = extract_gene_instances(tmp_path, {"codonGene"}, "gene")
    inst = instances["codonGene"][0]
    assert inst.codons[:3] == ["ATG", "GCT", "TAA"]
def test_extract_gene_instances_minus_strand(tmp_path):
    seq = Seq("ATGAAACCCGGG")
    record = SeqRecord(seq, id="TEST_MINUS", name="TEST_MINUS", description="")
    record.annotations["molecule_type"] = "DNA"
    feature = SeqFeature(
        FeatureLocation(0, len(seq), strand=-1),
        type="CDS",
        qualifiers={"gene": ["minusGene"]},
    )
    record.features = [feature]

    gbk_path = tmp_path / "minus.gbk"
    SeqIO.write(record, gbk_path, "genbank")

    instances = extract_gene_instances(tmp_path, {"minusGene"}, "gene")
    assert "minusGene" in instances
    assert len(instances["minusGene"]) == 1
    inst = instances["minusGene"][0]
    assert inst.gene_strand == "-"
    assert inst.gene_genomic == str(seq)
    assert inst.codons[:4] == ["CCC", "GGG", "TTT", "CAT"]
