import pathlib
import shutil

import cogent3
import pytest
from click.testing import CliRunner

from ensembl_tui import _cli_option as cli_opt
from ensembl_tui import _config as eti_config
from ensembl_tui import cli as eti_cli

RUNNER = CliRunner()


@pytest.mark.slow
@pytest.mark.internet
@pytest.mark.timeout(500)
def test_download(tmp_config_just_yeast):
    """runs download, install, drop according to a special test cfg"""
    tmp_dir = tmp_config_just_yeast.parent
    # now download

    r = RUNNER.invoke(
        eti_cli.download, [f"-c{tmp_config_just_yeast}"], catch_exceptions=False
    )
    assert r.exit_code == 0, r.output
    genome_dir = tmp_dir / "staging" / "genomes"
    dirnames = [dn.name for dn in genome_dir.iterdir() if dn.is_dir()]
    assert "saccharomyces_cerevisiae" in dirnames
    # make sure file sizes > 0
    paths = list((genome_dir / "saccharomyces_cerevisiae").glob("*"))
    size = sum(p.stat().st_size for p in paths)
    assert size > 0
    assert r.exit_code == 0, r.output


def test_download_no_config():
    r = RUNNER.invoke(eti_cli.download, ["-d"], catch_exceptions=False)
    assert r.exit_code != 0, r.output
    assert "No config" in r.output


@pytest.mark.internet
def test_demo_config(tmp_dir):
    """demo_config works correctly"""
    outdir = tmp_dir / "exported"
    r = RUNNER.invoke(eti_cli.demo_config, [f"-o{outdir}"])
    assert r.exit_code == 0, r.output
    fnames = {f.name for f in outdir.iterdir()}
    assert "species.tsv" in fnames
    assert len(fnames) == 3
    shutil.rmtree(tmp_dir)


@pytest.mark.internet
def test_demo_config_exists(tmp_dir):
    outdir = tmp_dir / "exported"
    outdir.mkdir(parents=True, exist_ok=True)
    r = RUNNER.invoke(eti_cli.demo_config, [f"-o{outdir}"])
    assert r.exit_code == 1
    r = RUNNER.invoke(eti_cli.demo_config, [f"-o{outdir}", "--force_overwrite"])
    assert r.exit_code == 0


@pytest.fixture(scope="module")
def installed(tmp_downloaded):
    # tmp_downloaded is a temp copy of the download folder
    # we add the verbose and force_overwrite flags to exercise
    # those conditional statements
    r = RUNNER.invoke(
        eti_cli.install,
        [f"-d{tmp_downloaded}", "-v", "-f"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    return tmp_downloaded.parent / "install"


def test_do_install(tmp_config_no_compara):
    r = RUNNER.invoke(
        eti_cli.install,
        [f"-d{tmp_config_no_compara}", "-v", "-f"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    install_dir = tmp_config_no_compara.parent / "install"
    assert install_dir.exists()
    config = eti_config.read_installed_cfg(install_dir)
    assert not config.homologies_path.exists()
    r = RUNNER.invoke(eti_cli.installed, [f"-i{install_dir}"], catch_exceptions=False)
    assert r.exit_code == 0, r.output
    assert "Ensembl release:" in r.output
    assert "Installed genomes" in r.output
    assert "saccharomyces_cerevisiae" in r.output


def test_installed(small_install_path):
    config = eti_config.read_installed_cfg(small_install_path)
    assert config.homologies_path.exists()
    assert sum(f.stat().st_size for f in config.homologies_path.iterdir()) > 8_000
    r = RUNNER.invoke(
        eti_cli.installed, [f"-i{small_install_path}"], catch_exceptions=False
    )
    assert r.exit_code == 0, r.output
    assert "Ensembl release:" in r.output
    assert "Installed genomes" in r.output
    assert "caenorhabditis_elegans" in r.output
    path = config.installed_genome("caenorhabditis_elegans")
    # should be 2 combined attr parquet files
    assert len(list(path.glob("*attr.parquet"))) == 2
    table = config.get_version_table()
    assert table.title in r.output


def test_installed_with_alignments(apes_install_path):
    r = RUNNER.invoke(
        eti_cli.installed,
        [f"-i{apes_install_path}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert "Installed alignments: ✅" in r.output


def test_installed_full_path(apes_install_path):
    r = RUNNER.invoke(
        eti_cli.installed,
        [f"-i{apes_install_path / eti_config.INSTALLED_CONFIG_NAME}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert "Installed alignments: ✅" in r.output


def test_installed_invalid_path():
    r = RUNNER.invoke(
        eti_cli.installed,
        ["-i", "not/a/valid/path"],
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


def test_check_one_cds_seq(small_install_path):
    # checking a single exon sequence with a rel_start > 0
    from ensembl_tui import _genome as eti_genome

    config = eti_config.read_installed_cfg(small_install_path)
    genome = eti_genome.load_genome(
        config=config,
        species="saccharomyces_cerevisiae",
    )
    cds = next(iter(genome.get_features(name="YMR242C", biotype="cds", canonical=True)))
    seq = cds.get_slice()
    expect = (
        "GCTCACTTTAAAGAATACCAAGTTATTGGCCGTCGTTTGCCAACTGAATCTGTTCCAGAA"
        "CCAAAGTTGTTCAGAATGAGAATCTTTGCTTCAAATGAAGTTATTGCCAAGTCTCGTTAC"
        "TGGTATTTCTTGCAAAAGTTGCACAAGGTTAAGAAGGCTTCTGGTGAAATTGTTTCCATC"
        "AACCAAATCAACGAAGCTCATCCAACCAAGGTCAAGAACTTCGGTGTCTGGGTTAGATAC"
        "GACTCCAGATCTGGTACTCACAATATGTACAAGGAAATCAGAGACGTCTCCAGAGTTGCT"
        "GCCGTCGAAACCTTATACCAAGACATGGCTGCCAGACACAGAGCTAGATTTAGATCTATT"
        "CACATCTTGAAGGTTGCTGAAATTGAAAAGACTGCTGACGTCAAGAGACAATACGTTAAG"
        "CAATTTTTGACCAAGGACTTGAAATTCCCATTGCCTCACAGAGTCCAAAAATCCACCAAG"
        "ACTTTCTCCTACAAGAGACCTTCCACTTTCTACTGA"
    )
    assert str(seq) == expect


def test_check_multi_exon_cds_seq_plus_strand(small_install_path):
    # checking a multi exon sequence with a rel_start > 0
    # and rel_end != exon length
    from ensembl_tui import _genome as eti_genome

    config = eti_config.read_installed_cfg(small_install_path)
    genome = eti_genome.load_genome(
        config=config,
        species="caenorhabditis_elegans",
    )
    cds = next(
        iter(genome.get_features(name="WBGene00185002", biotype="cds", canonical=True)),
    )
    aa = str(cds.get_slice().get_translation())
    # seq expected values from ensembl
    assert aa.startswith("MEMEDIDDDITVFYTDDRGTVQGPYGASTVLDWYQKGYFSDNHQMRFTDNGQRIGNLFTY")
    assert aa.endswith("IEKVKTNCRDAPSPLPPAMDPVAPYHVRDKCTQS")
    assert len(aa) == 274


def test_check_two_exon_cds_seq_rev_strand(small_install_path):
    # checking a two exon sequence with a rel_start > 0
    # and rel_end != exon length
    from ensembl_tui import _genome as eti_genome

    config = eti_config.read_installed_cfg(small_install_path)
    genome = eti_genome.load_genome(
        config=config,
        species="caenorhabditis_elegans",
    )
    cds = next(
        iter(genome.get_features(name="WBGene00184990", biotype="cds", canonical=True)),
    )
    aa = str(cds.get_slice().get_translation())
    # seq expected values from ensembl
    assert aa.startswith("MSGVYNNSGSRMRSKNFEKHQVPSDMAFFQKFRKQSHSNETVDCKKKQEE")
    assert aa.endswith("DGHYSDETVEEKHNREHRNKTKADNRTRRIAEIRRKHNINA")
    assert len(aa) == 161


def test_species_summary(small_install_path):
    r = RUNNER.invoke(
        eti_cli.species_summary,
        [f"-i{small_install_path}", "--species", "caenorhabditis_elegans"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert "Caenorhabditis elegans" in r.output
    assert "protein_coding" in r.output


@pytest.fixture(params=[1, 2, 3])
def bad_species(request, tmp_path):
    if request.param == 1:
        yield "caenorhabditis_elegans,saccharomyces_cerevisiae"
    if request.param == 2:
        outpath = tmp_path / "species.tsv"
        outpath.write_text("\n")
        yield str(outpath)
    if request.param == 3:
        yield "not-a-species"


def test_species_summary_invalid_species(apes_install_path, bad_species):
    r = RUNNER.invoke(
        eti_cli.species_summary,
        [f"-i{apes_install_path}", "--species", bad_species],
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


def test_dump_genes(small_install_path):
    species = "caenorhabditis_elegans"
    outdir = small_install_path.parent
    limit = 10
    args = [
        f"-i{small_install_path}",
        "--species",
        species,
        "--outdir",
        str(outdir),
        "--limit",
        f"{limit}",
    ]
    r = RUNNER.invoke(
        eti_cli.dump_genes,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    tsv_path = next(iter(outdir.glob("*.tsv")))
    assert tsv_path.name.startswith(species)
    table = cogent3.load_table(tsv_path)
    assert table.shape[0] == limit
    gene_biotype = str(table.columns["biotype"][0])
    transcript_biotype = set(table.columns["transcript_biotypes"][0].split(","))
    assert gene_biotype in transcript_biotype


def test_dump_genes_error(apes_install_path, bad_species):
    outdir = apes_install_path.parent
    args = [
        f"-i{apes_install_path}",
        "--species",
        bad_species,
        "--outdir",
        str(outdir),
    ]
    r = RUNNER.invoke(
        eti_cli.dump_genes,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


def test_homologs(small_install_path, tmp_dir):
    outdir = tmp_dir / "output"
    limit = 10
    args = [
        f"-i{small_install_path}",
        "--ref",
        "caenorhabditis_elegans",
        "--outdir",
        f"{outdir}",
        "--limit",
        str(limit),
        "-ht",
        "ortholog_one2one",
        "-v",
    ]

    r = RUNNER.invoke(
        eti_cli.homologs,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert len(dstore.completed) == limit


def test_homologs_error_no_ref(apes_install_path, tmp_dir):
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "-ht",
        "ortholog_one2one",
        "-v",
    ]

    r = RUNNER.invoke(
        eti_cli.homologs,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


def test_homologs_error_refgenes_coords(
    apes_install_path,
    tmp_dir,
    ref_genes,
):
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--ref",
        "homo_sapiens",
        "-ht",
        "ortholog_one2one",
        "--ref_genes",
        f"{ref_genes}",
        "--coord_names",
        "22",
    ]

    r = RUNNER.invoke(
        eti_cli.homologs,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


def test_homologs_coord_name(small_install_path, tmp_dir):
    outdir = tmp_dir / "output"
    limit = 10
    args = [
        f"-i{small_install_path}",
        "--ref",
        "saccharomyces_cerevisiae",
        "--outdir",
        f"{outdir}",
        "--limit",
        str(limit),
        "--coord_names",
        "I,XVI,II",
        "-ht",
        "ortholog_one2one",
        "-v",
    ]

    r = RUNNER.invoke(
        eti_cli.homologs,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert len(dstore.completed) == limit


def test_compara_summary(small_install_path):
    r = RUNNER.invoke(
        eti_cli.compara_summary,
        [f"-i{small_install_path}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert "homology_type" in r.output
    assert "ortholog_one2many" in r.output
    assert "alignments" not in r.output.lower()


def test_compara_summary_apes(apes_install_path):
    r = RUNNER.invoke(
        eti_cli.compara_summary,
        [f"-i{apes_install_path}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert "homology_type" in r.output
    assert "ortholog_one2many" in r.output
    assert "10_primates" in r.output


def test_compara_folder_not_created(tmp_config_no_compara):
    # ensure compara folder is not created if not specified in the config
    r = RUNNER.invoke(
        eti_cli.install,
        [f"-d{tmp_config_no_compara}"],
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    assert not (tmp_config_no_compara / "compara").exists()


def test_genome_coords_from_tsv(tmp_dir):
    species_tsv = tmp_dir / "genome_coords.tsv"
    with open(species_tsv, "w") as f:
        f.write("species\tseqid\tstart\tstop\tstrand\nhomo_sapiens\t1\t3000\t4000\t1\n")
    coords = cli_opt.genome_coords_from_tsv(None, None, species_tsv)
    assert len(coords) == 1
    got = coords[0]
    assert got.species == "homo_sapiens"
    assert got.seqid == "1"
    assert got.start == 3000
    assert got.stop == 4000
    assert got.strand == "1"


def test_genome_coords_from_tsv_noheader(tmp_dir, capsys):
    invalid = tmp_dir / "invalid.tsv"
    with open(invalid, "w") as f:
        f.write("homo_sapiens\t1\t3000\t4000\t1\n")
    with pytest.raises(SystemExit) as excinfo:
        cli_opt.genome_coords_from_tsv(None, None, invalid)

    captured = capsys.readouterr()
    assert "ERROR: failed to load file" in captured.out
    assert excinfo.value.code == 1


def test_genome_coords_from_tsv_missing_value(tmp_dir, capsys):
    invalid = tmp_dir / "invalid.tsv"
    with open(invalid, "w") as f:
        f.write("species\tseqid\tstart\tstop\tstrand\nhomo_sapiens\t1\t3000\t\t1\n")
    with pytest.raises(SystemExit) as excinfo:
        cli_opt.genome_coords_from_tsv(None, None, invalid)

    captured = capsys.readouterr()
    assert "ERROR: all values of 'stop' must be integers" in captured.out
    assert excinfo.value.code == 1


@pytest.fixture(params=[True, False])
def coord_name(request, tmp_dir):
    if request.param:
        outpath = pathlib.Path(tmp_dir) / "coord_names.tsv"
        outpath.write_text("22\n")
        return str(outpath)
    return "22"


def test_alignments_coord_names(apes_install_path, tmp_dir, coord_name):
    # coord_names as argument or file 22
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--coord_names",
        coord_name,
        "--limit",
        "2",
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert len(dstore.completed)
    logged = dstore.logs[0].read()
    assert isinstance(logged, str)
    # check the version of one dependency is present
    assert "version : cogent3_h5seqs" in logged


@pytest.fixture
def human_genes() -> list[str]:
    # aligned as plus, minus strand
    return ["ENSG00000100346", "ENSG00000100412"]


@pytest.fixture
def human_cds(apes, human_genes) -> list[str]:
    hsap = apes["homo_sapiens"]
    return [
        next(iter(hsap.get_features(seqid="22", biotype="gene", name=stable_id)))
        for stable_id in human_genes
    ]


@pytest.fixture
def ref_genes(human_genes, tmp_path) -> str:
    genes = cogent3.make_table(data={"stableid": human_genes})
    outpath = tmp_path / "ref_genes.tsv"
    genes.write(outpath)
    return str(outpath)


@pytest.fixture
def bad_ref_genes(human_genes, tmp_path) -> str:
    genes = cogent3.make_table(data={"name": human_genes})
    outpath = tmp_path / "ref_genes.tsv"
    genes.write(outpath)
    return str(outpath)


def _check_alignments(
    dstore,
    masked: bool = False,
    just_ref: bool = False,
    shadow: bool = False,
) -> bool:
    # checking sequences or whether they have masked, or not,
    # in just reference, or not
    loader = cogent3.get_app("load_aligned", moltype="dna")
    alns = [loader(m) for m in dstore.completed]
    if not masked:
        for aln in alns:
            all_seqs = "".join(aln.to_dict().values())
            if "?" in all_seqs:
                return False
    elif masked:
        for aln in alns:
            for s in aln.seqs:
                raw = str(s).replace("-", "")
                if "homo" in s.name:
                    if (
                        "?" not in raw
                        or raw.startswith("?")
                        or (shadow and not raw.startswith("?"))
                    ):
                        return False
                elif just_ref and "?" in str(s):
                    return False

    return True


def test_alignments_ref_genes(apes_install_path, tmp_dir, ref_genes, human_genes):
    # coord_names as argument or file 22
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_genes",
        ref_genes,
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert len(dstore.completed) >= 2
    record_ids = [m.unique_id.split(".")[0].split("-")[0] for m in dstore.completed]
    assert set(record_ids) == set(human_genes)
    assert _check_alignments(dstore, masked=False)


def test_alignments_bad_ref_genes(apes_install_path, tmp_dir, bad_ref_genes):
    # coord_names as argument or file 22
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_genes",
        bad_ref_genes,
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


# ref coords
@pytest.fixture
def ref_coords(human_cds, tmp_path) -> str:
    data = {"species": [], "seqid": [], "start": [], "stop": [], "strand": []}
    for cds in human_cds:
        data["species"].append("homo_sapiens")
        data["seqid"].append(cds.seqid)
        data["start"].append(cds.map.start)
        data["stop"].append(cds.map.end)
        data["strand"].append(-1 if cds.reversed else 1)
    table = cogent3.make_table(data=data)
    outpath = tmp_path / "ref_coords.tsv"
    table.write(outpath)
    return str(outpath)


def test_alignments_ref_coords(apes_install_path, tmp_dir, ref_coords):
    # coord_names as argument or file 22
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_coords",
        ref_coords,
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert len(dstore.completed) >= 2
    assert all(m.unique_id.startswith("homo_sapiens-22") for m in dstore.completed)


@pytest.fixture(params=[1, 2, 3, 4, 5, 6])
def bad_coords(ref_coords, request):
    path = pathlib.Path(ref_coords)
    new_out = path.parent / "bad_coords.txt"
    if request.param == 1:
        yield str(new_out)
    table = cogent3.load_table(ref_coords)
    if request.param == 2:
        table.write(new_out, sep=";")
        yield str(new_out)
    if request.param == 3:
        data = table.columns.to_dict()
        data.pop("strand")
        new_tab = cogent3.make_table(data=data)
        new_tab.write(new_out, sep="\t")
        yield str(new_out)
    if request.param == 4:
        data = table.columns.to_dict()
        data["a-strand"] = data.pop("strand")
        new_tab = cogent3.make_table(data=data)
        new_tab.write(new_out, sep="\t")
        yield str(new_out)
    if request.param == 5:
        header = "\t".join(table.header)
        new_out.write_text(f"{header}\n")
        yield str(new_out)
    if request.param == 6:
        data = table.columns.to_dict()
        data["start"] = list("a" * table.shape[0])
        new_tab = cogent3.make_table(data=data)
        new_tab.write(new_out, sep="\t")
        yield str(new_out)


def test_alignments_ref_coords_error(apes_install_path, tmp_dir, bad_coords):
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_coords",
        bad_coords,
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code != 0, r.output


# ref_genes, mask
@pytest.mark.parametrize("just_ref", [True, False])
def test_alignments_mask(apes_install_path, tmp_dir, ref_genes, just_ref):
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_genes",
        ref_genes,
        "--mask",
        "Simple_repeat",
    ]
    args += ["--mask_ref"] if just_ref else args

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert _check_alignments(dstore, masked=True, just_ref=just_ref)


# mask shadow
@pytest.mark.parametrize("shadow", [True, False][1:])
def test_alignments_mask_shadow(apes_install_path, tmp_dir, ref_genes, shadow):
    # coord_names as argument or file 22
    outdir = tmp_dir / "output"
    args = [
        f"-i{apes_install_path}",
        "--outdir",
        f"{outdir}",
        "--align_name",
        "*primate*",
        "--ref",
        "Human",
        "--ref_genes",
        ref_genes,
        "--mask_ref",
        "--mask_shadow" if shadow else "--mask",
        "Simple_repeat",
    ]

    r = RUNNER.invoke(
        eti_cli.alignments,
        args,
        catch_exceptions=False,
    )
    assert r.exit_code == 0, r.output
    dstore = cogent3.open_data_store(outdir, suffix="fa", mode="r")
    assert _check_alignments(dstore, masked=True, just_ref=True, shadow=shadow)
