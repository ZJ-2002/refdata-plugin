# 离线纯单测（审计 F11）：URL 构造 / network_type 校验 / 官方阈值档位。
# 不打网络：fetch 与主流程不在此覆盖（Egress 行为由部署侧 log 锚自证）。
import importlib.util
import sys

import pytest

SCRIPT = "scripts/string_ppi.py"
BASE = "https://stringdb-downloads.org/download"


def load_module():
    spec = importlib.util.spec_from_file_location("string_ppi", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["string_ppi"] = module  # 供 monkeypatch.setattr 按名替换 fetch
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return load_module()


def test_import_is_side_effect_free():
    # main 守卫：import 不得触发 env 读取/退出（引擎以 python3 <file> 运行，
    # __name__ == "__main__" 成立；见 container_command 的 script 插入契约）。
    load_module()


def test_functional_urls_match_deployed_shape(mod):
    links, info = mod.download_urls(BASE, "9606", "12.0", "functional")
    assert links == (
        "https://stringdb-downloads.org/download/protein.links.v12.0/"
        "9606.protein.links.v12.0.txt.gz"
    )
    assert info == (
        "https://stringdb-downloads.org/download/protein.info.v12.0/"
        "9606.protein.info.v12.0.txt.gz"
    )


def test_physical_urls_switch_links_file_only(mod):
    links, info = mod.download_urls(BASE, "10090", "12.0", "physical")
    assert links == (
        "https://stringdb-downloads.org/download/protein.physical.links.v12.0/"
        "10090.protein.physical.links.v12.0.txt.gz"
    )
    # protein.info 与网络类型无关：物理型也用同一份 ID→符号映射。
    assert info == (
        "https://stringdb-downloads.org/download/protein.info.v12.0/"
        "10090.protein.info.v12.0.txt.gz"
    )


def test_default_network_type_is_functional(monkeypatch):
    # 默认 functional：不传 STRING_NETWORK_TYPE 时不静默切换底图。
    import os

    monkeypatch.delenv("STRING_NETWORK_TYPE", raising=False)
    assert os.environ.get("STRING_NETWORK_TYPE", "functional") == "functional"


def test_invalid_network_type_fails_with_exit_2(mod, monkeypatch, capsys):
    monkeypatch.setenv("STRING_ORGANISM", "9606")
    monkeypatch.setenv("STRING_VERSION", "12.0")
    monkeypatch.setenv("STRING_NETWORK_TYPE", "physicalish")
    monkeypatch.setenv("STRING_MIN_SCORE", "150")
    with pytest.raises(SystemExit) as exc:
        mod.main()
    assert exc.value.code == 2
    assert "network_type" in capsys.readouterr().err


def test_links_stem(mod):
    assert mod.links_stem("functional") == "protein.links"
    assert mod.links_stem("physical") == "protein.physical.links"


def test_score_tiers_official_thresholds(mod):
    # STRING 官方档位（0-1000 尺度）：150=low / 400=medium / 700=high /
    # 900=highest。旧注释"150 = medium"是错的——此处锁定修正后的口径。
    assert dict(mod.SCORE_TIERS) == {
        150: "low", 400: "medium", 700: "high", 900: "highest",
    }


@pytest.mark.parametrize("score,label", [
    (0, "below_low"), (149, "below_low"), (150, "low"), (399, "low"),
    (400, "medium"), (699, "medium"), (700, "high"), (899, "high"),
    (900, "highest"), (999, "highest"), (1000, "highest"),
])
def test_score_label_boundaries(mod, score, label):
    assert mod.score_label(score) == label


def test_threshold_note_names_tier_mapping_and_source(mod):
    note = mod.threshold_note(150)
    assert "min_score=150" in note
    assert "official tier: low" in note  # 150 = low，不是 medium
    assert "150=low" in note and "400=medium" in note
    assert "700=high" in note and "900=highest" in note
    assert "0.400" in mod.SCORE_TIER_SOURCE  # 来源含官方 medium≥0.400 证据


# --- 离线全流程：fetch 打桩喂本地合成 gzip，走完 main() 的解析/过滤/落盘。---

INFO_TSV = (
    "#string_protein_id\tpreferred_name\tprotein_size\tannotation\n"
    "9606.ENSP00000000233\tARHGAP11B\t1037\tRho GTPase\n"
    "9606.ENSP00000000412\tARHGAP11A\t2013\tRho GTPase\n"
    "9606.ENSP00000001031\tFGR\t564\ttyrosine kinase\n"
    "9606.ENSP00000002146\tABCB11\t1321\ttransporter\n"
)
LINKS_TSV = (
    "protein1\tprotein2\tcombined_score\n"
    "9606.ENSP00000000233\t9606.ENSP00000000412\t999\n"   # 高分保留
    "9606.ENSP00000000233\t9606.ENSP00000001031\t400\n"   # 恰过 400（medium 下限）
    "9606.ENSP00000001031\t9606.ENSP00000002146\t150\n"   # 恰过 150（low 下限）
    "9606.ENSP00000000412\t9606.ENSP00000002146\t90\n"    # 低于 150 被滤
    "9606.ENSP00000000233\t9606.ENSP99999999999\t999\n"   # 未映射蛋白被剔除
)


def run_main(monkeypatch, tmp_path, network_type, min_score="150"):
    fake_fetch_setup(monkeypatch, tmp_path)
    outs = [str(tmp_path / name) for name in
            ("edges.tsv", "raw.tsv", "gene_info.tsv", "run.log")]
    monkeypatch.setenv("STRING_ORGANISM", "9606")
    monkeypatch.setenv("STRING_VERSION", "12.0")
    monkeypatch.setenv("STRING_NETWORK_TYPE", network_type)
    monkeypatch.setenv("STRING_MIN_SCORE", min_score)
    monkeypatch.setenv("AUTONOMICS_WORKDIR", str(tmp_path))
    for index, path in enumerate(outs):
        monkeypatch.setenv(f"AUTONOMICS_OUTPUT{index}", path)
    sys.modules["string_ppi"].main()
    return outs


def fake_fetch_setup(monkeypatch, tmp_path):
    import gzip as gzip_mod

    # fixture 放独立子目录：WORKDIR 也是 tmp_path，同名文件会自拷贝截断。
    fixture_dir = tmp_path / "fx"
    fixture_dir.mkdir(exist_ok=True)
    info_gz = fixture_dir / "info.txt.gz"
    links_gz = fixture_dir / "links.txt.gz"
    with open(info_gz, "wb") as handle:
        handle.write(gzip_mod.compress(INFO_TSV.encode()))
    with open(links_gz, "wb") as handle:
        handle.write(gzip_mod.compress(LINKS_TSV.encode()))

    def stub(url, destination):
        source = info_gz if "protein.info" in url else links_gz
        with open(source, "rb") as reader, open(destination, "wb") as writer:
            writer.write(reader.read())
        return 42, "deadbeef"

    load_module()  # 确保 sys.modules 里有 string_ppi
    monkeypatch.setattr(sys.modules["string_ppi"], "fetch", stub)


def test_main_offline_physical_flow_and_log(tmp_path, monkeypatch):
    outs = run_main(monkeypatch, tmp_path, "physical", min_score="400")
    edges, raw, info, log = (open(path, encoding="utf-8").read() for path in outs)

    # min_score=400：仅前两条边保留（150 与 90 被滤），未映射剔除。
    assert edges.splitlines() == [
        "gene_a\tgene_b",
        "ARHGAP11B\tARHGAP11A",
        "ARHGAP11B\tFGR",
    ]
    assert raw.count("\n") == 6  # 表头 + 5 条原始边（raw 不过滤）
    assert "ARHGAP11B" in info and "string_protein_id" in info

    # log 自锚：URL 含 physical 文件名；params 行带 release/网络类型/官方档位。
    assert "protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz" in log
    assert "network_type=physical" in log
    assert "release=12.0" in log
    assert "official tier: medium" in log  # min_score=400 = medium 下限
    assert "edges\tnetwork_type=physical" in log
    assert "kept=2" in log and "unmapped=1" in log
    assert "total=5" in log


def test_main_offline_functional_keeps_existing_url(tmp_path, monkeypatch):
    outs = run_main(monkeypatch, tmp_path, "functional")
    log = open(outs[3], encoding="utf-8").read()
    # 默认路径不静默变更：functional 仍下载 protein.links（非 physical）。
    assert "protein.links.v12.0/9606.protein.links.v12.0.txt.gz" in log
    assert "protein.physical.links" not in log
    assert "official tier: low" in log  # 默认 150 = low 下限（v3 修正口径）
    # min_score=150：三条边保留（150 恰过线），90 滤除。
    edges = open(outs[0], encoding="utf-8").read()
    assert edges.count("\n") == 4
