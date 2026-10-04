# 中文注释：STRING PPI 网络获取 + 符号映射。
# 输入构念：零输入源节点（Egress 下载）。规则来源：方案 §18.2——
# 网络邻近用 STRING 作底图。网络类型二选一（v3 修正，审计 F11）：
#   - functional（默认）：protein.links.v<V>，功能关联底图（含间接证据）；
#   - physical：protein.physical.links.v<V>，仅物理互作底图。
#     默认保持 functional 是为了不静默变更既有下游结果；方案要物理底图
#     必须显式 network_type=physical（README/manifest 同口径）。
# 分数口径（v3 修正）：combined_score 阈值按 STRING 官方档位——
# 150=low / 400=medium / 700=high / 900=highest（0-1000 整数尺度，
# 即 web 界面 0.150/0.400/0.700/0.900 档）。来源：string-db.org
# scores 页（"medium or better (score >= 0.400)"）+ web 界面档位下拉。
# 旧注释"150 = medium"是错的（150 实为 low 档下限）。
# 失败处理：参数格式错/下载失败 → 退出 2。解释边界：well-studied 基因
# 度偏倚是 STRING 固有混杂；解读在网络邻近层，本节点不做注释。
import gzip
import hashlib
import os
import re
import sys
import time
import urllib.request

BASE = "https://stringdb-downloads.org/download"

# STRING 官方置信档位（0-1000 整数尺度的下限分）；排序即从低到高。
SCORE_TIERS = (
    (150, "low"),
    (400, "medium"),
    (700, "high"),
    (900, "highest"),
)
SCORE_TIER_SOURCE = (
    "STRING official scores page (medium or better = score >= 0.400) "
    "+ web UI minimum-score presets 0.150/0.400/0.700/0.900"
)

NETWORK_TYPES = ("functional", "physical")


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


def links_stem(network_type):
    """functional → protein.links；physical → protein.physical.links。"""
    if network_type == "functional":
        return "protein.links"
    if network_type == "physical":
        return "protein.physical.links"
    fail(f"network_type must be functional|physical, got {network_type}")


def download_urls(base, organism, version, network_type):
    """返回 (links_url, info_url)。info 与网络类型无关，永远同版。"""
    stem = links_stem(network_type)
    links = f"{base}/{stem}.v{version}/{organism}.{stem}.v{version}.txt.gz"
    info = f"{base}/protein.info.v{version}/{organism}.protein.info.v{version}.txt.gz"
    return links, info


def score_label(score):
    """分数落进哪一档（档位下限 ≤ score 的最高档）；<150 无官方档位。"""
    label = "below_low"
    for threshold, name in SCORE_TIERS:
        if score >= threshold:
            label = name
    return label


def threshold_note(min_score):
    """min_score 的日志注记：官方档位名 + 全档位映射与来源。"""
    tiers = "/".join(f"{t}={n}" for t, n in SCORE_TIERS)
    return (
        f"min_score={min_score} (official tier: {score_label(min_score)}; "
        f"STRING v-scale tiers {tiers}; source: {SCORE_TIER_SOURCE})"
    )


def fetch(url, destination):
    """3 次指数退避下载，边下边算 sha256；非 2xx 视为失败。"""
    last = None
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "autonomics-refdata/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, open(destination, "wb") as handle:
                digest = hashlib.sha256()
                size = 0
                while True:
                    block = response.read(1 << 20)
                    if not block:
                        break
                    handle.write(block)
                    digest.update(block)
                    size += len(block)
            return size, digest.hexdigest()
        except Exception as error:  # noqa: BLE001 — 重试一切网络错误
            last = error
            time.sleep(2 ** attempt)
    fail(f"download failed after 3 attempts: {url}: {last}")


def main():
    organism = os.environ.get("STRING_ORGANISM", "9606")
    version = os.environ.get("STRING_VERSION", "12.0")
    min_score = int(os.environ.get("STRING_MIN_SCORE", "150"))
    network_type = os.environ.get("STRING_NETWORK_TYPE", "functional")

    if not re.fullmatch(r"[0-9A-Za-z.]+", organism):
        fail(f"organism must be alphanumeric: {organism}")
    if not re.fullmatch(r"[0-9A-Za-z.]+", version):
        fail(f"version must be alphanumeric: {version}")
    if not 0 <= min_score <= 1000:
        fail(f"min_score must be in [0,1000]: {min_score}")
    if network_type not in NETWORK_TYPES:
        fail(f"network_type must be {'|'.join(NETWORK_TYPES)}, got {network_type}")

    out_edges = os.environ["AUTONOMICS_OUTPUT0"]
    out_raw = os.environ["AUTONOMICS_OUTPUT1"]
    out_info = os.environ["AUTONOMICS_OUTPUT2"]
    out_log = os.environ["AUTONOMICS_OUTPUT3"]
    workdir = os.environ.get("AUTONOMICS_WORKDIR", "/work")

    log_lines = [
        "params\torganism={}\trelease={}\tnetwork_type={}\t{}".format(
            organism, version, network_type, threshold_note(min_score)
        )
    ]
    raw = {}

    links_url, info_url = download_urls(BASE, organism, version, network_type)
    links_path = os.path.join(workdir, "links.txt.gz")
    info_path = os.path.join(workdir, "info.txt.gz")

    for url, path in ((links_url, links_path), (info_url, info_path)):
        size, digest = fetch(url, path)
        raw[url] = (size, digest)
        log_lines.append(f"fetch\turl={url}\tbytes={size}\tsha256={digest}")

    # protein.info：#string_protein_id \t preferred_name \t protein_size \t annotation
    symbol = {}
    with gzip.open(info_path, "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if "#string_protein_id" not in header:
            fail("protein.info header missing #string_protein_id column")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1]:
                symbol[fields[0]] = fields[1]
    with open(out_info, "w", encoding="utf-8") as handle:
        handle.write("string_protein_id\tpreferred_name\n")
        for key in sorted(symbol):
            handle.write(f"{key}\t{symbol[key]}\n")
    log_lines.append(f"gene_info\tproteins={len(symbol)}")

    # protein.links（或 protein.physical.links）：空格分隔，末列 combined_score。
    n_total = n_kept = n_unmapped = 0
    with gzip.open(links_path, "rt", encoding="utf-8") as source, \
            open(out_raw, "w", encoding="utf-8") as raw_out, \
            open(out_edges, "w", encoding="utf-8") as edge_out:
        raw_out.write("protein_a\tprotein_b\tcombined_score\n")
        edge_out.write("gene_a\tgene_b\n")
        for line in source:
            fields = line.split()
            if len(fields) < 3 or fields[0] == "protein1":
                continue
            n_total += 1
            score = int(fields[-1])
            raw_out.write(f"{fields[0]}\t{fields[1]}\t{score}\n")
            if score < min_score:
                continue
            gene_a = symbol.get(fields[0])
            gene_b = symbol.get(fields[1])
            if gene_a is None or gene_b is None:
                n_unmapped += 1
                continue
            n_kept += 1
            edge_out.write(f"{gene_a}\t{gene_b}\n")

    log_lines.append(
        f"edges\tnetwork_type={network_type}\trelease=v{version}\t"
        f"total={n_total}\tkept={n_kept}\t{threshold_note(min_score)}\tunmapped={n_unmapped}"
    )
    log_lines.append("rule: edge table is symbol-mapped; degree bias caveat at interpretation layer")
    with open(out_log, "w", encoding="utf-8") as handle:
        handle.write("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
