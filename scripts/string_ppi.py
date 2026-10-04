# 中文注释：STRING PPI 网络获取 + 符号映射。
# 输入构念：零输入源节点（Egress 下载）。规则来源：方案 §18.2——
# 网络邻近用 STRING 作底图；边表口径 combined_score ≥ min_score
# （150 = STRING medium，v12 尺度 0-1000）。
# 失败处理：参数格式错/下载失败 → 退出 2。解释边界：well-studied 基因
# 度偏倚是 STRING 固有混杂；解读在网络邻近层，本节点不做注释。
import gzip
import hashlib
import os
import re
import sys
import time
import urllib.request

ORGANISM = os.environ.get("STRING_ORGANISM", "9606")
VERSION = os.environ.get("STRING_VERSION", "12.0")
MIN_SCORE = int(os.environ.get("STRING_MIN_SCORE", "150"))
WORKDIR = os.environ.get("AUTONOMICS_WORKDIR", "/work")
OUT_EDGES = os.environ["AUTONOMICS_OUTPUT0"]
OUT_RAW = os.environ["AUTONOMICS_OUTPUT1"]
OUT_INFO = os.environ["AUTONOMICS_OUTPUT2"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT3"]

BASE = "https://stringdb-downloads.org/download"


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


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


if not re.fullmatch(r"[0-9A-Za-z.]+", ORGANISM):
    fail(f"organism must be alphanumeric: {ORGANISM}")
if not re.fullmatch(r"[0-9A-Za-z.]+", VERSION):
    fail(f"version must be alphanumeric: {VERSION}")
if not 0 <= MIN_SCORE <= 1000:
    fail(f"min_score must be in [0,1000]: {MIN_SCORE}")

log_lines = []
raw = {}

links_url = f"{BASE}/protein.links.v{VERSION}/{ORGANISM}.protein.links.v{VERSION}.txt.gz"
info_url = f"{BASE}/protein.info.v{VERSION}/{ORGANISM}.protein.info.v{VERSION}.txt.gz"
links_path = os.path.join(WORKDIR, "links.txt.gz")
info_path = os.path.join(WORKDIR, "info.txt.gz")

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
with open(OUT_INFO, "w", encoding="utf-8") as handle:
    handle.write("string_protein_id\tpreferred_name\n")
    for key in sorted(symbol):
        handle.write(f"{key}\t{symbol[key]}\n")
log_lines.append(f"gene_info\tproteins={len(symbol)}")

# protein.links：空格分隔，末列 combined_score。
n_total = n_kept = n_unmapped = 0
with gzip.open(links_path, "rt", encoding="utf-8") as source, \
        open(OUT_RAW, "w", encoding="utf-8") as raw_out, \
        open(OUT_EDGES, "w", encoding="utf-8") as edge_out:
    raw_out.write("protein_a\tprotein_b\tcombined_score\n")
    edge_out.write("gene_a\tgene_b\n")
    for line in source:
        fields = line.split()
        if len(fields) < 3 or fields[0] == "protein1":
            continue
        n_total += 1
        score = int(fields[-1])
        raw_out.write(f"{fields[0]}\t{fields[1]}\t{score}\n")
        if score < MIN_SCORE:
            continue
        gene_a = symbol.get(fields[0])
        gene_b = symbol.get(fields[1])
        if gene_a is None or gene_b is None:
            n_unmapped += 1
            continue
        n_kept += 1
        edge_out.write(f"{gene_a}\t{gene_b}\n")

log_lines.append(
    f"edges\ttotal={n_total}\tkept={n_kept}\tmin_score={MIN_SCORE}\tunmapped={n_unmapped}"
)
log_lines.append("rule: edge table is symbol-mapped; degree bias caveat at interpretation layer")
with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write("\n".join(log_lines) + "\n")
