# 中文注释：MSigDB 基因集 GMT 下载 + 摊平。
# 输入构念：零输入源节点（Egress）。规则来源：引擎 bundle 表无 MSigDB，
# 方案的基因集对照（Hallmark/GO/KEGG）由此进 DAG。输出长表
# set_id/gene 直连 gene_set_score / netprox 集合输入。
# 失败处理：release/collection 格式错或下载失败 → 退出 2。
# 解释边界：MSigDB 版本间集合内容漂移，release 钉死在 log 里。
import hashlib
import os
import re
import sys
import time
import urllib.request

COLLECTION = os.environ.get("MSIGDB_COLLECTION", "").strip()
RELEASE = os.environ.get("MSIGDB_RELEASE", "2026.1.Hs").strip()
WORKDIR = os.environ.get("AUTONOMICS_WORKDIR", "/work")
OUT_SETS = os.environ["AUTONOMICS_OUTPUT0"]
OUT_INFO = os.environ["AUTONOMICS_OUTPUT1"]
OUT_LOG = os.environ["AUTONOMICS_OUTPUT2"]


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


if not re.fullmatch(r"[a-z0-9_.]+", COLLECTION or ""):
    fail(f"collection must be like h.all / c5.go_bp / c2.cp.kegg_legacy: {COLLECTION!r}")
if not re.fullmatch(r"[0-9]{4}\.[0-9]+\.(Hs|Mm)", RELEASE):
    fail(f"release must be YYYY.N.Hs or YYYY.N.Mm: {RELEASE!r}")

url = (
    f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    f"{RELEASE}/{COLLECTION}.v{RELEASE}.symbols.gmt"
)
destination = os.path.join(WORKDIR, "collection.gmt")

size = 0
last = None
for attempt in range(3):
    try:
        request = urllib.request.Request(url, headers={"User-Agent": "autonomics-refdata/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response, open(destination, "wb") as handle:
            digest = hashlib.sha256()
            while True:
                block = response.read(1 << 20)
                if not block:
                    break
                handle.write(block)
                digest.update(block)
                size += len(block)
        sha256 = digest.hexdigest()
        break
    except Exception as error:  # noqa: BLE001 — 重试一切网络错误
        last = error
        time.sleep(2 ** attempt)
else:
    fail(f"download failed after 3 attempts: {url}: {last}")

# GMT：set_id \t description \t gene1 \t gene2 ...（纯文本，MSigDB 不压缩）。
n_sets = n_rows = 0
with open(destination, "r", encoding="utf-8") as source, \
        open(OUT_SETS, "w", encoding="utf-8") as sets_out, \
        open(OUT_INFO, "w", encoding="utf-8") as info_out:
    sets_out.write("set_id\tgene\n")
    info_out.write("set_id\tdescription\tn_genes\n")
    for line in source:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        set_id, description, genes = fields[0], fields[1], fields[2:]
        n_sets += 1
        info_out.write(f"{set_id}\t{description}\t{len(genes)}\n")
        for gene in genes:
            if gene:
                sets_out.write(f"{set_id}\t{gene}\n")
                n_rows += 1

with open(OUT_LOG, "w", encoding="utf-8") as handle:
    handle.write(
        f"fetch\turl={url}\tbytes={size}\tsha256={sha256}\n"
        f"sets={n_sets}\trows={n_rows}\trelease={RELEASE}\n"
    )
