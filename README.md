# refdata

公共参考数据获取家族。补齐 2026-10-04 节点审计问题 3 的参考数据缺口：
geo 家族只盖 GEO/SRA 元数据，STRING PPI 网络与 MSigDB 基因集库
（引擎 bundle 表无 MSigDB）此前都得宿主侧备好后手递。LINCS/CMap 的
GCT/GCTx 矩阵走 geo 家族的 `geo_suppl` 节点（GSE92742 suppl 目录）。

## 节点

| kind | 输入 | 输出 | 说明 |
|---|---|---|---|
| `string_ppi` | 无（Egress 源节点） | `string_edges.tsv`、`string_edges_raw.tsv`、`string_gene_info.tsv`、`string_ppi.log` | STRING protein.links + protein.info，映射基因符号并按 combined_score 过滤 |
| `msigdb_gmt` | 无（Egress 源节点） | `msigdb_sets.tsv`、`msigdb_set_info.tsv`、`msigdb_gmt.log` | MSigDB collection GMT 下载并摊平成长表 |

参数：`string_ppi` 的 `organism`（默认 9606）/`version`（默认 12.0）/
`min_score`（默认 150 = STRING medium，v12 尺度 0-1000）；
`msigdb_gmt` 的 `collection`（必填，如 `h.all`、`c5.go_bp`）/`release`
（默认 2026.1.Hs）。

## 布局

```text
refdata/
├── manifest.toml
├── scripts/string_ppi.py
├── scripts/msigdb_gmt.py
├── Dockerfile          # 镜像 provenance（共享 biotools-py，见 _images/）
└── README.md
```

镜像：`localhost/autonomics/biotools-py@sha256:21e14c1582a9d4c261c0a23a11864293313ee2b49b9b7d85bcabc89210b48c7e`
（构建树在上级 `_images/biotools-py/`，本目录 Dockerfile 为副本）。

## 口径

- 下载 3 次指数退避，流式 sha256 边下边算；URL/字节数/摘要全部入 log
  自锚（§26.1），可直接进运行账本。
- STRING 边表两份：raw（STRING ID + combined_score）与符号映射版
  （network_proximity 的直接输入）；未映射计数如实记录。
- MSigDB release 钉死进 log；版本间基因集内容漂移，复现以 release 为准。
- 实测 URL 形态（2026-10-04 验证 200）：
  `https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz`
  （83MB）、`https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt`。

## 解释边界

well-studied 基因度偏倚是 STRING 固有混杂，解读在网络邻近层，本家族
不做注释；MSigDB 非商业学术许可（Broad），商用需另议。

## 冒烟记录

见各 log 锚（部署后在 /artifacts/string_ppi、/artifacts/msigdb_gmt）。
