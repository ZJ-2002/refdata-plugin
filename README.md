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
`network_type`（v3，审计 F11：`functional` 默认 = protein.links 功能关联
底图；`physical` = protein.physical.links 仅物理互作底图。默认 functional
保持既有行为不静默变更——**方案要物理底图必须显式 `network_type=physical`**）/
`min_score`（默认 150；官方档位 150=low / 400=medium / 700=high /
900=highest，0-1000 尺度——v3 修正，旧文档"150 = medium"是错的，默认
阈值实为 low 档下限）；`msigdb_gmt` 的 `collection`（必填，如 `h.all`、
`c5.go_bp`）/`release`（默认 2026.1.Hs）。

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
- 网络类型（v3，审计 F11）：functional 与 physical 是**不同的底图**——
  functional 含间接/共表达等功能关联证据，physical 只含物理互作。
  网络药理学语境（方案 §18.2）下游 network_proximity 吃哪种底图由方案
  决定；本节点不替方案选，默认 functional 仅是兼容既有运行。
- 分数阈值口径（v3 修正）：combined_score 官方档位 150=low /
  400=medium / 700=high / 900=highest（0-1000 整数尺度 = web 界面
  0.150/0.400/0.700/0.900 档；来源 string-db.org scores 页 "medium or
  better (score >= 0.400)"）。log 记录 release、网络类型、分数阈值、
  官方档位映射与来源。
- STRING 边表两份：raw（STRING ID + combined_score）与符号映射版
  （network_proximity 的直接输入）；未映射计数如实记录。
- MSigDB release 钉死进 log；版本间基因集内容漂移，复现以 release 为准。
- 实测 URL 形态（2026-10-04 验证 200）：
  `https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz`
  （83MB；physical 型为
  `protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz`）、
  `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt`。

## 解释边界

well-studied 基因度偏倚是 STRING 固有混杂，解读在网络邻近层，本家族
不做注释；MSigDB 非商业学术许可（Broad），商用需另议。

## 测试

`tests/test_string_ppi.py`（离线纯单测，不打网络）：URL 构造覆盖
functional/physical 两型与 info 不随类型变；network_type 非法值 →
退出 2；score_label 档位边界（150/400/700/900/0/999）；threshold_note
含官方档位映射与来源；`import string_ppi` 无副作用（main 守卫）。

## 冒烟记录

见各 log 锚（部署后在 /artifacts/string_ppi、/artifacts/msigdb_gmt）。
