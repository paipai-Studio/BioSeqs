# fq_qc_web — 浏览器端 FASTQ 质控仪表盘

FastQC 风格的纯前端 FASTQ 质控报告：MoonBit 计算（`src/fqqc.mbt`），HTML/SVG 渲染，**测序数据不出本机**。

## 使用

```bash
# 1. 生成报告(演示数据, 内含质量衰减/adapter 污染/GC 偏置特征)
cd examples/fq_qc_web
moon run . -- --demo > report.json 2>/dev/null

# 或分析真实数据
moon run . -- "$(cat sample.fq)" > report.json

# 2. 本地起服务并打开(file:// 协议下 fetch 受限)
python3 -m http.server 8000
# 浏览器访问 http://localhost:8000/examples/fq_qc_web/index.html
```

## 报告模块（对齐 FastQC）

| 模块 | 说明 |
|---|---|
| Basic Statistics | reads / 碱基 / 均长 / 平均 GC |
| Per base quality | 逐位质量箱线（低/四分/中位/四分/高），绿黄红区带 |
| Per sequence quality | 序列平均质量直方图 |
| Per sequence GC | 序列 GC 分布（混样/污染时出现双峰） |
| Per base GC | 逐位 GC 折线 |
| Length distribution | 读长分布 |
| Adapter content | Illumina adapter / polyG 前缀命中 |
| Overrepresented sequences | 高频重复序列 Top 20 |

## 计算侧说明

- 流式逐条处理（`fastq_each_record`），内存 O(单条 read)；
- 逐位质量按 0-50 直方图聚合（精确五数，内存与 reads 数无关）；
- 输出 JSON 亦可独立嵌入任意前端或流水线。
