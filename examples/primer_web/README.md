# primer_web — 浏览器端引物设计器

primer3 引擎驱动：MoonBit 计算（`src/primer3.mbt` + `src/primer3_report.mbt`），HTML 渲染候选表与引物对产物视图，**模板序列不出本机**。

## 使用

```bash
cd examples/primer_web
# 1. 生成报告(演示模板)
moon run . -- --demo > report.json 2>/dev/null

# 或对自定义模板
moon run . -- "$(cat template.txt)" > report.json

# 2. 本地起服务并打开
python3 -m http.server 8000
# 浏览器访问 http://localhost:8000/examples/primer_web/index.html
```

## 功能

| 视图 | 说明 |
|---|---|
| 推荐引物对 | 最优 F/R 组合 + 扩增子示意条（产物长度可视化） |
| 引物对组合 | 全部组合按组合得分排序（含产物 bp） |
| 正/反向候选 | 序列 / 起点 / 长度 / Tm / GC / 自互补 / 得分，**表头可点击排序** |

质量提示：Tm 52–68°C 绿色、GC 40–60% 绿色、自互补 0 为最佳（标签颜色分级）。

## 计算侧说明

- 复用库层 `primer3_design_primers`（Wallace Tm / GC / 自互补 / 发夹 / 交叉二聚体打分）；
- 引物对组合得分 = 两引物得分之和 − 交叉二聚体罚分（0.5×），产物长度 0–5000bp；
- 输出 JSON 亦可独立嵌入流水线。
