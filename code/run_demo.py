#!/usr/bin/env python3
"""
直接执行演示
"""

import pandas as pd
from pathlib import Path

# 基因位置数据
gene_positions = {
    'TP53': {'chromosome': 'chr17', 'start': 7661779, 'end': 7687550},
    'BRCA1': {'chromosome': 'chr17', 'start': 43044295, 'end': 43125483},
    'BRCA2': {'chromosome': 'chr13', 'start': 32889665, 'end': 32973808},
    'EGFR': {'chromosome': 'chr7', 'start': 55086724, 'end': 55275031},
    'KRAS': {'chromosome': 'chr12', 'start': 25205249, 'end': 25250929},
    'PIK3CA': {'chromosome': 'chr3', 'start': 178866882, 'end': 178952497},
    'AKT1': {'chromosome': 'chr14', 'start': 105235686, 'end': 105262119},
    'MYC': {'chromosome': 'chr8', 'start': 128748315, 'end': 128753680},
    'RB1': {'chromosome': 'chr13', 'start': 48877837, 'end': 49056148},
    'PTEN': {'chromosome': 'chr10', 'start': 87863113, 'end': 87971930}
}

# BED文件内容
bed_content = """chr17\t7661779\t7687550
chr13\t32889665\t32973808
chr7\t55086724\t55275031
chr12\t25205249\t25250929
chr10\t87863113\t87971930
chr17\t43044295\t43125483
chr8\t128748315\t128753680
chr14\t105235686\t105262119
chr3\t178866882\t178952497
chr13\t48877837\t49056148"""

# 解析BED文件
lines = [line.strip() for line in bed_content.strip().split('\n') if line.strip()]
data = []

for line in lines:
    fields = line.split('\t')
    if len(fields) >= 3:
        data.append({
            'chromosome': fields[0],
            'start': int(fields[1]),
            'end': int(fields[2])
        })

bed_df = pd.DataFrame(data)

# 基因注释
def find_overlapping_genes(chromosome, start, end):
    overlapping_genes = []
    for gene_name, gene_info in gene_positions.items():
        if gene_info['chromosome'] == chromosome:
            if not (end < gene_info['start'] or start > gene_info['end']):
                overlapping_genes.append(gene_name)
    return overlapping_genes

genes_list = []
for idx, row in bed_df.iterrows():
    chromosome = str(row['chromosome'])
    start = int(row['start'])
    end = int(row['end'])
    
    overlapping_genes = find_overlapping_genes(chromosome, start, end)
    genes_list.append(';'.join(overlapping_genes) if overlapping_genes else 'Unknown')

bed_df['genes'] = genes_list

# 基因功能分类
gene_functions = {
    'TP53': '肿瘤抑制基因',
    'BRCA1': 'DNA修复基因',
    'BRCA2': 'DNA修复基因',
    'EGFR': '生长因子受体',
    'KRAS': '原癌基因',
    'PIK3CA': 'PI3K通路',
    'AKT1': 'PI3K/AKT通路',
    'MYC': '转录因子',
    'RB1': '肿瘤抑制基因',
    'PTEN': '肿瘤抑制基因'
}

# 提取基因列表
all_genes = []
for genes_str in bed_df['genes']:
    if genes_str != 'Unknown':
        genes = genes_str.split(';')
        all_genes.extend(genes)

unique_genes = list(set(all_genes))

# 功能富集分析
functional_groups = {}
for gene in unique_genes:
    if gene in gene_functions:
        function = gene_functions[gene]
        if function not in functional_groups:
            functional_groups[function] = []
        functional_groups[function].append(gene)

# 创建输出目录
Path('/workspace/data').mkdir(exist_ok=True)

# 保存结果
bed_df.to_csv('/workspace/data/demo_annotated_genes.bed', sep='\t', index=False)

# 生成报告
report_content = f"""基因富集分析结果
{'=' * 50}

分析的基因总数: {len(unique_genes)}
基因列表: {', '.join(unique_genes)}

功能分组:
{'-' * 30}
"""

for function, genes in functional_groups.items():
    report_content += f"{function}: {', '.join(genes)}\n"

with open('/workspace/data/demo_enrichment_analysis.txt', 'w', encoding='utf-8') as f:
    f.write(report_content)

print("BED文件基因注释和富集分析完成!")
print(f"注释了 {len(bed_df)} 个基因组区域")
print(f"识别出 {len(unique_genes)} 个基因")
print(f"功能分组数量: {len(functional_groups)}")

print("\n注释结果:")
print(bed_df.to_string(index=False))

print("\n功能分组:")
for function, genes in functional_groups.items():
    print(f"  {function}: {', '.join(genes)}")

print("\n结果文件已保存:")
print("- /workspace/data/demo_annotated_genes.bed")
print("- /workspace/data/demo_enrichment_analysis.txt")