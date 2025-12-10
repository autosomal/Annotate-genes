#!/usr/bin/env python3
"""
实际执行BED文件基因注释和富集分析
"""

import pandas as pd
from pathlib import Path

# 基因位置数据库
GENE_POSITIONS = {
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

# 基因功能分类
GENE_FUNCTIONS = {
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

# 您的BED文件内容
bed_content = """chr17\t36143270\t36143271
chr4\t156041918\t156041919
chr8\t41004246\t41004247
chr18\t31929457\t31929458
chr17\t71965860\t71965861
chr10\t67103337\t67103338
chr10\t80123580\t80123581
chr6\t29374583\t29374584
chr19\t5477150\t5477151
chr8\t25512704\t25512705"""

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
    for gene_name, gene_info in GENE_POSITIONS.items():
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
    if gene in GENE_FUNCTIONS:
        function = GENE_FUNCTIONS[gene]
        if function not in functional_groups:
            functional_groups[function] = []
        functional_groups[function].append(gene)

# 创建输出目录
Path('/workspace/data').mkdir(exist_ok=True)

# 保存注释结果
bed_df.to_csv('/workspace/data/your_annotated_genes.bed', sep='\t', index=False)

# 生成详细报告
report_content = f"""BED文件基因注释和富集分析报告
{'=' * 60}

1. 输入数据概况
{'-' * 30}
基因组区域数量: {len(bed_df)}
染色体分布: {bed_df['chromosome'].value_counts().to_dict()}

2. 基因注释结果
{'-' * 30}
"""

for idx, row in bed_df.iterrows():
    report_content += f"区域 {idx+1}: {row['chromosome']}:{row['start']:,}-{row['end']:,} -> {row['genes']}\n"

report_content += f"""
3. 识别基因列表
{'-' * 30}
总基因数: {len(unique_genes)}
基因列表: {', '.join(unique_genes) if unique_genes else '无匹配基因'}

4. 功能富集分析
{'-' * 30}
"""

if functional_groups:
    for function, genes in functional_groups.items():
        report_content += f"{function}: {', '.join(genes)} ({len(genes)}个基因)\n"
else:
    report_content += "未识别到已知功能基因\n"

report_content += f"""
5. 分析建议
{'-' * 30}
• 当前示例使用的基因数据库较小，可能无法覆盖所有区域
• 建议使用更大的注释数据库，如GENCODE、Ensembl等
• 可以考虑使用在线工具如UCSC Genome Browser进行验证
• 对于未知区域，可能需要进一步的功能预测分析

6. 后续分析建议
{'-' * 30}
• 使用gseapy进行GO/KEGG富集分析
• 集成在线API进行更全面的注释
• 添加基因表达数据进行分析
• 进行通路分析和网络分析
"""

with open('/workspace/data/your_enrichment_analysis.txt', 'w', encoding='utf-8') as f:
    f.write(report_content)

# 生成总结
print("BED文件基因注释和富集分析完成!")
print(f"分析了 {len(bed_df)} 个基因组区域")
print(f"识别出 {len(unique_genes)} 个基因")

if unique_genes:
    print("\n识别基因:")
    for gene in unique_genes:
        print(f"  - {gene}")
    
    print("\n功能分组:")
    for function, genes in functional_groups.items():
        print(f"  - {function}: {', '.join(genes)}")
else:
    print("\n未识别到匹配基因，可能原因:")
    print("  - 基因组区域不在已知基因范围内")
    print("  - 基因数据库覆盖度不足")
    print("  - 染色体命名格式不匹配")

print("\n输出文件:")
print("- /workspace/data/your_annotated_genes.bed - 注释结果")
print("- /workspace/data/your_enrichment_analysis.txt - 详细分析报告")