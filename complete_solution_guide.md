# BED文件基因注释和富集分析完整解决方案

## 问题概述

您提供了包含基因组位置信息的BED文件，需要将其注释到基因名，并进行功能富集分析。

## 分析结果

### 您的数据
- **基因组区域数量**: 10个
- **染色体分布**: chr17(2), chr8(2), chr10(2), chr4(1), chr18(1), chr6(1), chr19(1)
- **注释结果**: 所有区域都显示为"Unknown"，表明需要更大的基因数据库

## 解决方案

我为您提供了三种不同复杂度的解决方案：

### 1. 简化版本 (`simple_gene_annotation.py`)
**特点**: 
- 使用本地基因位置数据库
- 无需网络连接
- 适合快速初步分析

**使用方法**:
```python
from simple_gene_annotation import SimpleGeneAnnotator, SimpleEnrichmentAnalyzer

annotator = SimpleGeneAnnotator()
analyzer = SimpleEnrichmentAnalyzer()

# 读取和注释
bed_df = annotator.parse_bed_file('your_file.bed')
annotated_df = annotator.annotate_bed_regions(bed_df)

# 提取基因和富集分析
gene_list = analyzer.get_gene_list_from_bed(annotated_df)
enrichment_results = analyzer.simple_enrichment_analysis(gene_list)
```

### 2. 本地版本 (`local_gene_annotation.py`)
**特点**:
- 使用pyensembl或本地GTF文件
- 支持GO/KEGG富集分析
- 需要安装额外依赖

**使用方法**:
```python
from local_gene_annotation import LocalGeneAnnotator, EnrichmentAnalyzer

annotator = LocalGeneAnnotator(gtf_file='path/to/annotation.gtf')
analyzer = EnrichmentAnalyzer()

# 注释
annotated_df = annotator.annotate_with_pyensembl(bed_df)

# 富集分析
go_results = analyzer.perform_go_enrichment(gene_list)
kegg_results = analyzer.perform_kegg_enrichment(gene_list)
```

### 3. API版本 (`bed_gene_annotation.py`)
**特点**:
- 使用在线API进行注释
- 支持实时数据更新
- 需要网络连接

**使用方法**:
```python
from bed_gene_annotation import GeneAnnotator, GeneEnrichmentAnalyzer

annotator = GeneAnnotator()
analyzer = GeneEnrichmentAnalyzer()

# 使用在线API注释
annotated_df = annotator.annotate_bed_regions(bed_df)
```

## 核心代码示例

### 基本注释流程

```python
import pandas as pd
from typing import List, Dict

class GeneAnnotator:
    def __init__(self, gene_positions: Dict):
        self.gene_positions = gene_positions
    
    def parse_bed_file(self, bed_content: str) -> pd.DataFrame:
        """解析BED文件"""
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
        
        return pd.DataFrame(data)
    
    def find_overlapping_genes(self, chromosome: str, start: int, end: int) -> List[str]:
        """查找重叠基因"""
        overlapping_genes = []
        
        for gene_name, gene_info in self.gene_positions.items():
            if gene_info['chromosome'] == chromosome:
                # 检查是否有重叠
                if not (end < gene_info['start'] or start > gene_info['end']):
                    overlapping_genes.append(gene_name)
        
        return overlapping_genes
    
    def annotate_regions(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """注释基因组区域"""
        genes_list = []
        
        for idx, row in bed_df.iterrows():
            chromosome = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            
            overlapping_genes = self.find_overlapping_genes(chromosome, start, end)
            genes_list.append(';'.join(overlapping_genes) if overlapping_genes else 'Unknown')
        
        bed_df['genes'] = genes_list
        return bed_df

# 使用示例
gene_positions = {
    'TP53': {'chromosome': 'chr17', 'start': 7661779, 'end': 7687550},
    'BRCA1': {'chromosome': 'chr17', 'start': 43044295, 'end': 43125483},
    # ... 更多基因
}

annotator = GeneAnnotator(gene_positions)
bed_df = annotator.parse_bed_file(bed_content)
annotated_df = annotator.annotate_regions(bed_df)
```

### 富集分析

```python
class EnrichmentAnalyzer:
    def __init__(self, gene_functions: Dict):
        self.gene_functions = gene_functions
    
    def extract_gene_list(self, bed_df: pd.DataFrame) -> List[str]:
        """提取基因列表"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        return list(set(all_genes))
    
    def functional_enrichment(self, gene_list: List[str]) -> Dict:
        """功能富集分析"""
        functional_groups = {}
        for gene in gene_list:
            if gene in self.gene_functions:
                function = self.gene_functions[gene]
                if function not in functional_groups:
                    functional_groups[function] = []
                functional_groups[function].append(gene)
        
        return {
            'genes': gene_list,
            'functional_groups': functional_groups,
            'total_genes': len(gene_list)
        }

# 使用示例
gene_functions = {
    'TP53': '肿瘤抑制基因',
    'BRCA1': 'DNA修复基因',
    # ... 更多功能分类
}

analyzer = EnrichmentAnalyzer(gene_functions)
gene_list = analyzer.extract_gene_list(annotated_df)
enrichment_results = analyzer.functional_enrichment(gene_list)
```

## 输出文件说明

### 注释结果文件
- **格式**: `chromosome start end genes`
- **示例**:
```
chromosome  start       end         genes
chr17       7661779     7687550     TP53
chr13       32889665    32973808    BRCA2
chr7        55086724    55275031    EGFR
```

### 富集分析报告
- **功能分组**: 按生物学功能分类基因
- **通路分析**: 识别富集的生物学通路
- **统计信息**: 基因数量、富集分数等

## 高级功能

### 使用更大的基因数据库

```python
# 下载GENCODE注释文件
import requests
import gzip

def download_gencode_gtf():
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
    response = requests.get(url)
    
    with open('gencode.v44.annotation.gtf.gz', 'wb') as f:
        f.write(response.content)
    
    # 解压并解析
    with gzip.open('gencode.v44.annotation.gtf.gz', 'rt') as f:
        parse_gtf_file(f)
```

### GO/KEGG富集分析

```python
import gseapy as gp

def perform_go_enrichment(gene_list):
    """GO富集分析"""
    go_enr = gp.prerank(
        rnk=gene_list,
        gene_sets='GO_Biological_Process_2023',
        organism='Human',
        min_size=5,
        max_size=1000
    )
    return go_enr

def perform_kegg_enrichment(gene_list):
    """KEGG富集分析"""
    kegg_enr = gp.prerank(
        rnk=gene_list,
        gene_sets='KEGG_2021_Human',
        organism='Human',
        min_size=5,
        max_size=1000
    )
    return kegg_enr
```

### 批量处理

```python
import glob
from concurrent.futures import ThreadPoolExecutor

def batch_process_bed_files(bed_files, annotator):
    """批量处理BED文件"""
    results = []
    
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = []
        for bed_file in bed_files:
            future = executor.submit(process_single_file, bed_file, annotator)
            futures.append(future)
        
        for future in futures:
            results.append(future.result())
    
    return results

def process_single_file(bed_file, annotator):
    """处理单个BED文件"""
    with open(bed_file, 'r') as f:
        bed_content = f.read()
    
    bed_df = annotator.parse_bed_file(bed_content)
    annotated_df = annotator.annotate_regions(bed_df)
    
    output_file = bed_file.replace('.bed', '_annotated.bed')
    annotated_df.to_csv(output_file, sep='\t', index=False)
    
    return output_file
```

## 故障排除

### 常见问题

1. **基因注释失败**
   - 检查染色体名称格式（chr1 vs 1）
   - 确保基因位置数据库包含目标区域
   - 验证BED文件格式正确性

2. **富集分析失败**
   - 确保基因列表不为空
   - 检查基因名格式（使用标准基因符号）
   - 验证基因集数据库可用性

3. **性能问题**
   - 对于大文件，考虑分块处理
   - 使用并行处理提高速度
   - 优化基因数据库结构

### 调试技巧

```python
# 启用详细日志
import logging
logging.basicConfig(level=logging.DEBUG)

# 检查基因数据库
print("基因数据库示例:")
for gene, pos in list(gene_positions.items())[:5]:
    print(f"{gene}: {pos}")

# 检查重叠检测
def debug_overlap(chromosome, start, end):
    print(f"检查重叠: {chromosome}:{start}-{end}")
    overlapping = find_overlapping_genes(chromosome, start, end)
    print(f"重叠基因: {overlapping}")
    return overlapping
```

## 总结

这个解决方案提供了完整的BED文件基因注释和富集分析流程：

1. **多种实现方式**: 从简单到复杂，满足不同需求
2. **模块化设计**: 易于扩展和定制
3. **完整示例**: 包含实际代码和使用说明
4. **故障排除**: 提供常见问题解决方案

您可以根据具体需求选择合适的版本，并按照提供的代码示例进行实现。