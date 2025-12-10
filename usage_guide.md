# BED文件基因注释和富集分析使用指南

## 概述

本工具包提供了三种不同复杂度的BED文件基因注释和富集分析解决方案：

1. **API版本** (`bed_gene_annotation.py`) - 使用在线API进行注释
2. **本地版本** (`local_gene_annotation.py`) - 使用本地数据库和pyensembl
3. **简化版本** (`simple_gene_annotation.py`) - 使用本地数据和简化分析

## 快速开始

### 环境要求

```bash
# 安装基本依赖
pip install pandas numpy requests

# 可选：安装高级功能依赖
pip install pyensembl gseapy
```

### 准备您的BED文件

BED文件格式要求：
```
chromosome    start    end
chr1          100000   101000
chr2          200000   201000
```

## 使用方法

### 方法1：简化版本（推荐入门）

```python
from simple_gene_annotation import SimpleGeneAnnotator, SimpleEnrichmentAnalyzer

# 初始化
annotator = SimpleGeneAnnotator()
analyzer = SimpleEnrichmentAnalyzer()

# 读取BED文件
bed_df = annotator.parse_bed_file('your_file.bed')

# 注释基因
annotated_df = annotator.annotate_bed_regions(bed_df)

# 提取基因列表
gene_list = analyzer.get_gene_list_from_bed(annotated_df)

# 执行富集分析
enrichment_results = analyzer.simple_enrichment_analysis(gene_list)

# 保存结果
analyzer.save_results(annotated_df, enrichment_results, 'output_dir')
```

### 方法2：本地版本（推荐高级用户）

```python
from local_gene_annotation import LocalGeneAnnotator, EnrichmentAnalyzer

# 初始化（可指定GTF文件）
annotator = LocalGeneAnnotator(gtf_file='path/to/annotation.gtf')
analyzer = EnrichmentAnalyzer()

# 执行分析
bed_df = annotator.parse_bed_file('your_file.bed')
annotated_df = annotator.annotate_with_pyensembl(bed_df)

# GO和KEGG富集分析
go_results = analyzer.perform_go_enrichment(gene_list)
kegg_results = analyzer.perform_kegg_enrichment(gene_list)

# 保存结果
analyzer.save_enrichment_results(go_results, kegg_results, 'output_dir')
```

### 方法3：API版本（需要网络连接）

```python
from bed_gene_annotation import GeneAnnotator, GeneEnrichmentAnalyzer

# 初始化
annotator = GeneAnnotator()
analyzer = GeneEnrichmentAnalyzer()

# 使用在线API进行注释和分析
annotated_df = annotator.annotate_bed_regions(bed_df)
go_results = analyzer.perform_go_enrichment(gene_list)
kegg_results = analyzer.perform_kegg_enrichment(gene_list)
```

## 输出文件说明

### 注释结果文件
- `annotated_genes.bed` - 包含原始BED信息和注释的基因名
- 格式：`chromosome    start    end    genes`

### 富集分析结果
- `go_enrichment_results.csv` - GO功能富集分析结果
- `kegg_enrichment_results.csv` - KEGG通路富集分析结果
- `enrichment_analysis_report.txt` - 分析报告

## 高级用法

### 使用自定义GTF文件

```python
# 下载人类基因组注释
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# 使用本地注释文件
annotator = LocalGeneAnnotator(gtf_file='gencode.v44.annotation.gtf')
```

### 批量处理多个BED文件

```python
import glob

bed_files = glob.glob('*.bed')
for bed_file in bed_files:
    print(f"处理文件: {bed_file}")
    
    # 处理每个文件
    bed_df = annotator.parse_bed_file(bed_file)
    annotated_df = annotator.annotate_bed_regions(bed_df)
    
    # 保存结果
    output_name = bed_file.replace('.bed', '_annotated.bed')
    annotated_df.to_csv(output_name, sep='\t', index=False)
```

## 故障排除

### 常见问题

1. **pyensembl安装失败**
   ```bash
   # 手动安装
   pip install --upgrade pyensembl
   
   # 或者使用conda
   conda install -c bioconda pyensembl
   ```

2. **网络连接问题**
   - 使用本地版本或简化版本
   - 检查防火墙设置

3. **基因注释失败**
   - 检查BED文件格式是否正确
   - 确保染色体名称格式一致（chr1 vs 1）

4. **富集分析失败**
   - 确保基因列表不为空
   - 检查基因名格式（使用标准基因符号）

### 调试模式

```python
# 启用详细输出
import logging
logging.basicConfig(level=logging.DEBUG)

# 检查数据
print("基因位置数据示例:")
for gene, pos in list(annotator.gene_positions.items())[:5]:
    print(f"{gene}: {pos}")
```

## 性能优化

### 大文件处理

```python
# 分块处理大文件
def process_large_bed(bed_file, chunk_size=1000):
    chunks = pd.read_csv(bed_file, sep='\t', header=None, 
                        names=['chromosome', 'start', 'end'], 
                        chunksize=chunk_size)
    
    for i, chunk in enumerate(chunks):
        print(f"处理第 {i+1} 块...")
        annotated_chunk = annotator.annotate_bed_regions(chunk)
        # 保存或合并结果
```

### 并行处理

```python
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

def parallel_annotation(bed_df, n_cores=None):
    if n_cores is None:
        n_cores = multiprocessing.cpu_count()
    
    # 将数据分成多个部分
    chunks = [bed_df[i::n_cores] for i in range(n_cores)]
    
    with ThreadPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(annotate_chunk, chunks))
    
    return pd.concat(results, ignore_index=True)
```

## 扩展功能

### 添加自定义基因数据库

```python
class CustomGeneAnnotator(SimpleGeneAnnotator):
    def load_custom_database(self, db_file):
        """加载自定义基因数据库"""
        custom_df = pd.read_csv(db_file)
        for _, row in custom_df.iterrows():
            gene_name = row['gene_name']
            self.gene_positions[gene_name] = {
                'chromosome': row['chromosome'],
                'start': row['start'],
                'end': row['end'],
                'description': row.get('description', '')
            }
```

### 自定义富集分析

```python
class CustomEnrichmentAnalyzer(SimpleEnrichmentAnalyzer):
    def custom_pathway_analysis(self, gene_list):
        """自定义通路分析"""
        # 添加您的分析逻辑
        pathway_results = {}
        
        # 示例：分析特定基因集
        your_gene_sets = {
            'Cell_Cycle': ['TP53', 'RB1', 'MYC'],
            'DNA_Repair': ['BRCA1', 'BRCA2'],
            'Growth_Factor': ['EGFR', 'AKT1']
        }
        
        for pathway, pathway_genes in your_gene_sets.items():
            overlap = set(gene_list) & set(pathway_genes)
            if overlap:
                pathway_results[pathway] = list(overlap)
        
        return pathway_results
```

## 联系和支持

如有问题或建议，请检查：
1. 文件格式是否正确
2. 依赖包是否正确安装
3. 网络连接是否正常
4. 查看错误日志进行调试