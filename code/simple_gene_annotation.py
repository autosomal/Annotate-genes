#!/usr/bin/env python3
"""
简化版BED文件基因注释工具
使用在线数据库和文件下载进行基因注释
"""

import pandas as pd
import requests
import gzip
import shutil
from pathlib import Path
import time
from typing import List, Dict, Optional
import subprocess
import sys

class SimpleGeneAnnotator:
    """简化的基因注释器"""
    
    def __init__(self):
        self.gene_positions = {}
        self.load_gene_positions()
    
    def load_gene_positions(self):
        """加载基因位置信息"""
        print("加载基因位置数据...")
        
        # 尝试下载GENCODE基因位置文件
        try:
            self.download_gencode_annotations()
        except Exception as e:
            print(f"下载GENCODE注释失败: {e}")
            print("使用示例数据进行演示")
            self.load_sample_data()
    
    def download_gencode_annotations(self):
        """下载GENCODE注释文件"""
        gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
        local_file = "/workspace/data/gencode.annotation.gtf.gz"
        
        print("下载GENCODE注释文件...")
        response = requests.get(gencode_url, stream=True)
        response.raise_for_status()
        
        with open(local_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print("解析GENCODE注释文件...")
        self.parse_gtf_file(local_file)
    
    def parse_gtf_file(self, gtf_file: str):
        """解析GTF文件"""
        try:
            # 解压文件
            with gzip.open(gtf_file, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        chromosome = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        attributes = fields[8]
                        
                        # 提取基因名
                        if 'gene_name' in attributes:
                            gene_name = attributes.split('gene_name "')[1].split('"')[0]
                            self.gene_positions[gene_name] = {
                                'chromosome': chromosome,
                                'start': start,
                                'end': end
                            }
        except Exception as e:
            print(f"解析GTF文件失败: {e}")
            self.load_sample_data()
    
    def load_sample_data(self):
        """加载示例基因位置数据"""
        # 一些常见基因的示例位置数据
        sample_genes = {
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
        self.gene_positions = sample_genes
        print(f"加载了 {len(self.gene_positions)} 个基因的位置信息")
    
    def parse_bed_file(self, bed_file: str) -> pd.DataFrame:
        """解析BED文件"""
        try:
            bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                               names=['chromosome', 'start', 'end'])
            print(f"成功读取BED文件，共{len(bed_df)}个区域")
            return bed_df
        except Exception as e:
            print(f"读取BED文件失败: {e}")
            return None
    
    def find_overlapping_genes(self, chromosome: str, start: int, end: int) -> List[str]:
        """查找重叠的基因"""
        overlapping_genes = []
        
        for gene_name, gene_info in self.gene_positions.items():
            if gene_info['chromosome'] == chromosome:
                # 检查是否有重叠
                if not (end < gene_info['start'] or start > gene_info['end']):
                    overlapping_genes.append(gene_name)
        
        return overlapping_genes
    
    def annotate_bed_regions(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """注释BED区域"""
        genes_list = []
        
        print("开始基因注释...")
        for idx, row in bed_df.iterrows():
            chromosome = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            
            # 查找重叠基因
            overlapping_genes = self.find_overlapping_genes(chromosome, start, end)
            
            if overlapping_genes:
                genes_list.append(';'.join(overlapping_genes))
            else:
                genes_list.append('Unknown')
            
            if (idx + 1) % 5 == 0:
                print(f"已处理 {idx + 1}/{len(bed_df)} 个区域")
        
        bed_df['genes'] = genes_list
        return bed_df

class SimpleEnrichmentAnalyzer:
    """简化的富集分析器"""
    
    def __init__(self):
        pass
    
    def get_gene_list_from_bed(self, bed_df: pd.DataFrame) -> List[str]:
        """从注释后的BED文件中提取基因列表"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        unique_genes = list(set(all_genes))
        print(f"提取到 {len(unique_genes)} 个唯一基因: {unique_genes}")
        return unique_genes
    
    def simple_enrichment_analysis(self, gene_list: List[str]) -> Dict:
        """简单的富集分析（基于基因功能分类）"""
        print("执行简单的富集分析...")
        
        # 癌症相关基因分类
        cancer_genes = {
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
        
        # 分析基因功能
        functional_analysis = {}
        for gene in gene_list:
            if gene in cancer_genes:
                function = cancer_genes[gene]
                if function not in functional_analysis:
                    functional_analysis[function] = []
                functional_analysis[function].append(gene)
        
        return {
            'genes': gene_list,
            'functional_groups': functional_analysis,
            'total_genes': len(gene_list)
        }
    
    def save_results(self, bed_df: pd.DataFrame, enrichment_results: Dict, output_dir: str):
        """保存分析结果"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存注释结果
        bed_df.to_csv(f"{output_dir}/annotated_genes.bed", sep='\t', index=False)
        
        # 保存富集分析结果
        with open(f"{output_dir}/enrichment_analysis.txt", 'w', encoding='utf-8') as f:
            f.write("基因富集分析结果\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"分析的基因总数: {enrichment_results['total_genes']}\n")
            f.write(f"基因列表: {', '.join(enrichment_results['genes'])}\n\n")
            
            f.write("功能分组:\n")
            f.write("-" * 30 + "\n")
            for function, genes in enrichment_results['functional_groups'].items():
                f.write(f"{function}: {', '.join(genes)}\n")
        
        print(f"结果已保存到: {output_dir}/")

def main():
    """主函数"""
    print("=== 简化版BED文件基因注释和富集分析 ===\n")
    
    # 创建数据目录
    Path('/workspace/data').mkdir(exist_ok=True)
    
    # 示例BED文件内容
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
    
    # 保存示例BED文件
    with open('/workspace/data/example.bed', 'w') as f:
        f.write(bed_content)
    
    # 初始化注释器和分析器
    annotator = SimpleGeneAnnotator()
    analyzer = SimpleEnrichmentAnalyzer()
    
    # 1. 解析BED文件
    print("1. 解析BED文件...")
    bed_df = annotator.parse_bed_file('/workspace/data/example.bed')
    if bed_df is None:
        return
    
    print("BED文件内容:")
    print(bed_df)
    print()
    
    # 2. 注释基因
    print("2. 注释基因...")
    annotated_df = annotator.annotate_bed_regions(bed_df)
    
    print("基因注释结果:")
    print(annotated_df)
    print()
    
    # 3. 提取基因列表
    print("3. 提取基因列表...")
    gene_list = analyzer.get_gene_list_from_bed(annotated_df)
    print()
    
    # 4. 富集分析
    print("4. 执行富集分析...")
    enrichment_results = analyzer.simple_enrichment_analysis(gene_list)
    print()
    
    # 5. 保存结果
    print("5. 保存结果...")
    analyzer.save_results(annotated_df, enrichment_results, '/workspace/data')
    
    print("\n=== 分析完成 ===")
    print("结果文件:")
    print("- 注释结果: /workspace/data/annotated_genes.bed")
    print("- 富集分析: /workspace/data/enrichment_analysis.txt")
    
    # 显示结果摘要
    print("\n结果摘要:")
    print(f"总共注释了 {len(annotated_df)} 个基因组区域")
    print(f"识别出 {len(gene_list)} 个基因")
    print(f"功能分组数量: {len(enrichment_results['functional_groups'])}")

if __name__ == "__main__":
    main()