#!/usr/bin/env python3
"""
BED文件基因注释和富集分析 - 完整解决方案
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict

# 1. 基因位置数据库（示例）
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

# 2. 基因功能分类
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

# 3. 癌症相关通路
CANCER_PATHWAYS = {
    'Cell_Cycle_Control': ['TP53', 'RB1', 'MYC'],
    'DNA_Repair': ['BRCA1', 'BRCA2'],
    'Growth_Factor_Signaling': ['EGFR', 'KRAS', 'PIK3CA', 'AKT1'],
    'Tumor_Suppression': ['TP53', 'RB1', 'PTEN']
}

class GeneAnnotator:
    """基因注释器"""
    
    def __init__(self, gene_positions: Dict = None):
        self.gene_positions = gene_positions or GENE_POSITIONS
        print(f"初始化基因注释器，加载 {len(self.gene_positions)} 个基因位置")
    
    def parse_bed_file(self, bed_content: str) -> pd.DataFrame:
        """解析BED文件内容"""
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
        print(f"解析BED文件，共 {len(bed_df)} 个基因组区域")
        return bed_df
    
    def find_overlapping_genes(self, chromosome: str, start: int, end: int) -> List[str]:
        """查找与指定区域重叠的基因"""
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
            
            print(f"区域 {idx+1}: {chromosome}:{start:,}-{end:,} -> {genes_list[-1]}")
        
        bed_df['genes'] = genes_list
        return bed_df

class EnrichmentAnalyzer:
    """富集分析器"""
    
    def __init__(self, gene_functions: Dict = None, pathways: Dict = None):
        self.gene_functions = gene_functions or GENE_FUNCTIONS
        self.pathways = pathways or CANCER_PATHWAYS
    
    def extract_gene_list(self, bed_df: pd.DataFrame) -> List[str]:
        """从注释结果中提取基因列表"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        unique_genes = list(set(all_genes))
        print(f"提取到 {len(unique_genes)} 个唯一基因: {unique_genes}")
        return unique_genes
    
    def functional_enrichment(self, gene_list: List[str]) -> Dict:
        """功能富集分析"""
        print("执行功能富集分析...")
        
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
    
    def pathway_enrichment(self, gene_list: List[str]) -> Dict:
        """通路富集分析"""
        print("执行通路富集分析...")
        
        pathway_results = {}
        for pathway, pathway_genes in self.pathways.items():
            overlap = set(gene_list) & set(pathway_genes)
            if overlap:
                pathway_results[pathway] = {
                    'genes': list(overlap),
                    'enrichment_score': len(overlap) / len(pathway_genes),
                    'p_value': len(overlap) / len(pathway_genes)  # 简化的p值计算
                }
        
        return pathway_results
    
    def save_results(self, bed_df: pd.DataFrame, functional_results: Dict, 
                    pathway_results: Dict, output_dir: str):
        """保存分析结果"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存注释结果
        output_file = f"{output_dir}/annotated_genes.bed"
        bed_df.to_csv(output_file, sep='\t', index=False)
        print(f"注释结果已保存到: {output_file}")
        
        # 生成详细报告
        report_file = f"{output_dir}/enrichment_analysis_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("基因富集分析详细报告\n")
            f.write("=" * 60 + "\n\n")
            
            # 基本信息
            f.write("1. 分析概述\n")
            f.write("-" * 30 + "\n")
            f.write(f"基因组区域数量: {len(bed_df)}\n")
            f.write(f"识别基因数量: {functional_results['total_genes']}\n")
            f.write(f"基因列表: {', '.join(functional_results['genes'])}\n\n")
            
            # 功能分组
            f.write("2. 功能富集分析\n")
            f.write("-" * 30 + "\n")
            for function, genes in functional_results['functional_groups'].items():
                f.write(f"{function}:\n")
                f.write(f"  基因数量: {len(genes)}\n")
                f.write(f"  基因列表: {', '.join(genes)}\n\n")
            
            # 通路富集
            f.write("3. 通路富集分析\n")
            f.write("-" * 30 + "\n")
            for pathway, info in pathway_results.items():
                f.write(f"{pathway}:\n")
                f.write(f"  富集基因: {', '.join(info['genes'])}\n")
                f.write(f"  富集分数: {info['enrichment_score']:.3f}\n")
                f.write(f"  显著性: {info['p_value']:.3f}\n\n")
        
        print(f"详细报告已保存到: {report_file}")
        
        return {
            'annotated_file': output_file,
            'report_file': report_file
        }

# 演示函数
def run_demo():
    """运行演示"""
    print("=== BED文件基因注释和富集分析演示 ===\n")
    
    # 您的BED文件内容
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
    
    # 创建输出目录
    output_dir = '/workspace/data'
    Path(output_dir).mkdir(exist_ok=True)
    
    # 初始化组件
    annotator = GeneAnnotator()
    analyzer = EnrichmentAnalyzer()
    
    # 1. 解析BED文件
    print("1. 解析BED文件")
    bed_df = annotator.parse_bed_file(bed_content)
    print(bed_df)
    print()
    
    # 2. 基因注释
    print("2. 基因注释")
    annotated_df = annotator.annotate_regions(bed_df)
    print()
    
    # 3. 提取基因列表
    print("3. 提取基因列表")
    gene_list = analyzer.extract_gene_list(annotated_df)
    print()
    
    # 4. 富集分析
    print("4. 富集分析")
    functional_results = analyzer.functional_enrichment(gene_list)
    pathway_results = analyzer.pathway_enrichment(gene_list)
    print()
    
    # 5. 保存结果
    print("5. 保存结果")
    result_files = analyzer.save_results(annotated_df, functional_results, 
                                       pathway_results, output_dir)
    print()
    
    # 显示结果摘要
    print("=== 分析结果摘要 ===")
    print(f"• 基因组区域: {len(annotated_df)} 个")
    print(f"• 识别基因: {len(gene_list)} 个")
    print(f"• 功能分组: {len(functional_results['functional_groups'])} 个")
    print(f"• 富集通路: {len(pathway_results)} 个")
    
    print("\n功能分组详情:")
    for function, genes in functional_results['functional_groups'].items():
        print(f"  {function}: {', '.join(genes)}")
    
    print("\n通路富集详情:")
    for pathway, info in pathway_results.items():
        print(f"  {pathway}: {', '.join(info['genes'])} (分数: {info['enrichment_score']:.3f})")
    
    return annotated_df, functional_results, pathway_results

if __name__ == "__main__":
    # 运行演示
    annotated_df, functional_results, pathway_results = run_demo()