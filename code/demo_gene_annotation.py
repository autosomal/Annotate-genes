#!/usr/bin/env python3
"""
演示版本：BED文件基因注释和富集分析
直接使用示例数据，无需下载
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict

class DemoGeneAnnotator:
    """演示版基因注释器"""
    
    def __init__(self):
        # 示例基因位置数据
        self.gene_positions = {
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
        print(f"加载了 {len(self.gene_positions)} 个基因的位置信息")
    
    def parse_bed_file(self, bed_content: str) -> pd.DataFrame:
        """解析BED内容"""
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
        print(f"解析BED文件，共 {len(bed_df)} 个区域")
        return bed_df
    
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
            
            print(f"区域 {idx+1}: {chromosome}:{start}-{end} -> {genes_list[-1]}")
        
        bed_df['genes'] = genes_list
        return bed_df

class DemoEnrichmentAnalyzer:
    """演示版富集分析器"""
    
    def __init__(self):
        # 基因功能分类
        self.gene_functions = {
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
    
    def functional_enrichment_analysis(self, gene_list: List[str]) -> Dict:
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
    
    def pathway_analysis(self, gene_list: List[str]) -> Dict:
        """通路分析"""
        print("执行通路分析...")
        
        # 癌症相关通路
        cancer_pathways = {
            'Cell_Cycle_Control': ['TP53', 'RB1', 'MYC'],
            'DNA_Repair': ['BRCA1', 'BRCA2'],
            'Growth_Factor_Signaling': ['EGFR', 'KRAS', 'PIK3CA', 'AKT1'],
            'Tumor_Suppression': ['TP53', 'RB1', 'PTEN']
        }
        
        pathway_enrichment = {}
        for pathway, pathway_genes in cancer_pathways.items():
            overlap = set(gene_list) & set(pathway_genes)
            if overlap:
                pathway_enrichment[pathway] = {
                    'genes': list(overlap),
                    'enrichment_score': len(overlap) / len(pathway_genes)
                }
        
        return pathway_enrichment
    
    def save_results(self, bed_df: pd.DataFrame, enrichment_results: Dict, 
                    pathway_results: Dict, output_dir: str):
        """保存分析结果"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存注释结果
        bed_df.to_csv(f"{output_dir}/annotated_genes.bed", sep='\t', index=False)
        print(f"注释结果已保存到: {output_dir}/annotated_genes.bed")
        
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
            
            f.write("\n\n通路富集分析:\n")
            f.write("-" * 30 + "\n")
            for pathway, info in pathway_results.items():
                f.write(f"{pathway}:\n")
                f.write(f"  富集基因: {', '.join(info['genes'])}\n")
                f.write(f"  富集分数: {info['enrichment_score']:.3f}\n\n")
        
        print(f"富集分析结果已保存到: {output_dir}/enrichment_analysis.txt")

def main():
    """主演示函数"""
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
    
    # 创建数据目录
    Path('/workspace/data').mkdir(exist_ok=True)
    
    # 保存示例BED文件
    with open('/workspace/data/demo_example.bed', 'w') as f:
        f.write(bed_content)
    
    # 初始化注释器和分析器
    annotator = DemoGeneAnnotator()
    analyzer = DemoEnrichmentAnalyzer()
    
    print("\n1. 解析BED文件...")
    bed_df = annotator.parse_bed_file(bed_content)
    print("BED文件内容:")
    print(bed_df)
    print()
    
    print("2. 注释基因...")
    annotated_df = annotator.annotate_bed_regions(bed_df)
    print()
    
    print("3. 提取基因列表...")
    gene_list = analyzer.get_gene_list_from_bed(annotated_df)
    print()
    
    print("4. 执行富集分析...")
    enrichment_results = analyzer.functional_enrichment_analysis(gene_list)
    pathway_results = analyzer.pathway_analysis(gene_list)
    print()
    
    print("5. 保存结果...")
    analyzer.save_results(annotated_df, enrichment_results, pathway_results, '/workspace/data')
    
    print("\n=== 演示完成 ===")
    print("\n结果摘要:")
    print(f"• 总共注释了 {len(annotated_df)} 个基因组区域")
    print(f"• 识别出 {len(gene_list)} 个基因")
    print(f"• 功能分组数量: {len(enrichment_results['functional_groups'])}")
    print(f"• 富集通路数量: {len(pathway_results)}")
    
    print("\n输出文件:")
    print("• /workspace/data/annotated_genes.bed - 注释结果")
    print("• /workspace/data/enrichment_analysis.txt - 富集分析报告")
    
    # 显示详细结果
    print("\n=== 详细结果 ===")
    print("\n注释结果:")
    print(annotated_df.to_string(index=False))
    
    print("\n功能分组:")
    for function, genes in enrichment_results['functional_groups'].items():
        print(f"  {function}: {', '.join(genes)}")
    
    print("\n通路富集:")
    for pathway, info in pathway_results.items():
        print(f"  {pathway}: {', '.join(info['genes'])} (富集分数: {info['enrichment_score']:.3f})")

if __name__ == "__main__":
    main()