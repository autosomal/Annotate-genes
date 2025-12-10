#!/usr/bin/env python3
"""
Demo Version: BED File Gene Annotation and Enrichment Analysis
直接使用示例数据，无需下载
"""

import pandas as pd
from pathlib import Path
from typing import List, Dict

class DemoGeneAnnotator:
    """Demo Gene Annotator"""
    
    def __init__(self):
        # Example gene position data
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
        print(f"Loaded {len(self.gene_positions)} gene position records")
    
    def parse_bed_file(self, bed_content: str) -> pd.DataFrame:
        """Parse BED content"""
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
        print(f"Parse BED file, total {len(bed_df)} regions")
        return bed_df
    
    def find_overlapping_genes(self, chromosome: str, start: int, end: int) -> List[str]:
        """Find overlapping genes"""
        overlapping_genes = []
        
        for gene_name, gene_info in self.gene_positions.items():
            if gene_info['chromosome'] == chromosome:
                # 检查是否有重叠
                if not (end < gene_info['start'] or start > gene_info['end']):
                    overlapping_genes.append(gene_name)
        
        return overlapping_genes
    
    def annotate_bed_regions(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """Annotate BED regions"""
        genes_list = []
        
        print("Start gene annotation...")
        for idx, row in bed_df.iterrows():
            chromosome = str(row['chromosome'])
            start = int(row['start'])
            end = int(row['end'])
            
            # Find overlapping genes
            overlapping_genes = self.find_overlapping_genes(chromosome, start, end)
            
            if overlapping_genes:
                genes_list.append(';'.join(overlapping_genes))
            else:
                genes_list.append('Unknown')
            
            print(f"Region {idx+1}: {chromosome}:{start}-{end} -> {genes_list[-1]}")
        
        bed_df['genes'] = genes_list
        return bed_df

class DemoEnrichmentAnalyzer:
    """Demo Enrichment Analyzer"""
    
    def __init__(self):
        # Gene function classification
        self.gene_functions = {
            'TP53': 'Tumor Suppressor Gene',
            'BRCA1': 'DNA Repair Gene',
            'BRCA2': 'DNA Repair Gene',
            'EGFR': 'Growth Factor Receptor',
            'KRAS': 'Oncogene',
            'PIK3CA': 'PI3K Pathway',
            'AKT1': 'PI3K/AKT Pathway',
            'MYC': '转录因子',
            'RB1': 'Tumor Suppressor Gene',
            'PTEN': 'Tumor Suppressor Gene'
        }
    
    def get_gene_list_from_bed(self, bed_df: pd.DataFrame) -> List[str]:
        """Extract gene list from annotated BED file"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        unique_genes = list(set(all_genes))
        print(f"Extracted {len(unique_genes)} unique genes: {unique_genes}")
        return unique_genes
    
    def functional_enrichment_analysis(self, gene_list: List[str]) -> Dict:
        """Functional enrichment analysis"""
        print("Perform functional enrichment analysis...")
        
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
        """Pathway analysis"""
        print("Perform pathway analysis...")
        
        # Cancer-related pathways
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
        """Save analysis results"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # Save annotation results
        bed_df.to_csv(f"{output_dir}/annotated_genes.bed", sep='\t', index=False)
        print(f"Annotation results saved to: {output_dir}/annotated_genes.bed")
        
        # Save enrichment analysis results
        with open(f"{output_dir}/enrichment_analysis.txt", 'w', encoding='utf-8') as f:
            f.write("Gene Enrichment Analysis Results\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total number of analyzed genes: {enrichment_results['total_genes']}\n")
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
        
        print(f"Enrichment analysis results saved to: {output_dir}/enrichment_analysis.txt")

def main():
    """Main demo function"""
    print("=== BED File Gene Annotation and Enrichment Analysis Demo ===\n")
    
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
    
    # Save example BED file
    with open('/workspace/data/demo_example.bed', 'w') as f:
        f.write(bed_content)
    
    # Initialize annotator and analyzer
    annotator = DemoGeneAnnotator()
    analyzer = DemoEnrichmentAnalyzer()
    
    print("\n1. Parse BED file...")
    bed_df = annotator.parse_bed_file(bed_content)
    print("BED file content:")
    print(bed_df)
    print()
    
    print("2. Gene annotation...")
    annotated_df = annotator.annotate_bed_regions(bed_df)
    print()
    
    print("3. Extract gene list...")
    gene_list = analyzer.get_gene_list_from_bed(annotated_df)
    print()
    
    print("4. Perform enrichment analysis...")
    enrichment_results = analyzer.functional_enrichment_analysis(gene_list)
    pathway_results = analyzer.pathway_analysis(gene_list)
    print()
    
    print("5. Save results...")
    analyzer.save_results(annotated_df, enrichment_results, pathway_results, '/workspace/data')
    
    print("\n=== Demo completed ===")
    print("\nResult summary:")
    print(f"• Annotated {len(annotated_df)} genomic regions")
    print(f"• Identified {len(gene_list)} genes")
    print(f"• Number of functional groups: {len(enrichment_results['functional_groups'])}")
    print(f"• Number of enriched pathways: {len(pathway_results)}")
    
    print("\nOutput files:")
    print("• /workspace/data/annotated_genes.bed - Annotation results")
    print("• /workspace/data/enrichment_analysis.txt - Enrichment analysis report")
    
    # 显示详细结果
    print("\n=== Detailed results ===")
    print("\nAnnotation results:")
    print(annotated_df.to_string(index=False))
    
    print("\nFunctional groups:")
    for function, genes in enrichment_results['functional_groups'].items():
        print(f"  {function}: {', '.join(genes)}")
    
    print("\nPathway enrichment:")
    for pathway, info in pathway_results.items():
        print(f"  {pathway}: {', '.join(info['genes'])} (enrichment score: {info['enrichment_score']:.3f})")

if __name__ == "__main__":
    main()