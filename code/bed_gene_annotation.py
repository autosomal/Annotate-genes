#!/usr/bin/env python3
"""
BED文件基因注释和富集分析工具
将BED格式的基因组位置信息注释到基因名，并进行功能富集分析
"""

import pandas as pd
import numpy as np
import requests
import json
from typing import List, Dict, Tuple, Optional
import time
from pathlib import Path

class GeneAnnotator:
    """基因组位置注释器"""
    
    def __init__(self):
        self.gencode_url = "https://rest.genenames.org/fetch/region/"
        self.header = {"Content-Type": "application/json"}
        
    def parse_bed_file(self, bed_file: str) -> pd.DataFrame:
        """解析BED文件"""
        try:
            # 读取BED文件，假设没有header
            bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                               names=['chromosome', 'start', 'end'])
            print(f"成功读取BED文件，共{len(bed_df)}个区域")
            return bed_df
        except Exception as e:
            print(f"读取BED文件失败: {e}")
            return None
    
    def annotate_with_gencode_api(self, chromosome: str, start: int, end: int) -> List[str]:
        """使用GENCODE API注释基因"""
        genes = []
        try:
            # 构建GENCODE API请求
            region = f"{chromosome}:{start}:{end}"
            url = f"{self.gencode_url}{region}"
            
            response = requests.get(url, headers=self.header)
            if response.status_code == 200:
                data = response.json()
                docs = data.get('response', {}).get('docs', [])
                for doc in docs:
                    if 'symbol' in doc:
                        genes.append(doc['symbol'])
            time.sleep(0.1)  # 避免请求过快
        except Exception as e:
            print(f"API请求失败 {chromosome}:{start}-{end}: {e}")
        
        return genes
    
    def annotate_bed_regions(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """注释BED区域到基因"""
        genes_list = []
        
        print("开始注释基因...")
        for idx, row in bed_df.iterrows():
            chromosome = f"chr{row['chromosome']}" if not str(row['chromosome']).startswith('chr') else str(row['chromosome'])
            genes = self.annotate_with_gencode_api(chromosome, row['start'], row['end'])
            genes_list.append(';'.join(genes) if genes else 'Unknown')
            
            if (idx + 1) % 10 == 0:
                print(f"已处理 {idx + 1}/{len(bed_df)} 个区域")
        
        bed_df['genes'] = genes_list
        return bed_df
    
    def save_annotated_results(self, bed_df: pd.DataFrame, output_file: str):
        """保存注释结果"""
        bed_df.to_csv(output_file, sep='\t', index=False)
        print(f"注释结果已保存到: {output_file}")

class GeneEnrichmentAnalyzer:
    """基因富集分析器"""
    
    def __init__(self):
        self.go_api = "https://api.geneontology.org/api/bioentity/function/"
        self.kegg_api = "https://rest.kegg.jp/link/pathway/hsa/"
        
    def get_gene_list_from_bed(self, bed_df: pd.DataFrame) -> List[str]:
        """从注释后的BED文件中提取基因列表"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        # 去重
        unique_genes = list(set(all_genes))
        print(f"提取到 {len(unique_genes)} 个唯一基因")
        return unique_genes
    
    def perform_go_enrichment(self, genes: List[str]) -> Dict:
        """执行GO功能富集分析"""
        enrichment_results = {}
        
        print("执行GO功能富集分析...")
        for gene in genes:
            try:
                # 获取基因的GO注释
                url = f"{self.go_api}{gene}"
                response = requests.get(url)
                if response.status_code == 200:
                    data = response.json()
                    enrichment_results[gene] = data
                time.sleep(0.1)
            except Exception as e:
                print(f"GO富集分析失败 {gene}: {e}")
        
        return enrichment_results
    
    def perform_kegg_enrichment(self, genes: List[str]) -> Dict:
        """执行KEGG通路富集分析"""
        kegg_results = {}
        
        print("执行KEGG通路富集分析...")
        for gene in genes:
            try:
                # 获取基因的KEGG通路注释
                url = f"https://rest.kegg.jp/link/pathway/hsa:{gene}"
                response = requests.get(url)
                if response.status_code == 200:
                    pathways = response.text.strip().split('\n')
                    kegg_results[gene] = [p.split('\t')[1] for p in pathways if p]
                time.sleep(0.1)
            except Exception as e:
                print(f"KEGG富集分析失败 {gene}: {e}")
        
        return kegg_results
    
    def save_enrichment_results(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """保存富集分析结果"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存GO结果
        if go_results:
            go_df = pd.DataFrame([
                {
                    'gene': gene,
                    'go_terms': json.dumps(data, indent=2)
                }
                for gene, data in go_results.items()
            ])
            go_df.to_csv(f"{output_dir}/go_enrichment_results.csv", index=False)
            print(f"GO富集分析结果已保存到: {output_dir}/go_enrichment_results.csv")
        
        # 保存KEGG结果
        if kegg_results:
            kegg_df = pd.DataFrame([
                {
                    'gene': gene,
                    'pathways': ';'.join(pathways)
                }
                for gene, pathways in kegg_results.items()
            ])
            kegg_df.to_csv(f"{output_dir}/kegg_enrichment_results.csv", index=False)
            print(f"KEGG富集分析结果已保存到: {output_dir}/kegg_enrichment_results.csv")

def main():
    """主函数"""
    # 示例BED文件内容（您可以替换为实际的BED文件路径）
    bed_content = """chr1\t36143270\t36143271
chr4\t156041918\t156041919
chr8\t41004246\t41004247
chr18\t31929457\t31929458
chr17\t71965860\t71965861
chr10\t67103337\t67103338
chr10\t80123580\t80123581
chr6\t29374583\t29374584
chr19\t5477150\t5477151
chr8\t25512704\t25512705"""
    
    # 保存示例BED文件
    with open('/workspace/data/example.bed', 'w') as f:
        f.write(bed_content)
    
    print("=== BED文件基因注释和富集分析 ===\n")
    
    # 初始化注释器和分析器
    annotator = GeneAnnotator()
    analyzer = GeneEnrichmentAnalyzer()
    
    # 1. 解析BED文件
    print("1. 解析BED文件...")
    bed_df = annotator.parse_bed_file('/workspace/data/example.bed')
    if bed_df is None:
        return
    
    print(f"BED文件内容预览:")
    print(bed_df.head())
    print()
    
    # 2. 注释基因
    print("2. 注释基因...")
    annotated_df = annotator.annotate_bed_regions(bed_df)
    
    # 保存注释结果
    annotator.save_annotated_results(annotated_df, '/workspace/data/annotated_genes.bed')
    print("注释结果:")
    print(annotated_df.head())
    print()
    
    # 3. 提取基因列表
    print("3. 提取基因列表...")
    gene_list = analyzer.get_gene_list_from_bed(annotated_df)
    print(f"提取的基因: {gene_list}")
    print()
    
    # 4. 执行富集分析（示例，仅分析前3个基因）
    print("4. 执行富集分析...")
    if len(gene_list) >= 3:
        test_genes = gene_list[:3]  # 只分析前3个基因作为示例
    else:
        test_genes = gene_list
    
    print(f"分析基因: {test_genes}")
    
    # GO富集分析
    go_results = analyzer.perform_go_enrichment(test_genes)
    
    # KEGG富集分析
    kegg_results = analyzer.perform_kegg_enrichment(test_genes)
    
    # 5. 保存结果
    print("5. 保存富集分析结果...")
    analyzer.save_enrichment_results(go_results, kegg_results, '/workspace/data')
    
    print("\n=== 分析完成 ===")
    print("结果文件:")
    print("- 注释结果: /workspace/data/annotated_genes.bed")
    print("- GO富集结果: /workspace/data/go_enrichment_results.csv")
    print("- KEGG富集结果: /workspace/data/kegg_enrichment_results.csv")

if __name__ == "__main__":
    main()