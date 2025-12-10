#!/usr/bin/env python3
"""
使用本地基因组注释数据库的BED文件基因注释工具
推荐使用pyensembl或本地GTF文件进行注释
"""

import pandas as pd
import numpy as np
import subprocess
import sys
from pathlib import Path
import time
from typing import List, Dict, Tuple, Optional

class LocalGeneAnnotator:
    """使用本地数据库的基因注释器"""
    
    def __init__(self, gtf_file: Optional[str] = None):
        self.gtf_file = gtf_file
        self.install_dependencies()
    
    def install_dependencies(self):
        """安装必要的依赖包"""
        packages = ['pyensembl', 'gseapy', 'pandas', 'numpy']
        for package in packages:
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
                print(f"成功安装 {package}")
            except subprocess.CalledProcessError:
                print(f"安装 {package} 失败，尝试继续...")
    
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
    
    def annotate_with_pyensembl(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """使用pyensembl进行基因注释"""
        try:
            from pyensembl import EnsemblRelease
            
            # 使用最新的人类基因组注释
            data = EnsemblRelease(75)  # 可以根据需要选择版本
            
            genes_list = []
            print("使用pyensembl进行基因注释...")
            
            for idx, row in bed_df.iterrows():
                chromosome = str(row['chromosome'])
                start = int(row['start'])
                end = int(row['end'])
                
                try:
                    # 查询重叠的基因
                    genes_at_locus = data.genes_at_locus(contig=chromosome, position=start)
                    
                    if genes_at_locus:
                        gene_names = [gene.gene_name for gene in genes_at_locus]
                        genes_list.append(';'.join(gene_names))
                    else:
                        genes_list.append('Unknown')
                        
                except Exception as e:
                    print(f"注释失败 {chromosome}:{start}-{end}: {e}")
                    genes_list.append('Unknown')
                
                if (idx + 1) % 10 == 0:
                    print(f"已处理 {idx + 1}/{len(bed_df)} 个区域")
                    time.sleep(0.1)  # 避免请求过快
            
            bed_df['genes'] = genes_list
            return bed_df
            
        except ImportError:
            print("pyensembl未安装或导入失败")
            return self.annotate_with_gtf(bed_df)
    
    def annotate_with_gtf(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """使用GTF文件进行注释（备用方法）"""
        print("使用GTF文件进行注释...")
        
        # 创建临时BED文件
        temp_bed = '/workspace/data/temp_regions.bed'
        bed_df.to_csv(temp_bed, sep='\t', index=False, header=False)
        
        genes_list = []
        
        try:
            # 使用bedtools进行注释（如果可用）
            result = subprocess.run([
                'bedtools', 'closest', '-a', temp_bed, 
                '-b', self.gtf_file, '-g', '/dev/null'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                # 解析bedtools输出
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if line:
                        fields = line.split('\t')
                        if len(fields) > 8:  # GTF文件通常有9列
                            gene_name = fields[8].split(';')[0].strip('"')
                            genes_list.append(gene_name)
                        else:
                            genes_list.append('Unknown')
            else:
                print("bedtools不可用，使用简单位置匹配")
                # 简单的位置匹配（仅作示例）
                for idx, row in bed_df.iterrows():
                    genes_list.append(f"Region_{row['chromosome']}_{row['start']}_{row['end']}")
                
        except FileNotFoundError:
            print("bedtools未安装，使用简单注释")
            for idx, row in bed_df.iterrows():
                genes_list.append(f"Region_{row['chromosome']}_{row['start']}_{row['end']}")
        
        bed_df['genes'] = genes_list
        return bed_df

class EnrichmentAnalyzer:
    """基因富集分析器"""
    
    def __init__(self):
        self.install_gseapy()
    
    def install_gseapy(self):
        """安装gseapy用于富集分析"""
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'gseapy'])
            print("成功安装gseapy")
        except subprocess.CalledProcessError:
            print("安装gseapy失败")
    
    def perform_go_enrichment(self, gene_list: List[str], organism: str = 'Human') -> Dict:
        """执行GO功能富集分析"""
        try:
            import gseapy as gp
            
            print("执行GO功能富集分析...")
            
            # GO富集分析
            go_enr = gp.prerank(rnk=gene_list, 
                              gene_sets='GO_Biological_Process_2023',
                              organism=organism,
                              min_size=5,
                              max_size=1000,
                              permutation_num=100)
            
            return {
                'go_bp': go_enr,
                'genes': gene_list
            }
            
        except ImportError:
            print("gseapy不可用，返回空结果")
            return {'error': 'gseapy not available'}
        except Exception as e:
            print(f"GO富集分析失败: {e}")
            return {'error': str(e)}
    
    def perform_kegg_enrichment(self, gene_list: List[str], organism: str = 'Human') -> Dict:
        """执行KEGG通路富集分析"""
        try:
            import gseapy as gp
            
            print("执行KEGG通路富集分析...")
            
            # KEGG富集分析
            kegg_enr = gp.prerank(rnk=gene_list,
                                 gene_sets='KEGG_2021_Human',
                                 organism=organism,
                                 min_size=5,
                                 max_size=1000,
                                 permutation_num=100)
            
            return {
                'kegg': kegg_enr,
                'genes': gene_list
            }
            
        except ImportError:
            print("gseapy不可用，返回空结果")
            return {'error': 'gseapy not available'}
        except Exception as e:
            print(f"KEGG富集分析失败: {e}")
            return {'error': str(e)}
    
    def save_enrichment_results(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """保存富集分析结果"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # 保存GO结果
        if 'go_bp' in go_results and not go_results.get('error'):
            try:
                go_results['go_bp'].to_csv(f"{output_dir}/go_enrichment_results.csv")
                print(f"GO富集分析结果已保存到: {output_dir}/go_enrichment_results.csv")
            except Exception as e:
                print(f"保存GO结果失败: {e}")
        
        # 保存KEGG结果
        if 'kegg' in kegg_results and not kegg_results.get('error'):
            try:
                kegg_results['kegg'].to_csv(f"{output_dir}/kegg_enrichment_results.csv")
                print(f"KEGG富集分析结果已保存到: {output_dir}/kegg_enrichment_results.csv")
            except Exception as e:
                print(f"保存KEGG结果失败: {e}")
        
        # 生成富集分析报告
        self.generate_enrichment_report(go_results, kegg_results, output_dir)
    
    def generate_enrichment_report(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """生成富集分析报告"""
        report_path = f"{output_dir}/enrichment_analysis_report.txt"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("基因富集分析报告\n")
            f.write("=" * 50 + "\n\n")
            
            # GO结果
            f.write("1. GO功能富集分析\n")
            f.write("-" * 30 + "\n")
            if 'go_bp' in go_results and not go_results.get('error'):
                try:
                    f.write("前10个显著富集的GO term:\n")
                    # 这里需要根据实际的gseapy输出格式调整
                    f.write("具体结果请查看 go_enrichment_results.csv\n")
                except Exception as e:
                    f.write(f"GO分析结果解析错误: {e}\n")
            else:
                f.write("GO富集分析失败\n")
            
            f.write("\n")
            
            # KEGG结果
            f.write("2. KEGG通路富集分析\n")
            f.write("-" * 30 + "\n")
            if 'kegg' in kegg_results and not kegg_results.get('error'):
                try:
                    f.write("前10个显著富集的KEGG通路:\n")
                    # 这里需要根据实际的gseapy输出格式调整
                    f.write("具体结果请查看 kegg_enrichment_results.csv\n")
                except Exception as e:
                    f.write(f"KEGG分析结果解析错误: {e}\n")
            else:
                f.write("KEGG富集分析失败\n")
        
        print(f"富集分析报告已保存到: {report_path}")

def main():
    """主函数"""
    print("=== BED文件基因注释和富集分析（本地版本）===\n")
    
    # 创建数据目录
    Path('/workspace/data').mkdir(exist_ok=True)
    
    # 示例BED文件内容
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
    
    # 初始化注释器和分析器
    annotator = LocalGeneAnnotator()
    analyzer = EnrichmentAnalyzer()
    
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
    try:
        annotated_df = annotator.annotate_with_pyensembl(bed_df)
        
        # 保存注释结果
        annotated_df.to_csv('/workspace/data/annotated_genes.bed', sep='\t', index=False)
        print("基因注释结果:")
        print(annotated_df)
        print()
        
        # 3. 提取基因列表
        print("3. 提取基因列表...")
        all_genes = []
        for genes_str in annotated_df['genes']:
            if genes_str != 'Unknown' and not genes_str.startswith('Region_'):
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        unique_genes = list(set(all_genes))
        print(f"提取的唯一基因: {unique_genes}")
        print()
        
        # 4. 执行富集分析
        print("4. 执行富集分析...")
        if len(unique_genes) >= 3:
            # GO富集分析
            go_results = analyzer.perform_go_enrichment(unique_genes)
            
            # KEGG富集分析
            kegg_results = analyzer.perform_kegg_enrichment(unique_genes)
            
            # 5. 保存结果
            print("5. 保存富集分析结果...")
            analyzer.save_enrichment_results(go_results, kegg_results, '/workspace/data')
        else:
            print("基因数量不足，跳过富集分析")
        
    except Exception as e:
        print(f"注释过程出错: {e}")
        print("可能需要安装额外的数据或依赖")
    
    print("\n=== 分析完成 ===")
    print("结果文件:")
    print("- 注释结果: /workspace/data/annotated_genes.bed")
    print("- 富集分析报告: /workspace/data/enrichment_analysis_report.txt")

if __name__ == "__main__":
    main()