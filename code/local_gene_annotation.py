#!/usr/bin/env python3
"""
BED File Gene Annotation Tool Using Local Genome Annotation Database
Recommended to use pyensembl or local GTF file for annotation
"""

import pandas as pd
import numpy as np
import subprocess
import sys
from pathlib import Path
import time
from typing import List, Dict, Tuple, Optional

class LocalGeneAnnotator:
    """Gene annotator using local database"""
    
    def __init__(self, gtf_file: Optional[str] = None):
        self.gtf_file = gtf_file
        self.install_dependencies()
    
    def install_dependencies(self):
        """Install necessary dependencies"""
        packages = ['pyensembl', 'gseapy', 'pandas', 'numpy']
        for package in packages:
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
                print(f"Successfully installed {package}")
            except subprocess.CalledProcessError:
                print(f"Failed to install {package}, trying to continue...")
    
    def parse_bed_file(self, bed_file: str) -> pd.DataFrame:
        """Parse BED file"""
        try:
            bed_df = pd.read_csv(bed_file, sep='\t', header=None, 
                               names=['chromosome', 'start', 'end'])
            print(f"Successfully read BED file, total {len(bed_df)} regions")
            return bed_df
        except Exception as e:
            print(f"Failed to read BED file: {e}")
            return None
    
    def annotate_with_pyensembl(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """Perform gene annotation using pyensembl"""
        try:
            from pyensembl import EnsemblRelease
            
            # Use latest human genome annotation
            data = EnsemblRelease(75)  # Can select version as needed
            
            genes_list = []
            print("Performing gene annotation using pyensembl...")
            
            for idx, row in bed_df.iterrows():
                chromosome = str(row['chromosome'])
                start = int(row['start'])
                end = int(row['end'])
                
                try:
                    # Query overlapping genes
                    genes_at_locus = data.genes_at_locus(contig=chromosome, position=start)
                    
                    if genes_at_locus:
                        gene_names = [gene.gene_name for gene in genes_at_locus]
                        genes_list.append(';'.join(gene_names))
                    else:
                        genes_list.append('Unknown')
                        
                except Exception as e:
                    print(f"Annotation failed {chromosome}:{start}-{end}: {e}")
                    genes_list.append('Unknown')
                
                if (idx + 1) % 10 == 0:
                    print(f"Processed {idx + 1}/{len(bed_df)} regions")
                    time.sleep(0.1)  # Avoid too fast requests
            
            bed_df['genes'] = genes_list
            return bed_df
            
        except ImportError:
            print("pyensembl not installed or import failed")
            return self.annotate_with_gtf(bed_df)
    
    def annotate_with_gtf(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """Annotate using GTF file (backup method)"""
        print("Annotating using GTF file...")
        
        # Create temporary BED file
        temp_bed = '/workspace/data/temp_regions.bed'
        bed_df.to_csv(temp_bed, sep='\t', index=False, header=False)
        
        genes_list = []
        
        try:
            # Use bedtools for annotation (if available)
            result = subprocess.run([
                'bedtools', 'closest', '-a', temp_bed, 
                '-b', self.gtf_file, '-g', '/dev/null'
            ], capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse bedtools output
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if line:
                        fields = line.split('\t')
                        if len(fields) > 8:  # GTF file usually has 9 columns
                            gene_name = fields[8].split(';')[0].strip('"')
                            genes_list.append(gene_name)
                        else:
                            genes_list.append('Unknown')
            else:
                print("bedtools not available, using simple position matching")
                # Simple position matching (for example only)
                for idx, row in bed_df.iterrows():
                    genes_list.append(f"Region_{row['chromosome']}_{row['start']}_{row['end']}")
                
        except FileNotFoundError:
            print("bedtools not installed, using simple annotation")
            for idx, row in bed_df.iterrows():
                genes_list.append(f"Region_{row['chromosome']}_{row['start']}_{row['end']}")
        
        bed_df['genes'] = genes_list
        return bed_df

class EnrichmentAnalyzer:
    """Gene enrichment analyzer"""
    
    def __init__(self):
        self.install_gseapy()
    
    def install_gseapy(self):
        """Install gseapy for enrichment analysis"""
        try:
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'gseapy'])
            print("Successfully installed gseapy")
        except subprocess.CalledProcessError:
            print("Failed to install gseapy")
    
    def perform_go_enrichment(self, gene_list: List[str], organism: str = 'Human') -> Dict:
        """Perform GO functional enrichment analysis"""
        try:
            import gseapy as gp
            
            print("Performing GO functional enrichment analysis...")
            
            # GO enrichment analysis
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
            print("gseapy not available, returning empty result")
            return {'error': 'gseapy not available'}
        except Exception as e:
            print(f"GO enrichment analysis failed: {e}")
            return {'error': str(e)}
    
    def perform_kegg_enrichment(self, gene_list: List[str], organism: str = 'Human') -> Dict:
        """Perform KEGG pathway enrichment analysis"""
        try:
            import gseapy as gp
            
            print("Performing KEGG pathway enrichment analysis...")
            
            # KEGG enrichment analysis
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
            print("gseapy not available, returning empty result")
            return {'error': 'gseapy not available'}
        except Exception as e:
            print(f"KEGG enrichment analysis failed: {e}")
            return {'error': str(e)}
    
    def save_enrichment_results(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """Save enrichment analysis results"""
        Path(output_dir).mkdir(exist_ok=True)
        
        # Save GO results
        if 'go_bp' in go_results and not go_results.get('error'):
            try:
                go_results['go_bp'].to_csv(f"{output_dir}/go_enrichment_results.csv")
                print(f"GO enrichment analysis results saved to: {output_dir}/go_enrichment_results.csv")
            except Exception as e:
                print(f"Failed to save GO results: {e}")
        
        # Save KEGG results
        if 'kegg' in kegg_results and not kegg_results.get('error'):
            try:
                kegg_results['kegg'].to_csv(f"{output_dir}/kegg_enrichment_results.csv")
                print(f"KEGG enrichment analysis results saved to: {output_dir}/kegg_enrichment_results.csv")
            except Exception as e:
                print(f"Failed to save KEGG results: {e}")
        
        # Generate enrichment analysis report
        self.generate_enrichment_report(go_results, kegg_results, output_dir)
    
    def generate_enrichment_report(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """Generate enrichment analysis report"""
        report_path = f"{output_dir}/enrichment_analysis_report.txt"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("Gene Enrichment Analysis Report\n")
            f.write("=" * 50 + "\n\n")
            
            # GO results
            f.write("1. GO Functional Enrichment Analysis\n")
            f.write("-" * 30 + "\n")
            if 'go_bp' in go_results and not go_results.get('error'):
                try:
                    f.write("Top 10 significantly enriched GO terms:\n")
                    # This needs to be adjusted according to the actual gseapy output format
                    f.write("See go_enrichment_results.csv for details\n")
                except Exception as e:
                    f.write(f"GO analysis result parsing error: {e}\n")
            else:
                f.write("GO enrichment analysis failed\n")
            
            f.write("\n")
            
            # KEGG results
            f.write("2. KEGG Pathway Enrichment Analysis\n")
            f.write("-" * 30 + "\n")
            if 'kegg' in kegg_results and not kegg_results.get('error'):
                try:
                    f.write("Top 10 significantly enriched KEGG pathways:\n")
                    # This needs to be adjusted according to the actual gseapy output format
                    f.write("See kegg_enrichment_results.csv for details\n")
                except Exception as e:
                    f.write(f"KEGG analysis result parsing error: {e}\n")
            else:
                f.write("KEGG enrichment analysis failed\n")
        
        print(f"Enrichment analysis report saved to: {report_path}")

def main():
    """Main function"""
    print("=== BED File Gene Annotation and Enrichment Analysis (Local Version)===\n")
    
    # Create data directory
    Path('/workspace/data').mkdir(exist_ok=True)
    
    # Example BED file content
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
    
    # Save example BED file
    with open('/workspace/data/example.bed', 'w') as f:
        f.write(bed_content)
    
    # Initialize annotator and analyzer
    annotator = LocalGeneAnnotator()
    analyzer = EnrichmentAnalyzer()
    
    # 1. Parse BED file
    print("1. Parse BED file...")
    bed_df = annotator.parse_bed_file('/workspace/data/example.bed')
    if bed_df is None:
        return
    
    print("BED文件内容:")
    print(bed_df)
    print()
    
    # 2. Annotate genes
    print("2. Annotate genes...")
    try:
        annotated_df = annotator.annotate_with_pyensembl(bed_df)
        
        # Save annotation results
        annotated_df.to_csv('/workspace/data/annotated_genes.bed', sep='\t', index=False)
        print("Gene annotation results:")
        print(annotated_df)
        print()
        
        # 3. Extract gene list
        print("3. Extract gene list...")
        all_genes = []
        for genes_str in annotated_df['genes']:
            if genes_str != 'Unknown' and not genes_str.startswith('Region_'):
                genes = genes_str.split(';')
                all_genes.extend(genes)
        
        unique_genes = list(set(all_genes))
        print(f"Extracted unique genes: {unique_genes}")
        print()
        
        # 4. Perform enrichment analysis
        print("4. Perform enrichment analysis...")
        if len(unique_genes) >= 3:
            # GO enrichment analysis
            go_results = analyzer.perform_go_enrichment(unique_genes)
            
            # KEGG enrichment analysis
            kegg_results = analyzer.perform_kegg_enrichment(unique_genes)
            
            # 5. Save results
            print("5. Save enrichment analysis results...")
            analyzer.save_enrichment_results(go_results, kegg_results, '/workspace/data')
        else:
            print("Insufficient gene count, skipping enrichment analysis")
        
    except Exception as e:
        print(f"Annotation process error: {e}")
        print("Additional data or dependencies may need to be installed")
    
    print("\n=== Analysis completed ===")
    print("Result files:")
    print("- Annotation results: /workspace/data/annotated_genes.bed")
    print("- Enrichment analysis report: /workspace/data/enrichment_analysis_report.txt")

if __name__ == "__main__":
    main()