#!/usr/bin/env python3
"""
BED File Gene Annotation and Enrichment Analysis Tool
Annotate genomic position information in BED format to gene names and
perform functional enrichment analysis
"""

import json
import time
from pathlib import Path
from typing import List, Dict

import pandas as pd
import requests

class GeneAnnotator:
    """Genomic Position Annotator"""

    def __init__(self):
        self.gencode_url = "https://rest.genenames.org/fetch/region/"
        self.header = {"Content-Type": "application/json"}

    def parse_bed_file(self, bed_file: str) -> pd.DataFrame:
        """Parse BED file"""
        try:
            # Read BED file, assuming no header
            bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                               names=['chromosome', 'start', 'end'])
            print(f"Successfully read BED file, {len(bed_df)} regions in total")
            return bed_df
        except Exception as e:
            print(f"Failed to read BED file: {e}")
            return None

    def annotate_with_gencode_api(self, chromosome: str, start: int, end: int) -> List[str]:
        """Annotate genes using GENCODE API"""
        genes = []
        try:
            # Build GENCODE API request
            region = f"{chromosome}:{start}:{end}"
            url = f"{self.gencode_url}{region}"

            response = requests.get(url, headers=self.header, timeout=30)
            if response.status_code == 200:
                data = response.json()
                docs = data.get('response', {}).get('docs', [])
                for doc in docs:
                    if 'symbol' in doc:
                        genes.append(doc['symbol'])
            time.sleep(0.1)  # Avoid requests being too fast
        except requests.exceptions.RequestException as e:
            print(f"API request failed {chromosome}:{start}-{end}: {e}")

        return genes

    def annotate_bed_regions(self, bed_df: pd.DataFrame) -> pd.DataFrame:
        """Annotate BED regions to genes"""
        genes_list = []

        print("Starting gene annotation...")
        for idx, row in bed_df.iterrows():
            chromosome = str(row['chromosome'])
            if not chromosome.startswith('chr'):
                chromosome = f"chr{chromosome}"
            genes = self.annotate_with_gencode_api(chromosome, row['start'], row['end'])
            genes_list.append(';'.join(genes) if genes else 'Unknown')

            if (idx + 1) % 10 == 0:
                print(f"Processed {idx + 1}/{len(bed_df)} regions")

        bed_df['genes'] = genes_list
        return bed_df

    def save_annotated_results(self, bed_df: pd.DataFrame, output_file: str):
        """Save annotation results"""
        bed_df.to_csv(output_file, sep='\t', index=False)
        print(f"Annotation results saved to: {output_file}")

class GeneEnrichmentAnalyzer:
    """Gene Enrichment Analyzer"""

    def __init__(self):
        self.go_api = "https://api.geneontology.org/api/bioentity/function/"
        self.kegg_api = "https://rest.kegg.jp/link/pathway/hsa/"

    def get_gene_list_from_bed(self, bed_df: pd.DataFrame) -> List[str]:
        """Extract gene list from annotated BED file"""
        all_genes = []
        for genes_str in bed_df['genes']:
            if genes_str != 'Unknown':
                genes = genes_str.split(';')
                all_genes.extend(genes)

        # Remove duplicates
        unique_genes = list(set(all_genes))
        print(f"Extracted {len(unique_genes)} unique genes")
        return unique_genes

    def perform_go_enrichment(self, genes: List[str]) -> Dict:
        """Perform GO functional enrichment analysis"""
        enrichment_results = {}

        print("Performing GO functional enrichment analysis...")
        for gene in genes:
            try:
                # Get GO annotations for the gene
                url = f"{self.go_api}{gene}"
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    enrichment_results[gene] = data
                time.sleep(0.1)
            except requests.exceptions.RequestException as e:
                print(f"GO enrichment analysis failed for {gene}: {e}")

        return enrichment_results

    def perform_kegg_enrichment(self, genes: List[str]) -> Dict:
        """Perform KEGG pathway enrichment analysis"""
        kegg_results = {}

        print("Performing KEGG pathway enrichment analysis...")
        for gene in genes:
            try:
                # Get KEGG pathway annotations for the gene
                url = f"https://rest.kegg.jp/link/pathway/hsa:{gene}"
                response = requests.get(url)
                if response.status_code == 200:
                    pathways = response.text.strip().split('\n')
                    kegg_results[gene] = [p.split('\t')[1] for p in pathways if p]
                time.sleep(0.1)
            except Exception as e:
                print(f"KEGG enrichment analysis failed for {gene}: {e}")

        return kegg_results

    def save_enrichment_results(self, go_results: Dict, kegg_results: Dict, output_dir: str):
        """Save enrichment analysis results"""
        Path(output_dir).mkdir(exist_ok=True)

        # Save GO results
        if go_results:
            go_df = pd.DataFrame([
                {
                    'gene': gene,
                    'go_terms': json.dumps(data, indent=2)
                }
                for gene, data in go_results.items()
            ])
            output_path = f"{output_dir}/go_enrichment_results.csv"
            go_df.to_csv(output_path, index=False)
            print(f"GO enrichment analysis results saved to: {output_path}")

        # Save KEGG results
        if kegg_results:
            kegg_df = pd.DataFrame([
                {
                    'gene': gene,
                    'pathways': ';'.join(pathways)
                }
                for gene, pathways in kegg_results.items()
            ])
            output_path = f"{output_dir}/kegg_enrichment_results.csv"
            kegg_df.to_csv(output_path, index=False)
            print(f"KEGG enrichment analysis results saved to: {output_path}")

def main():
    """Main function"""
    # Example BED file content (you can replace with actual BED file path)
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

    print("=== BED File Gene Annotation and Enrichment Analysis ===\n")

    # Initialize annotator and analyzer
    annotator = GeneAnnotator()
    analyzer = GeneEnrichmentAnalyzer()

    # 1. Parse BED file
    print("1. Parsing BED file...")
    bed_df = annotator.parse_bed_file('/workspace/data/example.bed')
    if bed_df is None:
        return

    print(f"BED file content preview:")
    print(bed_df.head())
    print()

    # 2. Annotate genes
    print("2. Annotating genes...")
    annotated_df = annotator.annotate_bed_regions(bed_df)

    # Save annotation results
    annotator.save_annotated_results(annotated_df, '/workspace/data/annotated_genes.bed')
    print("Annotation results:")
    print(annotated_df.head())
    print()

    # 3. Extract gene list
    print("3. Extracting gene list...")
    gene_list = analyzer.get_gene_list_from_bed(annotated_df)
    print(f"Extracted genes: {gene_list}")
    print()

    # 4. Perform enrichment analysis (example, only analyze first 3 genes)
    print("4. Performing enrichment analysis...")
    if len(gene_list) >= 3:
        test_genes = gene_list[:3]  # Only analyze first 3 genes as example
    else:
        test_genes = gene_list

    print(f"Analyzing genes: {test_genes}")

    # GO enrichment analysis
    go_results = analyzer.perform_go_enrichment(test_genes)

    # KEGG enrichment analysis
    kegg_results = analyzer.perform_kegg_enrichment(test_genes)

    # 5. Save results
    print("5. Saving enrichment analysis results...")
    analyzer.save_enrichment_results(go_results, kegg_results, '/workspace/data')

    print("\n=== Analysis Complete ===")
    print("Result files:")
    print("- Annotation results: /workspace/data/annotated_genes.bed")
    print("- GO enrichment results: /workspace/data/go_enrichment_results.csv")
    print("- KEGG enrichment results: /workspace/data/kegg_enrichment_results.csv")

# Add a final newline to the file
if __name__ == "__main__":
    main()