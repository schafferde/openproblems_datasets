import mygene
import scanpy as sc
import sys

mg = mygene.MyGeneInfo()

def extract_ensembl(ens_data):
    """Helper to pull the first gene ID from single dict or list of dicts."""
    if isinstance(ens_data, list):
        # Return first dict in list that contains a 'gene' key
        return next((item['gene'] for item in ens_data if 'gene' in item), None)
    elif isinstance(ens_data, dict):
        return ens_data.get('gene')
    return None

def bulk_convert_robust(genes, species='human'):
    # Initial batch query
    results = mg.querymany(genes, scopes='symbol,alias', 
                           fields='ensembl.gene,symbol', species=species)
    
    final_mapping = {}
    to_retry = {} # Map official_symbol -> original_query

    for res in results:
        query = res['query']
        
        # Priority 1: Direct Ensembl ID found
        ens_id = extract_ensembl(res.get('ensembl'))
        
        if ens_id:
            # Keep the first/highest score match
            if query not in final_mapping or final_mapping[query] is None:
                final_mapping[query] = ens_id
        else:
            # Priority 2: No Ensembl ID, but we found an official symbol (potential rename)
            official_symbol = res.get('symbol')
            if official_symbol and official_symbol != query:
                to_retry[official_symbol] = query
            
            # Ensure query exists in dict even if currently None
            if query not in final_mapping:
                final_mapping[query] = None

    # Secondary Pass: Lookup the official symbols for the "misses"
    if to_retry:
        retry_symbols = list(to_retry.keys())
        retry_results = mg.querymany(retry_symbols, scopes='symbol', 
                                     fields='ensembl.gene', species=species)
        
        for r_res in retry_results:
            r_ens_id = extract_ensembl(r_res.get('ensembl'))
            if r_ens_id:
                # Map it back to the original input name
                original_query = to_retry[r_res['query']]
                if final_mapping.get(original_query) is None:
                    final_mapping[original_query] = r_ens_id

    return final_mapping

# Execute
adata = sc.read_h5ad(sys.argv[1])
genes = list(adata.var.index)
if "muris" in sys.argv[1]:
    species = "mouse"
else:
    species = "human"
mapping = bulk_convert_robust(genes, species=species)

ids = []
missing = []
for gene in genes:
    id = mapping[gene]
    if id is None:
        ids.append("")
        missing.append(gene)
    else:
        ids.append(id)

print(f"Missing {len(missing)} genes: {missing}")

adata.var["gene_ids"] = ids
adata.write_h5ad(sys.argv[1].replace(".h5ad", "_ensid.h5ad"))

