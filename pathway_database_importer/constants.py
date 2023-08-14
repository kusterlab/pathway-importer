ORGANISMS = [
    {'identifier': 'hsa', 'taxcode': 9606},
    {'identifier': 'eco', 'taxcode': 562},
    {'identifier': 'sce', 'taxcode': 4932},
    {'identifier': 'cel', 'taxcode': 6239},
    {'identifier': 'dme', 'taxcode': 7227},
    {'identifier': 'ath', 'taxcode': 3702},
    {'identifier': 'dre', 'taxcode': 7955},
    {'identifier': 'mmu', 'taxcode': 10090},
    {'identifier': 'mtu', 'taxcode': 1773},
    {'identifier': 'osa', 'taxcode': 4530},
]

KEGG_PATHWAY_LIST_ENDPOINT = 'http://rest.kegg.jp/list/pathway/{}'

KEGG_PATHWAY_KGML_GENERIC_ENDPOINT = 'http://rest.kegg.jp/get/{}/kgml'

KEGG_GENERIC_IDENTIFIER_LIST_ENDPOINT = 'http://rest.kegg.jp/list/{}'

UNIPROT_MAPPING_ENDPOINT = 'https://rest.uniprot.org/idmapping/run'

UNIPROT_RESULT_ENDPOINT = 'https://rest.uniprot.org/idmapping/stream/'

UNIPROT_EC_ENDPOINT = 'https://rest.uniprot.org/uniprotkb/search?fields=accession&query=((organism_id={}) AND (reviewed=true) AND (ec:{}))'

WIKIPATHWAYS_ALL_PATHWAYS_ZIPPED_ENDPOINTS_MAP = {
    'hsa': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Homo_sapiens.zip',
    'eco': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Escherichia_coli.zip',
    'sce': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Saccharomyces_cerevisiae.zip',
    'cel': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Caenorhabditis_elegans.zip',
    'dme': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Drosophila_melanogaster.zip',
    'ath': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Arabidopsis_thaliana.zip',
    'dre': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Danio_rerio.zip',
    'mmu': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Mus_musculus.zip',
    'mtu': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Mycobacterium_tuberculosis.zip',
    'osa': 'http://data.wikipathways.org/current/gpml/wikipathways-{}-gpml-Oryza_sativa.zip',
}

GENE_ID_DATABASES = ['EC', 'Ensembl', 'Ensembl_Transcript', 'GeneID', 'HGNC', 'GeneCards', 'KEGG', 'RefSeq_Protein',
                     'Uniprot']

LINK_TYPE_DICTIONARY = {
    # All link types that do not appear here are mapped to themselves
    'Arrow': 'activation',
    'LigandRound': 'binding/association',
    'LigandSquare': 'binding/association',
    'Receptor': 'binding/association',
    'ReceptorRound': 'binding/association',
    'ReceptorSquare': 'binding/association',
    'TBar': 'inhibition',
    'mim-binding': 'binding/association',
    'mim-catalysis': 'catalysis',
    'mim-cleavage': 'cleavage',
    'mim-conversion': 'conversion',
    'mim-covalent-bond': 'binding/association',
    'mim-gap': 'gap',
    'mim-inhibition': 'inhibition',
    'mim-modification': 'modification',
    'mim-necessary-stimulation': 'stimulation',
    'mim-stimulation': 'stimulation',
    'mim-transcription-translation': 'transcription/translation',
    'mim-translocation': 'translocation',
    'mim-branching-right': 'branching-right',
    'SBGN-Production': 'production',
    'SBGN-Catalysis': 'catalyis',
    'SBGN-Inhibition': 'inhibition',
    'indirect effect': 'indirect',
    'repression': 'inhibition',
    None: 'undefined',
}

LINK_LABEL_DICTIONARY = {
    'phosphorylation': '+p',
    'dephosphorylation': '-p',
    'ubiquitination': '+u',
    'glycosylation': '+g',
    'methylation': '+m',
}

NODE_TYPE_DICTIONARY = {
    # KEGG
    'compound': 'compound',
    'gene': 'gene_protein',
    'ortholog': 'gene_protein',
    'map': 'pathway',
    'group': 'group',
    # WIKIPATHWAYS
    'Complex': 'gene_protein',
    'Event': 'misc',
    'GeneProduct': 'gene_protein',
    'Label': 'misc',
    'Metabolite': 'compound',
    None: 'misc',
    'Pathway': 'pathway',
    'Pathways': 'pathway',
    'Protein': 'gene_protein',
    'RNA': 'gene_protein',
    'Rna': 'gene_protein'
}

GENE_NODETYPE = 'gene_protein'
