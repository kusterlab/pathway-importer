# Pathway Importer

This is a Python package that converts pathways from KEGG and WikiPathways into JSON format. These files can then be used with [`biowc-pathwaygraph`](https://www.npmjs.com/package/biowc-pathwaygraph).

## Installation
- Clone this repository: 
```sh
git clone https://github.com/jmueller95/pathway-importer.git
cd pathway-importer
```
- Install and activate the conda environment: 
```sh
conda env create -f ./environment.yml
conda activate pathway_importer
```
(or install the dependencies yourself)

## Usage
- Run the `main` script as follows:
```sh
cd pathway_database_importer
python main.py --download
```
This will create a folder `raw` next to `pathway_database_importer` and download the latest versions of all pathway diagrams from KEGG and WikiPathways.
The raw pathways are in KGML (KEGG) and GPML (WikiPathways) format. _Pathway Importer_ downloads automatically downloads diagrams of ten organisms:

- Homo sapiens    `(hsa)`  
- Escherichia coli    `(eco)`  
- Saccharomyces cerevisiae    `(sce)`  
- Caenorhabditis elegans    `(cel)`  
- Drosophila melanogaster    `(dme)`  
- Arabidopsis thaliana    `(ath)`  
- Danio rerio    `(dre)`  
- Mus musculus    `(mmu)`  
- Mycobacterium tuberculosis   `(mtu)`  
- Oryza sativa    `(osa)`  

If you don't want to download all of them, delete the respective line from `constants.py`.

_Pathway Importer_ then converts the pathways in the `raw` directory into JSON and stores the output in a folder `json`.  
If you have already obtained the raw pathway diagrams before and don't want to run the download, omit the `--download` parameter.
