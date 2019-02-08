PubChem_Actives_Inactives
=========================

Protein targets are listed on PubChem BioAssay record's as they are
provided by submitters. There are several work-arounds that need to be
performed to ensure all data is extracted.

Using Eutils, we link from UniProt Accessions to Genes (or start with
EGIDs) and identical protein products for that Uniprot (even the ones
not from UniProt). This way, regardless of whether the submitter used
UniProt IDs or GenPept IDs or RefSeq IDs or GeneIDs, we are able to find
BioAssays with a defined target the same as indicated by the UniProt
accession. We then link to the active or inactive compounds from the
BioAssays. We retrive the SMILES for the unique CIDs and map those back
to the links between the compounds and Uniprot, and finally calculate
the fingerprints for the SMILES.

This method uses the following resources:-
'Eutils: Entrez Programming Utilities' (http://www.ncbi.nlm.nih.gov/books/NBK25501/)
'PubChem PUG'(https://pubchem.ncbi.nlm.nih.gov/pug/pughelp.html)


This script is uses threading in order to achieve an acceptable
completion time

Dependencies : 'lxml', 'rdkit', 'Chemaxon JChem', 'OrderedDict'

--------------

Usage: PubChem_Inactives_Actives.py [options]
---------------------------------------------

Options:
  -h, --help            show this help message and exit
  -i FILE, --input=FILE
                        Input file (Uniprots [default] or EGIDs [must supply
                        the --geneids option])
  --geneids             Toggle Entrez Gene ID input format
  --actives             Toggle extract actives
  --inactives           Toggle extract inactives
