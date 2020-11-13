# Annotation of protease structures

* Rodrigo Ochoa, Mikhail Magnitov, Roman A Laskowski, Pilar Cossio, Janet M Thornton
* *An automated protocol for modelling peptide substrates to proteases*
* **Manuscript under the review**

## Purpose

The goal of these scripts is to provide an example workflow of the protease-peptide structure annotation prior to the modelling step described in [https://github.com/rochoa85/Modelling-Protease-Substrates](https://github.com/rochoa85/Modelling-Protease-Substrates). The first script is focused on the annotation of the proteases with necessary information from publicly available resources. The second script allows the assignment of the binding pockets to each substrate residue by using the annotation of the reference substrate and protease.

## Extra files required

To perform the annotation of the ligand and protein chains, two files from the PDBsum database are required:

1. rasmol.dfn file that contains the information on ligand chain and residues. This file is accessible for each PDB entry via https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/PDB_code/rasmol.dfn

2. grow.out file that contains the information on protein chains and residues bound to ligand. This file is accessible for each PDB entry via https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/PDB_code/grow.out

**NOTE:** Put these files into the annotation folder (see repository for an example).

To perform the annotation of the IDs and catalytic residues, two files from the PDBsum and M-CSA databases are required:

1. Annotation of PDB structures with EC numbers and UniProt IDs. This file can be downloaded from: [http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/pdb_chain_sp_ec](http://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/pdb_chain_sp_ec) (36.7 Mb)

2. Information about the catalytic residues, including the residues found by homology. This file can be downloaded from: [https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json](https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json) (295.4 Mb)

## Implementation

#### Step 1: Annotating the dataset
First, proteases are annotated with ligand and protein PDB chains, EC number, UniProt and MEROPS IDs. The search for EC number and UniProt ID is performed using the obtained PDBsum database file, and the search for MEROPS ID is performed with UniProt REST API. Second, the protein chains bound to the ligand are annotated with a list of the catalytic residues. This step is performed using the obtained M-CSA database file.

**NOTE:** If any protein chains are not annotated with catalytic residues, one can perform a homology search using the protease sequence in order to infer them.

The basic command line to run the annotation script is:

```
python 01_run_annotation.py --input INPUT [--output OUTPUT]

arguments:
  --input INPUT    Dataset with protease-ligand pairs
  --output OUTPUT  Path to save the annotated dataset
```

#### Step 2: Identification of pockets
To identify the pocket residues we use PDB coordinates of the reference substrate residues. The reference substrate should be stretched across pockets S4-S4’, or, failing that, the longest available. Protease residues with at least one atom within a defined distance from the substrate are considered to form a binding pocket. We then project this annotation of pockets to other substrates within the protease family.

**NOTE:** This analysis requires two additional steps. The first is selecting a reference structure. The second is superimposing the proteases of interest with a selected reference structure (we suggest using PDBeFold) and obtaining new superimposed coordinates of the structures.

## Examples

#### Step 1: Annotating the dataset

For this step we provided a sample data of 7 proteases from major protease families as an example in **sample_data_annotation.csv**. It includes a list of PDB codes followed by ligand sequence:

```
3tjv,PRO-THR-SER-TYR-ALA-GLY-ASP-ASP-SER
```

The following command runs the annotation: `python 01_run_annotation.py --input sample_data_annotation.csv`

After running this comand, the annotated dataset will be saved to **annotated_dataset.csv**. It will include to following fields:

```
PDB: 3tjv
Ligand: PRO-THR-SER-TYR-ALA-GLY-ASP-ASP-SER
Ligand_size: 9
Ligand_chain: B
Ligand_residues_number: 1, 2, 3, 4, 5, 6, 7, 8, 9
Protein_chains: A
Protein_residues: A100PHE, A144TYR, A150LEU, A174ASN, A192GLY, A193PHE, A194LYS, A195GLY, A196ASP, A197SER, A212SER, A213TYR, A214GLY, A215ASN, A216LYS, A218GLY, A32PHE, A34GLN, A41ARG, A42LYS, A43ARG, A44CYS, A59HIS, A75LYS
EC_number: 3.4.21.-
UniProt_ID: P20718
MEROPS_ID: S01.147
Catalytic_residues: A103ASN, A191THR, A194LYS, A195GLY, A196ASP, A197SER, A198GLY, A59HIS, A99ASN
```

#### Step 2: Identification of pockets

