from goatools.obo_parser import GODag
from goatools import obo_parser
from goatools.base import download_ncbi_associations
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
import mygene
import pandas as pd



def Return_GOA_stuff(): 
    #parse GO stuff
    gene2go_fname = download_ncbi_associations()
    objanno = Gene2GoReader(gene2go_fname ,taxids=[559292])#, evidence_set='IGI')  # Read data (only for Saccharomyces cerevisiae S288C)
    ns2assoc = objanno.get_ns2assc()

    #'gene2go.obo' was downloaded from http://current.geneontology.org/ontology/go-basic.obo on the 15.10.21 at 2:30pm
    #'go.obo' was downloaded from  http://current.geneontology.org/ontology/go.oboon the 23.11.21 at 10:00am

    # we us the optional_attrs option to load 'part_of' relationships, see:
    #https://github.com/tanghaibao/goatools/blob/main/notebooks/get_isa_and_partof.ipynb
    obo = obo_parser.GODag('Databases/go.obo', optional_attrs={'relationship'})
    return(ns2assoc,obo)



def Gene_Name_Convert(gene_list):

    ## fetch gene name conversion ################################################
    mg = mygene.MyGeneInfo()
    #create 


    ID_list  = mg.querymany(gene_list, scopes='all', fields='entrezgene', species=559292)  # S. cerevisiae S288C
    not_found=[]
    gene2entrez ={} # a dictionary storing gene name conversion
    entrez2gene ={}


    for D in ID_list:
        if 'entrezgene' in D.keys(): #check that a uniprot ID was found
            entrez2gene[D[ 'entrezgene']] =D['query']
        else:
            not_found.append(D['query'])

    print(len(not_found ),'ID were not found for: ', not_found )

    gene2entrez ={}
    #manually add wild type with made up Entrez
    gene2entrez['WT']='000000'
    entrez2gene['000000']='WT'
    for D in ID_list:
        if 'entrezgene' in D.keys(): #check that a uniprot ID was found
            gene2entrez[D['query']]=  D[ 'entrezgene']
            
    #convert to int
    
    return gene2entrez, entrez2gene
        
    

def return_GEA(gene_group, population, ns2assoc, obo, p_thresh=.05):
    '''Return a DF with the GO terms which p-values are below 0.05 
    also returns those terms as a set (go_set)'''

    goeaobj = GOEnrichmentStudyNS(
    population,  # Population (all mapped yeast genes in Entrez format)
    ns2assoc,  # geneid/GO associations
    obo,  # Ontologies
    propagate_counts = False,  # ???
    alpha = 0.05,  # default significance cut-off
    methods = ['fdr_bh']) 

    goea_results_all = goeaobj.run_study(gene_group, prt=None)
    top=0
    GOs=[]
    for term in goea_results_all:
        if term.p_fdr_bh <1:
            top+=1
        if term.p_fdr_bh < p_thresh: #if term.p_uncorrected <0.012:#
            top+=1

            GOs.append([term.GO,term.name, term.p_fdr_bh])

    GOs=pd.DataFrame(GOs, columns =('go','description','p'))
    GEA= GOs.iloc[GOs['p'].argsort()]

    go_set=GEA['go'].values
    GEA=GEA.drop_duplicates()
    #for some reason some terms are duplicated, we keep the first occurence
    i2keep = GEA['go'].duplicated(keep='first')==False
    GEA=GEA[i2keep ]

    return GEA,  go_set
