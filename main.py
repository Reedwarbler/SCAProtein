import os
import sys
from Bio import Entrez,SeqIO
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

def fetch_from_NCBI(term):
    Entrez.email = "chongchu.cs@gmail.com"     # Always tell NCBI who you are
    batchSize = 5000
    retmax=10000
    db="protein"

    #first get GI for query accesions
    fields=term.split("-")
    search_term=""
    if len(fields)==1:
        search_term=term
    else:
        search_term=" ".join(fields)
    handle = Entrez.esearch( db=db,term=search_term, retmax=retmax)
    giList = Entrez.read(handle)['IdList']

    #post NCBI query
    search_handle     = Entrez.epost(db=db, id=",".join(giList))
    search_results    = Entrez.read(search_handle)
    webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]

    #fecth all results in batch of batchSize entries at once
    filename = term+".gbk"
    if os.path.exists(filename):
        cmd="rm {0}".format(filename)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    out_handle = open(filename, "a")
    for start in range( 0 ,len(giList), batchSize):
        #fetch entries in batch
        print "Fetching start from ", start, "queries!!"
        net_handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
        out_handle.write(net_handle.read())
        print "Saved ", batchSize, "results!!!"
    out_handle.close()


def cvt_genebank_fq(term):
    gbk_file = term+".gbk"
    fa_file=term+".fa"

    m_fa={}
    if os.path.exists(gbk_file):
        records=SeqIO.parse(gbk_file,"genbank")
        for record in records:
            name,seq=rtn_species_seq(record)
            name_fields=name.split()
            name_id="_".join(name_fields)

            if m_fa.has_key(name_id)==False:
                m_fa[name_id]=seq

    with open(fa_file,"w") as fout_fa:
        for id in m_fa:
            fout_fa.write(">"+id+"\n")
            fout_fa.write(m_fa[id]+"\n")


def rtn_species_seq(record):
    for feature in record.features:
        if feature.type=="source":
            return str(feature.qualifiers["organism"][0]), str(record.seq)


def prepare_database(l_st, l_su):
    for term in l_st:
        print "Fetching sequences for ", term
        fetch_from_NCBI(term)
        cvt_genebank_fq(term)

    for term in l_su:
        print "Fetching sequences for ", term
        fetch_from_NCBI(term)
        cvt_genebank_fq(term)

#def targetP_cleave():

def concatenate(l_st, l_su):
    for term_su in l_su:
        m_sunit={}
        sf_su=term_su+".fa"

        if os.path.exists(sf_su):
            for record in SeqIO.parse(sf_su, "fasta"):
                m_sunit[str(record.id)]=str(record.seq)

        for term_st in l_st:
            m_stract={}
            sf_st=term_st+".fa"

            if os.path.exists(sf_st):
                for record in SeqIO.parse(sf_st, "fasta"):
                    m_stract[str(record.id)]=str(record.seq)

            sf_cat="concate_"+term_su+"_"+term_st+".fa"
            with open(sf_cat,"w") as fout_cat:
                for species in m_sunit:
                    if m_stract.has_key(species):
                        fout_cat.write(">"+species+"\n")
                        fout_cat.write(m_sunit[species]+m_stract[species]+"\n")

def run_msa(sf_id):
    sf_cat="concate_"+sf_id
    sf_algn="align_"+sf_id
    cmd="muscle3.8.31_i86linux64 -in {0} -out {1}".format(sf_cat, sf_algn)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def multi_seq_alignment(l_st, l_su, n_threads):
    l_align=[]
    for term_su in l_su:
        for term_st in l_st:
            sf_cat="concate_"+term_su+"_"+term_st+".fa"
            if os.path.exists(sf_cat):
                l_align.append(term_su+"_"+term_st+".fa")

    pool = Pool(n_threads)
    pool.map(run_msa, l_align)
    pool.close()
    pool.join()


def clear():
    cmd="rm *.fa"
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def get_terms(sf_subtract, sf_subunit):
    l_subtract=[]
    l_subunit=[]
    with open(sf_subtract) as fin_st:
        for line in fin_st:
            fields=line.split()
            l_subtract.append("-".join(fields))

    with open(sf_subunit) as fin_su:
        for line in fin_su:
            fields=line.split()
            l_subunit.append("-".join(fields))

    return l_subtract, l_subunit


def usage():
    print ""
    print 'Usage: python {0} Options subtrate_term_file subunit_term_file\n'.format(sys.argv[0]),
    print 'Options:\n',
    print '    fetch      Fetch from NCBI by give list\n',
    print '    cat            Concatenate sequences\n',
    print '    align          Run multiple sequence alignment\n',
    print '    all            Run the whole pipeline\n',
    print '    clear          Clear all old files\n'


if __name__ == "__main__":
    if len(sys.argv) <= 3:
        usage()
        raise SystemExit

    option=sys.argv[1]
    sf_subtract=sys.argv[2]
    sf_subunit=sys.argv[3]
    n_threads=4

    if option=="fetch":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        prepare_database(l_st, l_su)
    elif option=="cat":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        concatenate(l_st, l_su)
    elif option=="align":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        multi_seq_alignment(l_st, l_su, n_threads)
    elif option=="all":
        clear()
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        prepare_database(l_st, l_su)
        concatenate(l_st, l_su)
        multi_seq_alignment(l_st, l_su, n_threads)
    elif option=="clear":
        clear()
    else:
        print "Wrong option, please check!!!"
        usage()
