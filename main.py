import os
import sys
from Bio import Entrez,SeqIO
from subprocess import *
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

batchSize=5000
retmax=10000
n_threads=4
window_size=15
n_species_gap=1


def read_config(sf_config):
    global batchSize
    global retmax
    global n_threads
    global window_size
    global n_species_gap

    with open(sf_config,"rU") as fin_config:
        for line in fin_config:
            fields=line.split()
            if fields[0]=="BatchSize":
                batchSize=int(fields[1])
            elif fields[0]=="RetMax":
                retmax=int(fields[1])
            elif fields[0]=="N_Threads":
                n_threads=int(n_threads)
            elif fields[0]=="WindowSize":
                window_size=int(fields[1])
            elif fields[0]=="N_Species_Gap":
                n_species_gap=int(fields[1])


def fetch_from_NCBI(term):
    global batchSize
    global retmax

    Entrez.email = "sally.chamberland@uconn.edu"     #Always tell NCBI who you are
    #batchSize = 5000
    #retmax=10000
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


def convert_2_fa(term):
    with open(term+".txt","rU") as fin_raw:
        with open("cleaved_"+term+".fa","w") as fout_fa:
            for line in fin_raw:
                fields=line.split()
                if len(fields)<2:
                    continue
                seq=fields[-1]
                name_id="_".join(fields[:-1])
                fout_fa.write(">"+name_id+"\n")
                fout_fa.write(seq+"\n")


def convert_list_2_fa(l_st, l_su):
    for term in l_st:
        convert_2_fa(term)
    for term in l_su:
        convert_2_fa(term)

def cleave(l_st, l_su, m_st_cleave, m_su_cleave):
    for term in l_st:
        sf_original=term+".fa"
        if os.path.exists(sf_original)==False:
            continue

        if m_st_cleave.has_key(term)==True:
            run_targetP(term)
            b_su=False
            sf_targetp_out=term+".targetP.out"
            parse_targetP_out(sf_targetp_out, b_su, sf_original)
        else:
            cmd="cp {0} {1}".format(sf_original, "cleaved_"+sf_original)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

    for term in l_su:
        sf_original=term+".fa"
        if os.path.exists(sf_original)==False:
            continue

        if m_su_cleave.has_key(term)==True:
            run_targetP(term)
            b_su=True
            sf_targetp_out=term+".targetP.out"
            parse_targetP_out(sf_targetp_out, b_su, sf_original)
        else:
            cmd="cp {0} {1}".format(sf_original, "cleaved_"+sf_original)
            Popen(cmd, shell = True, stdout = PIPE).communicate()


def run_targetP(term):
    cmd="/scratch/scratch2/chc12015/SCAProtein/targetp-1.1/targetp -N -c {0}.fa > {1}.targetP.out".format(term, term)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def parse_targetP_out(sf_targetp_out, b_su, sf_original):
    l_fa_name=[]
    l_fa_seq=[]
    for record in SeqIO.parse(sf_original, "fasta"):
        l_fa_name.append(str(record.id))
        l_fa_seq.append(str(record.seq))

    m_new_fa={}
    with open(sf_targetp_out,"rU") as fin_tgt:
        index=0
        for line in fin_tgt:
            fields=line.split()
            if len(fields)<8 or fields[0]=="Name":
                continue

            #Name                  Len            mTP     SP  other  Loc  RC  TPlen
            name=l_fa_name[index]
            tplen=fields[7]

            if tplen!="-":
                tlen=int(tplen)
                if b_su==True:
                     m_new_fa[name]=l_fa_seq[index][tlen:]
                else:
                    m_new_fa[name]=l_fa_seq[index][:tlen]
            index=index+1

    sf_new="cleaved_"+sf_original
    with open(sf_new,"w") as fout:
        for name in m_new_fa:
            fout.write(">"+name+"\n")
            fout.write(m_new_fa[name]+"\n")



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


def concatenate(l_st, l_su):
    for term_su in l_su:
        m_sunit={}
        sf_su="cleaved_"+term_su+".fa"

        if os.path.exists(sf_su):
            for record in SeqIO.parse(sf_su, "fasta"):
                m_sunit[str(record.id)]=str(record.seq)

        for term_st in l_st:
            m_stract={}
            sf_st="cleaved_"+term_st+".fa"

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
    #print cmd
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def multi_seq_alignment(l_st, l_su):
    global n_threads
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


def kick_alignment(l_st, l_su):
    for term_su in l_su:
        for term_st in l_st:
            sf_align="align_"+term_su+"_"+term_st+".fa"
            m_kick_out=pick_alignments(sf_align)
            sf_cat="concate_"+term_su+"_"+term_st+".fa"
            sf_cat_after_kick="concate_af_kick_"+term_su+"_"+term_st+".fa"

            with open(sf_cat_after_kick,"w") as fout_af_kick:
                cnt=0
                for record in SeqIO.parse(sf_cat, "fasta"):
                    if m_kick_out.has_key(cnt)==False:
                        fout_af_kick.write(">"+str(record.id)+"\n")
                        fout_af_kick.write(str(record.seq)+"\n")
                    cnt=cnt+1

def pick_alignments(sf_align):
    global window_size
    global n_species_gap

    l_name=[]
    l_seq=[]
    for record in SeqIO.parse(sf_align, "fasta"):
        l_name.append(str(record.id))
        l_seq.append(str(record.seq))

    n_species=len(l_name)
    if n_species<1:
        return
    seq_len=len(l_seq[0])

    label=""
    for i in range(window_size):
        label=label+"-"

    m_kick_out={}
    for i in range(seq_len-window_size):
        non_gap=0
        temp={}
        for j in range(n_species):
            if l_seq[j][i:i+window_size]!=label:
                temp[j]=1
                non_gap=non_gap+1
                if non_gap>n_species_gap:
                    break
        if non_gap<=n_species_gap:
            for key in temp:
                m_kick_out[key]=1

    return m_kick_out

def run_msa_round2(sf_id):
    sf_cat="concate_af_kick_"+sf_id
    sf_algn="align_final_"+sf_id
    cmd="./muscle3.8.31_i86linux64 -in {0} -out {1}".format(sf_cat, sf_algn)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def multi_seq_alignment_round2(l_st, l_su):
    global n_threads
    l_align=[]
    for term_su in l_su:
        for term_st in l_st:
            sf_cat="concate_af_kick_"+term_su+"_"+term_st+".fa"
            if os.path.exists(sf_cat):
                l_align.append(term_su+"_"+term_st+".fa")

    pool = Pool(n_threads)
    pool.map(run_msa_round2, l_align)
    pool.close()
    pool.join()

def clear():
    cmd="rm *.fa"
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def get_terms(sf_subtract, sf_subunit):
    l_subtract=[]
    l_subunit=[]
    with open(sf_subtract,"rU") as fin_st:
        for line in fin_st:
            line=line.rstrip()
            fields=line.split()
            l_subtract.append("-".join(fields))

    with open(sf_subunit,"rU") as fin_su:
        for line in fin_su:
            line=line.rstrip()
            fields=line.split()
            l_subunit.append("-".join(fields))

    return l_subtract, l_subunit

def get_terms_map(sf_subtract, sf_subunit):
    m_subtract={}
    m_subunit={}
    with open(sf_subtract,"rU") as fin_st:
        for line in fin_st:
            fields=line.split()
            m_subtract["-".join(fields)]=1

    with open(sf_subunit,"rU") as fin_su:
        for line in fin_su:
            fields=line.split()
            m_subunit["-".join(fields)]=1

    return m_subtract, m_subunit

def usage():
    print ""
    print 'Usage: python {0} Options subtrate_term_file subunit_term_file subtrate_cleave_file subunit_cleave_file' \
          '\n'.format(sys.argv[0]),
    print 'Options:\n',
    print '    fetch          Fetch from NCBI by give list\n',
    print '    convert        Convert plain to fasta\n',
    print '    cleave         Cleave pre-seq\n'
    print '    cat            Concatenate sequences\n',
    print '    align          Run multiple sequence alignment\n',
    print '    kick           Kick out the badly aligned sequences\n'
    print '    realign        Re-run multiple sequence alignment after kicking'
    print '    all            Run the whole pipeline\n',
    print '    clear          Clear all old files\n'


if __name__ == "__main__":
    if len(sys.argv) <= 5:
        usage()
        raise SystemExit

    option=sys.argv[1]
    sf_subtract=sys.argv[2]
    sf_subunit=sys.argv[3]
    sf_subtract_cleave=sys.argv[4]
    sf_subunit_cleave=sys.argv[5]

    read_config("./config.txt") ##read in config file

    if option=="fetch":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        prepare_database(l_st, l_su)
    elif option=="convert":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        #print l_st,l_su###########################################################
        convert_list_2_fa(l_st, l_su)
    elif option=="cleave":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        m_st_cleave, m_su_cleave=get_terms_map(sf_subtract_cleave, sf_subunit_cleave)
        cleave(l_st, l_su, m_st_cleave, m_su_cleave)
    elif option=="cat":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        concatenate(l_st, l_su)
    elif option=="align":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        multi_seq_alignment(l_st, l_su)
    elif option=="kick":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        kick_alignment(l_st, l_su)
    elif option=="realign":
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        multi_seq_alignment_round2(l_st, l_su)
    elif option=="all":
        clear()
        l_st, l_su=get_terms(sf_subtract, sf_subunit)
        prepare_database(l_st, l_su)
        cleave(l_st, l_su)
        concatenate(l_st, l_su)
        multi_seq_alignment(l_st, l_su)
        kick_alignment(l_st, l_su)
        multi_seq_alignment_round2(l_st, l_su)
    elif option=="clear":
        clear()
    else:
        print "Wrong option, please check!!!"
        usage()
