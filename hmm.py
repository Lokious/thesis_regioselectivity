"""
Run hmmscan
"""
import os
def main():
    seeds = ["PF05175_seed","PF08241_seed","PF08242_seed","PF13489_seed","PF13649_seed","PF13847_seed"]
    # for seed in seeds:
    #     buildfile = "data/hmm_out/{}_build".format(seed)
    #     os.system("hmmbuild {} data/hmm_out/top_ten_hits_exclude_nomethylrelated/{}.txt".format(buildfile,seed))
    # for seed in seeds:
    #     buildfile = "data/hmm_out/{}_build".format(seed)
    #     os.system("hmmpress {}".format(buildfile))
    for seed in seeds:
        inputfile="data/uniprot_ec2.1.1.fasta"
        buildfile = "data/hmm_out/{}_build".format(seed)
        o = "data/hmm_out/{}_hit.out".format(seed)
        domtblout="data/hmm_out/{}_hit.txt".format(seed)
        tblout="data/hmm_out/{}_hit.tsv".format(seed)
        os.system("nohup hmmscan -o {} --domtblout {} --tblout {} --cpu 2 {} {} &\n".format(o,domtblout,tblout,buildfile,inputfile))
if __name__ == "__main__":
    main()