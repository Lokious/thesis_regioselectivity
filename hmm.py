"""
Run hmmscan
"""
import os
def main():
    seeds = ["PF05175","PF08241","PF08242","PF13489","PF13649","PF13847"]
    # for seed in seeds:
    #     buildfile = "data/hmm_out/{}_build".format(seed)
    #     os.system("hmmbuild {} data/hmm_out/top_ten_hits_exclude_nomethylrelated/{}_seed.txt".format(buildfile,seed))
    # for seed in seeds:
    #     buildfile = "data/hmm_out/{}_seed_build".format(seed)
    #     os.system("hmmpress {}".format(buildfile))
    # for seed in seeds:
    #     inputfile="data/uniprot_ec2.1.1.fasta"
    #     buildfile = "data/hmm_out/{}_seed_build".format(seed)
    #     o = "data/hmm_out/{}_hit.out".format(seed)
    #     domtblout="data/hmm_out/{}_seed_hit.txt".format(seed)
    #     tblout="data/hmm_out/{}_seed_hit.tsv".format(seed)
    #     os.system("nohup hmmscan -o {} --domtblout {} --tblout {} --cpu 2 {} {} &\n".format(o,domtblout,tblout,buildfile,inputfile))
    for seed in seeds:
        os.system("nohup hmmalign --amino --outformat A2M data/hmm_out/top_ten_hits_exclude_nomethylrelated/{0}.hmm data/{0}_rm.fasta > {0}.a2m&\n".format(seed))

if __name__ == "__main__":
    main()