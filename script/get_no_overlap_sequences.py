# coding=utf-8
import os
import pandas as pd
import parse_data


def hhalign(domian_list):
    """This function is to run hhalign to merge hmm and build new hmm for searching sequences"""

    template = domian_list.pop(0)
    log_file=open("../autodata/align/separate_by_domain/no_overlap_sequences/log_file_hhalign.txt",'a')
    while domian_list:
        domain_new = domian_list.pop(0)
        out_align = template + domain_new
        hhalign_cmd = "hhalign -i ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.a2m -t ../autodata/align/separate_by_domain/no_overlap_sequences/{1}_hmmalign_out_trim.a2m -oa3m ../autodata/align/separate_by_domain/no_overlap_sequences/{2}_hmmalign_out_trim.a2m".format(domain_new,
                                                               template,
                                                               out_align)
        os.system("wait")
        os.system(hhalign_cmd)
        log_file.write("{}\n".format(hhalign_cmd))
        log_file.write("\n")
        hmmsearch_cmd = "hmmsearch --domtblout {0}_hmmalign_out_trim_domtblout.tsv {0}_hmmalign_out_trim.hmm ../autodata/sequences/uniprot_ec2.1.1.fasta".format(out_align)
        os.system("wait")
        os.system(hmmsearch_cmd)
        log_file.write("{}\n".format(hmmsearch_cmd))
        log_file.write("\n")
        hmmsearch_pdb_cmd = "hmmsearch {0}_hmmalign_out_trim.hmm pdbaa >{0}_pdb.tsv".format(out_align)
        os.system("wait")
        os.system(hmmsearch_pdb_cmd)
        log_file.write("{}\n".format(hmmsearch_pdb_cmd))
        log_file.write("\n")
        template = out_align
    print(out_align)
    log_file.close()
    return out_align

def hmmsearch_for_no_overlap_sequence(domains):
    for domain in domains:
        print(domain)
        os.system("wait")
        os.system(
            "hmmfetch ../autodata/align/Pfam35.0/Pfam-A.hmm {0} > ../autodata/align/{0}.hmm".format(
                domain))
        os.system("wait")
        cmd1 = 'hmmalign --amino  --outformat A2M ../autodata/align/{0}.hmm ../autodata/sequences/no_overlap_sequences/{0}.fasta > ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out.a2m'.format(
            domain)
        cmd2 = "hmmbuild ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.hmm ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.a2m".format(
            domain)

        cmd3 = "hmmsearch --noali ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.hmm pdbaa > ../autodata/align/separate_by_domain/no_overlap_sequences/pdb_{0}_trim.hmmer".format(
            domain)
        print("runing: {}".format(cmd1))
        os.system(cmd1)
        # os.system("wait")
        # print("runing: {}".format(cmd2))
        # os.system(cmd2)
        # os.system("wait")
        # print("runing: {}".format(cmd3))
        # os.system(cmd3)


def main():
    #count number of sequences for most frequennt domains
    seq_domain_df = parse_data.read_hmmsearch_out(
       "../autodata/align/different_version_pfam/Pfam35.0uniprot_2_1_1_domout.tsv")
    domains = parse_data.sepreate_sequence_based_on_domain_without_overlap(seq_domain_df)
    # for domain in domains:
    #
    #     cmd_hmmfetch ="hmmfetch ../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm {0} > ../autodata/align/{0}.hmm".format(domain)
    #     os.system("wait")
    #     os.system(cmd_hmmfetch)
    #     name_description="".join(open("{0}.hmm".format(domain)).readlines()[:4])
    #     if ("DNA" or "RNA") in name_description.upper():
    #         print("rm {}.hmm".format(domain))
    #         os.system("wait")
    #         os.system("rm {0}.hmm".format(domain))
    #         domains.remove(domain)
    # print(domains)
    #use hmmsearch for closest pdb structure
    hmmsearch_for_no_overlap_sequence(domains.copy())
    #hhalign(domains.copy())
    #parse_data.save_sequences_from_hmmscan_result()
if __name__ == "__main__":
    main()