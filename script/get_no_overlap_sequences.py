# coding=utf-8
import os
import pandas as pd
import parse_data
from sequence import Sequences
import time
from Model_class import Model_class
def hhalign(domian_list):
    """This function is to run hhalign to merge hmm and build new hmm for searching sequences"""

    template = domian_list.pop(0)
    log_file=open("../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/log_file_hhalign.txt",'a')
    while domian_list:
        domain_new = domian_list.pop(0)
        out_align = template.split(".")[0] + domain_new.split(".")[0]
        #use hhalign merge hmm model to MSA
        hhalign_cmd = "hhalign -i ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.a2m -t ../autodata/align/separate_by_domain/no_overlap_sequences/{1}_hmmalign_out_trim.a2m -oa3m ../autodata/align/separate_by_domain/no_overlap_sequences/{2}_hmmalign_out_trim.a2m".format(domain_new,
                                                               template,
                                                               out_align)
        os.system("wait")
        print("running command line:\n{}".format(hhalign_cmd))
        os.system(hhalign_cmd)
        log_file.write("{}\n".format(hhalign_cmd))
        log_file.write("\n")

        #build hmm model with MSA
        hmmbuild_cmd="hmmbuild ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm  ../autodata/align/separate_by_domain/no_overlap_sequences/{0}_hmmalign_out_trim.a2m".format(out_align)
        os.system("wait")
        print("running command line:\n{}".format(hmmbuild_cmd))
        os.system(hmmbuild_cmd)
        log_file.write("{}\n".format(hmmbuild_cmd))
        log_file.write("\n")

        #use hmmsearch (1) to search sequences with built hmm model
        hmmsearch_cmd = "hmmsearch --domT 15 -T 15 --domtblout ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim_domtblout.tsv ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm ../autodata/sequences/uniprot_ec2.1.1.fasta".format(out_align)
        os.system("wait")
        print("######running command line:\n{}######".format(hmmsearch_cmd))
        os.system(hmmsearch_cmd)
        log_file.write("{}\n".format(hmmsearch_cmd))
        log_file.write("\n")

        #use hmmsearch (2
        # 0for pdb structure wfrom pdbaa
        hmmsearch_pdb_cmd = "hmmsearch ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm pdbaa > ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_pdb.tsv".format(out_align)
        os.system("wait")
        print("######running command line:\n{}######".format(hmmsearch_pdb_cmd))
        os.system(hmmsearch_pdb_cmd)
        log_file.write("{}\n".format(hmmsearch_pdb_cmd))
        log_file.write("\n")

        #save sequences got from hmmsearch (1)
        parse_data.get_fasta_file_from_hmmsearch_hit("../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim_domtblout.tsv".format(out_align),out_align)
        seq = Sequences()

        #remove high similarity sequences
        ##run mmseqs
        mmseqs_build_database="mmseqs createdb ../autodata/sequences/{0}.fasta ../autodata/sequences/{0}_db --createdb-mode 1".format(out_align)
        mmseqs_map = "mmseqs map ../autodata/sequences/{0}_db ../autodata/sequences/{0}_db ../autodata/sequences/{0}_map_result tmp".format(out_align)
        mmseqs_covert_to_tab = "mmseqs convertalis ../autodata/sequences/{0}_db ../autodata/sequences/{0}_db ../autodata/sequences/{0}_map_result ../autodata/sequences/{0}.tab".format(out_align)
        os.system(mmseqs_build_database)
        os.system("wait")
        os.system(mmseqs_map)
        os.system("wait")
        os.system(mmseqs_covert_to_tab)
        os.system("wait")
        ## remove sequences
        seq.remove_sequences_from_result_of_mmseqs(tab_file="../autodata/sequences/{0}.tab".format(out_align), seq_file="../autodata/sequences/{}.fasta".format(out_align))


        #drop too long and too short sequences
        seq.drop_sequences(
            sequences_file="../autodata/sequences/{}.fasta".format(out_align))
        file = open("../autodata/sequences/{}.fasta".format(out_align), "a")

        #add structure sequences to seq file

        structure_seq_entry=open("../autodata/sequences/rcsb_pdb_5WP4_{}.fasta".format(out_align)).readlines()[0].strip(">")
        structure_seq=open("../autodata/sequences/rcsb_pdb_5WP4_{}.fasta".format(out_align)).readlines()[1]
        print(structure_seq_entry)
        print(structure_seq)
        file.write(">{}".format(structure_seq_entry))
        file.write("{}".format(structure_seq))
        file.close()

        # use hmmalign align searched sequences to get MSA
        hmmalign_cmd = "hmmalign --amino --outformat clustal ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm ../autodata/sequences/{0}.fasta > ../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/{0}_hmmalign_out_{1}.aln".format(
            out_align,"pdb_5WP4")
        os.system("wait")
        print("######running command line:\n{}######".format(hmmalign_cmd))
        os.system(hmmalign_cmd)
        log_file.write("{}\n".format(hmmalign_cmd))
        log_file.write("\n")
        template = out_align
        #create protein encoding dataframe
        parse_data.use_atom_properties_for_sequences_encoding(file_name="../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/{}_hmmalign_out_pdb_5WP4.aln".format(out_align),
            group=out_align, file_format="clustal", start=0,
            structure_chain="5WP4_1|Chain", pdb_name="5wp4.pdb")
        #create input dataframe
        X = pd.read_csv(
            "../autodata/fingerprint/fingerprint_bit128_radius3_all_data_drop_atom_19_09.csv",
            header=0, index_col=0)
        add_dataframe = pd.read_csv(
            "../autodata/protein_encoding/active_site/{}_AA_properties_encoding.csv".format(out_align),
            header=0, index_col=0)
        add_dataframe["Entry"] = add_dataframe.index
        add_dataframe.reset_index(drop=True, inplace=True)
        print(add_dataframe)
        input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
        print(input_dataframe)
        print("remove NA")
        input_dataframe = input_dataframe.dropna(axis=0, how="any")
        print(input_dataframe)
        input_dataframe.to_csv(
            "../autodata/input_data/active_site/{}_ACS_bit128_3_remove_redundant.csv".format(out_align))

    #     #train model
    #     mo_del = Model_class()
    #     print(input_dataframe)
    #     X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
    #         input_dataframe)
    #
    #     X_train = X_train.drop(columns=["methyl_type"])
    #     X_test = X_test.drop(columns=["methyl_type"])
    #     y_train = y_train.drop(columns=["methyl_type"])
    #     y_test = y_test.drop(columns=["methyl_type"])
    #     # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
    #     #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
    #     model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
    #                              "active_site_128fg_bi_type_bond3_rf_{}_remove_redundant".format(
    #                                  "PF08241_PF03602_ACS"), i=0)
    # print(out_align)
    # log_file.close()
    # return out_align

def hmmsearch_for_no_overlap_sequence(domains):
    """
    This function is to use hmmsearch to search pdb structure against hmm model
    built by align sequences only hits from one domain to hmm

    :param domains: domain name list["PF08241.15","PF03602.18"...]
    :return:None
    """
    for domain in domains:
        print(domain)
        os.system("wait")
        os.system(
            "hmmfetch ../autodata/align/separate_by_domain/Pfam35.0/Pfam-A.hmm {0} > ../autodata/align/{0}.hmm".format(
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

def hmmsearch_for_all_domains() ->list :
    """
    Run hmmsearch for all sequences from ec2.1.1

    :return: list of domains
    """
    #hmmsearch uniprot 2.1.1 against all pfamA domains
    hmmsearch_cmd="bash hmmsearch_different_bit_score.sh"
    os.system(hmmsearch_cmd)
    seq_domain_df = parse_data.read_hmmsearch_out(
       "../autodata/align/different_version_pfam/Pfam35.0uniprot_2_1_1_domout.tsv")
    domains = parse_data.sepreate_sequence_based_on_domain_without_overlap(seq_domain_df)

    return domains
def main():
    #count number of sequences for most frequennt domains
    domains=hmmsearch_for_all_domains()
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
    #hmmsearch_for_no_overlap_sequence(domains.copy())
    # domains=["PF08241.15","PF03602.18"]
    # hhalign(domains.copy())
    #parse_data.save_sequences_from_hmmscan_result()
if __name__ == "__main__":
    main()