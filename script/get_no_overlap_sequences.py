# coding=utf-8
import os
import pandas as pd
import parse_data
from sequence import Sequences
import time
import seaborn as sns
import matplotlib.pyplot as plt
from Model_class import Model_class
def hhalign(domian_list, bit_score:int=0, coverage:float=0):
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
        hmmsearch_cmd = "hmmsearch --domT {1} -T {1} --domtblout ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim_domtblout_{1}.tsv ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm ../autodata/sequences/uniprot_ec2.1.1.fasta".format(out_align,bit_score)
        os.system("wait")
        print("######running command line:\n{}######".format(hmmsearch_cmd))
        os.system(hmmsearch_cmd)
        log_file.write("{}\n".format(hmmsearch_cmd))
        log_file.write("\n")

        #use hmmsearch (2)
        # 0for pdb structure wfrom pdbaa
        hmmsearch_pdb_cmd = "hmmsearch ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm pdbaa > ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_pdb.tsv".format(out_align)
        os.system("wait")
        print("######running command line:\n{}######".format(hmmsearch_pdb_cmd))
        os.system(hmmsearch_pdb_cmd)
        log_file.write("{}\n".format(hmmsearch_pdb_cmd))
        log_file.write("\n")

        #save sequences got from hmmsearch (1) after remove redundant sequences
        parse_data.remove_not_fully_aligned_domain(hmmsearch_file="../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim_domtblout_{1}.tsv".format(out_align,bit_score), domain=out_align,
                                                   coverage=coverage, bit_score=bit_score)
        seq = Sequences()

        #remove high similarity sequences
        ##run mmseqs
        mmseqs_build_database="mmseqs createdb ../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db --createdb-mode 1".format(out_align, coverage, bit_score)
        mmseqs_map = "mmseqs map ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_map_result tmp".format(out_align, coverage, bit_score)
        mmseqs_covert_to_tab = "mmseqs convertalis ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_map_result ../autodata/sequences/{0}_coverage{1}_bit_score{2}.tab".format(out_align, coverage, bit_score)
        os.system(mmseqs_build_database)
        os.system("wait")
        os.system(mmseqs_map)
        os.system("wait")
        os.system(mmseqs_covert_to_tab)
        os.system("wait")
        ## remove sequences
        seq.remove_sequences_from_result_of_mmseqs(tab_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.tab".format(out_align, coverage, bit_score), seq_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(out_align, coverage, bit_score))

        #drop too long and too short sequences
        seq_dictionary=seq.drop_sequences(
            sequences_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(out_align, coverage, bit_score))
        print("seq number: {}".format(len(seq_dictionary.keys())))

        return len(seq_dictionary.keys()),out_align
    # print(out_align)
    log_file.close()
    # return out_align

def hmmalign_for_combination_of_domains(out_align, bit_score, coverage, sturcture="pdb_5WP4"):
    """
    Use the choosen bitscore result to build MSA  and for protein encoding
    :param out_align:
    :param bit_score:
    :return:
    """
    file = open("../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(out_align, coverage, bit_score), "a")

    # add structure sequences to seq file

    structure_seq_entry = open(
        "../autodata/sequences/rcsb_{}_{}.fasta".format(sturcture,
            out_align)).readlines()[0].strip(">")
    structure_seq =  open(
        "../autodata/sequences/rcsb_{}_{}.fasta".format(sturcture,
            out_align)).readlines()[1]
    print(structure_seq_entry)
    print(structure_seq)
    file.write(">{}".format(structure_seq_entry))
    file.write("{}".format(structure_seq))
    file.close()
    # use hmmalign align searched sequences to get MSA
    print("#########hmmalign################")
    hmmalign_cmd = "hmmalign --amino --outformat clustal ../autodata/align/separate_by_domain/no_overlap_sequences/hhalign/{0}_hmmalign_out_trim.hmm ../autodata/sequences/{0}_coverage{3}_bit_score{2}.fasta > ../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/{0}_hmmalign_out_{1}_bitscore{2}_coverage{3}.aln".format(
        out_align, sturcture,bit_score,coverage)
    os.system("wait")
    print("######running command line:\n{}######".format(hmmalign_cmd))
    os.system(hmmalign_cmd)

def protein_encoding_and_trainning(out_align,sturcture,structure_chain="5WP4_1|Chain",pdb_name="5wp4.pdb",start_point:int=1,bit_score=0, coverage=0.0,):

    os.system("wait")
    # create protein encoding dataframe
    add_dataframe=parse_data.use_atom_properties_for_sequences_encoding(
        file_name="../autodata/align/separate_by_domain/no_overlap_sequences/hmmalign/{0}_hmmalign_out_{1}_bitscore{2}_coverage{3}.aln".format(
            out_align,sturcture,bit_score,coverage),
        group="{}_bit_score{}_coverage{}".format(out_align,bit_score,coverage), file_format="clustal", start=start_point,
        structure_chain=structure_chain, pdb_name=pdb_name)
    # create input dataframe
    X = pd.read_csv(
        "../autodata/fingerprint/fingerprint_bit128_radius3_all_data_drop_atom_19_09.csv",
        header=0, index_col=0)
    # add_dataframe = pd.read_csv(
    #     "../autodata/protein_encoding/active_site/{}_AA_properties_encoding.csv".format(
    #         out_align),
    #     header=0, index_col=0)
    add_dataframe["Entry"] = add_dataframe.index
    add_dataframe.reset_index(drop=True, inplace=True)
    print(add_dataframe)
    input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
    print(input_dataframe)
    print("remove NA")
    input_dataframe = input_dataframe.dropna(axis=0, how="any")
    print(input_dataframe)
    input_dataframe.to_csv(
        "../autodata/input_data/active_site/{}_ACS_bit128_3_remove_redundant.csv".format(
            "{}_bit_score{}_coverage{}".format(out_align,bit_score,coverage)))
    #train model
    mo_del = Model_class()
    print(input_dataframe)
    X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
        input_dataframe)

    X_train = X_train.drop(columns=["methyl_type"])
    X_test = X_test.drop(columns=["methyl_type"])
    y_train = y_train.drop(columns=["methyl_type"])
    y_test = y_test.drop(columns=["methyl_type"])
    # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
    #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
    model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                             "active_site_128fg_bi_type_bond3_rf_{}_ACS_remove_redundant_{}_{}".format(
                                 out_align,bit_score, str(int(coverage*100)),i=0))


def hmmsearch_for_sequence_and_structure(domains, coverage, bit_score,sturcture,structure_chain="5WP4_1|Chain",pdb_name="5wp4.pdb",start_point:int=1,num_bit=128):
    """
    This function is to use hmmsearch to search pdb structure against hmm model
    built by align sequences from the domain to hmm

    :param domains: domain name list["PF08241.15","PF03602.18"...]
    :return:None
    """
    domian_dict={}
    for domain in domains:
        print(domain)
        domain_full=domain
        domain = domain.split(".")[0]
        #hmmfetch for domain's hmm
        os.system(
            "hmmfetch ../autodata/align/different_version_pfam/Pfam35.0/Pfam-A.hmm {0} > ../autodata/align/{1}.hmm".format(
                domain_full,domain))
        os.system("wait")

        #hmmsearch for separate domain
        hmmsearch_cmd = "hmmsearch --domT {1} -T {1} --domtblout ../autodata/align/separate_by_domain/{0}_hmmalign_out_trim_domtblout_{1}.tsv ../autodata/align/{0}.hmm ../autodata/sequences/uniprot_ec2.1.1.fasta".format(
            domain, bit_score)
        os.system(hmmsearch_cmd)
        os.system("wait")

        #save sequences got from hmmsearch (1) after remove redundant sequences
        parse_data.remove_not_fully_aligned_domain(hmmsearch_file="../autodata/align/separate_by_domain/{0}_hmmalign_out_trim_domtblout_{1}.tsv".format(domain,bit_score), domain=domain,
                                                   coverage=coverage, bit_score=bit_score)
        seq = Sequences()

        ## remove high similarity sequences save them under the directory sequences/
        ### run mmseqs
        mmseqs_build_database="mmseqs createdb ../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db --createdb-mode 1".format(domain, coverage, bit_score)
        mmseqs_map = "mmseqs map ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_map_result tmp".format(domain, coverage, bit_score)
        mmseqs_covert_to_tab = "mmseqs convertalis ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_db ../autodata/sequences/{0}_coverage{1}_bit_score{2}_map_result ../autodata/sequences/{0}_coverage{1}_bit_score{2}.tab".format(domain, coverage, bit_score)
        os.system(mmseqs_build_database)
        os.system("wait")
        os.system(mmseqs_map)
        os.system("wait")
        os.system(mmseqs_covert_to_tab)
        os.system("wait")

        ### remove sequences
        seq.remove_sequences_from_result_of_mmseqs(tab_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.tab".format(domain, coverage, bit_score), seq_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(domain, coverage, bit_score))

        ## drop too long and too short sequences
        seq_dictionary=seq.drop_sequences(
            sequences_file="../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(domain, coverage, bit_score))
        print("{} seq number: {}".format(domain,len(seq_dictionary.keys())))
        domian_dict[domain]=len(seq_dictionary.keys())

        hmmalign_sequence_tohmm = 'hmmalign --amino  --outformat clustal ../autodata/align/{0}.hmm ../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta > ../autodata/align/separate_by_domain/{0}_hmmalign_out_bitscore{2}_coverage{1}.aln'.format(
            domain, coverage, bit_score)
        hmmbuild_with_sequences_from_MSA = "hmmbuild ../autodata/align/separate_by_domain/{0}_hmmalign_out_trim_bit_score{2}_coverage{1}.hmm ../autodata/align/separate_by_domain/{0}_hmmalign_out_bitscore{2}_coverage{1}.aln".format(
            domain, coverage, bit_score)

        use_built_hmm_searching_for_pdb_structure = "hmmsearch --noali ../autodata/align/separate_by_domain/{0}_hmmalign_out_trim_bit_score{2}_coverage{1}.hmm pdbaa > ../autodata/align/separate_by_domain/pdb_{0}_trim_bitscore{2}_coverage{1}.tsv".format(
            domain, coverage, bit_score)

        print("runing: {}".format(hmmalign_sequence_tohmm))
        os.system(hmmalign_sequence_tohmm)
        os.system("wait")

        print("runing: {}".format(hmmbuild_with_sequences_from_MSA))
        os.system(hmmbuild_with_sequences_from_MSA)
        os.system("wait")

        print("runing: {}".format(use_built_hmm_searching_for_pdb_structure))
        os.system(use_built_hmm_searching_for_pdb_structure)
        os.system("wait")

        ###add structure sequences to alignment
        file = open(
            "../autodata/sequences/{0}_coverage{1}_bit_score{2}.fasta".format(
                domain, coverage, bit_score), "a")

        # add structure sequences to seq file

        structure_seq_entry = open(
            "../autodata/sequences/rcsb_{}_{}.fasta".format(sturcture,
                                                            domain)).readlines()[
            0].strip(">")
        structure_seq = open(
            "../autodata/sequences/rcsb_{}_{}.fasta".format(sturcture,
                                                            domain)).readlines()[
            1]
        print(structure_seq_entry)
        print(structure_seq)
        file.write(">{}".format(structure_seq_entry))
        file.write("{}".format(structure_seq))
        file.close()
        # use hmmalign align searched sequences to get MSA
        print("#########hmmalign################")
        hmmalign_cmd = "hmmalign --amino --outformat clustal ../autodata/align/separate_by_domain/{0}_hmmalign_out_trim.hmm ../autodata/sequences/{0}_coverage{3}_bit_score{2}.fasta > ../autodata/align/separate_by_domain/{0}_hmmalign_out_{1}_bitscore{2}_coverage{3}.aln".format(
            domain, sturcture, bit_score, coverage)
        os.system("wait")
        print("######running command line:\n{}######".format(hmmalign_cmd))
        os.system(hmmalign_cmd)

        os.system("wait")
        # create protein encoding dataframe

        add_dataframe = parse_data.use_atom_properties_for_sequences_encoding(
            file_name="../autodata/align/separate_by_domain/{0}_hmmalign_out_{1}_bitscore{2}_coverage{3}.aln".format(
                domain, sturcture, bit_score, coverage),
            group="{}_bit_score{}_coverage{}".format(domain, bit_score,
                                                     coverage),
            file_format="clustal", start=start_point,
            structure_chain=structure_chain, pdb_name=pdb_name)
        # create input dataframe
        # X = pd.read_csv(
        #     "../autodata/fingerprint/fingerprint_bit128_radius3_all_data_drop_atom_19_09.csv",
        #     header=0, index_col=0)
        X = pd.read_csv(
            "../autodata/fingerprint/MACCS_fingerprint_bit167_radius3_all_data.csv",
            header=0, index_col=0)


        add_dataframe["Entry"] = add_dataframe.index
        add_dataframe.reset_index(drop=True, inplace=True)
        print(add_dataframe)
        input_dataframe = X.merge(add_dataframe, on="Entry", how="left")
        print(input_dataframe)
        print("remove NA")
        input_dataframe = input_dataframe.dropna(axis=0, how="any")
        print("before drop duplicate")
        print(input_dataframe)
        drop_duplicate_columns=list(input_dataframe.columns)
        drop_duplicate_columns.remove("molecular_id")
        drop_duplicate_columns.remove("atom_index")
        input_dataframe.drop_duplicates(subset=drop_duplicate_columns)
        print("after drop duplicate")
        print(input_dataframe)
        input_dataframe.to_csv(
            "../autodata/input_data/active_site/{}_ACS_bit128_3_remove_redundant_MACCS.csv".format(
                "{}_bit_score{}_coverage{}".format(domain, bit_score,
                                                   coverage)))
        # train model
        mo_del = Model_class()
        print(input_dataframe)
        X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
            input_dataframe)

        X_train = X_train.drop(columns=["methyl_type","molecular_id","atom_index"])

        #save x test for further analysis result
        y_test.to_csv(
            "../autodata/model/active_site_128fg_bi_type_bond3_rf_{}_ACS_remove_redundant_{}_{}_MACCS_y_test.csv".format(
                domain, bit_score,
                str(int(coverage * 100)), i=0))
        X_test.to_csv("../autodata/model/active_site_128fg_bi_type_bond3_rf_{}_ACS_remove_redundant_{}_{}_MACCS_X_test.csv".format(domain, bit_score,
                                     str(int(coverage * 100)), i=0))
        X_test = X_test.drop(columns=["methyl_type", "molecular_id","atom_index"])
        # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
        #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
        print(X_train)
        model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                                 "active_site_128fg_bi_type_bond3_rf_{}_ACS_remove_redundant_{}_{}_MACCS".format(
                                     domain, bit_score,
                                     str(int(coverage * 100)), i=0))

    return domian_dict

def hmmsearch_for_all_domains() -> list:
    """
    Run hmmsearch for all sequences from ec2.1.1

    :return: list of domains
    """

    #hmmsearch uniprot 2.1.1 against all pfamA domains
    hmmsearch_cmd="bash hmmsearch_different_bit_score.sh"
    os.system(hmmsearch_cmd)
    seq_domain_df = parse_data.read_hmmsearch_out(
       "../autodata/align/different_version_pfam/Pfam35.0uniprot_2_1_1_domout.tsv")
    domains,seq_df = parse_data.sepreate_sequence_based_on_domain_without_overlap(seq_domain_df)

    return domains


def separate_model_for_different_domain():
    """
    Build separate model for sequences from different domain

    """
    # seq_domain_df = parse_data.read_hmmsearch_out(
    #    "../autodata/align/different_version_pfam/Pfam35.0uniprot_2_1_1_domout.tsv")
    # domains, seq_df = parse_data.sepreate_sequence_based_on_domain_without_overlap(
    #     seq_domain_df)
    bit_score=[11,15]
    #coverage=[round(x*0.1,2) for x in list(range(1,11,2))]
    coverage=[0.5,0.6,0.7,0.8,0.9]
    '''
    seq_num_dict={}
    for score in bit_score:
        #seq_num_dict[score]={}
        for j in coverage:
            seq_num_dict[score][j]={}
            #search for sequences
            domains_dict=hmmsearch_for_sequence_and_structure(domains, coverage=j, bit_score=score,sturcture="pdb_5WP4")
            #seq_num_dict[score][j]=domains_dict

    print(seq_num_dict)
    '''
    for score in bit_score:
        for cover in coverage:
            hmmsearch_for_sequence_and_structure(["PF08241.15"], bit_score=score,
                                                coverage=cover,
                                                sturcture="pdb_5WP4")
    # for domain in domains:
    #     domain=domain.split(".")[0]
    #     seq_number_df = pd.DataFrame(index=coverage,columns=bit_score)
    #     for score in bit_score:
    #         for j in coverage:
    #             length=seq_num_dict[score][j][domain]
    #             seq_number_df.loc[j,score] = length
    #             print(domain)
    #             print(seq_number_df)
    #     seq_number_df=seq_number_df.astype(int)
    #     f, ax = plt.subplots(figsize=(20,20))
    #     ax=sns.heatmap(seq_number_df, annot=True, fmt='d')
    #     ax.set(xlabel='Bit score', ylabel='coverage')
    #     plt.savefig('heatmap_{}_.pdf'.format(domain), dpi=800)


def align_PF08241_PF01795_domain():
    domains=["PF08241.15","PF01795.22"]

    bit_score=[5,7,9,11,13,15,17,19,21]
    coverage=[round(x*0.1,2) for x in list(range(1,11,2))]
    #seq_number_df = pd.DataFrame(index=coverage,columns=bit_score)
    for score in bit_score:
        for j in coverage:
            # print("####bit_score = {}#####".format(score))
            # seq_number, combined_domains = hhalign(domains.copy(),
            #                                        bit_score=score, coverage=j)
            # print(seq_number)
            # hmmalign_for_combination_of_domains("PF08241PF01795", bit_score=score,
            #                                     coverage=j,
            #                                     sturcture="pdb_3TKA")
            # protein_encoding_and_trainning("PF08241PF01795",
            #                                sturcture="pdb_3TKA",
            #                                structure_chain="3TKA_1|Chain",
            #                                pdb_name="3tka.pdb",
            #                                start_point=-33, bit_score=score,
            #                                coverage=j)
            input_dataframe=pd.read_csv(
                "../autodata/input_data/active_site/{}_ACS_bit128_3_remove_redundant.csv".format(
                    "{}_bit_score{}_coverage{}".format('PF08241PF01795', score,
                                                       j)),header=0,index_col=0)
            # train model
            mo_del = Model_class()
            print(input_dataframe)
            X_train, X_test, y_train, y_test = mo_del.prepare_train_teat_data(
                input_dataframe)

            X_train = X_train.drop(columns=["methyl_type"])
            X_test = X_test.drop(columns=["methyl_type"])
            y_train = y_train.drop(columns=["methyl_type"])
            y_test = y_test.drop(columns=["methyl_type"])
            # model1 = mo_del.SVM(X_train, X_test, y_train, y_test,
            #                         "_input128fg_bi_type_bond2_svm{}".format(d1),i=0)
            model2 = mo_del.RF_model(X_train, X_test, y_train, y_test,
                                     "active_site_128fg_bi_type_bond3_rf_{}_ACS_remove_redundant_{}_{}".format(
                                         'PF08241_PF01795', score,
                                         str(int(j * 100)), i=0))
            # seq_number_df.loc[j,score]=seq_number
            # print(seq_number_df)
    # seq_number_df=seq_number_df.astype(int)
    # f, ax = plt.subplots(figsize=(20,20))
    # ax=sns.heatmap(seq_number_df, annot=True, fmt='d')
    # ax.set(xlabel='Bit score', ylabel='coverage')
    # plt.savefig('heatmap_PF08241_PF01795_on_filtered_uniprot.png', dpi=800)
    #

def main():
    #align_PF08241_PF01795_domain()
    separate_model_for_different_domain()
    #count number of sequences for most frequennt domains
    #domains=hmmsearch_for_all_domains()
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
    # domains=["PF08241.15"," PF01795.22"]
    # bit_score=[5,7,9,11,13,15,17,19,21]
    # coverage=[round(x*0.1,2) for x in list(range(1,11,1))]
    # seq_number_df = pd.DataFrame(index=bit_score,columns=coverage)
    # for score in bit_score:
    #     for j in coverage:
    #         print("####bit_score = {}#####".format(score))
    #         seq_number,combined_domains=hhalign(domains.copy(), score,
    #                                             coverage=j)
    #         seq_number_df.loc[score,j]=seq_number
    #         print(seq_number_df)
    # seq_number_df=seq_number_df.astype(int)
    # f, ax = plt.subplots(figsize=(20,20))
    # ax=sns.heatmap(seq_number_df, annot=True, fmt='d')
    # plt.savefig('heatmap.png', dpi=800)

if __name__ == "__main__":
    main()