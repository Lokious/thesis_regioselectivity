import os

import pandas as pd
import parse_data
def main():
    seq_domain_df = parse_data.read_hmmsearch_out(
        "../autodata/align/Pfam35.0uniprot_2_1_1_domout.tsv")
    domains = parse_data.sepreate_sequence_based_on_domain_without_overlap(seq_domain_df)
    for domain in domains:
        print(domain)
        os.system("hmmfetch ../autodata/align/Pfam35.0/Pfam-A.hmm {0} > {0}.hmm".format(domain))
if __name__ == "__main__":
    main()