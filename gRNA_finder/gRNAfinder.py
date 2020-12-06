# coding utf-8
import os
import subprocess

import pandas as pd
from Bio import SeqIO

from gRNA_finder.DeepCpf1_py3 import main as efficiency_estimator

cas_offinder = "./cas-offinder"
genomes_db = "genomes_db"
PAM = "TTKGT"
gRNAlength = 20
PAM_grna_pattern = "{pam}{N}".format(pam=PAM, N="N" * gRNAlength)
column_names = {0: "pattern", 1: "fasta_ID", 2: "coord", 3: "PAM_gRNA", 4: "strand", 5: "mismatch"}


class gRNAfinder(object):

    @staticmethod
    def precount_gf(gf, output_text, PAM):

        input_text = os.path.join(genomes_db, "db_queries",
                                  "{}_{}.precount.txt".format(genome_file.replace(os.sep, "_"), PAM))
        grna_list = [PAM_grna_pattern, ]
        return run_cas_offinder(gf, grna_list, input_text, output_text, mismatch=0)

    def __init__(self, PAM, genome_file=None, genome_file_precount=None):

        self.genome_file = genome_file
        self.genome_file_precount = genome_file_precount

        if not os.path.isfile(genome_file_precount):
            try:
                f = os.path.dirname(genome_file_precount)
                os.makedirs(f)
            except Exception as ex:
                pass
            genome_file_df = self.precount_gf(
                genome_file, genome_file_precount, PAM
            )

        else:
            genome_file_df = pd.read_csv(
                self.genome_file_precount,
                sep="\t", header=None
            )
            genome_file_df = genome_file_df.rename(
                columns=column_names)

        self.genome_file_df = genome_file_df
        self.PAM = PAM

    def findPAM(self, location):
        """
        по координатам в геноме найти вхождения PAM
        вернуть список guideRNA
        :return:
        """
        pams = self.genome_file_df[
            (self.genome_file_df["fasta_ID"] == location["fasta_ID"]) &
            (self.genome_file_df["coord"] < location["end"]) &
            (self.genome_file_df["coord"] >= location["start"])
            ]
        return pams

    def findOffTagets(self, location, gRNAdf):
        """
        поискать offtarget вхождения
        :param location:
        :param gRNAlist:
        :return:
        """
        genome_name = genome_file.replace(os.path.sep, "_").split(".")[0]
        file_pref = "{genome_name}_{PAM}".format(
            genome_name=genome_name,
            PAM=PAM
        )
        input_text = os.path.join("genomes_db", "db_queries", file_pref + ".input.txt")

        output_text = os.path.join("genomes_db", file_pref + ".txt")

        all_gRNAs = run_cas_offinder(self.genome_file, gRNAdf["PAM_gRNA"], input_text, output_text)
        off_gRNAs = all_gRNAs[
            ((all_gRNAs["fasta_ID"] != location["fasta_ID"]) |
             (all_gRNAs["coord"] > location["end"]) |
             (all_gRNAs["coord"] <= location["start"]))
        ]

        return off_gRNAs

    def estimate_efficiency(self, guide_df, genome_fasta, location):

        input_data = []

        for seq_record in genome_fasta:
            if seq_record.id == location["fasta_ID"].split()[0]:
                i = 1
                for index, grna in guide_df.iterrows():
                    start = grna["coord"]
                    end = start + len(grna["PAM_gRNA"])
                    if grna["strand"] == "+":
                        start = start - 3 if start > 3 else 0
                        end = end + 6
                        input_data.append(
                            [i, seq_record.seq[start:end], 0]
                        )
                    else:
                        start = start - 6 if start > 6 else 0
                        end = end + 3
                        input_data.append(
                            [i, seq_record.seq[start:end].reverse_complement(), 0]
                        )
                    i += 1

                break

        # chromatin_data_available = False
        # if chromatin_data_available:
        #     chromatin_dnase_data = get_chromatin_dnase_data()
        #     input_data = input_data

        scores_out_list = efficiency_estimator(input_data)

        scores = [score[0] for score in scores_out_list]

        guide_df["score"] = scores

        return guide_df


def run_cas_offinder(genome_file, grna_list, input_text, output_text, mismatch=1):
    """
    ./cas-offinder test.input.txt C test.output.txt
    :return:
    """

    with open(input_text, "w") as input:
        input.write("{}\n".format(
            os.path.join(genomes_db, genome_file))
        )
        input.write("{}\n".format(PAM_grna_pattern))

        for grna in grna_list:
            input.write("{grna} {mm}\n".format(grna=grna, mm=mismatch))

    subprocess.check_call(
        [
            cas_offinder,
            input_text,
            "C",
            output_text
        ]
    )

    if os.path.isfile(output_text):
        all_gRNAs = pd.read_csv(output_text, sep="\t", header=None)
        all_gRNAs = all_gRNAs.rename(columns=column_names)
        return all_gRNAs
    else:
        return pd.DataFrame(columns=column_names)


def get_chromatin_dnase_data():
    raise NotImplementedError()


if __name__ == "__main__":
    genome_file = os.path.join("hg38", "chr1.fna")
    genome_file_precount = "PAMS/{PAM}/hg38_chr1.txt".format(PAM=PAM)
    genome_file_precount_path = os.path.join(genomes_db, genome_file_precount)
    gf = gRNAfinder(PAM, genome_file=genome_file, genome_file_precount=genome_file_precount_path)

    # ASPM chr1:197,084,127-197,146,669
    location = {
        "fasta_ID": "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p13 Primary Assembly",
        "start": 197084127,
        "end": 197096669
        # "end": 197146669
    }

    our_gRNA = gf.findPAM(location)
    off_gRNA = gf.findOffTagets(location, our_gRNA)

    genome_fasta = SeqIO.parse(os.path.join(genomes_db, genome_file), "fasta")

    our_gRNA = gf.estimate_efficiency(our_gRNA, genome_fasta, location)
    # off_gRNA = gf.estimate_efficiency(off_gRNA, genome_fasta)


    # RESPONSE
    response_fasta_peace = "response/peace.fasta"
    response_gRNA = "response/gRNA.tsv"
    response_offtarget = "response/off_gRNA.tsv"

    with open(response_fasta_peace, "w") as f:
        genome_fasta = SeqIO.parse(os.path.join(genomes_db, genome_file), "fasta")
        for seq_record in genome_fasta:
            if seq_record.id == location["fasta_ID"].split()[0]:
                f.write(">{}:{}-{}\n".format(location["fasta_ID"], location["start"], location["end"]))
                f.write(str(seq_record.seq[location["start"]:location["end"]]))
                break
    our_gRNA.to_csv(response_gRNA, sep="\t", index=False, header=True)
    off_gRNA.to_csv(response_offtarget, sep="\t", index=False, header=True)
