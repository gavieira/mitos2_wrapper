#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# mitofree.py

"""Mitochondrial genome assembly using public data"""

__author__ = "Gabriel Alves Vieira"
__contact__ = "gabrieldeusdeth@gmail.com"

import subprocess, os, shutil, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class mitoannotation():
    def __init__(self, fasta_file, refdir, gencode=2):
        self.fasta = fasta_file
        self.gbk = "{}.gbk".format(os.path.splitext(self.fasta)[0])
        self.results = "{}_mitos".format(os.path.splitext(self.fasta)[0])
        self.gencode = gencode
        self.refdir = refdir
        self.bed_file = "{}/result.bed".format(self.results)
        self.feat_dict = {
        "cox1": {"name" : "COX1", "product" : "cytochrome c oxidase subunit I"},
        "cox2": {"name": "COX2", "product": "cytochrome c oxidase subunit II"},
        "atp8": {"name": "ATP8", "product": "ATP synthase F0 subunit 8"},
        "atp6": {"name": "ATP6", "product": "ATP synthase F0 subunit 6"},
        "cox3": {"name": "COX3", "product": "cytochrome c oxidase subunit III"},
        "nad3": {"name":"ND3", "product": "NADH dehydrogenase subunit 3"},
        "nad5": {"name": "ND5", "product": "NADH dehydrogenase subunit 5"},
        "nad4l": {"name": "ND4L", "product": "NADH dehydrogenase subunit 4L"},
        "nad4": {"name":"ND4", "product": "NADH dehydrogenase subunit 4"},
        "nad6": {"name":"ND6", "product":"NADH dehydrogenase subunit 6"},
        "cob": {"name":"CYTB", "product":"cytochrome b"},
        "nad1": {"name":"ND1", "product":"NADH dehydrogenase subunit 1"},
        "nad2": {"name":"ND2", "product":"NADH dehydrogenase subunit 2"},
        "trnL2": {"name":"trnL2", "product":"tRNA-Leu"},
        "trnK": {"name":"trnK", "product":"tRNA-Lys"},
        "trnK": {"name":"trnK", "product":"tRNA-Lys"},
        "trnD": {"name":"trnD", "product":"tRNA-Asp"},
        "trnG": {"name":"trnG", "product":"tRNA-Gly"},
        "trnA": {"name":"trnA", "product":"tRNA-Ala"},
        "trnR": {"name":"trnR", "product":"tRNA-Arg"},
        "trnN": {"name":"trnN", "product":"tRNA-Asn"},
        "trnS1": {"name":"trnS1", "product":"tRNA-Ser"},
        "trnE": {"name":"trnE", "product":"tRNA-Glu"},
        "trnF": {"name":"trnF", "product":"tRNA-Phe"},
        "trnH": {"name":"trnH", "product":"tRNA-His"},
        "trnT": {"name":"trnT", "product":"tRNA-Thr"},
        "trnP": {"name":"trnP", "product":"tRNA-Pro"},
        "trnS2": {"name":"trnS2", "product":"tRNA-Ser"},
        "trnL1": {"name":"trnL1", "product":"tRNA-Leu"},
        "trnV": {"name":"trnV", "product":"tRNA-Val"},
        "trnM": {"name":"trnM" , "product":"tRNA-Met"},
        "trnI": {"name":"trnI", "product":"tRNA-Ile"},
        "trnQ": {"name":"trnQ", "product":"tRNA-Gln"},
        "trnW": {"name":"trnW", "product":"tRNA-Trp"},
        "trnC": {"name":"trnC", "product":"tRNA-Cys"},
        "trnY": {"name":"trnY", "product":"tRNA-Tyr"},
        "rrnL": {"name":"rrnL", "product":"large subunit ribosomal RNA"},
        "rrnS": {"name":"rrnS", "product":"small subunit ribosomal RNA"}
        }
        self.cds_list = [x for x in self.feat_dict.keys() if not x.startswith("trn") and not x.startswith("rrn")]
        self.trna_list = [x for x in self.feat_dict.keys() if x.startswith("trn")]
        self.rrna_list = [x for x in self.feat_dict.keys() if x.startswith("rrn")]
        
    def check_mitos_results(self):
        if os.path.isfile(self.bed_file):
            return True
        else:
            if os.path.isdir(self.results):
                shutil.rmtree(self.results)
            os.mkdir(self.results)
            return False


    def run_mitos(self):
        if not self.check_mitos_results():
            print("Running MITOS for file {}...".format(self.fasta))
            with open("mitos.out", "w+") as output, open('mitos.err', 'w+') as error:
                mitos = subprocess.run(["runmitos.py", "-i", self.fasta, "-c", str(self.gencode), "-o", self.results, "-r", self.refdir, "--linear", "--ncbicode", "--noplots", "--best", "--alarab", "--intron", "0", "--oril", "0", "--orih", "0"], stdout=output, stderr=error)                
                print("Finished MITOS with exit status {}".format(str(mitos.returncode)))
                output.seek(0)
                missing_feat = False
                for line in output:
                    if line.startswith("missing:"):
                        missing_feat = True
                        print('WARNING - {}'.format(line))
                if not missing_feat:
                    print("All features were succesfully annotated")
            os.remove("mitos.out")
            os.remove("mitos.err")
        else:
            print("Annotation process already finished. Skipping to generation of genbank file (if any)...")


    def generate_gbk(self):
        print("Generating file {}...".format(self.gbk))
        with open(self.gbk, "w") as gbk:
            gbk.write(self.format_features() + self.format_sequence())
        print("{} successfully created!".format(self.gbk))

    def format_features(self):
        formatted_feats = ''
        with open(self.bed_file) as bed:
            for feature in bed:
                (feature_name, feature_type, product, anticodon, inipos, endpos, strand) = self.parse_bed(feature)
                inipos = int(inipos) + 1
                if strand.strip() == "+":
                    formatted_feats += "{}{:<16}{:<}\n".format(5*" ", feature_type, "{}..{}".format(str(inipos), endpos))
                if strand.strip() == "-":
                    formatted_feats += "{}{:<16}{:<}\n".format(5*" ", feature_type, "complement({}..{})".format(str(inipos), endpos))
                if anticodon:
                    formatted_feats += "{}{:<}\n".format(21*" ", '/note="anticodon:{}"'.format(anticodon))
                formatted_feats += "{}{:<}\n".format(21*" ", '/product="{}"'.format(product))
                formatted_feats += "{}{:<}\n".format(21*" ", '/gene="{}"'.format(feature_name))
        return formatted_feats

    def format_sequence(self):
        formatted_seq = 'ORIGIN\n'
        mitoseq = SeqIO.read(self.fasta, "fasta").seq.lower()
        subsequences = ((mitoseq[0+i:60+i], i) for i in range(0, len(mitoseq), 60))
        for subseq, index in subsequences:
            formatted_seq += "{:>10} {}\n".format(index+1, subseq)
        formatted_seq += "//"
        return formatted_seq
    
    def parse_bed(self, feature):
        line = feature.split("\t")
        (feature_name, feature_type, product, anticodon) = self.parse_feat_name(line[3])
        inipos = line[1]
        endpos = line[2]
        strand = line[5]
        return (feature_name, feature_type, product, anticodon, inipos, endpos, strand)
        
    def parse_feat_name(self, feat_name):
        anticodon = ''
        if feat_name.strip().startswith("trn"): #'(' index used to separate codon/name - e.g. "trnL1(tag)" (MITOS output in .bed file) to "trnL1" and "tag". The codon is then reverse translated to obtain the anticodon sequence ('cua', in this case).
            separator_index = feat_name.find("(")
            codon = feat_name[separator_index+1:-1]
            anticodon = str(Seq(codon, generic_dna).reverse_complement().transcribe())
            feat_base_name = feat_name[:separator_index]
            final_feat_name = "{}-{}".format(feat_base_name, anticodon)
            product = self.feat_dict.get(feat_base_name).get("product")
            feat_type = 'tRNA' 
        elif feat_name.startswith("rrn"):
            final_feat_name = feat_name
            product = self.feat_dict.get(feat_name).get("product")
            feat_type = 'rRNA'
        elif feat_name in self.cds_list:
            final_feat_name = self.feat_dict.get(feat_name).get("name")
            product = self.feat_dict.get(feat_name).get("product")
            feat_type = 'CDS'
        return (final_feat_name, feat_type, product, anticodon)

import argparse, traceback

def getArgs():
    default_refdir = "{}/refseq63m".format(os.path.dirname(os.path.realpath(__file__)))
    parser = argparse.ArgumentParser(description="A wrapper around MITOS (help on how to install and run here: https://gitlab.com/Bernt/MITOS) to annotate metazoan mitochondrial contigs using hmm (--alarab option)")
    parser.add_argument("-k", "--keep", action="store_true", default=False, help="Keeps MITOS annotation files. Default: False")
    parser.add_argument("-c", "--gencode", type=int, metavar="GENETIC CODE", default=2, help="NCBI's codon table. Default: 2 (Vertebrate Mitochondrial)")
    parser.add_argument("-r", "--refdir", type=str, metavar="REFERENCE DATA", default=default_refdir, help="Custom path to RefSeq63m directory (downloadable from https://zenodo.org/record/2672835). Default: {}".format(default_refdir))
    parser.add_argument("fasta", type=str, nargs="*", metavar="FASTA", help="Fasta file(s) to be annotated")
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    try:
        for fasta in args.fasta:
            annotation = mitoannotation(fasta, args.refdir, gencode=args.gencode)
            annotation.run_mitos()
            annotation.generate_gbk()
            if not args.keep:
                shutil.rmtree(annotation.results)
    except Exception as error:
        fullerror = traceback.format_exc()
        print("An error has occurred for this annotation/gbk conversion process: {}\n\nFULL ERROR:\n\n{}".format(error, fullerror))
