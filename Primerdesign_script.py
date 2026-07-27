import primer3
import subprocess
from Bio import SeqIO

gff_file = "Brassica_juncea_Varuna_pk.gff"
fasta_file = "Brassica_juncea_Varuna.genome_pk.fasta"
blast_db = "Brassica_juncea_varuna"
output_file = "unique_primers.txt"
target_gene_id = "BjuVaA01g004457"

def extract_gene_sequence(gff_file, fasta_file, target_gene_id):
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue

            chrom = fields[0]
            start, end = int(fields[3]), int(fields[4])
            attributes = fields[8]
            gene_id = attributes.split(";")[0].split("=")[-1]

            if gene_id == target_gene_id and chrom in genome_dict:
                return chrom, str(genome_dict[chrom].seq[start:end])

    return None, None

def check_primer_uniqueness(primer_seq, target_chrom, blast_db):
    with open("temp_primer.fasta", "w") as temp_fasta:
        temp_fasta.write(f">primer\n{primer_seq}\n")

    blast_cmd = [
        "blastn",
        "-query", "temp_primer.fasta",
        "-db", blast_db,
        "-task", "blastn-short",
        "-word_size", "20",
        "-evalue", "1e-3",
        "-outfmt", "6"
    ]

    result = subprocess.run(blast_cmd, capture_output=True, text=True)
    hits = result.stdout.strip().split("\n")

    unique_hit = False
    for hit in hits:
        cols = hit.split("\t")
        if len(cols) < 2:
            continue
        chrom = cols[1]

        if chrom == target_chrom:  
            unique_hit = True
        else:
            return False

    return unique_hit

def design_unique_primers(sequence, target_chrom, blast_db):
    primer_candidates = primer3.design_primers(
        {"SEQUENCE_TEMPLATE": sequence},
        {
            "PRIMER_OPT_SIZE": 22,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 23,
            "PRIMER_PRODUCT_SIZE_RANGE": [[800, 2000]],
            "PRIMER_MIN_GC": 50,
            "PRIMER_MAX_GC": 60,
            "PRIMER_MAX_SELF_ANY": 4,
            "PRIMER_PAIR_MAX_COMPL_ANY": 4,
            "PRIMER_NUM_RETURN": 100,
        }
    )

    for i in range(100):
        fwd_primer = primer_candidates.get(f"PRIMER_LEFT_{i}_SEQUENCE")
        rev_primer = primer_candidates.get(f"PRIMER_RIGHT_{i}_SEQUENCE")

        if fwd_primer and rev_primer:
            if fwd_primer[0] in "GC" and rev_primer[0] in "GC":
                fwd_unique = check_primer_uniqueness(fwd_primer, target_chrom, blast_db)
                rev_unique = check_primer_uniqueness(rev_primer, target_chrom, blast_db)

                if fwd_unique and rev_unique:
                    return fwd_primer, rev_primer  

    return None, None

def main():
    target_chrom, sequence = extract_gene_sequence(gff_file, fasta_file, target_gene_id)

    if sequence:
        fwd_primer, rev_primer = design_unique_primers(sequence, target_chrom, blast_db)

        if fwd_primer and rev_primer:
            with open(output_file, "w") as out:
                out.write("GeneID\tForward_Primer\tReverse_Primer\n")
                out.write(f"{target_gene_id}\t{fwd_primer}\t{rev_primer}\n")
            
            print(f"✔ Unique GC-start primers found for {target_gene_id}")
            print(">fp")
            print(fwd_primer)
            print(">rp")
            print(rev_primer)
        else:
            print(f"✖ No unique GC-start primers found for {target_gene_id}")
    else:
        print(f"✖ Gene {target_gene_id} not found in the GFF file.")

if __name__ == "__main__":
    main()

