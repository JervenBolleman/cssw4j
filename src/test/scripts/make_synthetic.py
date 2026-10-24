#!/usr/bin/env python3
"""Generate the deterministic synthetic part of the protein test dataset.

Writes synthetic_queries.fasta and synthetic_targets.fasta next to the
UniProt-derived files. The sequences are chosen to exercise edge cases of the
striped kernels: lengths around SIMD lane boundaries (8/16/32/64 lanes),
very short peptides, low complexity repeats, ambiguity codes (B Z X J U O *),
lower case input, and mutated / indel-carrying copies of real proteins so
alignments contain gaps in both sequences.
"""
import random
import sys
from pathlib import Path

AA = "ARNDCQEGHILKMFPSTWYV"


def read_fasta(path):
    name, seq, out = None, [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        out.append((name, "".join(seq)))
    return out


def mutate(rng, seq, sub_rate, indel_rate, max_indel=6):
    out = []
    i = 0
    while i < len(seq):
        r = rng.random()
        if r < indel_rate / 2:  # deletion
            i += rng.randint(1, max_indel)
            continue
        if r < indel_rate:  # insertion
            out.append("".join(rng.choice(AA) for _ in range(rng.randint(1, max_indel))))
        c = seq[i]
        if rng.random() < sub_rate:
            c = rng.choice(AA)
        out.append(c)
        i += 1
    return "".join(out)


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width)) or ""


def main(dataset_dir):
    d = Path(dataset_dir)
    rng = random.Random(20261024)
    human = dict(read_fasta(d / "uniprot_human.fasta"))
    homologs = dict(read_fasta(d / "uniprot_homologs.fasta"))
    hbb = human["sp|P68871|HBB_HUMAN"]
    p53 = human["sp|P04637|P53_HUMAN"]
    cdk1 = human["sp|P06493|CDK1_HUMAN"]
    egfr = human["sp|P00533|EGFR_HUMAN"]

    queries = []
    # Lengths around lane boundaries, taken from a real protein.
    for n in (1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129):
        queries.append((f"syn_p53_prefix_{n}", p53[:n]))
    # Random sequences (unrelated, mostly low scores).
    for n in (20, 100, 300):
        queries.append((f"syn_random_{n}", "".join(rng.choice(AA) for _ in range(n))))
    # Mutated copies with indels: alignments with gaps on both sides.
    queries.append(("syn_hbb_mut10", mutate(rng, hbb, 0.10, 0.03)))
    queries.append(("syn_hbb_mut30", mutate(rng, hbb, 0.30, 0.05)))
    queries.append(("syn_cdk1_mut20", mutate(rng, cdk1, 0.20, 0.04)))
    queries.append(("syn_egfr_kinase_mut15", mutate(rng, egfr[700:1000], 0.15, 0.03, 12)))
    # Low complexity / repeats.
    queries.append(("syn_polyQ", "M" + "Q" * 40 + "PPPPLPPPP" + "Q" * 20))
    queries.append(("syn_repeat_WCWC", "WC" * 30))
    # Ambiguity codes, stop, rare letters and lower case.
    queries.append(("syn_ambiguous", "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAF"
                    .replace("G", "B", 3).replace("E", "Z", 3).replace("K", "X", 4) + "*UOJ*"))
    queries.append(("syn_lowercase_hbb", hbb[10:90].lower()))
    # Chimera: halves of two unrelated proteins (two separate local hits).
    queries.append(("syn_chimera_hbb_cdk1", hbb[:70] + cdk1[50:150]))
    # Long, highly similar query: score far above 255 to force the 16-bit kernel.
    queries.append(("syn_egfr_mut5", mutate(rng, egfr, 0.05, 0.01)))

    targets = []
    targets.append(("syn_target_random_500", "".join(rng.choice(AA) for _ in range(500))))
    targets.append(("syn_target_short_5", p53[100:105]))
    targets.append(("syn_target_hbb_embedded",
                    "".join(rng.choice(AA) for _ in range(200)) + hbb + "".join(rng.choice(AA) for _ in range(200))))
    targets.append(("syn_target_hbb_tandem", hbb + "GSGSGSGS" + mutate(rng, hbb, 0.25, 0.02)))
    targets.append(("syn_target_polyQ", "Q" * 100))
    targets.append(("syn_target_p53_mut", mutate(rng, homologs["sp|P02340|P53_MOUSE"], 0.1, 0.02)))
    targets.append(("syn_target_stop", "MSTOP*" * 10))

    for fname, recs in (("synthetic_queries.fasta", queries), ("synthetic_targets.fasta", targets)):
        with open(d / fname, "w") as f:
            for name, seq in recs:
                f.write(f">{name}\n{wrap(seq)}\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "src/test/resources/dataset")
