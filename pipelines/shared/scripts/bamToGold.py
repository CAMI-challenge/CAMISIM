#!/usr/bin/env python
"""
Convert BAM file to gold standard assembly.
"""

import argparse
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate gold standard assembly from BAM file using samtools mpileup"
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Reference FASTA file (indexed with faidx)"
    )
    parser.add_argument(
        "-b", "--bam", required=True,
        help="Sorted BAM file"
    )
    parser.add_argument(
        "-l", "--length", type=int, required=True,
        help="Length cutoff for output sequences"
    )
    parser.add_argument(
        "-c", "--coverage", type=int, required=True,
        help="Minimum coverage threshold"
    )
    parser.add_argument(
        "-st", "--samtools", default="samtools",
        help="Path to samtools executable (default: samtools)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Build samtools mpileup command
    cmd = [
        args.samtools, "mpileup",
        "-B", "-Q", "0",
        "-f", args.reference,
        args.bam
    ]

    pos = -42
    seq = ""
    start = 0
    stop = 0
    seqname = ""
    previous_sequence_name = ""

    # Run samtools mpileup and process output
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True
        )

        for line in process.stdout:
            # mpileup output format: seqname, pos, ref_base, coverage, ...
            # We only need the first 4 fields (equivalent to cut -f1,2,3,4)
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            seqname = fields[0]
            current_pos = int(fields[1])
            ref_base = fields[2]
            coverage = int(fields[3])

            # Check if contig is done (new sequence or non-consecutive position)
            if seqname != previous_sequence_name or (pos + 1) != current_pos:
                # Output previous sequence if it meets length cutoff
                if len(seq) >= args.length and previous_sequence_name != "":
                    stop = pos
                    header = f">{previous_sequence_name}_from_{start}_to_{stop}_total_{len(seq)}"
                    print(header)
                    print(seq)

                previous_sequence_name = seqname
                start = current_pos
                seq = ""

            pos = current_pos

            # Add base to sequence if coverage threshold is met
            if coverage >= args.coverage:
                seq += ref_base

        process.wait()

        # Catch the last sequence
        if len(seq) >= args.length:
            stop = pos
            header = f">{previous_sequence_name}_from_{start}_to_{stop}_total_{len(seq)}"
            print(header)
            print(seq)

    except FileNotFoundError:
        print(f"Error: samtools not found at '{args.samtools}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
