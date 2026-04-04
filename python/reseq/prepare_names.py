#!/usr/bin/env python3
# Changes read names, so that pairs have identical names and tiles are not stripped of during the mapping

import argparse
import gzip
import sys
from pathlib import Path

# Ignore broken pipe error
from signal import SIG_DFL, SIGPIPE, signal

signal(SIGPIPE, SIG_DFL)


def modName(name, num_spaces, num_colons):
    split_line = name[:-1].split(" ")
    split_part = split_line[num_spaces].split(":")
    part = ":".join(split_part[: num_colons + 1])
    start = "_".join(split_line[:num_spaces])
    if len(start):
        return start + "_" + part
    else:
        return part


def open_text(path):
    file_path = Path(path)
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt", encoding="utf-8")
    return file_path.open("r", encoding="utf-8")


def prepareNames(file1, file2):
    with open_text(file2) as f2:
        first_line2 = f2.readline()

    with open_text(file1) as f1:
        # Find the space and the semicolon where to separate
        first_line1 = f1.readline()

        num_spaces = 0
        num_colons = 0
        last_num_colons = 0
        read_pos = 0
        while (
            read_pos < min(len(first_line1), len(first_line2))
            and first_line1[read_pos] == first_line2[read_pos]
            and (num_colons < 2 or first_line1[read_pos] != " ")
        ):
            if first_line1[read_pos] == " ":
                num_spaces += 1
                last_num_colons = num_colons
                num_colons = 0
            if first_line1[read_pos] == ":":
                num_colons += 1

            read_pos += 1

        # In case we ended with a difference, reduce it by one so the difference is not included
        if read_pos < min(len(first_line1), len(first_line2)) and first_line1[read_pos] != first_line2[read_pos]:
            if num_colons:
                num_colons -= 1
            elif num_spaces:
                num_spaces -= 1
                num_colons = last_num_colons

        # Pass file1 to stdout
        print(modName(first_line1, num_spaces, num_colons))
        n_lines = 1  # We already have the first line
        for line in f1:
            if n_lines % 4 == 0:
                print(modName(line, num_spaces, num_colons))
            else:
                print(line[:-1])  # Remove the line break from line before printing

            n_lines += 1

    return


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=(
            "Return File1 on stdout with changed read names, so paired reads keep identical names "
            "without stripping tile information during mapping."
        )
    )
    parser.add_argument("file1")
    parser.add_argument("file2")
    return parser.parse_args(argv)


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = parse_args(argv)
    prepareNames(args.file1, args.file2)


if __name__ == "__main__":
    main()
