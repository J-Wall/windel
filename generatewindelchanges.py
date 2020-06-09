#!/usr/bin/env python

import csv
import sys

header = sys.stdin.readline().strip().split()
reader = csv.DictReader(sys.stdin, fieldnames=header, dialect="excel-tab")

ref_name = None
diff = 0
for line in reader:
    if line["ref_name"] != ref_name:
        ref_name = line["ref_name"]
        diff = 0
        print("new contig", file=sys.stderr)

    original_coord = f"{ref_name}:{line['position']}"
    if line["deletion"] == "True" and int(line["indel_length"]) > 1:
        original_coord += f"-{int(line['position']) + int(line['indel_length']) - 1}"

    new_coord = f"{ref_name}:{int(line['position']) + diff}"
    if line["insertion"] == "True" and int(line["indel_length"]) > 1:
        new_coord += f"-{int(line['position']) + diff + int(line['indel_length']) - 1}"

    print(f"{original_coord} {new_coord}")

    print(diff, line["length_mode"], file=sys.stderr)
    diff += int(line["length_mode"])
