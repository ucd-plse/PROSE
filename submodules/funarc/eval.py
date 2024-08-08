import os
import sys
import subprocess
import numpy as np
from decimal import Decimal

REQ_MATCHING_DIGITS = 2
EPS = 3e-4

def read_stdout(stdout_filepath):
    with open(stdout_filepath, "r") as f:
        for line in f.readlines():
            line = line.strip()
            
            if line.startswith("out:"):
                result = float(line.split(":")[-1].strip())

            elif line.startswith("time:"):
                time = float(line.split(":")[-1].strip())

    return result, time

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()


def check_matching_digits(x1, x2, n_digits):

    if fexp(x1) != fexp(x2):
        return False
    else:

        displace = 10**n_digits
        x1_dig = int(fman(x1) * displace) / displace
        x2_dig = int(fman(x2) * displace) / displace

        return x1_dig == x2_dig


if __name__ == "__main__":
    new_logPath = sys.argv[1]
    original_logPath = sys.argv[2]

    if "prose_logs/0000" in new_logPath:
        outLog = sorted([ os.path.join(new_logPath, x) for x in os.listdir(new_logPath) if x == "outlog.txt"])[-1]
        new_result, new_time = read_stdout(outLog)

        print(new_time)
        exit()
    else:
        original_outLog = sorted([ os.path.join(original_logPath, x) for x in os.listdir(original_logPath) if x == "outlog.txt" ])[-1]
        original_result, original_time = read_stdout(original_outLog)

        outLog = sorted([ os.path.join(new_logPath, x) for x in os.listdir(new_logPath) if x == "outlog.txt" ])[-1]

        new_result, new_time = read_stdout(outLog)

        if abs((new_result-original_result)/(original_result)) > EPS:
            print(-1 * new_time)
        else:
            print(new_time)