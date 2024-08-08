#!/usr/bin/env python3
import os

configDirs = sorted([x for x in os.listdir(os.getcwd()) if x.isnumeric() and len(x) == 3])
configDirs = configDirs + sorted([x for x in os.listdir(os.getcwd()) if x.isnumeric() and len(x) == 4])

for dir in configDirs:

    flags = [x for x in os.listdir(dir) if x.startswith("FLAG_")]
    cost = [x[len("COST_"):] for x in os.listdir(dir) if x.startswith("COST_")]
    configFile = [x for x in os.listdir(dir) if "config" in x][0]

    if configFile == "config":
        break

    # effRatioFile = [x for x in os.listdir(dir) if x.startswith("eff_ratio")][0]
    # runoutFile = [x for x in os.listdir(dir) if x.startswith("000_runout")][0]

    config = []

    with open(os.path.join(dir, configFile), "r") as f:
        buffer = f.readline()
        while buffer:
            buffer = buffer.split(",")
            if len(buffer) == 2 and buffer[1].strip().isnumeric():
                config.append(int(buffer[1]))
            buffer = f.readline()

    # effRatioText = ""

    # with open(os.path.join(dir, effRatioFile), "r") as f:
    #     effRatioText = f.readline().strip()
    
    # runoutText = ""
    # with open(os.path.join(dir, runoutFile), "r") as f:
    #     runoutText = f.read()

    outString = "["
    for i in range(0, len(config), max(1, len(config)//100)):#(int(os.get_terminal_size()[0]*0.7)))):
        if config[i] == 4:
            outString += '_'
        elif config[i] == 8:
            outString += '-'
        elif config[i] == 10:
            outString += '\''
        elif config[i] == 16:
            outString += '^'
    outString += "]"

    print("{}_{:<7}: {}{}".format(dir, configFile[configFile.find("_") + 1:], outString, cost))
    # print(effRatioText)
    # print(runoutText)
    if flags:
        for flag in flags:
            print("\t {}".format(flag[len("FLAG_"):]))
    print()