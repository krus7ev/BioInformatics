#!/usr/bin/python

import sys, getopt
from os import path
from subprocess import check_call, check_output

# Verifies MUSCLE path has been configured
def hasMuscle():
    global musclePath
    if not musclePath or musclePath is '':
        return False
    return True

# Verifies PhyML path has been configured
def hasPhyml():
    global phymlPath
    if not phymlPath or phymlPath is '':
        return False
    return True

# Checks if the BLAST database exists yet
def checkDbCreated(dbName):
    if path.isfile(dbName + ".phr") and path.isfile(dbName + ".pin") and path.isfile(dbName + ".psq"):
        return True
    return False

# Checks there exists a ".pal" file in the folder listing UniProt subparts
def checkIndexListExists(dbFileName):
    if path.isfile(dbFileName + ".pal"):
        return True
    return False

# Gets the contents of a file as a list of lines
def getFileContents(fileName):
    if not path.isfile(fileName):
        raise IOError("Unable to find file \"" + fileName + "\"")
    with open(fileName) as f:
        content = f.readlines()
    if not content or len(content) is 0:
        raise IOError("File \"" + fileName + "\" appears to be empty")
    return content

# Creates the BLAST database from a given FASTA file
def createDb(dbName):
    dbFile = dbName + ".fa"
    print("Trying to create database from \"" + dbFile + "\"")
    if not path.isfile(dbFile):
        raise IOError("Unable to find FASTA file: \"" + dbFile + "\" when making database.")
    check_call(["makeblastdb", "-in", dbName + ".fa", "-dbtype", "prot", "-out", dbName])

# Ensures the BLAST database exists and creates it if it does not
def verifyDb(dbName):
    if not checkDbCreated(dbName):
        print("Can't find database for \"" + dbName + "\"")
        createDb(dbName)
    else:
        print("Found database for \"" + dbName+ "\"")

# Makes sure we're in UniProt folder
def verifyFastaDb(dbFileName):
    if not checkIndexListExists(dbFileName):
        print("Can't find Uniprot Index list for " + dbFileName)
    return

# Performs a query from a file against a database
def queryDb(dbName, queryFile):
    print("Querying \"" + dbName + "\" with \"" + queryFile + "\"")
    if not path.isfile(queryFile):
        raise IOError("Unable to find query file \"" + queryFile + "\"")
    outFile =  queryFile + ".out"
    check_call(["blastp", "-db", dbName, "-query", queryFile, "-out", outFile])
    print("Wrote output to \"" + outFile + "\"")
    return outFile

# Recursively processes lines of content until a given character is found
def processLinesUntil(content, index, until, stripLines, skipFirst):
    # Quit if we reach the end sequence
    if index >= len(content) or (content[index].startswith(until) and not skipFirst):
        return "";
    line = content[index]
    if stripLines:
        line = line.replace("\n", '').replace("\r\n", '')
    # Keep processing the subsequent lines
    return line + processLinesUntil(content, index + 1, until, stripLines, False)

# Processes lines of a match file
def processMatchLine(content, index):
    # Length denotes the end of the output
    line = processLinesUntil(content, index, "Length", True, False)
    if line.startswith("> "):
        line = line[2:] # split off the "> " character
    return line

# Gets the full description of significant hits in the file
def findQueryMatches(outFile):
    print("Parsing query match results")
    content = getFileContents(outFile)
    matches = []
    for i in range(0, len(content)):
        line = content[i]
        if line.startswith(">"):
            match = processMatchLine(content, i)
            matches.append(match)
    return matches

# Gets the identifier of a match
def getMatchIdentifier(match):
    return match.split(' ', 1)[0]

# Extracts the matched sequences from the source db
# The extracted sequences will be over multiple lines
def extractMatchedSequences(dbName, queryFile, matches):
    dbFile = dbName + ".fa"
    print("Extracting matching sequences from FASTA file \"" + dbFile + "\"")
    content = getFileContents(dbFile)
    sequences = []
    for i in range(0, len(content)):
        if content[i].startswith(">"):
            for match in matches:
                identifier = getMatchIdentifier(match)
                if content[i].startswith(">" + identifier):
                    sequence = content[i] + processLinesUntil(content, i + 1, ">", False, False)
                    sequences.append(sequence)
    return sequences

# Extract matched sequences reading source db file
# line by line to optimise the process for large databases
def extractReadMatchedSequences(dbFileName, matches):
    if not dbFileName.endswith(".fasta"):
        print("Warning: Db file not in \".fasta\" format!")
    if not path.isfile(dbFileName):
        raise IOError("Unable to find file \"" + dbFileName + "\"")
    else:
        print("Extracting matching sequences from \".fasta\" file \"" + dbFileName + "\"")

    matchIds = []
    for match in matches:
        matchIds += [getMatchIdentifier(match)]
    sequences = []
    sequence = ""

    f = open(dbFileName)
    print ("Start Matching")
    line = f.readline()
    j = 0
    if not line: raise Exception("Problem reading file during matching!")
    while matchIds:
        if line.startswith(">"):
            j+= 1
            m = 0
            for i in matchIds:
                if line.startswith(">" + i):
                    matchIds.remove(i)
                    m=1
                    print ("Found hit matching uniprot entry #" +str(j))
                    sequence += line
                    line = f.readline()
                    if not line:
                        print ("EOF!")
                        break
                    while not line.startswith(">"):
                        print ("...adding line to sequence match at #" +str(j))
                        sequence += line
                        line = f.readline()
                        if not line:
                            print ("EOF!")
                            break
                    sequences.append(sequence)
                    sequence = ""
                    print("Sequence match .fasta data added.")
                    break
            if not m:
                line = f.readline()
                if not line:
                    print ("EOF!")
                    break
        else:
            line = f.readline()
            if not line:
                print ("EOF!")
                break
        if not line:
            print ("EOF!")
            break
    f.close()
    return sequences

# Outputs the entire match description
def printMatchesVerbose(matches):
    if len(matches) > 0:
        print("Significant matches:")
        for match in matches:
            print(match)
    else:
        print("Unable to find any significant matches")

# Outputs just the identifier of a match
def printMatchesIdentifiers(matches):
    if len(matches) > 0:
        print("Identifiers of significant matches:")
        for match in matches:
            print(getMatchIdentifier(match))
    else:
        print("Unable to find any significant matches")

# Outputs a series of sequences to a file
def writeSequencesToFile(queryFile, sequences):
    sequenceFile = queryFile + ".matched_sequences.fa"
    print("Writing matched sequences to \"" + sequenceFile + "\"")
    fh = open(sequenceFile,"w")
    for sequence in sequences:
        fh.write(sequence)
    fh.close()
    return sequenceFile

# Creates a multiple sequence alignment file using MUSCLE
def sequenceAlignment(queryFile, sequenceFile):
    global musclePath
    print("Using \"" + musclePath + "\" to calculate alignment of matched sequences")
    alignmentFile = queryFile + ".alignment.fa"
    check_call([musclePath, "-in", sequenceFile, "-out", alignmentFile])
    print("Wrote alignment to \"" + alignmentFile + "\"")
    return alignmentFile

# Creates a phylogenetic tree from an alignment file
def phylogeneticTree(queryFile, alignmentFile):
    print("Converting \"" + alignmentFile + "\" to PhyML format")
    phymlFile = alignmentFile + ".phylip"
    check_call(["py", "fastatophyml.py", alignmentFile, phymlFile])
    global phymlPath
    print("Using \"" + phymlPath + "\" to generated phylogenetic tree")
    treeFile = phymlFile + "_phyml_tree.txt"
    check_call([phymlPath, "-i", phymlFile, "-d", "aa"])
    print("Wrote tree to \"" + treeFile + "\"")
    return treeFile

# Parses a given query file and splits it up into multiple query files if required
def extractQueryFiles(queryFile):
    print("Finding queries in \"" + queryFile + "\"")
    content = getFileContents(queryFile)
    queries = []
    for i in range(0, len(content)):
        if content[i].startswith(">"):
            sequence = processLinesUntil(content, i, ">", False, True)
            queries.append(sequence)
    if len(queries) is 0:
        raise IOError("Unable to find any queries in \"" + queryFile + "\"")
    # If the file only contained one query then just use the file itself
    if len(queries) is 1:
        return [queryFile]
    queryFiles = []
    count = 0
    for query in queries:
        filename = queryFile + "." + str(count)
        queryFiles.append(filename)
        fh = open(filename,"w")
        for sequencePart in query:
            fh.write(sequencePart)
        fh.close()
        count = count + 1
    return queryFiles

# Pipes the entire set of operations together
def blastpipe(dbName, queryFile):
    global printVerbose
    global printIds
    global firstStageOnly

    # We process .fasta databases differently to .fa dbs
    if dbName.endswith(".fasta"):
        verifyFastaDb(dbName)
    else:
        if dbName.endswith(".fa"):
            dbName = dbName[:-3]
        verifyDb(dbName)
    queryFiles = extractQueryFiles(queryFile)
    treeFiles = []
    for query in queryFiles:
        outFile = queryDb(dbName, query)
        matches = findQueryMatches(outFile)
        if printVerbose:
            printMatchesVerbose(matches)
        if printIds:
            printMatchesIdentifiers(matches)
        if not firstStageOnly:
            if dbName.endswith(".fasta"):
                sequences = extractReadMatchedSequences(dbName, matches)
            else:
                sequences = extractMatchedSequences(dbName, query, matches)
            sequenceFile = writeSequencesToFile(query, sequences)
            # Process muscle and phyml if they are configured
            if hasMuscle():
                alignmentFile = sequenceAlignment(query, sequenceFile)
                if hasPhyml():
                    treeFile = phylogeneticTree(query, alignmentFile)
                    treeFiles.append(treeFile)
    # If there is more than one tree file they need to be combined
    if len(treeFiles > 1):
        finalTree = queryFile + "_phyml_tree.txt"
        print("Combining trees to \"" + finalTree + "\"")
        fh = open(filename,"w")
        fh.close()
        for treeFile in treeFiles:
            contents = getFileContents(treeFile)
            fh = open(filename,"a")
            for line in contents:
                fh.write(line)
            fh.close()

# Generate output for specific question
def printQuestion(number):
    print("Question " + number)
    print("--------------------------")

def question17(dbName, queryFile):
    printQuestion("17")
    global printVerbose
    global printIds
    global firstStageOnly
    firstStageOnly = True
    printVerbose = False
    printIds = True
    blastpipe(dbName, queryFile)

def question18(dbName, queryFile):
    printQuestion("18")
    global printVerbose
    global printIds
    printVerbose = False
    printIds = False
    blastpipe(dbName, queryFile)

def question19(dbName, queryFile):
    printQuestion("19")
    global printVerbose
    global printIds
    printVerbose = False
    printIds = False
    if not hasMuscle():
        raise ValueError("Specify a MUSCLE path or command with \"-m <muscle command>\"")
    blastpipe(dbName, queryFile)

def question20(dbName, queryFile):
    printQuestion("20")
    global printVerbose
    global printIds
    printVerbose = False
    printIds = False
    if not hasMuscle():
        raise ValueError("Specify a MUSCLE path or command with \"-m <muscle command>\"")
    if not hasPhyml():
        raise ValueError("Specify a PhyML path or command with \"-p <phyml command>\"")
    blastpipe(dbName, queryFile)

def useage():
    print("blastpipe.py -q <queryfile> [-d <database>] [-v (verbose)] [-i (output identifiers)] [-m <muscle command>]")

# Program start
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"viq:d:?:m:p:")
        dbName = "astral.1.75.fa"
        queryFile = ''
        questionNumber = ''
        global musclePath
        musclePath = ''
        global phymlPath
        phymlPath = ''
        global printVerbose
        global printIds
        global firstStageOnly
        firstStageOnly = False
        printVerbose = False
        printIds = False
        for opt, arg in opts:
            if opt == '-v':
                printVerbose = True
            elif opt == '-i':
                printIds = True
            elif opt == '-d':
                dbName = arg
            elif opt == '-q':
                queryFile = arg
            elif opt == '-?':
                questionNumber = arg
            elif opt == '-m':
                musclePath = arg
            elif opt == '-p':
                phymlPath = arg
            else:
                useage()
                sys.exit(1)
        if len(queryFile) is 0:
            raise IOError("Specify a query file with -q <filename>")
        if questionNumber == '17':
            question17(dbName, queryFile)
        elif questionNumber == '18':
            question18(dbName, queryFile)
        elif questionNumber == '19':
            question19(dbName, queryFile)
        elif questionNumber == '20':
            question20(dbName, queryFile)
        else:
            blastpipe(dbName, queryFile)
        print ("Done")
    except getopt.GetoptError:
        useage()
        sys.exit(1)
    except Exception as e:
        print("Program Error:")
        print(e)
        sys.exit(1)
    sys.exit(0)

if __name__ == "__main__":
    main()
