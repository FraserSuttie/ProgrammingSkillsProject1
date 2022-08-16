# importing the "argparse" library for command line interaction
import argparse

# creating a parser to contain the arguments
parser = argparse.ArgumentParser(description = "Input a fasta file to output a file of its open reading frames")
# parsing an argument for the input file to be processed
parser.add_argument("input", type = str, metavar = "input", help = "file to be processed")
# parsing an argument for the minimum open reading frame length, the default value is 50 if the argument is not specified in the command line
parser.add_argument("-s", "--size", type = int, metavar = "", default = 50, help = "minimum open reading frame length, defaults to 50 if not specified")
# parsing an argument for the pathway and name of the output file, the default name is "group1_task1_output.fasta" if the argument is not specified in the command line
parser.add_argument("-o", "--output", type = str, metavar = "", default = "group1_task1_output.fasta", help = "pathway and name of output file, defaults to 'group1_task1_output.fasta'
if not specified")
# parsing the arguments
args = parser.parse_args()

# defining a function to input a nucleotide fasta file and output the sequence as a string
def fileToSequence(file):
    # open the file
    openFile = open(file)
    # counting the number of lines in the file
    numLines = len(open(file).readlines())
    # creating a blank string for the sequence variable
    sequence = ""
    # loop for adding each line of the file to our sequence variable, excluding the first line
    for i in range(numLines):
        line = openFile.readline()
        if line.startswith(">"):
            continue
        else:
            sequence += line
    # removing the line breaks from the string
    sequence = sequence.replace("\n", "")
    # return the new "sequence" string variable
    return sequence

# defining a function to input a sequence string and output the same string with the first character removed
def frame2(string):
    frame2 = string[1:]
    return frame2

# defining a function to input a sequence string and output the same string with the first and second character removed
def frame3(string):
    frame3 = string[2:]
    return frame3

# defining a function to input a sequence string and output the reverse complementary sequence string
def revComp(string):
    # convert each character to its complementary base
    transTab = str.maketrans("atgc", "tacg")
    # reverse the order that the sequence is written
    comp = string.translate(transTab)
    revComp = comp[::-1]
    # return the reverse complementary sequence string
    return revComp

# defining a function to input a sequence as a string and output a list of codons
def sequenceToCodons(string):
    # insert a " " every three characters
    codonsString = ' '.join(string[i:i + 3] for i in range(0, len(string), 3))
    # convert this new string into a list, where elements are deliniated by " "
    codonsList = codonsString.split()
    # return the new list of codons
    return codonsList

# defining a function to input a list of codons and output a list of open reading frames (orfs)
def codonsToOrfs(codonsList):
    # the "read" variable lets us know if the current codon is part of a reading frame or not, it begins at 0
    read = 0
    # creating blank lists for "currentOrf" and "orfsList"
    currentOrf = []
    orfsList = []
    for n in range(len(codonsList)):
        # if we are not currently transcribing an orf
        if (read == 0):
            # if the current codon is a START codon, add the current codon to the "currentOrf" list and set read to 1
            if ((codonsList[n]) == "atg"):
                currentOrf.append(codonsList[n])
                read = 1
            else:
                # if the current codon is not a START codon, move on to the next
                continue
        # if we are currently transcribing an orf
        elif (read == 1):
            # if the current codon contains an ambiguous read, "n", then wipe the "currentOrfs" list and set read to 0
            if ("n" in codonsList[n]):
                currentOrf = []
                read = 0
            else:
                # if the current codon is a STOP codon, add the current codon to the "currentOrf" list, add the completed "currentOrf" list to the "orfsList", wipe the "currentOrfs" list
		  and set read to 0
                if (((codonsList[n]) == "tga") or ((codonsList[n]) == "taa") or ((codonsList[n]) == "tag")):
                    currentOrf.append(codonsList[n])
                    orfsList.append(currentOrf)
                    currentOrf = []
                    read = 0
                # if the current orf is not a STOP codon, add the current codon to the "currentOrf" list
                else:
                    currentOrf.append(codonsList[n])
    # when all codons have been read, orfs are added to "newOrfsList" if their length is greater than or equal to amount defined by the user in the command line (-s/--size)
    newOrfsList = [x for x in orfsList if len(x) >= args.size]
    return newOrfsList

# defining a function to input a list of orfs and output a list of peptides
def orfsToPeptides(orfsList):
    # "for" loop to cycle through the differnt orfs in "orfsList"
    for n in range(len(orfsList)):
        # "for" loop to cycle through the different codons in each orf
        for m in range(len(orfsList[n])):
            # if the current codon matches with the codons in the genetic code, it will be converted to its corresponding amino acid
            if ((orfsList[n])[m] == "aaa"):
                (orfsList[n])[m] = "K"
            elif ((orfsList[n])[m] == "aca"):
                (orfsList[n])[m] = "T"
            elif ((orfsList[n])[m] == "aga"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "ata"):
                (orfsList[n])[m] = "I"
            elif ((orfsList[n])[m] == "caa"):
                (orfsList[n])[m] = "Q"
            elif ((orfsList[n])[m] == "cca"):
                (orfsList[n])[m] = "P"
            elif ((orfsList[n])[m] == "cga"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "cta"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "gaa"):
                (orfsList[n])[m] = "E"
            elif ((orfsList[n])[m] == "gca"):
                (orfsList[n])[m] = "A"
            elif ((orfsList[n])[m] == "gga"):
                (orfsList[n])[m] = "G"
            elif ((orfsList[n])[m] == "gta"):
                (orfsList[n])[m] = "V"
            elif ((orfsList[n])[m] == "taa"):
                (orfsList[n])[m] = ""
            elif ((orfsList[n])[m] == "tca"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "tga"):
                (orfsList[n])[m] = ""
            elif ((orfsList[n])[m] == "tta"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "aac"):
                (orfsList[n])[m] = "N"
            elif ((orfsList[n])[m] == "acc"):
                (orfsList[n])[m] = "T"
            elif ((orfsList[n])[m] == "agc"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "atc"):
                (orfsList[n])[m] = "I"
            elif ((orfsList[n])[m] == "cac"):
                (orfsList[n])[m] = "H"
            elif ((orfsList[n])[m] == "ccc"):
                (orfsList[n])[m] = "P"
            elif ((orfsList[n])[m] == "cgc"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "ctc"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "gac"):
                (orfsList[n])[m] = "D"
            elif ((orfsList[n])[m] == "gcc"):
                (orfsList[n])[m] = "A"
            elif ((orfsList[n])[m] == "ggc"):
                (orfsList[n])[m] = "G"
            elif ((orfsList[n])[m] == "gtc"):
                (orfsList[n])[m] = "V"
            elif ((orfsList[n])[m] == "tac"):
                (orfsList[n])[m] = "Y"
            elif ((orfsList[n])[m] == "tcc"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "tgc"):
                (orfsList[n])[m] = "C"
            elif ((orfsList[n])[m] == "ttc"):
                (orfsList[n])[m] = "F"
            elif ((orfsList[n])[m] == "aag"):
                (orfsList[n])[m] = "K"
            elif ((orfsList[n])[m] == "acg"):
                (orfsList[n])[m] = "T"
            elif ((orfsList[n])[m] == "agg"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "atg"):
                (orfsList[n])[m] = "M"
            elif ((orfsList[n])[m] == "cag"):
                (orfsList[n])[m] = "Q"
            elif ((orfsList[n])[m] == "ccg"):
                (orfsList[n])[m] = "P"
            elif ((orfsList[n])[m] == "cgg"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "ctg"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "gag"):
                (orfsList[n])[m] = "E"
            elif ((orfsList[n])[m] == "gcg"):
                (orfsList[n])[m] = "A"
            elif ((orfsList[n])[m] == "ggg"):
                (orfsList[n])[m] = "G"
            elif ((orfsList[n])[m] == "gtg"):
                (orfsList[n])[m] = "V"
            elif ((orfsList[n])[m] == "tag"):
                (orfsList[n])[m] = ""
            elif ((orfsList[n])[m] == "tcg"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "tgg"):
                (orfsList[n])[m] = "W"
            elif ((orfsList[n])[m] == "ttg"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "aat"):
                (orfsList[n])[m] = "N"
            elif ((orfsList[n])[m] == "act"):
                (orfsList[n])[m] = "T"
            elif ((orfsList[n])[m] == "agt"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "att"):
                (orfsList[n])[m] = "I"
            elif ((orfsList[n])[m] == "cat"):
                (orfsList[n])[m] = "H"
            elif ((orfsList[n])[m] == "cct"):
                (orfsList[n])[m] = "P"
            elif ((orfsList[n])[m] == "cgt"):
                (orfsList[n])[m] = "R"
            elif ((orfsList[n])[m] == "ctt"):
                (orfsList[n])[m] = "L"
            elif ((orfsList[n])[m] == "gat"):
                (orfsList[n])[m] = "D"
            elif ((orfsList[n])[m] == "gct"):
                (orfsList[n])[m] = "A"
            elif ((orfsList[n])[m] == "ggt"):
                (orfsList[n])[m] = "G"
            elif ((orfsList[n])[m] == "gtt"):
                (orfsList[n])[m] = "V"
            elif ((orfsList[n])[m] == "tat"):
                (orfsList[n])[m] = "Y"
            elif ((orfsList[n])[m] == "tct"):
                (orfsList[n])[m] = "S"
            elif ((orfsList[n])[m] == "tgt"):
                (orfsList[n])[m] = "C"
            elif ((orfsList[n])[m] == "ttt"):
                (orfsList[n])[m] = "F"
            # if the current codon does not match with a codon from the genetic code, it will be registered as an "X"
            else:
                (orfsList[n])[m] = "X"
    # return the "orfsList", which will now consist of peptides
    return orfsList

# defining a function to input a list of peptides from 5'3 frame 1 and output the peptides as a string in fasta format
def peptidesToFasta1(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 5'3' frame 1 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-5'3'-1-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a function to input a list of peptides from 5'3 frame 2 and output the peptides as a string in fasta format
def peptidesToFasta2(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 5'3' frame 2 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-5'3'-2-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a function to input a list of peptides from 5'3 frame 3 and output the peptides as a string in fasta format
def peptidesToFasta3(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 5'3' frame 3 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-5'3'-3-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a function to input a list of peptides from 3'5 frame 1 and output the peptides as a string in fasta format
def peptidesToFasta4(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 3'5' frame 1 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-3'5'-1-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a function to input a list of peptides from 3'5 frame 2 and output the peptides as a string in fasta format
def peptidesToFasta5(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 3'5' frame 2 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-3'5'-2-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a function to input a list of peptides from 3'5 frame 3 and output the peptides as a string in fasta format
def peptidesToFasta6(peptidesList):
# create a blank string for variable "fastaString"
    fastaString = ""
    # format the string with 3'5' frame 3 heading
    for n in range(len(peptidesList)):
        fastaString += ">CLAUD-3'5'-3-" + str(n + 1)
        fastaString += "\n"
        fastaString += "".join(peptidesList[n])
        fastaString += "\n"
    # return the printed "fastaString" to include the line breaks
    return fastaString

# defining a compound function to input a nucleotide fasta file and output the peptides as a string in fasta format
def task1(file):
    # creating a blank string for the "task1" variable
    task1 = ""
    # setting the "frame" variable to 1
    frame = 1
    # beginning a "while" loop to iterate through all the different frames
    while (frame <= 6):
        # convert the fasta file to a sequence string
        sequence = fileToSequence(file)
        # if we want the sequence read 3'5', then convert the sequence to the reverse complement
        if (frame >= 4):
            sequence = revComp(sequence)
        # if we want the second frame of the 5'3' or 3'5', then remove the first nucleotide from the sequence
        if ((frame == 2) or (frame == 5)):
            sequence = frame2(sequence)
        # if we want the third frame of the 5'3' or 3'5', then remove the first and second nucleotides from the sequence
        if ((frame == 3) or (frame == 6)):
            sequence = frame3(sequence)
        # convert the sequence string to a list of codons
        codons = sequenceToCodons(sequence)
        # convert the list of codons to a list of orfs
        orfs = codonsToOrfs(codons)
        # convert the list of orfs to a list of peptides
        peptides = orfsToPeptides(orfs)
        # convert the list of peptides to a string with appropriate labelling based on which frame it is from
        if (frame == 1):
            fasta = peptidesToFasta1(peptides)
        elif (frame == 2):
            fasta = peptidesToFasta2(peptides)
        elif (frame == 3):
            fasta = peptidesToFasta3(peptides)
        elif (frame == 4):
            fasta = peptidesToFasta4(peptides)
        elif (frame == 5):
            fasta = peptidesToFasta5(peptides)
        elif (frame == 6):
            fasta = peptidesToFasta6(peptides)
        # add the complete frame read to the "task1" variable
        task1 += str(fasta)
        # move on to the next frame
        frame += 1
    # when all frames have been read and added to the "task1" variable, return "task1"
    return task1

# calling the "task1" function for the file defined by the user in the command line (-i/--input)
task1Output = task1(args.input)
# open a new file with pathway and name defined by the user in the command line (-o/--output)
task1OutputFile = open(args.output, "w")
# write the contents of "task1Output" to the file
task1OutputFile.write(task1Output)
# close the output file
task1OutputFile.close()