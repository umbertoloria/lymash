import os

from factorization import get_fingers_after_split, subdivide, get_factors
from format import read_from_fasta
from lyndon.utils import _complement


def main():
    dir = input("Tell me the directory with the file to preprocess> ")
    factorization = "cfl_icfl_comb"
    # input("Tell me the algorithm for factoring\n   (cfl,icfl,cfl_icfl,cfl_comb,icfl_comb,cfl_icfl_comb)> ")

    if not os.path.isdir(dir + "_preprocessed"):
        os.mkdir(dir + "_preprocessed")
    for root, dirs, filepaths in os.walk(dir):
        for filepath in filepaths:
            if filepath.endswith(".fasta"):
                title, read = read_from_fasta(dir + "/" + filepath)
                comp_read = _complement(read)
                print("read : {} reversed : {}".format(read, comp_read))
                comp_fingerprint = get_fingers_after_split(factorization, comp_read,100)
                fingerprint = get_fingers_after_split(factorization, read,100)
                file = open(os.path.join(dir + "_preprocessed", "preprocessed_" + filepath), "w")
                file.write( title + " preprocessed encoding lyndon fact lengths with ascii code\n")
                alf = " 0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZàèìòù%^-+:_çé!?£$§°*"
                print(filepath)
                print(fingerprint)
                i=0
                for finger in fingerprint:
                    if i%70==0 and i>0:
                        file.write("\n")
                    file.write(alf[finger])
                    i+=1


                file.write("#")
                i=0
                print(comp_fingerprint)
                for finger in comp_fingerprint:
                    if i % 70 == 0:
                        file.write("\n")
                    file.write(alf[finger])
                    i += 1

                file.close()


main()