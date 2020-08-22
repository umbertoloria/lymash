import os
from factorization import get_fingers_after_text_subdividing
from format import read_from_fasta


def preprocess(dir, factorization_method, subdivision):
    if not os.path.isdir(dir + "_preprocessed"):
        os.mkdir(dir + "_preprocessed")
    for root, dirs, filepaths in os.walk(dir):
        for filepath in filepaths:
            if filepath.endswith(".fasta"):
                title, read = read_from_fasta(dir + "/" + filepath)
                fingers = get_fingers_after_text_subdividing(factorization_method, read, subdivision)
                file = open(os.path.join(dir + "_preprocessed", "preprocessed_" + filepath), "w")
                file.write(title +
                           " preprocessed using lyndon fac and encode the " +
                           "lenght of the factors with ascii code\n")
                i = 1
                for finger in fingers:
                    file.write(chr(ord("0") + finger))
                    if i == 70:
                        file.write("\n")
                        i = 1
                    else:
                        i += 1
                file.close()
