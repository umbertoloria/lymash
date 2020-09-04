import os
import random
import Sequence


def create_dataset_fromread(seq: Sequence.FastaSequence, n):
    #create a dataset of n short reads from one long read
    start=1
    lenght = 600
    read = seq.get_data()
    print(len(read))
    name = seq.get_name()
    title = seq.get_title()

    dirname = "dataset_" + name.split(".")[0]
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    for i in range(1,n+1):
        file = open(os.path.join(dirname + "/"+str(i)+".fasta"), "w")
        file.write(title + " #short read generated\n")

        for x in range(lenght):
            file.write(read[start+x])
        file.write("\n")
        rand = random.randrange(250,lenght)
        start = start+lenght - rand
    print(start)
    return dirname

z=0
while 1:
    z+=1
    seq1 = Sequence.FastaSequence("in/300.fasta")
    dirname = create_dataset_fromread(seq1,100)


    import os
    from Sequence import FastaSequence

    from main_estimate_jaccard import estimate_jaccard_difference

    # CLI inputs
    files = []
    path = dirname
    while path != "":
        if not os.path.exists(path):
            print("File o directory non esistente")
        elif os.path.isdir(path):
            new_files = os.listdir(path)
            new_files.sort()
            for new_file in new_files:
                new_file = path + "/" + new_file
                if new_file not in files:
                    files.append(new_file)
        else:
            if path not in files:
                files.append(path)
        path = ""

    kmer_size = 16
    kmer_size = 21 if kmer_size == "" else int(kmer_size)
    print()
    factorization = "cfl_icfl_comb"
    if factorization == "":
        factorization = "cfl_icfl_comb"
    kfinger_size = 3
    kfinger_size = 3 if kfinger_size == "" else int(kfinger_size)
    files = sorted(files, key=lambda x: int(x.split(".")[0].split("/")[1]))

    max = 0

    for i in range(len(files)):
        file = files[i]
        seq1 = FastaSequence(file)
        read1 = seq1.get_data()
        from factorization import get_factors
        from lyndon.utils import _complement

        factors = [len(f) for f in get_factors(factorization, read1)]

        comp_fingerprint = [len(f) for f in get_factors(factorization, _complement(read1))]
        fingerprint = [len(f) for f in get_factors(factorization, read1)]
        for x in fingerprint:
            if int(x) > max:
                max = int(x)

        for w in comp_fingerprint:
            if int(w)> max:
                max = int(w)



    y = 0
    for i in range(len(files) - 1):
        file1 = files[i]
        file2 = files[i + 1]
        #print("Comparing", file1, "-", file2)
        seq1 = FastaSequence(file1)
        seq2 = FastaSequence(file2)
        calc, estim = estimate_jaccard_difference(seq1, seq2, kmer_size, factorization, kfinger_size)
        #print(calc, estim)
        if abs(calc - estim) > 0.30:
            y += 1
    print(y)
    print("il max è {} ".format(max))
    if y <= 900 and max <500:
        break

print(y)