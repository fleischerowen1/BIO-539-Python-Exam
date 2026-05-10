import sys


def validate_sequence(sequence, k):
    """
    Check whether a DNA sequence is valid for k-mer analysis.

    A valid sequence must be at least as long as k and contain only
    the DNA bases A, C, G, and T.

    Parameters
    ----------
    sequence : str
        DNA sequence to validate.
    k : int
        Length of the k-mer.

    Returns
    -------
    bool
        True if the sequence is valid, False otherwise.
    """
    if len(sequence) < k:
        return False

    for nucleotide in sequence:
        if nucleotide not in "ACGT":
            return False

    return True


def update_kmer_count(kmer_data, kmer, next_char):
    """
    Update the count for one k-mer and the character that follows it.

    Parameters
    ----------
    kmer_data : dict
        Dictionary storing k-mer counts and following character counts.
    kmer : str
        The k-mer being counted.
    next_char : str or None
        Character immediately after the k-mer. None is used when the
        k-mer occurs at the end of a sequence.

    Returns
    -------
    dict
        Updated k-mer data dictionary.
    """
    # Create a new entry for this k-mer if it has not appeared before.
    if kmer not in kmer_data:
        kmer_data[kmer] = {"count": 0, "next_chars": {}}

    # Increase the total frequency of this k-mer.
    kmer_data[kmer]["count"] += 1

    # Only count a following character if one exists.
    if next_char is not None:
        if next_char not in kmer_data[kmer]["next_chars"]:
            kmer_data[kmer]["next_chars"][next_char] = 0

        kmer_data[kmer]["next_chars"][next_char] += 1

    return kmer_data


def count_kmers_with_context(sequence, k):
    """
    Count all k-mers in a DNA sequence and record following characters.

    Parameters
    ----------
    sequence : str
        DNA sequence to analyze.
    k : int
        Length of each k-mer.

    Returns
    -------
    dict
        Dictionary containing total counts for each k-mer and counts of
        characters that immediately follow each k-mer.
    """
    kmer_data = {}

    # Include the final k-mer, even though it has no following character.
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]

        if i + k < len(sequence):
            next_char = sequence[i + k]
        else:
            next_char = None

        update_kmer_count(kmer_data, kmer, next_char)

    return kmer_data


def write_results_to_file(kmer_data, output_filename):
    """
    Write k-mer counts and following character frequencies to a file.

    Parameters
    ----------
    kmer_data : dict
        Dictionary of k-mer counts and following character counts.
    output_filename : str
        Name of the output file to write.
    """
    sorted_kmers = sorted(kmer_data.keys())

    with open(output_filename, "w") as f:
        for kmer in sorted_kmers:
            count = kmer_data[kmer]["count"]
            next_chars = kmer_data[kmer]["next_chars"]

            next_char_str = " ".join(
                f"{char}:{freq}"
                for char, freq in sorted(next_chars.items())
            )

            f.write(f"{kmer} count:{count} {next_char_str}\n")


def main():
    """
    Run the k-mer analysis from the command line.

    Expected command:
    python kmer_analysis.py input_file k output_file
    """
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]

    kmer_data = {}

    print(f"Reading sequences from {sequence_file}...")

    with open(sequence_file, "r") as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print("  Warning: Skipping invalid sequence")
                continue

            sequence_kmers = count_kmers_with_context(sequence, k)

            # Combine results from all sequences in the input file.
            for kmer in sequence_kmers:
                for _ in range(sequence_kmers[kmer]["count"]):
                    update_kmer_count(kmer_data, kmer, None)

                for next_char, freq in sequence_kmers[kmer]["next_chars"].items():
                    for _ in range(freq):
                        kmer_data[kmer]["next_chars"][next_char] = (
                            kmer_data[kmer]["next_chars"].get(next_char, 0) + 1
                        )

    write_results_to_file(kmer_data, output_file)


if __name__ == "__main__":
    main()
    
    
