import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    joint_prob = 1.0
    
    for person in people:
        # Determine number of copies of gene for this person
        gene_copies = 2 if person in two_genes else 1 if person in one_gene else 0
        
        # Determine if person has the trait
        has_trait = person in have_trait
        
        # Calculate probability based on gene copies and trait
        person_prob = 1.0
        mother = people[person]["mother"]
        father = people[person]["father"]
        
        # Case 1: Person has no parents in dataset
        if mother is None and father is None:
            person_prob = PROBS["gene"][gene_copies]
            
        # Case 2: Person has parents in dataset
        else:
            # Calculate probability of gene inheritance from mother
            mother_has_gene = 2 if mother in two_genes else 1 if mother in one_gene else 0
            father_has_gene = 2 if father in two_genes else 1 if father in one_gene else 0
            
            # Calculate probability of gene inheritance based on parents' genes
            if gene_copies == 0:
                # Need to NOT inherit gene from either parent
                prob_no_gene_from_mother = 1 - (mother_has_gene / 2) if mother_has_gene > 0 else 0.99
                prob_no_gene_from_father = 1 - (father_has_gene / 2) if father_has_gene > 0 else 0.99
                
                # Apply mutation probability
                prob_no_gene_from_mother = (1 - PROBS["mutation"]) if prob_no_gene_from_mother > 0.5 else PROBS["mutation"]
                prob_no_gene_from_father = (1 - PROBS["mutation"]) if prob_no_gene_from_father > 0.5 else PROBS["mutation"]
                
                person_prob = prob_no_gene_from_mother * prob_no_gene_from_father
                
            elif gene_copies == 1:
                # Need to inherit gene from exactly one parent
                # Either from mother and not father, or from father and not mother
                
                # Probability of getting gene from mother
                prob_gene_from_mother = mother_has_gene / 2 if mother_has_gene > 0 else 0.01
                # Probability of not getting gene from father
                prob_no_gene_from_father = 1 - (father_has_gene / 2) if father_has_gene > 0 else 0.99
                
                # Apply mutation probability
                prob_gene_from_mother = (1 - PROBS["mutation"]) if prob_gene_from_mother > 0.5 else PROBS["mutation"]
                prob_no_gene_from_father = (1 - PROBS["mutation"]) if prob_no_gene_from_father > 0.5 else PROBS["mutation"]
                
                # Probability of getting gene from father
                prob_gene_from_father = father_has_gene / 2 if father_has_gene > 0 else 0.01
                # Probability of not getting gene from mother
                prob_no_gene_from_mother = 1 - (mother_has_gene / 2) if mother_has_gene > 0 else 0.99
                
                # Apply mutation probability
                prob_gene_from_father = (1 - PROBS["mutation"]) if prob_gene_from_father > 0.5 else PROBS["mutation"]
                prob_no_gene_from_mother = (1 - PROBS["mutation"]) if prob_no_gene_from_mother > 0.5 else PROBS["mutation"]
                
                # Probability of getting exactly one gene
                person_prob = (prob_gene_from_mother * prob_no_gene_from_father) + (prob_gene_from_father * prob_no_gene_from_mother)
                
            else:  # gene_copies == 2
                # Need to inherit gene from both parents
                prob_gene_from_mother = mother_has_gene / 2 if mother_has_gene > 0 else 0.01
                prob_gene_from_father = father_has_gene / 2 if father_has_gene > 0 else 0.01
                
                # Apply mutation probability
                prob_gene_from_mother = (1 - PROBS["mutation"]) if prob_gene_from_mother > 0.5 else PROBS["mutation"]
                prob_gene_from_father = (1 - PROBS["mutation"]) if prob_gene_from_father > 0.5 else PROBS["mutation"]
                
                person_prob = prob_gene_from_mother * prob_gene_from_father
                
        # Multiply by probability of having the trait given gene copies
        person_prob *= PROBS["trait"][gene_copies][has_trait]
        
        # Multiply joint probability
        joint_prob *= person_prob
    
    return joint_prob


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        # Update gene distribution
        if person in two_genes:
            probabilities[person]["gene"][2] += p
        elif person in one_gene:
            probabilities[person]["gene"][1] += p
        else:
            probabilities[person]["gene"][0] += p
        
        # Update trait distribution
        probabilities[person]["trait"][person in have_trait] += p


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        # Normalize gene distribution
        gene_sum = sum(probabilities[person]["gene"].values())
        if gene_sum > 0:
            for gene_count in probabilities[person]["gene"]:
                probabilities[person]["gene"][gene_count] /= gene_sum
        
        # Normalize trait distribution
        trait_sum = sum(probabilities[person]["trait"].values())
        if trait_sum > 0:
            for trait_value in probabilities[person]["trait"]:
                probabilities[person]["trait"][trait_value] /= trait_sum


if __name__ == "__main__":
    main()
