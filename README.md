## Problem Overview
This program calculates the probability that individuals in a family have a specific genetic trait based on observed trait data and family relationships. The implementation uses a Bayesian Network to model the relationships between genes and traits across generations.

## What I Did
I had to implement three key functions to make the heredity inference system work:

### 1. Joint Probability Calculator
This function calculates the probability of a specific combination of genes and traits across all people:
- For people without parents in the dataset, it uses population statistics
- For people with parents, it calculates inheritance probabilities based on:
  - How many copies of the gene each parent has
  - The probability of gene mutation during inheritance
  - The probability of exhibiting the trait based on gene count

```python
# My approach for joint_probability:
# 1. Process each person one at a time
# 2. Handle people with vs. without parents differently
# 3. Calculate exact probabilities of gene inheritance from each parent
# 4. Apply mutation probabilities where appropriate
# 5. Multiply everything together for the final joint probability
```

### 2. Probability Distribution Updater
This function adds new evidence to our probability distributions:
- Keeps track of probabilities for different numbers of genes (0, 1, or 2)
- Keeps track of probabilities for having vs. not having the trait
- Updates these probabilities for each person based on new joint probabilities

### 3. Distribution Normalizer
This function ensures our probability distributions are valid:
- Makes sure gene probabilities (0, 1, 2 copies) sum to 1 for each person
- Makes sure trait probabilities (True, False) sum to 1 for each person
- Preserves the relative proportions between values

## Implementation Insights
- The code models a complex Bayesian Network without needing explicit graph structures
- We compute probabilities for all possible combinations of genes and traits
- The final distribution emerges from accumulating all these possibilities
- Probability normalization ensures mathematically valid distributions

