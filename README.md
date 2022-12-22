# pyMissingAHP

## Introduction
The multi-criteria technique Analytic Hierarchy Process (AHP) needs a complete Pairwise Comparison Matrix (PCM) to generate results. With an incomplete PCM, our algorithm can infer the best (continuous or discrete) values to complete the missing data. The values can be calculated based on the minimum inconsistency (f0), target rank preservation (f1), or both (f0_f1). The target rank preservation can be total (when all criteria are ranked) or partial (when only a set of criteria are ranked). We also allow ties in ranks (criteria with the same rank). For small problems with discrete scale, we offer a brute force method that can find all available solutions.

It's worth noting that our implementation can deal with AHP and Fuzzy AHP. The Fuzzy AHP needs a fuzzy triangular scale to work, and although the user can define his scale, we have implemented a default fuzzy triangular scale that can be used in most problems: 


| Crisp Number |   Fuzzy Number  | 
|--------------|-----------------|
|     1/9      | (1/9, 1/9, 1/9) |
|     1/8      | (1/9, 1/8, 1/7) |
|     1/7      | (1/8, 1/7, 1/6) |
|     1/6      | (1/7, 1/6, 1/5) |
|     1/5      | (1/6, 1/5, 1/4) |
|     1/4      | (1/5, 1/4, 1/3) |
|     1/3      | (1/4, 1/3, 1/2) |
|     1/2      | (1/3, 1/2,   1) |
|       1      | (  1,   1,   1) |
|       2      | (  1,   2,   3) |
|       3      | (  2,   3,   4) |
|       4      | (  3,   4,   5) |
|       5      | (  4,   5,   6) |
|       6      | (  5,   6,   7) |
|       7      | (  6,   7,   8) |
|       8      | (  7,   8,   9) |
|       9      | (  9,   9,   9) |


## Usage
1. Install

```bash
pip install pyMissingAHP
```

2. Try it in **Colab**:

Single Objective AHP: 
- Example 1a - (AHP; f0; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/11FoDq0i5WGY7IH1Kxf7FBboWGAk6Mw9A?usp=sharing))
- Example 1b - (AHP; f0; Continuous): ([ Colab Demo ](https://colab.research.google.com/drive/1Jebj8Dqzm96DAmabF_i1RrS-d_-Au_YI?usp=sharing))
- Example 1c - (AHP; f1; Discrete; Different Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1n9hcYCW85bK5qU_LpNyZcaTalSnvT-de?usp=sharing))
- Example 1d - (AHP; f1; Continuous; Different Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1kB3nJl4jlSWUoviKZXblqMgIJk8iz_VA?usp=sharing))
- Example 1e - (AHP; f1; Discrete; Same Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1D6ae7wgcZg-yNFr_gj5pmxEriL-oG09X?usp=sharing))
- Example 1f - (AHP; f1; Continuous; Same Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1-wMDIPN4ZRgWX3JpyltUjpI8xiR-BKlh?usp=sharing))
- Example 1g - (AHP; f1; Discrete; Partial Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1LScLnOoSFI4FMR5qMRuyykwIcnj_S2lU?usp=sharing))
- Example 1h - (AHP; f1; Continuous; Partial Rank Positions): ([ Colab Demo ](https://colab.research.google.com/drive/1QjqU3uo0pnW4CuyTTmnaEyElpRdfDiE6?usp=sharing))
- Example 1  -  Brute Force - (AHP; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/1y1tycNbDFxFYiSb3_BrHmP2dUnOOHIqG?usp=sharing))

Single Objective Fuzzy AHP: 
- Example 2a - (Fuzzy AHP; f0; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/1aBEP7lYbSvpHJxJGxYrg4QS4na8Jk49f?usp=sharing))
- Example 2b - (Fuzzy AHP; f1; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/18aeD00Q2jmc_P6QSHGjuEKDUKeEoIiq4?usp=sharing))
- Example 2c - (Fuzzy AHP; f0; Discrete; Custom Fuzzy Scale): ([ Colab Demo ](https://colab.research.google.com/drive/1vPBq4CzNXS503W-ANdW8-WYacdDgSOVr?usp=sharing))
- Example 2d - (Fuzzy AHP; f1; Discrete; Custom Fuzzy Scale): ([ Colab Demo ](https://colab.research.google.com/drive/1sfpmhM7U3xvSKfGlbRlVNKhSszN4vAA4?usp=sharing))
- Example 2  -  Brute Force - (Fuzzy AHP; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/1FmhWnZw3SA7sCGxLYLK6-ISsrYX9kEZU?usp=sharing))

Multiobjective AHP:
- Example 3a - (AHP; f0 & f1; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/1kDo5Ur0_xK2LzGmDOPd0kLjwQwuOKezE?usp=sharing))
- Example 3b - (AHP; f0 & f1; Continuous): ([ Colab Demo ](https://colab.research.google.com/drive/1IwRxyxHXMEAAdDPSTr6yy8otEkv7l3kW?usp=sharing))

Multiobjective Fuzzy AHP:
- Example 4a - (Fuzzy AHP; f0 & f1; Discrete): ([ Colab Demo ](https://colab.research.google.com/drive/1_zRxMOmGgoEoiddF94383OHYq_ztF0nT?usp=sharing))
- Example 4b - (Fuzzy AHP; f0 & f1; Discrete; Custom Fuzzy Scale): ([ Colab Demo ](https://colab.research.google.com/drive/1Jn6KElsYwN6W9IXR4XbBDy2AW6JYoh9t?usp=sharing))

3. Others

- [pyDecision](https://github.com/Valdecy/pyDecision) - A library for many MCDA methods
- [3MOAHP](https://github.com/Valdecy/Method_3MOAHP) - Inconsistency Reduction Technique for AHP and Fuzzy-AHP Methods
- [ELECTRE-Tree](https://github.com/Valdecy/ELECTRE-Tree) - Algorithm to infer the ELECTRE Tri-B method parameters
- [Ranking-Trees](https://github.com/Valdecy/Ranking-Trees) - Algorithm to infer the ELECTRE II, III, IV and PROMETHEE I, II, III, IV method parameters
