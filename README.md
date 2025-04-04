# CATAM-9.3-Protein-Comparison-in-Bioinformatics

A rewrite of my Cambridge CATAM undergraduate project, "Protein Comparison in Bioinformatics".

## Go

This project is written in [Go](https://go.dev/learn/), a performant language that combines the readability of languages like python with the speed and safety of languages like C++. Since this project primarily focuses on string manipulation, readability is important for ensuring code functions as expected.

## The project

The aim of this project is to study protein & DNA comparison. Efficient identification and comparison of biomolecular sequences has become a central tool of modern molecular biology, allowing biologists to better understand the function and prevalence of particular sequences.

### The edit distance

We will work with two strings $S$ and $T$ of lengths $m$ and $n$ respectively, composed of characters from some finite alphabet. Write $S_i$ for the $i^{th}$ character of $S$, and $S[i, j]$ for the substring $S_i, ..., S_j$. If $i > j$ then $S[i, j]$ is the empty string. A prefix of $S$ is a substring $S[1, k]$ for some $k \leq m$ (possibly the empty prefix). Similarly, a suffix of $S$ is a substring $S[k, m]$ with $k > 1$.

We can perform a sequence of 'edit' operations on $S$ to obtain $T$. Suppose, for instance, that $S = fruit$ and $T = berry$. Then we can:

- **R**eplace f with b;
- **I**nsert e;
- **M**atch r with r;
- **R**eplace u with r;
- **R**eplace i with y;
- **D**elete t.

We call RIMRRD the _editing transcript_. We can use this to obtain the edit sequence:
```
RIMMRD
f ruit
berry
```
The optimal edit transcripts are those that involve the least number of edit operations (Replace, Insert, and Delete). We define the edit distance $d(S, T)$ to be the minimal number of edits to go from $S$ to $T$. We also write $D(i, j) = d(S[1, i], T[1, j])$, and observe that $d(S, T) = D(m, n)$.

### Scoring matrix

A protein is essentially a long sequence of amino acids. Approximately twenty types of amino acid (the exact number is species dependent) are involved in the construction of each protein. A gene is a sequence of DNA which can be translated into a sequence of amino acids, i.e., a protein. Mutations in DNA will lead to changes in the sequence of amino acids, and some mutations are more likely than others.

In [Problem 1](#problem-1), we derive a "scoring function" $s(a, b)$, which is used to measure the cost of replacing $a$ with $b$. However, there are alternative scoring functions we can use, which better model biological processes. The two dominant schemes for assessing the probability of a mutation from amino acid $a$ to amino acid $b$ are the PAM matrices introduced by Dayhoff, and the BLOSUM matrices of Henikoff and Henikoff.

### Scoring for gaps

Some mechanisms for DNA mutations involve the deletion or insertion of large chunks of DNA. Proteins are often composed of combinations of different domains from a relatively small repertoire; so two protein sequences might be relatively similar over several regions, but differ in other regions where one protein contains a certain domain but the other does not.

At some computational cost, we can still align two protein strings taking gaps into account. Let $w(l) < 0$, $l \geq 1$, be the score of deleting (or inserting) a sequence of amino acids of length $l$ from (or into) a protein. Let $v_{gap}(S, T)$ be the gap-weighted score between $S$ and $T$, and write $V_{gap}(i, j)$ for $v_{gap}(S[1, i], T[1, j])$. Then 

$$V_{gap}(i, j) = \max \{E(i, j), F(i, j), G(i, j)\},$$

$$E(i, j) = \max_{0 \leq k \leq j-1} \{V_{gap}(i, k) + w(j - k)\},$$

$$F(i, j) = \max_{0 \leq k \leq i-1} \{V_{gap}(k, j) + w(i - k)\},$$

$$G(i, j) = V_{gap}(i - 1, j - 1) + s(S_i, T_j).$$

Iterating the above equations on the $n \times m$ grid has complexity of $O(mn^2 + nm^2)$. Happily, if $w(l)$ takes some fixed value $u$ for all $l \geq 1$, then there exists an algorithm for finding $v_{gap}$ which has complexity $O(mn)$.

### Statistical significance

We may now ask at what threshold a score $v_{gap}(S, T)$ should be declared to have biological significance.

Let us simplify the problem slightly. Suppose there are only two letters in our alphabet, $a$ and $b$, corresponding, say, to hydrophobic and hydrophilic amino acids. Let $s(a, a) = s(b, b) = 1$ and $s(a, b) = s(b, a) = -1$. Let $U^n$ be a random protein of length $n$: all the amino acids $U^n_1 , ..., U^n_n$ are independent and identically distributed, with $\mathbb{P}(U^n_i = a) = p$ and $\mathbb{P}(U^n_i = b) = 1 - p$. What should we expect $\mathbb{E}(v_{gap}(U^n, V^n))$ to look like? This question is explored in [Problem 7](#problem-7).

### Local Alignment

Full alignment of proteins is meaningful when the two strings are members of the same family. For example, the full sequences of the oxygen-binding proteins myoglobin and haemoglobin are very similar. Often, though, only a small region of the protein is critical to its function and only this region will be conserved throughout the evolutionary process. When we identify two proteins which perform similar functions but look superficially different, it is useful to identify these highly conserved regions.

We aim to find a pair of substrings $S_0$ and $T_0$ of $S$ and $T$ with the highest alignment score, namely, 

$$v_{sub}(S, T) = \max\{v(S_0, T_0) | S_0 \text{ a substring of }S, T_0 \text{ a substring of }T\}.$$

(For simplicity, we will use the ungapped BLOSOM matrix score. We will also write $s(-, a) = s(a, -) < 0$ for the score of an insertion or deletion.)

Finding $v_{sub}(S, T)$ seems to be of much higher complexity than solving the global alignment problem, as there are $\Theta(n^2m^2)$ combinations of substrings of $S$ and $T$. Amazingly, we can solve it using an algorithm whose complexity is still only $O(mn)$.

This question is tackled in [Problem 9](#problem-9).

## Problems

The original CATAM project involved certain explicit questions and problems, which are reproduced (and solved) here.

### Problem 1

Prove that for all $i, j > 0$, 

$$D(i,j) = \min\{D(i-1,j) + 1, D(i,j-1) + 1, D(i-1,j-1) + s(S_i,T_j)\},$$

where $s(a,b)$ is some suitable function which you should determine. What boundary conditions $D(0,0)$, $D(i,0)$, and $D(0,j)$ did you use?

#### Solution

For the boundary conditions, it is immediately clear that $D(i,0) = D(0,i) = i$: it is both necessary and sufficient to make $i$ insertions.

In general, any substring of an optimal edit transcript must also be optimal. From this, the given formula is straightforward to derive. Consider the final character $c$ in the optimal edit transcript for $S[1,i]$ to $T[1,j]$.

- If $c=D$, then the rest of the edit transcript describes an optimal edit from $S[1,i-1]$ to $T[1,j]$;
- If $c=I$, then the rest of the edit transcript describes an optimal edit from $S[1,i]$ to $T[1,j-1]$;
- If $c=R$ or $c=M$, then the rest of the edit transcript describes an optimal edit from $S[1,i-1]$ to $T[1,j-1]$.

From this, we can derive the given formula:

$$D(i,j) = \min\{D(i-1,j) + 1, D(i,j-1) + 1, D(i-1,j-1) + 1 - \delta_{S_i,T_j}\},$$

### Problem 2

Write a program to find the edit distance between two strings. Use it to find the edit distance between 'shesells' and 'seashells'. What is the complexity of your algorithm?

#### Solution

This program is implemented as `EditDistance`.

To find the edit distance between 'shesells' and 'seashells', we use the following code:

```go
package main
import (
	"fmt"
)

func main() {
	empty_map := make(map[edit_key]edit_info)

	str1 := "shesells"
	str2 := "seashells" 

	fmt.Println(EditDistance(str1, str2, empty_map)[edit_key{len(str1), len(str2)}])
}
```

The edit distance between the given strings is 3, with an optimal edit transcript given below:

```
MDMIMIMMMM
she s ells
s eashells
```
Since we use memoization, each value $D(i,j)$ (along with the associated optimal edit transcript) is computed exactly once, and can be computed in constant time; so the total complexity is $O(nm)$, where $n$ and $m$ are the lengths of our two strings.

### Problem 3

Take proteins $A$ and $B$ from the file `proteins.txt` on the CATAM website. (Both proteins are myoglobin, protein $A$ is for the duckbill platypus and protein $B$ for yellowfin tuna.) 

Find the edit distance between them, and give the first 50 steps of an optimal alignment.

#### Solution

We have included the file `proteins.json`, which contains the same data as `proteins.txt`.

We use the following code:

```go
import (
	"encoding/json"
	"fmt"
	"os"
)

type Proteins struct {
	A string
	B string
	C string
	D string
}

func main() {
	var proteins Proteins
	content, _ := os.ReadFile("./assets/proteins.json")
	json.Unmarshal(content, &proteins)
	str1 := proteins.A
	str2 := proteins.B

	empty_map := make(map[edit_key]edit_info)

	edit_details := EditDistance(str1, str2, empty_map)[edit_key{len(str1), len(str2)}]
	fmt.Println(edit_details.distance)
	fmt.Println(edit_details.transcript[0:50])
}
```
This results in an edit distance of 83. The first 50 steps of an optimal alignment are given below:

```
MRRRMRDDDD MMMRMMRMMR MRRRRMRRMM RMMMMRMMMM RDMIMRMMRM
------------------------------------------------------
MGLSDGEWQL VLKVWGKVEG DLPGHGQEVL IRLFKTHPET LEK FDKFKG
MADFDA     VLKCWGPVEA DYTTMGGLVL TRLFKEHPET Q KLFPKFAG
```

### Problem 4

For historical reasons, we now talk about maximising a score rather than minimizing a distance. Let $v(S,T)$ be the maximum score of all edit transcripts from S to T.

Using the BLOSUM matrix `blosum.txt` from the CATAM website for the scoring function $s$, and scoring $-8$ for each Insert or Delete, find the score $v$ between proteins $A$ and $B$ and give the first 50 steps of the optimal alignment.

#### Solution

We have included the json file `BLOSOM.json`, which contains the same data as as `BLOSOM.txt`.

This optimisation algorithm is implemented as `EditDistanceBLOSUM`. Running it on proteins $A$ & $B$, we obtain a score of 290. The first 50 steps of an optimal alignment are given below.

```
MDDRMDDRRR MMMRMMRMMR MRRRRMRRMM RMMMMRMMMM RRRMRMMRMR
------------------------------------------------------
MGLSDGEWQL VLKVWGKVEG DLPGHGQEVL IRLFKTHPET LEKFDKFKGL
M  AD  FDA VLKCWGPVEA DYTTMGGLVL TRLFKEHPET QKLFPKFAGI
```

### Problem 5

Find and implement an algorithm for computing $v_{gap}$ when $w(l) = u$ is fixed for $l \geq 1$. Explain how your algorithm works, and why it has complexity $O(mn)$. What boundary conditions do you use?

#### Solution

Since $w(l)$ is constant for $l \geq 1$, and since 

$$\max_{0\leq k\leq n} \{x_k\} = \max\{x_n, \max_{0\leq k \leq n-1} \{x_k\}\},$$

we can simplify the computation of $V_{gap}(i,j)$:

$$V_{gap}(i, j) = \max \{E(i, j), F(i, j), G(i, j)\},$$

$$E(i, j) = \max \{V_{gap}(i, j-1) + u, E(i,j-1)\},$$

$$F(i, j) = \max \{V_{gap}(i-1, j) + u, F(i-1,j)\},$$

$$G(i, j) = V_{gap}(i - 1, j - 1) + s(S_i, T_j),$$

At the boundary, we have $V_{gap}(0,k) = V_{gap}(k,0) = u$ for $k > 0$, and $V_{gap}(0,0) = 0$.

We also have $E(0,k) = F(k,0) = u$, $k>0$, and $E(0,0) = F(0,0) = 0$.

We set $E(k,0) = F(0,k) = 3 \times u$ to ensure they are never the maximum in any calculations.

We implement this algorithm as `GappedEditDistanceBLOSOM`.

Again using memoization, we compute each of the values $V_{gap}(i,j), E(i,j), F(i,j), G(i,j)$ exactly once in linear time, giving a complexity of $O(nm)$.

### Problem 6

Take proteins $C$ and $D$ from the file `proteins.txt` on the CATAM website. Both proteins are keratin structures in humans. Using the BLOSUM matrix from Section 2 for the scoring function $s$, and $u = −12$ as the fixed score of insertion/deletion, find the gap-weighted score $v_{gap}(C, D)$ and give the first 50 steps of the optimal alignment.

#### Solution

We compute $v_{gap}(C,D) = 1236$, and the first 50 steps of the optimal alignment are

```
MRIIIIIIII IIIRRMMMRR MRMRMMRRDD DDDDDDDDDD DDDDDDDDDD
------------------------------------------------------
MT            SDCSSTH CSPESCGTAS GCAPASSCSV ETACLPGTCA
MPYNFCLPSL SCRTSCSSRP CVPPSCHS 
```

### Problem 7

Consider two random proteins $U^n$ and $V^n$, independent and identically distributed as in the section on [statistical significance](#statistical-significance). Let the score of inserting/deleting a sequence of length $l$ be fixed: $w(l) = u$ for all $l > 1$. Prove that for all $0 \leq p \leq 1$, and for all $u \leq 0$, 

$$\varliminf_{n\rightarrow\infty} \frac{\mathbb{E}(v_{gap}(U^n, V^n))}{n} > 0.$$

(Note: if $x_n$ is a sequence of real numbers, then $\varliminf_{n\rightarrow\infty} x_n$ is defined by $\varliminf_{n\rightarrow\infty} x_n = \lim_{n\rightarrow\infty} \inf_{m>n} x_m$. Note that the limit is always guaranteed to exist, though it may be $+\infty$ or $-\infty$. This is discussed further in most textbooks on analysis.)

In fact, $\lim_{n\rightarrow\infty} n^{-1}\mathbb{E}(v_{gap}(U^n, V^n))$ exists and is strictly positive. One way to obtain an excellence mark in this project, though not the only way, is to show that this limit exists.

#### Solution

Let $w$ be a fixed string of length $k$. Noting that we can delete or insert arbitrarily long sequences for a fixed cost $u$, we can compute $$v_{gap}(U^awU'^b, V^awV'^b) \geq v_{gap}(U^n,V^n) + v_{gap}(w,w) + v_{gap}(U'^n,V'^n) \geq k + 4u.$$ 

Consider the probability that a consecutive $k$ character string of $U^n$ and $V^n$ match. Setting $q = p^2 + (1-p)^2$, this probability is $(n-k+1)\times q^k$. So $$\mathbb{E}(v_{gap}(U^n,V^n)) \geq (n-k+1)\times q^k\times k + 4u$$

Therefore, $\inf_{n}n^{-1}\mathbb{E}(v_{gap}(U^n, V^n)) \geq q^k \times k$. On the other hand, $v_{gap}(U^n,V^n) \leq v_{gap}(U^n,U^n) = n$, so $n^{-1}\mathbb{E}(v_{gap}(U^n, V^n))$ is a bounded sequence, and therefore $\varliminf_{n\rightarrow\infty} \frac{\mathbb{E}(v_{gap}(U^n, V^n))}{n} \geq q^k\times k > 0$; in particular, the limit exists.

To see that $\lim_{n\rightarrow\infty} n^{-1}\mathbb{E}(v_{gap}(U^n, V^n))$ exists and is positive, we use [Fekete's lemma](https://proofwiki.org/wiki/Fekete%27s_Subadditive_Lemma). Write $\alpha_n = \mathbb{E}(v_{gap}(U^n, V^n))$. Then 

$$\begin{equation*} \begin{split}
\alpha_{n+m} & \leq \mathbb{E}(v_{gap}(U^{n+m}[1,n], V^{n+m}[1,n]) \\
 & + v_{gap}(U^{n+m}[n+1,n+m], V^{n+m}[n+1,n+m])) \\
 & = \mathbb{E}(v_{gap}(U^{n}, V^{n})) + \mathbb{E}(v_{gap}(U^{m}, V^{m})) \\
 & = \alpha_n + \alpha_m
\end{split} \end{equation*}$$

using the fact that expectation is linear. So the $\alpha_n$ are subadditive, so $\lim_{n\rightarrow\infty} n^{-1}\alpha_n$ exists and is equal to $\varliminf_{n\rightarrow\infty} n^{-1}\alpha_n > 0$.

### Problem 8

Let $u = −3$, $p = \frac{1}{2}$. Write a program to estimate $n^{-1}\mathbb{E}(v_{gap}(U^n, V^n))$.

Now vary $n$ and estimate the limit of $n^{-1}\mathbb{E}(v_{gap}(U^n, V^n))$.

#### Solution

This program is implemented as `AverageDistanceEstimator`.

Below, we tabulate some values of $\hat{\alpha_n} := \frac{1}{k}\sum_{i=1}^kv_{gap}(U_i^n, V_i^n)$ for i.i.d. $U_i^n, V_i^n$ as above; we use $k=1000$.

| $n$ | $\alpha_n$ | $\alpha_n/n$ | | $n$ | $\alpha_n$ | $\alpha_n/n$ |   | $n$ | $\alpha_n$ | $\alpha_n/n$ |
|--|------|------|-|--|-----|-----|-|--|-----|-----|
|1 |-0.006|-0.006| |11|1.083|0.098| |21|4.478|0.213|
|2 | 0.120| 0.060| |12|1.417|0.118| |22|4.917|0.224|
|3 |-0.076|-0.025| |13|1.633|0.126| |23|5.400|0.235|
|4 | 0.074| 0.019| |14|2.051|0.147| |24|5.592|0.233|
|5 | 0.081| 0.016| |15|2.416|0.161| |25|6.183|0.247|
|6 | 0.073| 0.012| |16|2.586|0.162| |26|6.497|0.250|
|7 | 0.213| 0.030| |17|2.979|0.175| |27|6.876|0.255|
|8 | 0.596| 0.075| |18|3.296|0.183| |28|7.190|0.257|
|9 | 0.628| 0.070| |19|3.961|0.208| |29|7.734|0.267|
|10| 0.709| 0.071| |20|4.139|0.206| |30|8.029|0.268|

Based on these data, we estimate that $\alpha_n/n \rightarrow 0.27$.

### Problem 9

Suppose we restrict ourselves to suffixes of $S$ and $T$:

$$v_{sfx}(S, T) = \max\{v(S_0, T_0) | S_0 \text{ a suffix of }S, T_0 \text{ a suffix of }T\}.$$

Prove carefully that $$v_{sub}(S, T) = \max\{v_{sfx}(S_0, T_0) | S_0 \text{ a prefix of } S, T_0 \text{ a prefix of } T\}.$$

#### Solution

Suppose $S'$, $T'$ are substrings of $S$ and $T$ respectively such that $v_{sub}(S,T) = v(S',T')$. Let $S' = S[i_1,j_1]$, $T'=T[i_2,j_2]$. Then $v(S',T') = v_{sfx}(S[1,j_1], T[1,j_2])$, by our assumption on $(S', T')$. Moreover, for any prefixes $S_0 = S[1,k_1], T_0 = T[1,k_2]$, $v_{sfx}(S_0,T_0) = v(S_0',T_0') \leq v(S',T')$ for some $S_0',T_0'$ suffixes of $S_0, T_0$, by definition. So 

$$\begin{equation*} \begin{split} 
\max\{v_{sfx}(S_0, T_0) | S_0 \text{ a prefix of } S, T_0 \text{ a prefix of } T\}  & = v_{sfx}(S[1,j_1], T[1,j_2]) \\
 & = v(S',T') \\
 & = v_{sub}(S,T),
\end{split} \end{equation*}$$

as required.

### Problem 10

Write $V_{sfx}(i, j)$ for $v_{sfx}(S[1, i], T[1, j])$. Prove that
$$V_{sfx}(i, j) = \max \begin{cases} 0, \\ V_{sfx}(i-1, j-1) + s(S_i, T_j), \\ V_{sfx}(i-1, j) + s(S_i, -), \\ V_{sfx}(i, j-1) + s(-, T_j), \end{cases}$$

with boundary conditions $V_{sfx}(i, 0) = V_{sfx}(0, j) = 0$.

#### Solution

First, note that $v(-,-) = 0$, so we can assume $V_{sfx}(i,j) > 0$. Suppose $S', T'$ are such that $V_{sfx}(i,j) = v(S',T')$. Then consider the last character of the optimal edit transcript from $S'$ to $T'$, which must exist because $S'$ and $T'$ are not both zero. If it is a $D$, then necessarily $V_{sfx}(i,j) = V_{sfx}(i-1,j) + s(S_i,-)$. If it is an $I$, then necessarily $V_{sfx}(i,j) = V_{sfx}(i, j-1) + s(-, T_j)$. If it is an $M$ or $R$, then $V_{sfx}(i,j) = V_{sfx}(i-1, j-1) + s(S_i, T_j)$. At the boundary, it is clearly cheaper to take the empty suffix than to use any number of deletions or insertions. So we are done.

### Problem 11

Find $v_{sub}$ for proteins $C$ and $D$, using the BLOSUM matrix from Section 2 and $s(a, -) = s(-, a) = -2$ for all amino acids $a$.

#### Solution

An algorithm for computing $v_{sub}$ is given as `SubstringAlignment`. 

We find that $v_{sub}(C,D) = 1312$.