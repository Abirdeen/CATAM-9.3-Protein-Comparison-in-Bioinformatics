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
	content, _ := os.ReadFile("./proteins.json")
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
