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
