# SRHE

This is an experimental implementation of the decomposition ring based homomorphic encryption scheme, which is proposed in [^AH2017] and [^AH2020].


## Prerequesites

- GNU Make (make)
- GNU C++ Compiler (g++)
- GNU Multi-Precision Library (GMP)
- Number Theory Library (NTL)


## Instructions for building and testing programs

1. Get these source files and into the directory.

```
cd SRHE
```

2. Compile test programs (Note: make sure your library path including GMP and NTL).

```
make
```

3. Run the test programs.

```
./nsgenTest
```

or

```
./ringTest
```

or

```
./sheTest
```



## References

[^AH2017]: [S. Arita and S. Handa, "Subring Homomorphic Encryption," Information Security and Cryptology – ICISC 2017](https://link.springer.com/chapter/10.1007/978-3-319-78556-1_7)

[^AH2020]: [S. Arita and S. Handa, "Fully Homomorphic Encryption Scheme Based on Decomposition Ring," IEICE Transactions on Fundamentals of Electronics, Communications and Computer Sciences, 2020, Volume E103.A, Issue 1, Pages 195-211](https://search.ieice.org/bin/summary.php?id=e103-a_1_195) 