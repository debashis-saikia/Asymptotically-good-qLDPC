## Tensor Product Codes

Let $C_1 \subseteq \mathbb{F}_q^{n_1}$ and $C_2 \subseteq \mathbb{F}_q^{n_2}$ be linear codes.  

The **tensor product code** $C_1 \otimes C_2$ is defined as the set of all matrices  
$c \in \mathbb{F}_q^{n_1 \times n_2}$ such that every row of $c$ is a codeword of $C_2$ and every column of $c$ is a codeword of $C_1$.

Equivalently, $C_1 \otimes C_2$ consists of all codewords that satisfy local constraints along both coordinate directions. This construction increases blocklength and distance while preserving algebraic structure, making tensor product codes a fundamental object in coding theory.
