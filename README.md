## Formulas for Determinant and Characteristic Polynomial of Seven-Like Matrices
The experiments conducted to verify the correctness of my algorithms for what I called seven-like matrices on certain data.

To see the definition of seven-like matrices and the derivation of the formulas, please refer to the following links:

[My paper](https://bphm.knu.ua/index.php/bphm/issue/view/79/79)

[My conference theses](https://probability.knu.ua/shv2024/ShV_2024.pdf)

## Original Motivation
Faster computation of Usher's model [10] target values (see the Acknowledgements section below).

For a matrix from Usher's paper [10], our algorithms indeed compute the maximum eigenvalue faster than by using SymPy's out-of-the-box means to calculate a characteristic polynomial.

## Acknowledgments
I would like to express my deepest gratitude
* to my family for their support and specifically to my father, Andriy Fisunenko, for suggesting improvements to my paper;
* to my academic supervisor, Vasyl Tereschenko, for his guidance throughout this research;
* and to assistant lecturer Tetiana Kolianova for giving me an inspiration for this paper through her academic assignment to find an optimal way of calculating the goal values of Usher's forest management model, and for providing me with literature from which I discovered Usher's paper [10].

## References
1. Alman, J., Duan, R., Williams, V. V., Xu, Y., Xu, Z., & Zhou, R. (2024). More Asymmetry Yields Faster Matrix Multiplication. _arXiv preprint_ arXiv:2404.16349.
2. Chiantini, L., Hauenstein, J. D., Ikenmeyer, C., Landsberg, J. M., & Ottaviani, G. (2018). Polynomials and the exponent of matrix multiplication. _Bulletin of the London Mathematical Society_, 50(3), 369-389.
3. Dumas, J. G., Pernet, C., & Wan, Z. (2005, July). Efficient computation of the characteristic polynomial. _In Proceedings of the 2005 international symposium on Symbolic and algebraic computation_ (pp. 140-147).
4. Fawzi, A., Balog, M., Huang, A. et al. "Discovering faster matrix multiplication algorithms with reinforcement learning." _Nature_ 610.7930 (2022): 47-53. https://doi.org/10.1038/s41586-022-05172-4
5. Harris, C. R., Millman, K. J., Van Der Walt, S. J. et al. (2020). Array programming with NumPy. _Nature_, 585(7825), 357-362. https://doi.org/10.1038/s41586-020-2649-2
6. Keller-Gehrig, W. (1985). Fast algorithms for the characteristics polynomial. _Theoretical computer science_, 36, 309–317.
7. Leslie, P. H. (1945). On the use of matrices in certain population mathematics. _Biometrika_, 33(3), 183–212.
8. Meurer, A., Smith, C. P., Paprocki, M. et al. (2017). SymPy: symbolic computing in Python. _PeerJ Computer Science_, 3, e103. https://doi.org/10.7717/peerj-cs.103
9. Usher, Michael B. "A matrix approach to the management of renewable resources, with special reference to selection forests." _Journal of Applied Ecology_ (1966): 355-367.
