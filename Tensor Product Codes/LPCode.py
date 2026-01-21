import code
import numpy as np
from algebra import (
    RingElement, RingMatrix,
    GroupAlgebraElement, RingLifter,
    GroupAlgebraMatrix,
)
from HGPCode import HGP


class LPC:
    """
    Full Lifted Product Code constructor.

    Pipeline:
    Abstract Ring Matrices (A, B)
        → Group Algebra Lift
        → Binary Permutation Expansion
        → Hypergraph Product (CSS Code)

    Outputs:
    - H_X, H_Z
    - n (number of physical qubits)
    - k (number of logical qubits)
    """

    def __init__(self, A_ring: RingMatrix, B_ring: RingMatrix,
                 generator_images: list[GroupAlgebraElement]):
        """
        Parameters
        ----------
        A_ring, B_ring : RingMatrix
            Abstract ring matrices satisfying AB = 0

        generator_images : list[GroupAlgebraElement]
            Images of ring generators (e.g. [g, g^{-1}])
        """

        if not (A_ring @ B_ring).is_zero():
            raise ValueError("Abstract condition AB = 0 is not satisfied.")

        self.A_ring = A_ring
        self.B_ring = B_ring

        self.lifter = RingLifter(generator_images)

        self.GA = GroupAlgebraMatrix([
            [self.lifter.lift(e) for e in row] for row in A_ring.data
        ])

        self.GB = GroupAlgebraMatrix([
            [self.lifter.lift(e) for e in row] for row in B_ring.data
        ])

        if not (self.GA @ self.GB).is_zero():
            raise ValueError("Lifted condition AB = 0 failed in group algebra.")


        self.HA = lift_to_binary(self.GA)
        self.HB = lift_to_binary(self.GB)

        if not np.all((self.HA @ self.HB) % 2 == 0):
            raise ValueError("Binary CSS condition HA HB = 0 failed.")

        self.hgp = HGP(A=self.HA, B=self.HB)

        self.H_X = self.hgp.H_X
        self.H_Z = self.hgp.H_Z

        self.n = self.hgp.n
        self.k = self.hgp.k

    def parameters(self):
        """Return code parameters [[n, k]]."""
        return self.n, self.k

    def parity_checks(self):
        """Return H_X and H_Z."""
        return self.H_X, self.H_Z


'''# Example usage:

x = RingElement([1, 0])
y = RingElement([0, 1])

A = RingMatrix([[x, y]])
B = RingMatrix([[y], [x]])

L = 4
g = GroupAlgebraElement([0,1,0,0], L)
g_inv = GroupAlgebraElement([0,0,0,1], L)

code = LPC(A, B, [g, g_inv])

print("n, k =", code.parameters())

H_X, H_Z = code.parity_checks()
print("H_X shape:", H_X.shape)
print("H_Z shape:", H_Z.shape)'''


