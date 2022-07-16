# -*- coding: utf-8 -*-
r"""
Algebraic Tori
==============


Tori are implemented using the equivalence of categories with  character lattices
with an action of the Galois group (at least as large as the Galois group of
a splitting field). For now this Galois group will just be an abstract group,
either a permutation group or a finite matrix group in ``GL(n,ZZ)``.


To define a torus we use ``AlgebraicTorus(character_lattice)``

::

    sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
    sage: L = GLattice([], 1)
    sage: AlgebraicTorus(L)
    Split algebraic torus of rank 1

This is the split torus ``\mathbb{G}_m``, with action of the trivial Galois group.

::

    sage: LL = GLattice(SymmetricGroup(3), 1)
    sage: AlgebraicTorus(LL)
    Split algebraic torus of rank 1

This is still ``\mathbb{G}_m``, with trivial action of a Galois group isomorphic to `S_3`. Note that
this Galois group is not necessarily the one of a minimal splitting extension.

::

    sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
    sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
    sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
    sage: T3 = AlgebraicTorus(LLL); T3
    Algebraic torus of rank 3 split by a degree 6 extension

This is a non-split anisotropic torus with Galois group of splitting field isomorphic to ``S_3``.

::

    sage: SL = L.sublattice([2*L.basis()[0]])
    sage: AlgebraicTorus(SL)
    Split algebraic torus of rank 1

This torus is obtained from the sublattice of the first lattice `L`. The torus obtained is isomorphic to `L`.

One can define an algebraic torus form a number field. We can compute restriction of scalars.

::

    sage: from sage.schemes.group_schemes.tori import RestrictionOfScalars
    sage: x = polygen(QQ);  K.<a> = NumberField(x^8 + 68*x^6 + 986*x^4 + 4624*x^2 + 4624)
    sage: RestrictionOfScalars(K)
    Algebraic torus of rank 8 over Rational Field split by a degree 8 extension

We can similarly define norm-one tori.

::

    sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
    sage: F.<b> = QuadraticField([2])
    sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
    sage: NormOneRestrictionOfScalars(F)
    Algebraic torus of rank 1 over Rational Field split by a degree 2 extension


Attributes of a Torus
=====================

- ``torus._lattice`` -- the character lattice of the torus

- ``torus._base_field`` -- lazy attribute, the field over which the torus is defined.

- ``torus._top_field`` -- lazy attribute, a field over which the torus splits

- ``torus._splitting_field`` -- lazy attribute, splitting field of the torus

Methods of a Torus
==================

- :meth:`AlgebraicTorus.rank` -- the rank of the torus.

- :meth:`AlgebraicTorus.galois_group` -- the Galois group (as abstract group) of a splitting field of the torus.

- :meth:`AlgebraicTorus.character_lattice` -- the character lattice of the torus.

- :meth:`AlgebraicTorus.cocharacter_lattice` -- the cocharacter lattice of the torus.

- :meth:`AlgebraicTorus.is_rational` -- tests if a point is rational.

- :meth:`AlgebraicTorus.tate_cohomology` -- the isomorphism type of Tate Cohomology groups of the Torus over a local field.

- :meth:`AlgebraicTorus.tamagawa_number` -- the Tamagawa number.

- :meth:`AlgebraicTorus.product` -- the product of the torus with another specified torus with same base and splitting field.

- :meth:`AlgebraicTorus.restriction_of_scalars` -- returns the torus obtained by restriction of scalars.

- :meth:`AlgebraicTorus.norm_one_restriction` -- the torus of norm 1 elements in the restriction of scalars.
"""

###########################################################################
#       Copyright (C) 2018-2019 Thomas Rüd <tompa.rud@gmail.com>
#                                David Roe <roed.math@gmail.com>

#
#  Distributed under the terms of the GNU General Public License (GPL)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
###########################################################################


from __future__ import print_function, absolute_import

from sage.schemes.generic.scheme import Scheme
from sage.matrix.constructor import matrix
from sage.modules.glattice import GLattice
from sage.misc.lazy_attribute import lazy_attribute


def RestrictionOfScalars(field, rk=1):
    r"""
    Create the torus defined by restriction of scalars of a split torus to the base field (rationals, rational function field, `\QQ_p` or ``GF(p)``).

    INPUT:

    - ``field`` -- Field extension over which the restriction of scalars is taken.  Usually a number field, function field, p-adic extension field or finite field.

    - ``rk`` -- rank (default ``1``) of the split torus of which we take the restriction of scalars.

    EXAMPLES::

        sage: from sage.schemes.group_schemes.tori import RestrictionOfScalars
        sage: F.<a> = QuadraticField([2])
        sage: T1 = RestrictionOfScalars(F); T1
        Algebraic torus of rank 2 over Rational Field split by a degree 2 extension
        sage: K.<t> = NumberField(x^4+x^3+x^2+x+1)
        sage: T2 = RestrictionOfScalars(K); T2
        Algebraic torus of rank 4 over Rational Field split by a degree 4 extension
        sage: T2.galois_group()
        Galois group 4T1 (4) with order 4 of x^4 + x^3 + x^2 + x + 1
        sage: T3 = RestrictionOfScalars(K,3); T3
        Algebraic torus of rank 12 over Rational Field split by a degree 4 extension

    The construction also works for non-Galois fields::

        sage: x = polygen(QQ);  K.<a> = NumberField(x^4 - x^3 + 3*x^2 - 2*x + 4)
        sage: RestrictionOfScalars(K)
        Algebraic torus of rank 4 over Rational Field split by a degree 8 extension

    and over finite fields::

        sage: RestrictionOfScalars(GF(9))
        Algebraic torus of rank 2 over Finite Field of size 3 split by a degree 2 extension
    """
    if field.is_galois():
        L = GLattice(rk).induced_lattice(field.galois_group())
        return AlgebraicTorus(L)
    else:
        G = field.galois_group()
        oG = G.order()
        oF = field.degree()
        subG = [h for h in G.conjugacy_classes_subgroups() if oG/h.order() == oF and not h.is_normal()]
        for H in subG:
            if any(h(e) == G.identity()(e) for e in field.gens() for h in H.gens()):
                L = GLattice(H, rk).induced_lattice(G)
                return AlgebraicTorus(L)


def NormOneRestrictionOfScalars(field, rk=1):
    r"""
    Create the torus of norm one elements inside the restriction of scalars of a split torus to the base field (rationals, rational function field, `\QQ_p` or ``GF(p)``).

    INPUT:

    - ``field`` -- Field extension over which the restriction of scalars is taken.  Usually a number field, function field, p-adic extension field or finite field.

    - ``rk`` -- rank (default ``1``) of the split torus of which we take the restriction of scalars.

    EXAMPLES::

        sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
        sage: F.<a> = QuadraticField([2])
        sage: T1 = NormOneRestrictionOfScalars(F); T1
        Algebraic torus of rank 1 over Rational Field split by a degree 2 extension
        sage: T1.character_lattice()._action_matrices
        [[-1]]
        sage: K.<b> = NumberField([x^4+1])
        sage: T2 = NormOneRestrictionOfScalars(K); T2
        Algebraic torus of rank 3 over Rational Field split by a degree 4 extension
        sage: T3 = NormOneRestrictionOfScalars(K,3); T3
        Algebraic torus of rank 9 over Rational Field split by a degree 4 extension

    The construction also works for non-Galois fields::

        sage: x = polygen(QQ);  K.<a> = NumberField(x^4 - x^3 + 3*x^2 - 2*x + 4)
        sage: NormOneRestrictionOfScalars(K)
        Algebraic torus of rank 3 over Rational Field split by a degree 8 extension

    and for finite fields::

        sage: NormOneRestrictionOfScalars(GF(9))
        Algebraic torus of rank 1 over Finite Field of size 3 split by a degree 2 extension
    """
    if field.is_galois():
        L = GLattice(rk).norm_one_restriction_of_scalars(field.galois_group())
    else:
        G = field.galois_group()
        H = G.stabilizer()
        L = GLattice(H, rk).norm_one_restriction_of_scalars(G)
    return AlgebraicTorus(L)


class AlgebraicTorus(Scheme):
    r"""
    Create an algebraic torus through its equivalence of categories with the action
    of a Galois group on an integral lattice.

    INPUT:

    - ``lattice`` -- the character lattice with Galois action
      defining the torus.
    """
    def __init__(self, lattice):
        r"""
        Construct an object of the algebraic torus class.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(DihedralGroup(4), 4)
            sage: AlgebraicTorus(L)
            Split algebraic torus of rank 4
        """
        Scheme.__init__(self)
        self._galois_group = lattice.group()
        self._lattice = lattice

    def _repr_(self):
        r"""
        The print representation of an algebraic torus.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: AlgebraicTorus(GLattice(CyclicPermutationGroup(3), 2))
            Split algebraic torus of rank 2

        ::

            sage: L = GLattice(CyclicPermutationGroup(3), 2).norm_one_restriction_of_scalars(SymmetricGroup(3
            ....: ))
            sage: AlgebraicTorus(L)
            Algebraic torus of rank 2 split by a degree 2 extension
        """
        split = self._lattice.action_is_trivial()
        if split:
            s = "Split algebraic"
        else:
            s = "Algebraic"
        s += " torus of rank %s" % self.rank()
        if self._base_field is not None:
            s += " over %s" % (self._base_field)
        if not split:
            s += " split by a degree %s extension" % self.splitting_degree()
        return s

    def rank(self):
        r"""
        The rank of the torus.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice([])
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0, 1, 0, 0, 0, 1, 1, 0, 0])
            sage: act2 = matrix(3, [0, 1, 0, 1, 0, 0, 0, 0, 1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.rank()
            1
            sage: T2.rank()
            1
            sage: T3.rank()
            3
        """
        return self._lattice._rank

    @lazy_attribute
    def _top_field(self):
        """
        The top field of the Galois group acting on the torus.

        EXAMPLES::

            sage: L.<a , b , c , d> = NumberField([x^2-5, x^2-29 , x^2-109 , x^2-281])
            sage: K = L.absolute_field('e')
            sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
            sage: T = NormOneRestrictionOfScalars(K)
            sage: T._top_field == K
            True
        """
        return self._galois_group.top_field()

    @lazy_attribute
    def _splitting_field(self):
        """
        The splitting field of the algebraic torus.

        This is the minimal Galois extension of the base field so that the induced action
        on the character lattice is trivial.

        EXAMPLES::

            sage: L.<a , b> = NumberField([x^2-5, x^2-3])
            sage: K = L.absolute_field('e')
            sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
            sage: T1 = NormOneRestrictionOfScalars(K)
            sage: T1._splitting_field
            Number Field in e with defining polynomial x^4 - 16*x^2 + 4
            sage: G = T1.galois_group()
            sage: H = G.subgroup([G.gens()[0]])
            sage: lat1 = GLattice(G, 1)
            sage: lat2 = GLattice(H, 1).induced_lattice(G)
            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: T2 = AlgebraicTorus(lat1); T2
            Split algebraic torus of rank 1 over Rational Field
            sage: T3 = AlgebraicTorus(lat2); T3
            Algebraic torus of rank 2 over Rational Field split by a degree 2 extension
            sage: T2._splitting_field
            Rational Field
            sage: T3._top_field
            Number Field in e with defining polynomial x^4 - 16*x^2 + 4
            sage: T3._splitting_field
            Number Field in e0 with defining polynomial x^2 - 3 with e0 = 1/4*e^3 - 7/2*e
        """
        ker = self.character_lattice().action_kernel()
        if ker.order() == 1:
            # In this case the splitting field is the full splitting field for the group.
            return self._galois_group.splitting_field()
        return ker.fixed_field()[0]

    @lazy_attribute
    def _base_field(self):
        """
        The field of definition of the torus.

        EXAMPLES::

            sage: L.<a , b , c , d> = NumberField([x^2-5, x^2-29 , x^2-109 , x^2-281])
            sage: K = L.absolute_field('e')
            sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
            sage: T = NormOneRestrictionOfScalars(K)
            sage: T._base_field
            Rational Field
        """
        return self._top_field.base_field()

    def splitting_degree(self):
        """
        The degree of the splitting field of the torus. It is used in :meth:`_repr_`.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: G = QuaternionGroup()
            sage: Z = G.center()
            sage: L = GLattice(Z, 1)
            sage: IL = L.induced_lattice(G)
            sage: T = AlgebraicTorus(IL)
            sage: T # indirect doctest
            Algebraic torus of rank 4 split by a degree 4 extension
            sage: T2 = AlgebraicTorus(GLattice(G, 1)); T2
            Split algebraic torus of rank 1
            sage: T3 = AlgebraicTorus(GLattice(G)); T3
            Algebraic torus of rank 8 split by a degree 8 extension
        """
        return self._lattice._group.order() // self._lattice.action_kernel().cardinality()

    def galois_group(self):
        r"""
        The abstract Galois group of a splitting field of the torus.

        .. NOTE::

            It doesn't have to be a minimal splitting field, therefore we
            allow trivial Galois action.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1,act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.galois_group()
            Permutation Group with generators [()]
            sage: T2.galois_group()
            Permutation Group with generators [(1,2,3), (1,2)]
            sage: T3.galois_group()
            Permutation Group with generators [(1,2,3), (1,2)]
        """
        return self._lattice._group

    def is_rational(self, ruple, subgroup=None):
        r"""
        Detect if a point is rational over the base field.

        INPUT:

        - ``ruple`` -- a point of the torus given as an r-tuple of points over the
          splitting field, where r is the rank of the torus

        - ``subgroup`` -- optional, subgroup of the Galois group corresponding to an
          intermediate extension; if specified, the algorithm will test the rationality
          over the intermediate extension

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
            sage: K.<w> = QuadraticField(5)
            sage: T = NormOneRestrictionOfScalars(K)
            sage: T.is_rational([(w-1)/2])
            True
            sage: T.is_rational([-1])
            True
            sage: T.is_rational([w])
            False
            sage: T.is_rational([2])
            False
        """
        lat = self.cocharacter_lattice()
        if subgroup is None:
            subgroup = self.galois_group()
        else:
            lat = lat.subgroup_lattice(subgroup)
        elt = matrix(self.rank(), ruple)
        gens = subgroup.gens()
        for g in range(len(gens)):
            res = lat._action_matrices[g] * matrix(self.rank(), [gens[g](x) for x in ruple])
            if res != elt:
                return False
        return True

    def character_lattice(self):
        r"""
        The character lattice of the torus.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.character_lattice()
            Ambient lattice of rank 1 with an action by a group of order 1
            sage: T2.character_lattice()
            Ambient lattice of rank 1 with the trivial action of a group of order 6
            sage: T3.character_lattice()
            Ambient lattice of rank 3 with a faithful action by a group of order 6
        """
        return self._lattice

    def cocharacter_lattice(self):
        r"""
        The cocharacter lattice of the torus.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T3.character_lattice()._action_matrices
            [
            [0 1 0]  [0 1 0]
            [0 0 1]  [1 0 0]
            [1 0 0], [0 0 1]
            ]
            sage: T3.cocharacter_lattice()._action_matrices
            [
            [0 1 0]  [0 1 0]
            [0 0 1]  [1 0 0]
            [1 0 0], [0 0 1]
            ]

        The matrices are all orthogonal so we get the same lattice. The action will be
        different in the next example.

        ::

            sage: L = GLattice(PermutationGroup([(1,2,3,4,5,6)]), [matrix(2, [0,1,-1,-1])]); L
            Ambient lattice of rank 2 with an action by a group of order 6
            sage: T = AlgebraicTorus(L); T
            Algebraic torus of rank 2 split by a degree 3 extension
            sage: T.character_lattice()._action_matrices
            [
            [ 0  1]
            [-1 -1]
            ]
            sage: T.cocharacter_lattice()._action_matrices
            [
            [-1  1]
            [-1  0]
            ]
        """
        return self._lattice.colattice()

    def tate_cohomology(self, n):
        r"""
        Gives the isomorphism type of the nth cohomology group using Tate-Nakayama duality.
        Only works when the base field is local.

        INPUT:

        - ``n`` -- which cohomology group to compute

        .. NOTE::

            This currently only works for tori over local p-adic fields. For global fields,
            Tate-Nakayama gives the cohomology of the class group, not the torus itself.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+" : ", T1. tate_cohomology (i)) # optional - gap_packages
            H^-5 :  []
            H^-4 :  []
            H^-3 :  []
            H^-2 :  []
            H^-1 :  []
            H^0 :  []
            H^1 :  []
            H^2 :  []
            H^3 :  []
            H^4 :  []
            H^5 :  []

        The Galois group is trivial and has obviously trivial cohomology

        ::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+" : ", T2. tate_cohomology (i)) # optional - gap_packages
            H^-5 :  []
            H^-4 :  [2]
            H^-3 :  []
            H^-2 :  [6]
            H^-1 :  []
            H^0 :  [2]
            H^1 :  []
            H^2 :  [6]
            H^3 :  []
            H^4 :  [2]
            H^5 :  []

        We can recognize from class field theory that `H^2` (the Brauer group of
        our extension) is isomorphic to the cyclic group Cn where n is the order of the
        Galois group. Here the group is `S_3`, which has order 6 so we get `C_6`.

        Another way to see it is seeing this `H^2` as `H^0` with coefficients in its character lattice. Since the group
        acts trivially, the fixed elements are the whole lattice, and the trace map is multiplication
        by the order of the group, which is 6, so we get `C_6` ^(rank of ``T1``)

        Also, `H^0` can be seen as the abelianization of the Galois group (indeed, by Tate-Nakayama it is
        `H^2(\ZZ)`), which here has order 2 (it is the group of signatures)::

            sage: for i in range(-5, 6):
            ....:     print("H^"+str(i)+" : ", T3. tate_cohomology (i)) # optional - gap_packages
            H^-5 :  []
            H^-4 :  [2]
            H^-3 :  []
            H^-2 :  [2]
            H^-1 :  []
            H^0 :  [2]
            H^1 :  []
            H^2 :  [2]
            H^3 :  []
            H^4 :  [2]
            H^5 :  []

        In this example, we can see the 2-periodicity of the cohomology groups, consequence
        of the group being cyclic.
        """
        return self._lattice.tate_cohomology(2 - n)

    def product(self, torus):
        r"""
        Return the product of two tori defined and (for now) also splitting over the same field.

        INPUT:

        - ``torus`` -- another algebraic torus

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: T1 = AlgebraicTorus(GLattice([], 1));
            sage: G = SymmetricGroup(3)
            sage: T2 = T1.norm_one_restriction(G)
            sage: T3 = T1.restriction_of_scalars(G)
            sage: PT = T2.product(T3); PT; PT.character_lattice()._action_matrices
            Algebraic torus of rank 11 split by a degree 6 extension
            [
            [ 0 -1  0  0  1| 0  0  0  0  0  0]  [ 0  0  1 -1  0| 0  0  0  0  0  0]
            [ 0 -1  1  0  0| 0  0  0  0  0  0]  [ 0  0  0 -1  1| 0  0  0  0  0  0]
            [ 0 -1  0  0  0| 0  0  0  0  0  0]  [ 1  0  0 -1  0| 0  0  0  0  0  0]
            [ 1 -1  0  0  0| 0  0  0  0  0  0]  [ 0  0  0 -1  0| 0  0  0  0  0  0]
            [ 0 -1  0  1  0| 0  0  0  0  0  0]  [ 0  1  0 -1  0| 0  0  0  0  0  0]
            [--------------+-----------------]  [--------------+-----------------]
            [ 0  0  0  0  0| 0  0  0  0  1  0]  [ 0  0  0  0  0| 0  0  1  0  0  0]
            [ 0  0  0  0  0| 0  0  1  0  0  0]  [ 0  0  0  0  0| 0  0  0  0  1  0]
            [ 0  0  0  0  0| 0  0  0  0  0  1]  [ 0  0  0  0  0| 1  0  0  0  0  0]
            [ 0  0  0  0  0| 1  0  0  0  0  0]  [ 0  0  0  0  0| 0  0  0  0  0  1]
            [ 0  0  0  0  0| 0  0  0  1  0  0]  [ 0  0  0  0  0| 0  1  0  0  0  0]
            [ 0  0  0  0  0| 0  1  0  0  0  0], [ 0  0  0  0  0| 0  0  0  1  0  0]
            ]
        """
        L = self.character_lattice()
        M = torus.character_lattice()
        G = L.group()
        H = M.group()
        if G is H:
            return AlgebraicTorus(L.direct_sum(M))
        if L._base_field is not M._base_field:
            raise ValueError("Tori must be defined over the same field")
        raise NotImplementedError("Product of tori only defined when they have the same splitting field")

    def restriction_of_scalars(self, group):
        r"""
        The torus obtained through restriction of scalars.

        INPUT:

        - ``group`` -- the bigger group corresponding the the Galois group
          of the splitting field over the subfield one wishes to restrict scalars.

        .. NOTE::

            The user has to input the (larger) Galois group for this extension.
            More concretely, if the torus is defined over K, splits over L, and
            is defined by the action of Gal(L/K) on its character lattice, then
            if one wants the restriction of scalars to a smaller field k,
            one has to enter Gal(L/k) as argument of this method. This is because
            galois groups of relative extensions are not supported so far. For
            now the user has to input the Galois group of the relative extension
            as abstract group.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)
            sage: act1 = matrix(3, [0,1,0,0,0,1,1,0,0])
            sage: act2 = matrix(3, [0,1,0,1,0,0,0,0,1])
            sage: LLL = GLattice(SymmetricGroup(3), [act1, act2])
            sage: T3 = AlgebraicTorus(LLL)

        ::

            sage: T1.restriction_of_scalars(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            Algebraic torus of rank 16 split by a degree 16 extension

        ::

            sage: T2.restriction_of_scalars(SymmetricGroup(4))
            Algebraic torus of rank 4 split by a degree 24 extension
            sage: _.character_lattice()._action_matrices
            [
            [0|0|0|1]  [1|0|0|0]
            [-+-+-+-]  [-+-+-+-]
            [1|0|0|0]  [0|1|0|0]
            [-+-+-+-]  [-+-+-+-]
            [0|1|0|0]  [0|0|0|1]
            [-+-+-+-]  [-+-+-+-]
            [0|0|1|0], [0|0|1|0]
            ]

        ::

            sage: T3.restriction_of_scalars(SymmetricGroup(4))
            Algebraic torus of rank 12 split by a degree 24 extension
            sage: _.character_lattice()._action_matrices
            [
            [0 0 0|0 0 0|0 0 0|1 0 0]  [0 1 0|0 0 0|0 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 0|0 1 0]  [1 0 0|0 0 0|0 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 0|0 0 1]  [0 0 1|0 0 0|0 0 0|0 0 0]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 1 0|0 0 0|0 0 0|0 0 0]  [0 0 0|0 1 0|0 0 0|0 0 0]
            [0 0 1|0 0 0|0 0 0|0 0 0]  [0 0 0|1 0 0|0 0 0|0 0 0]
            [1 0 0|0 0 0|0 0 0|0 0 0]  [0 0 0|0 0 1|0 0 0|0 0 0]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 0 0|0 1 0|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|1 0 0]
            [0 0 0|0 0 1|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|0 1 0]
            [0 0 0|1 0 0|0 0 0|0 0 0]  [0 0 0|0 0 0|0 0 0|0 0 1]
            [-----+-----+-----+-----]  [-----+-----+-----+-----]
            [0 0 0|0 0 0|0 1 0|0 0 0]  [0 0 0|0 0 0|1 0 0|0 0 0]
            [0 0 0|0 0 0|0 0 1|0 0 0]  [0 0 0|0 0 0|0 1 0|0 0 0]
            [0 0 0|0 0 0|1 0 0|0 0 0], [0 0 0|0 0 0|0 0 1|0 0 0]
            ]
        """
        return AlgebraicTorus(self._lattice.induced_lattice(group))

    def tamagawa_number(self, subgroups = []):
        """
        Computes the Tamagwa number of an algebraic torus.

        INPUT:

        - ``subgroups`` -- To be used when the group is not the Galois group
          of some number field, but just a permutation group. In that case
          this argument should be the list of ramified decomposition groups,
          which have to be manually added. If left blank, the method will assume
          that all decomposition groups are cyclic.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import NormOneRestrictionOfScalars
            sage: F.<a> = QuadraticField([2])
            sage: T = NormOneRestrictionOfScalars(F)
            sage: T.tamagawa_number() # optional - gap_packages
            2
            sage: from sage.schemes.group_schemes.tori import RestrictionOfScalars
            sage: T2 = RestrictionOfScalars(F)
            sage: T2.tamagawa_number() # optional - gap_packages
            1

        ::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: G = QuaternionGroup()
            sage: Z = G.center()
            sage: IL = GLattice(1).induced_lattice(G)
            sage: SL = IL.fixed_sublattice(Z).zero_sum_sublattice()
            sage: L = IL.quotient_lattice(SL); L
            Ambient lattice of rank 5 with a faithful action by a group of order 8
            sage: T2 = AlgebraicTorus(L)
            sage: T2.tamagawa_number() # optional - gap_packages
            1/2

        ::

            sage: G = GLattice([2, 2, 2, 2]).group()
            sage: L = GLattice(1).norm_one_restriction_of_scalars(G)
            sage: T = AlgebraicTorus(L)
            sage: T.tamagawa_number() # optional - gap_packages
            1/4

        The latter example is the Tamagawa number computed by Ono, the first example of non-integral Tamagawa number.
        """
        from sage.misc.misc_c import prod
        if self._base_field is None:
            ram_decomp = subgroups
            lat2 = self.character_lattice()
        else:
            from sage.groups.perm_gps.permgroup import PermutationGroup
            perm_group = PermutationGroup(self.galois_group().gens())
            ram_primes = self._top_field.discriminant().prime_divisors()
            ram_decomp = []
            for p in ram_primes:
                SG = self._top_field.prime_above(p).decomposition_group()
                SGperm = PermutationGroup(SG.gens())
                if not(SGperm.is_cyclic()):
                    ram_decomp += [SGperm]
            #unram_decomp = [h for h in perm_group.conjugacy_classes_subgroups() if h.is_cyclic()]
            lat = self.character_lattice()
            lat2 = GLattice(perm_group, lat._action_matrices)
        num = prod(self.tate_cohomology(1))
        denom = prod(lat2.Tate_Shafarevich_lattice(2))
        if ram_decomp == [] or denom == 1:
            return num/denom
        else:
            denom = lat2.character_lattice().Tate_Shafarevich_lattice(2, ram_decomp)[0]
            return num/denom



    def norm_one_restriction(self,group):
        r"""
        Torus of norm one elements in the restriction of scalars.

        INPUT:

        - ``group`` -- the abstract Galois group of the restriction of scalars.

        EXAMPLES::

            sage: from sage.schemes.group_schemes.tori import AlgebraicTorus
            sage: L = GLattice(PermutationGroup([()]), 1)
            sage: T1 = AlgebraicTorus(L)
            sage: LL = GLattice(SymmetricGroup(3), 1)
            sage: T2 = AlgebraicTorus(LL)

        ::

            sage: T1.norm_one_restriction(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            Algebraic torus of rank 15 split by a degree 16 extension
            sage: _.character_lattice()._action_matrices[0]
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1]
            [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0 -1]
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0 -1]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0 -1]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0 -1]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0  0 -1]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0  0 -1]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1]
            sage: T2.norm_one_restriction(SymmetricGroup(4))
            Algebraic torus of rank 3 split by a degree 24 extension
            sage: _.character_lattice()._action_matrices
            [
            [-1  1  0]  [0 1 0]
            [-1  0  1]  [1 0 0]
            [-1  0  0], [0 0 1]
            ]
            

        We now compute the cohomologies of all those tori.
        For the two latter examples, we check that we get different cohomologies with the
        same group, and the norm one restriction of the split torus has nontrivial cohomology::

            sage: ROS = T1.norm_one_restriction(PermutationGroup([(1,2), (3,4), (5,6), (7,8)]))
            sage: for i in range(-4, 6):
            ....:     print("H^"+str(i)+" : ", ROS.tate_cohomology (i)) # optional - gap_packages
            H^-4 :  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^-3 :  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^-2 :  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^-1 :  [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            H^0 :  [2, 2, 2, 2, 2, 2]
            H^1 :  [2, 2, 2, 2]
            H^2 :  []
            H^3 :  [2, 2, 2, 2]
            H^4 :  []
            H^5 :  [2, 2, 2, 2]

        This torus ``ROS`` is the example of Ono where he applies his formula for the Tamagawa
        number of a Torus. See his paper 'On the Tamagawa Number of Algebraic Tori'.

        ::

            sage: ROS2 = T2.norm_one_restriction(SymmetricGroup(4))
            sage: for i in range(-4, 6):
            ....:     print("H^"+str(i)+" : ", T2.tate_cohomology (i)) # optional - gap_packages
            H^-4 :  [2]
            H^-3 :  []
            H^-2 :  [6]
            H^-1 :  []
            H^0 :  [2]
            H^1 :  []
            H^2 :  [6]
            H^3 :  []
            H^4 :  [2]
            H^5 :  []
        """
        return AlgebraicTorus(self._lattice.norm_one_restriction_of_scalars(group))
