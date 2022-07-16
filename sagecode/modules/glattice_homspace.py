import inspect
import sage.matrix.all as matrix
import sage.modules.free_module_homspace




def is_GLatticeHomspace(x):
    return isinstance(x, GLatticeHomspace)

class GLatticeHomspace(sage.modules.free_module_homspace.FreeModuleHomspace):

    def __call__(self, A, check=True, **kwds):

        from .glattice_morphism import is_GLatticeMorphism, GLatticeMorphism
        D = self.domain()
        C = self.codomain()
        side = kwds.get("side", "left")
        from sage.structure.element import is_Matrix
        if is_Matrix(A):
            pass
        elif is_GLatticeMorphism(A):
            A = A.matrix()
        elif inspect.isfunction(A):
            try:
                images = [A(g) for g in D.basis()]
            except (ValueError, TypeError, IndexError) as e:
                msg = 'function cannot be applied properly to some basis element because\n' + e.args[0]
                raise ValueError(msg)
            try:
                A = matrix.matrix(D.dimension(), C.dimension(), [C.coordinates(C(a)) for a in images])
            except (ArithmeticError, TypeError) as e:
                msg = 'some image of the function is not in the codomain, because\n' + e.args[0]
                raise ArithmeticError(msg)
            if side == "right":
                A = A.transpose()
        elif isinstance(A, (list, tuple)):
            if len(A) != len(D.basis()):
                msg = "number of images should equal the size of the domain's basis (={0}), not {1}"
                raise ValueError(msg.format(len(D.basis()), len(A)))
            try:
                v = [C(a) for a in A]
                A = matrix.matrix(D.dimension(), C.dimension(), [C.coordinates(a) for a in v])
            except (ArithmeticError, TypeError) as e:
                msg = 'some proposed image is not in the codomain, because\n' + e.args[0]
                raise ArithmeticError(msg)
            if side == "right":
                A = A.transpose()
        else:
            msg = 'lattice homspace can only coerce matrices, lattice morphisms, functions or lists, not {0}'
            raise TypeError(msg.format(A))
        return GLatticeMorphism(self, A, side=side)

    def _repr_(self):

        msg = 'Set of Morphisms from {0} to {1}'
        return msg.format(self.domain(), self.codomain())
