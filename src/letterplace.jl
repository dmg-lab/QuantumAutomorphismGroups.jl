
using Oscar


# freeAlgebra and makeLetterplaceRing are the same?

strng = """
LIB "freegb.lib";
ring r = 0,(x,y,z),dp;
def R = freeAlgebra(r, 4);  // degree (length) bound 4; the ordering will be degree right lex
setring R;
ideal I = x*y + y*z, x*x + x*y - z; // a non-graded ideal
ideal J = letplaceGBasis(I);
"""



Singular.call_interpreter(strng);

ring, dct = Singular.lookup_library_symbol("Top", "R")
SI = dct[:J]



R, (x, y, z) = free_associative_algebra(QQ, ["x", "y","z"])

f1 = x*y + y*z
f2 = x*x + x*y - z




