#include "types.h"

int main() {
    using real_ty = quantized<15>;
    using real_ty_bl = double;
    using Matrix = SquareComplexMatrix<real_ty>;
    using Matrix_bl = SquareComplexMatrix<real_ty_bl>;

    Statevector<real_ty> sv(12, true);
    Statevector<real_ty_bl> sv_bl(12, true);

    sv_bl.randomize();
    sv.copyValuesFrom(sv_bl);

    sv.applyH(0);
    sv_bl.applyH(0);

    sv.applyCX(1, 0);
    sv_bl.applyCX(1, 0);

    sv_bl.print(std::cerr);
    sv.print(std::cerr);


    std::cerr << "norm = " << sv.norm() << "\n";


    return 0;
}