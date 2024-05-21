#include "types.h"

int main() {
    // using real_ty = float;
    using real_ty = quantized<14>;
    using Matrix = SquareComplexMatrix<real_ty>;

    auto hGate = Matrix::H();
    hGate.print(std::cerr);

    auto sv = Statevector<real_ty>(3, true);
    sv.randomize();
    sv.print(std::cerr);
    
    std::cerr << "norm=" << sv.norm() << "\n";

    sv.applySingleQubit(hGate, 0);
    // sv.applySingleQubit(hGate, 1);
    // sv.applySingleQubit(hGate, 2);

    sv.print(std::cerr);
    std::cerr << "\n";

    sv.applyTwoQubit(Matrix::CX(), 1, 0);
    sv.applyTwoQubit(Matrix::CZ(), 1, 2);

    sv.print(std::cerr);





    return 0;
}