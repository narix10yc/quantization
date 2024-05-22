#ifndef QUANTIZATION_TYPES_H
#define QUANTIZATION_TYPES_H

#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>

template<unsigned width>
class quantized {
    int64_t data;
public:
    quantized() : data() {}
    quantized(double value) : data(value * (1 << width)) {}
    quantized(uint64_t) = delete;

    double getValue() const {
        return static_cast<double>(data) / (1 << width);
    }

    quantized operator+(const quantized& other) const {
        quantized q;
        q.data = data + other.data;
        return q;
    }

    quantized operator-(const quantized& other) const {
        quantized q;
        q.data = data - other.data;
        return q;
    }

    quantized operator*(const quantized& other) const {
        quantized q;
        q.data = (data * other.data) >> width;
        return q;
    }

    quantized operator/(const quantized& other) const {
        quantized q;
        q.data = data / (other.data >> width);
        return q;
    }

    bool operator>(double n) const {
        return getValue() > n;
    }

    bool operator>=(double n) const {
        return getValue() >= n;
    }

    bool operator<=(double n) const {
        return getValue() <= n;
    }

    quantized& operator+=(const quantized& other) {
        data += other.data;
        return *this;
    }

    quantized& operator-=(const quantized& other) {
        data -= other.data;
        return *this;
    }

    quantized& operator*=(const quantized& other) {
        data *= other.data;
        data >>= width;
        return *this;
    }

    quantized& operator/=(const quantized& other) {
        data <<= width;
        data /= other.data;
        return *this;
    }

    operator double() const {
        return getValue();
    }

    friend std::ostream& operator<<(std::ostream& os, const quantized& q) {
        os << q.getValue();
        return os;
    }

};

template<typename real_ty>
class Complex {
public:
    real_ty real, imag;
    Complex() : real(0.0), imag(0.0) {}
    Complex(real_ty real, real_ty imag) : real(real), imag(imag) {}

    Complex operator+(const Complex& other) const {
        return Complex(real + other.real, imag + other.imag);
    }

    Complex operator-(const Complex& other) const {
        return Complex(real - other.real, imag - other.imag);
    }

    Complex operator*(const Complex& other) const {
        return Complex(real * other.real - imag * other.imag,
                       real * other.imag + imag * other.real);
    }

    Complex& operator+=(const Complex& other) {
        real += other.real;
        imag += other.imag;
        return *this;
    }

    Complex& operator-=(const Complex& other) {
        real -= other.real;
        imag -= other.imag;
        return *this;
    }

    Complex& operator*=(const Complex& other) {
        real_ty r = real * other.real - imag * other.imag;
        imag = real * other.imag + imag * other.real;
        real = r;
        return *this;
    }

    double normSquared() const {
        return static_cast<double>(real) * real + static_cast<double>(imag) * imag;
    }

    double norm() const {
        return sqrt(normSquared());
    }

    friend std::ostream& operator<<(std::ostream& os, const Complex& c) {
        if (c.real >= 0.0)
            os << " ";
        os << c.real;
        if (c.imag >= 0.0)
            os << "+";
        os << c.imag << "i";
        return os;
    }
};


template<typename real_ty>
class SquareComplexMatrix {
    size_t size;
public:
    std::vector<Complex<real_ty>> data;

    SquareComplexMatrix(size_t size) : size(size), data(size * size) {}
    SquareComplexMatrix(size_t size, std::initializer_list<Complex<real_ty>> data) 
            : size(size), data(data) {}

    size_t getSize() const { return size; }

    static SquareComplexMatrix Identity(size_t size) {
        SquareComplexMatrix m(size);
        for (size_t r = 0; r < size; r++)
            m.data[r*size + r].real = 1.0;
        return m;
    }

    static SquareComplexMatrix H() {
        return { 2, {{ 0.7071067811865476, 0.0}, { 0.7071067811865476, 0.0},
                     { 0.7071067811865476, 0.0}, {-0.7071067811865476, 0.0}}};
    }

    static SquareComplexMatrix X() {
        return { 2, {{0.0, 0.0}, {1.0, 0.0},
                     {1.0, 0.0}, {0.0, 0.0}}};
    }

    static SquareComplexMatrix Y() {
        return { 2, {{0.0, 0.0}, {0.0, -1.0},
                     {0.0, 1.0}, {0.0,  0.0}}};
    }

    static SquareComplexMatrix Z() {
        return { 2, {{1.0, 0.0}, {-1.0, 0.0},
                     {0.0, 0.0}, { 0.0, 0.0}}};
    }

    static SquareComplexMatrix P(double phi) {
        return { 2, {{1.0, 0.0}, {0.0, 0.0},
                     {0.0, 0.0}, {cos(phi), sin(phi)}}};    
    }

    static SquareComplexMatrix U3(double theta, double phi, double lambd) {
        return { 2, {
            {cos(0.5*theta), 0.0},
            {-cos(lambd) * sin(0.5*theta), -sin(lambd) * sin(0.5*theta)},
            {cos(phi) * sin(0.5*theta), sin(phi) * sin(0.5*theta)},
            {cos(phi+lambd) * cos(0.5*theta), sin(phi+lambd) * cos(0.5*theta)}}
        };    
    }

    static SquareComplexMatrix CX() {
        return { 4, {
            {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0},
            {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}
        }};
    }

    static SquareComplexMatrix CZ() {
        return { 4, {
            {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, { 0.0, 0.0},
            {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, { 0.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, { 0.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {-1.0, 0.0}
        }};
    }

    static SquareComplexMatrix CP(double phi) {
        return { 4, {
            {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
            {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {0.0, 0.0},
            {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {cos(phi), sin(phi)}
        }};
    }

    SquareComplexMatrix matmul(const SquareComplexMatrix& other) {
        if (size != other.size)
            throw std::runtime_error("matrix size mismatch!");

        SquareComplexMatrix m(size);
        for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
        for (size_t k = 0; k < size; k++) {
            // C_{ij} = A_{ik} B_{kj}
            m.data[i*size + j] += data[i*size + k] * other.data[k*size + j];
        } } }
        return m;
    }

    SquareComplexMatrix kron(const SquareComplexMatrix& other) const {
        size_t lsize = size;
        size_t rsize = other.size;
        size_t msize = lsize * rsize;
        SquareComplexMatrix m(msize);
        for (size_t lr = 0; lr < lsize; lr++) {
        for (size_t lc = 0; lc < lsize; lc++) {
        for (size_t rr = 0; rr < rsize; rr++) {
        for (size_t rc = 0; rc < rsize; rc++) {
            size_t r = lr * rsize + rr;
            size_t c = lc * rsize + rc;
            m.data[r*msize + c] = data[lr*lsize + lc] * other.data[rr*rsize + rc];
        } } } }
        return m;
    }

    SquareComplexMatrix leftKronI() const {
        SquareComplexMatrix m(size * size);
        for (size_t i = 0; i < size; i++) {
        for (size_t r = 0; r < size; r++) {
        for (size_t c = 0; c < size; c++) {
            m.data[(i*size + r) * size * size + (i*size + c)] = data[r*size + c];
        } } }
        return m;
    }

    SquareComplexMatrix rightKronI() const {
        SquareComplexMatrix m(size * size);
        for (size_t i = 0; i < size; i++) {
        for (size_t r = 0; r < size; r++) {
        for (size_t c = 0; c < size; c++) {
            m.data[(r*size + i) * size * size + (c*size + i)] = data[r*size + c];
        } } }
        return m;
    }

    SquareComplexMatrix swapTargetQubits() const {
        if (size != 4)
            throw std::runtime_error("swapTargetQubits only valid to 4x4 matrices");
        return { data[ 0], data[ 2], data[ 1], data[ 3],
                 data[ 8], data[10], data[ 9], data[11],
                 data[ 4], data[ 6], data[ 5], data[ 7],
                 data[12], data[14], data[13], data[15] };
    }

    void print(std::ostream& os) const {
        for (size_t r = 0; r < size; r++) {
            for (size_t c = 0; c < size; c++) {
                os << data[r*size + c] << " ";
            }
            os << "\n";
        }
    }

};


template<typename real_ty>
class Statevector {
public:
    unsigned nqubits;
    uint64_t N;
    Complex<real_ty>* data;

    Statevector(unsigned nqubits, bool initialize=false)
            : nqubits(nqubits), N(1 << nqubits) {
        data = new Complex<real_ty>[N];
        if (initialize) {
            for (size_t i = 0; i < N; i++)
                data[i] = { 0.0, 0.0 };
            data[0].real = 1.0;
        }
    }

    Statevector(const Statevector& that) : nqubits(that.nqubits), N(that.N) {
        data = new Complex<real_ty>[N];
        for (size_t i = 0; i < that.N; i++)
            data[i] = that.data[i];
    }

    Statevector(Statevector&&) = delete;

    ~Statevector() { delete[](data); }

    Statevector& operator=(const Statevector& that) {
        if (this != &that) {
            for (size_t i = 0; i < N; i++)
                data[i] = that.data[i];
        }
        return *this;
    }

    Statevector& operator=(Statevector&&) = delete;

    template<typename another_real_ty>
    void copyValuesFrom(Statevector<another_real_ty> other) {
        for (size_t i = 0; i < N; i++) {
            data[i].real = static_cast<double>(other.data[i].real);
            data[i].imag = static_cast<double>(other.data[i].imag);
        }
    }

    double normSquared() const {
        double s = 0;
        for (size_t i = 0; i < N; i++)
            s += data[i].normSquared();
        return s;
    }

    double norm() const { 
        return sqrt(normSquared());
    }

    void normalize() {
        double n = norm();
        for (size_t i = 0; i < N; i++) {
            // data[i].real = static_cast<double>(data[i].real) / n;
            // data[i].imag = static_cast<double>(data[i].imag) / n;
            data[i].real /= n;
            data[i].imag /= n;
        } 
    }

    void randomize() {
        std::random_device rd;
        std::mt19937 gen { rd() };
        std::normal_distribution<double> d { 0, 1 };
        for (size_t i = 0; i < N; i++) {
            data[i].real = d(gen);
            data[i].imag = d(gen);
        }
        normalize();
    }

    void 
    applySingleQubit(const SquareComplexMatrix<real_ty>& m, unsigned k) {
        if (m.getSize() != 2)
            throw std::runtime_error("single-qubit gates are 2x2 matrices");

        size_t K = 1 << k;
        Complex<real_ty> lo, hi;
        for (size_t t = 0; t < N; t += 2*K) {
        for (size_t tt = 0; tt < K; tt++) {
            lo = m.data[0] * data[t+tt] + m.data[1] * data[t+tt+K];
            hi = m.data[2] * data[t+tt] + m.data[3] * data[t+tt+K];
            data[t+tt] = lo;
            data[t+tt+K] = hi;
        } }
    }

    void 
    applyTwoQubit(const SquareComplexMatrix<real_ty>& m, unsigned k, unsigned l) {
        if (m.getSize() != 4)
            throw std::runtime_error("two-qubit gates are 4x4 matrices");

        // K is more significant than L (i.e. K > L)
        size_t K = 1 << k;
        size_t L = 1 << l;

        Complex<real_ty> v00, v01, v10, v11;
        Complex<real_ty> *p00, *p01, *p10, *p11;
        if (k > l) {
            for (size_t t = 0; t < N; t += 2*K) {
            for (size_t tt = 0; tt < K; tt += 2*L) {
            for (size_t ttt = 0; ttt < L; ttt++) {
                p00 = data + (t+tt+ttt);   p01 = data + (t+tt+ttt+L);
                p10 = data + (t+tt+ttt+K); p11 = data + (t+tt+ttt+L+K);

                v00 = m.data[ 0] * *p00 + m.data[ 1] * *p01 
                     +m.data[ 2] * *p10 + m.data[ 3] * *p11;
                v01 = m.data[ 4] * *p00 + m.data[ 5] * *p01
                     +m.data[ 6] * *p10 + m.data[ 7] * *p11;
                v10 = m.data[ 8] * *p00 + m.data[ 9] * *p01
                     +m.data[10] * *p10 + m.data[11] * *p11;
                v11 = m.data[12] * *p00 + m.data[13] * *p01
                     +m.data[14] * *p10 + m.data[15] * *p11;

                *p00 = v00; *p01 = v01; *p10 = v10; *p11 = v11;
            } } }
        } else {
            for (size_t t = 0; t < N; t += 2*L) {
            for (size_t tt = 0; tt < L; tt += 2*K) {
            for (size_t ttt = 0; ttt < K; ttt++) {
                p00 = data + (t+tt+ttt);   p01 = data + (t+tt+ttt+L);
                p10 = data + (t+tt+ttt+K); p11 = data + (t+tt+ttt+L+K);

                v00 = m.data[ 0] * *p00 + m.data[ 2] * *p01
                     +m.data[ 1] * *p10 + m.data[ 3] * *p11;
                v01 = m.data[ 8] * *p00 + m.data[10] * *p01
                     +m.data[ 9] * *p10 + m.data[11] * *p11;
                v10 = m.data[ 4] * *p00 + m.data[ 6] * *p01
                     +m.data[ 5] * *p10 + m.data[ 7] * *p11;
                v11 = m.data[12] * *p00 + m.data[14] * *p01
                     +m.data[13] * *p10 + m.data[15] * *p11;
                
                *p00 = v00; *p01 = v01; *p10 = v10; *p11 = v11;
            } } }
        }
    }

    void applyX(unsigned k) {
        return applySingleQubit(SquareComplexMatrix<real_ty>::X(), k);
    }

    void applyY(unsigned k) {
        return applySingleQubit(SquareComplexMatrix<real_ty>::Y(), k);
    }

    void applyZ(unsigned k) {
        return applySingleQubit(SquareComplexMatrix<real_ty>::Z(), k);
    }

    void applyH(unsigned k) {
        return applySingleQubit(SquareComplexMatrix<real_ty>::H(), k);
    }

    void applyP(unsigned k, double phi) {
        return applySingleQubit(SquareComplexMatrix<real_ty>::P(phi), k);
    }

    void applyU3(unsigned k, double theta, double phi, double lambd) {
        return applySingleQubit(
            SquareComplexMatrix<real_ty>::U3(theta, phi, lambd), k);
    }

    void applyCX(unsigned k, unsigned l) {
        return applyTwoQubit(SquareComplexMatrix<real_ty>::CX(), k, l);
    }

    void applyCZ(unsigned k, unsigned l) {
        return applyTwoQubit(SquareComplexMatrix<real_ty>::CZ(), k, l);
    }

    void applyCP(unsigned k, unsigned l, double phi) {
        return applyTwoQubit(SquareComplexMatrix<real_ty>::CP(phi), k, l);
    }

    void print(std::ostream& os) const {
        if (N > 32) {
            const char* cyan = "\033[36m";
            const char* bold = "\033[1m";
            const char* reset = "\033[0m";
            os << bold << cyan << "Warning: " << reset << "statevector has "
                "more than 5 qubits, only the first 32 amplitudes are shown.\n";
        }
        for (size_t i = 0; i < ((N > 32) ? 32 : N); i++) {
            if (N >= 10 && i < 10)
                os << std::setw(2) << std::setfill(' ');
            os << i << ": " << data[i] << "\n";
        }
    }

};


#endif // QUANTIZATION_TYPES_H