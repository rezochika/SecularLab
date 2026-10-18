#include "secular/exact.hpp"
#include "secular/solver.hpp"
#include <iostream>
#include <random>
#include <sstream>

using namespace secular;

int main() {
    unsigned checks = 0;
    auto check = [&](bool ok, const char* message) {
        ++checks;
        if (!ok) throw std::runtime_error(message);
    };
    try {
        std::mt19937 generator(20470311);
        auto small = [&]() { return static_cast<long>(generator() % 7) - 3; };
        for (std::size_t n = 1; n <= 8; ++n) {
            for (unsigned trial = 0; trial < 20; ++trial) {
                // A unit upper-triangular congruence preserves known inertia.
                auto p = identity<Z>(n);
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t j = i + 1; j < n; ++j) p(i, j) = small();
                ZMatrix d(n), transpose(n);
                Inertia expected;
                Z expected_determinant = 1;
                for (std::size_t i = 0; i < n; ++i) {
                    long value = small();
                    d(i, i) = value;
                    expected_determinant *= value;
                    if (value > 0) ++expected.positive;
                    else if (value < 0) ++expected.negative;
                    else ++expected.zero;
                    for (std::size_t j = 0; j < n; ++j) transpose(i, j) = p(j, i);
                }
                auto a = multiply(multiply(transpose, d), p);
                auto actual = inertia(a);
                check(actual.positive == expected.positive &&
                      actual.negative == expected.negative &&
                      actual.zero == expected.zero, "congruence inertia mismatch");
                // The diagonal sign count must agree with elimination on the congruent matrix.
                auto direct = inertia(d);
                check(direct.positive == actual.positive && direct.negative == actual.negative &&
                      direct.zero == actual.zero, "diagonal inertia differs from elimination");
                check(determinant(a) == expected_determinant, "congruence determinant mismatch");
                if (expected_determinant != 0) {
                    auto inverse_times_determinant = adjugate(a);
                    auto product = multiply(a, inverse_times_determinant);
                    auto wanted = identity<Z>(n);
                    for (auto& v : wanted.data) v *= expected_determinant;
                    check(product.data == wanted.data, "adjugate identity mismatch");
                }
            }
        }
        for (std::size_t n = 1; n <= 6; ++n) {
            for (unsigned trial = 0; trial < 12; ++trial) {
                auto e = identity<Z>(n);
                ZMatrix f(n), f1(n);
                for (std::size_t i = 0; i < n; ++i) {
                    e(i, i) = (trial % 2 ? -2L : 3L);
                    for (std::size_t j = 0; j < n; ++j) {
                        if (j > i) e(i, j) = small();
                        f(i, j) = small();
                        f1(i, j) = small();
                    }
                }
                const auto polynomial = characteristic(e, f);
                check(polynomial.terminal_zero && verify_coefficients(e, f, polynomial.coefficients),
                      "general nonsymmetric pencil coefficient mismatch");
                if (n <= 3) {
                    const auto affine = affine_characteristic(e, f, f1);
                    for (long t : {-2L, 0L, 3L}) {
                        auto specialized_matrix = f;
                        for (std::size_t i = 0; i < f.data.size(); ++i)
                            specialized_matrix.data[i] += t * f1.data[i];
                        std::vector<Z> specialized;
                        for (auto row : affine) {
                            std::reverse(row.begin(), row.end());
                            specialized.push_back(evaluate(row, Z(t)));
                        }
                        check(verify_coefficients(e, specialized_matrix, specialized),
                              "affine coefficient specialization mismatch");
                    }
                }
            }
        }
        // Diagonal sign count: long and zero entries, and rational entries cleared to integers.
        ZMatrix wide_diagonal(4);
        wide_diagonal(0, 0) = pow10(900);
        wide_diagonal(1, 1) = -pow10(700) - 1;
        wide_diagonal(3, 3) = 1;
        const auto wide_inertia = inertia(wide_diagonal);
        check(wide_inertia.positive == 2 && wide_inertia.negative == 1 && wide_inertia.zero == 1,
              "long diagonal inertia");
        QMatrix rational_diagonal(3);
        rational_diagonal(0, 0) = Q(-1, 3);
        rational_diagonal(1, 1) = Q(5, 7);
        rational_diagonal(2, 2) = Q(2, 9);
        const auto rational_diagonal_inertia = inertia(rational_diagonal);
        check(rational_diagonal_inertia.positive == 2 && rational_diagonal_inertia.negative == 1 &&
              rational_diagonal_inertia.zero == 0, "rational diagonal inertia");
        // One off-diagonal pair leaves the sign count: [[1,2],[2,1]] has eigenvalues 3 and -1.
        ZMatrix coupled(2);
        coupled(0, 0) = coupled(1, 1) = 1;
        coupled(0, 1) = coupled(1, 0) = 2;
        const auto coupled_inertia = inertia(coupled);
        check(coupled_inertia.positive == 1 && coupled_inertia.negative == 1 && coupled_inertia.zero == 0,
              "coupled matrix inertia");
        // A valid file must survive roundtrip; independent fields must reject tampering.
        QMatrix k(2); k(0, 0) = 2; k(1, 1) = 3; k(0, 1) = k(1, 0) = 1;
        Options options; options.digits = 40;
        auto result = solve(k, identity<Q>(2), options);
        check(result.certified, "irrational example did not certify");
        std::stringstream serialized;
        write_certificates(serialized, result.roots, options.digits);
        unsigned digits;
        auto roots = read_certificates(serialized, digits);
        check(verify_root_certificates(k, identity<Q>(2), roots, digits), "roundtrip failed");
        roots[0].radius = 0;
        check(!verify_root_certificates(k, identity<Q>(2), roots, digits), "zero-radius mutation accepted");
        roots = result.roots; roots[0].index = 1;
        check(!verify_root_certificates(k, identity<Q>(2), roots, digits), "wrong-index mutation accepted");
        roots = result.roots; roots[0].lower += 1; roots[0].upper += 1;
        check(!verify_root_certificates(k, identity<Q>(2), roots, digits), "interval-shift mutation accepted");
        QMatrix scalar(1); scalar(0, 0) = Q(7, 3);
        auto one = solve(scalar, identity<Q>(1), options);
        check(one.certified && verify_root_certificates(scalar, identity<Q>(1), one.roots, 40),
              "order-one case failed");
        std::cout << "PASS " << checks << " extended exact-algebra assertions\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL after " << checks << " assertions: " << error.what() << '\n';
        return 1;
    }
}
