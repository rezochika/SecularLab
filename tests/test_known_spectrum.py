#!/usr/bin/env python3
"""Independent exact-oracle CLI regressions for rational generalized pencils.

The fixture identity is det(z*M-K) = det(C)^2 * product(z-lambda), where
K=C^T*diag(lambda)*C and M=C^T*C.  Only standard-library Fraction arithmetic
is used here; no solver code, saved coefficient file, or numeric eigensolver
supplies an expected value.  Normal runs leave no artifacts.  --report and
--artifacts-dir explicitly retain timing data and generated runs respectively.
"""

import argparse
from contextlib import ExitStack
from dataclasses import dataclass, replace
from fractions import Fraction
import json
from math import gcd, prod
from pathlib import Path
import re
import subprocess
import tempfile
import time


FAMILIES = ("separated", "cluster", "ties")
METHODS = {"auto": "polynomial", "inertia": "exact_inertia"}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rounded_decimal(value, digits):
    """Round a rational directly by quotient/remainder, half away from zero."""
    scaled = abs(value) * 10**digits
    whole, remainder = divmod(scaled.numerator, scaled.denominator)
    if 2 * remainder >= scaled.denominator:
        whole += 1
    sign = "-" if value < 0 and whole else ""
    if not digits:
        return sign + str(whole)
    integral, fractional = divmod(whole, 10**digits)
    return f"{sign}{integral}.{fractional:0{digits}d}"


def polynomial(roots):
    """Multiply linear factors, keeping every descending-power coefficient."""
    coefficients = [Fraction(1)]
    for root in roots:
        product = [Fraction(0)] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            product[index] += coefficient
            product[index + 1] -= root * coefficient
        coefficients = product
    return coefficients


def spectrum(family, order):
    if family == "separated":
        roots = [Fraction(4 * i + 1, 3) for i in range(order)]
        digits = 50
    elif family == "cluster":
        # The gap is 100 times smaller than one printed unit.  These must
        # remain distinct, disjoint certified roots despite identical text.
        roots = [Fraction(1), Fraction(1) + Fraction(1, 10**12)]
        roots += [Fraction(i) for i in range(2, order)]
        digits = 10
        require(roots[1] - roots[0] == Fraction(1, 10**12), "cluster gap")
        require(rounded_decimal(roots[0], digits) ==
                rounded_decimal(roots[1], digits), "cluster must share display")
    else:
        roots = [Fraction(i) + Fraction(1, 4) for i in range(order)]
        digits = 1
        require(all((root * 10).denominator == 2 for root in roots),
                "every tie must be exactly halfway between printed decimals")
    require(len(set(roots)) == order and all(root > 0 for root in roots),
            "spectrum must be positive and simple")
    return sorted(roots, reverse=True), digits


def congruence(order, variant):
    """Two fixed nonsingular rational upper-triangular changes of basis."""
    change = [[Fraction(int(i == j)) for j in range(order)] for i in range(order)]
    if variant == "a":
        change[0][0] = Fraction(2)
        for i in range(order - 1):
            change[i][i + 1] = Fraction((-1)**i, 2)
    else:
        change[0][0] = Fraction(3, 2)
        change[-1][-1] = Fraction(2)
        for i in range(order - 1):
            change[i][i + 1] = Fraction(1, 3)
        for i in range(order - 2):
            change[i][i + 2] = Fraction(-1, 4)
    # Triangularity makes this determinant independent of any determinant
    # implementation, and nonzero diagonal entries prove nonsingularity.
    determinant = prod(change[i][i] for i in range(order))
    require(determinant == (2 if variant == "a" else 3), "fixture determinant")
    require(all(change[i][j] == 0 for i in range(order) for j in range(i)),
            "fixture must be upper triangular")
    return change, determinant


def pencil(roots, change):
    order = len(roots)
    stiffness = [[sum((change[k][i] * roots[k] * change[k][j]
                       for k in range(order)), Fraction(0))
                  for j in range(order)] for i in range(order)]
    mass = [[sum((change[k][i] * change[k][j] for k in range(order)), Fraction(0))
             for j in range(order)] for i in range(order)]
    require(any(mass[i][j] for i in range(order) for j in range(i)),
            "generalized fixture mass must be non-diagonal")
    multiplier = 1
    for matrix in (stiffness, mass):
        for row in matrix:
            for value in row:
                denominator = value.denominator
                multiplier = multiplier // gcd(multiplier, denominator) * denominator
    require(multiplier > 1, "fixture must exercise denominator clearing")
    return stiffness, mass, multiplier


def write_matrix(path, matrix):
    order = len(matrix)
    entries = [(i + 1, j + 1, matrix[i][j])
               for i in range(order) for j in range(i + 1) if matrix[i][j]]
    with path.open("x", encoding="utf-8") as output:
        output.write("%%MatrixMarket matrix coordinate real symmetric\n")
        output.write("% Exact rational known-spectrum regression fixture\n")
        output.write(f"{order} {order} {len(entries)}\n")
        for row, column, value in entries:
            output.write(f"{row} {column} {value}\n")


def check_coefficients(text, expected, multiplier):
    scale_match = re.search(r"\bL=(\d+)(?:;|\s)", text)
    require(scale_match is not None and int(scale_match.group(1)) == multiplier,
            "coefficient header denominator-clearing multiplier")
    rows = [line.split("\t") for line in text.splitlines()
            if line.strip() and not line.startswith("#")]
    require(len(rows) == len(expected), "complete coefficient row count")
    require(all(len(row) == 2 for row in rows), "coefficient column count")
    require([int(row[0]) for row in rows] == list(range(len(expected) - 1, -1, -1)),
            "complete descending coefficient powers")
    actual = [Fraction(row[1]) for row in rows]
    require(all(value.denominator == 1 for value in actual),
            "cleared determinant coefficients must be integers")
    for power, actual_value, expected_value in zip(range(len(expected) - 1, -1, -1),
                                                  actual, expected):
        require(actual_value == expected_value,
                f"coefficient z^{power}: got {actual_value}, expected {expected_value}")


@dataclass(frozen=True)
class RootRow:
    index: int
    decimal: str
    lower: Fraction
    upper: Fraction
    radius: Fraction


def read_certificate(text, digits):
    lines = text.splitlines()
    require(f"# decimal_places={digits}" in lines, "certificate decimal places")
    data = [line for line in lines if line.strip() and not line.startswith("#")]
    require(data and data[0] == "index\tdecimal\tlower\tupper\tradius",
            "certificate column header")
    rows = []
    for line in data[1:]:
        fields = line.split("\t")
        require(len(fields) == 5, "certificate column count")
        rows.append(RootRow(int(fields[0]), fields[1],
                            *(Fraction(field) for field in fields[2:])))
    return rows


def check_roots(rows, roots, digits):
    require(len(rows) == len(roots), "complete certificate row count")
    for index, (row, root) in enumerate(zip(rows, roots)):
        require(row.index == index, f"root {index}: index/order")
        require(row.lower <= root <= row.upper, f"root {index}: exact known root containment")
        if index:
            require(row.upper < rows[index - 1].lower, f"root {index}: disjoint ordered bounds")
        require(row.decimal == rounded_decimal(root, digits),
                f"root {index}: canonical rounded decimal")
        require(rounded_decimal(row.lower, digits) == row.decimal and
                rounded_decimal(row.upper, digits) == row.decimal,
                f"root {index}: endpoint rounding")
        printed = Fraction(row.decimal)
        require(row.radius >= max(abs(printed - row.lower), abs(printed - row.upper),
                                  abs(printed - root)),
                f"root {index}: radius covers bounds and known root")


def expect_rejected(label, check):
    try:
        check()
    except AssertionError:
        return
    raise AssertionError(f"oracle mutation escaped detection: {label}")


def check_oracle_sensitivity():
    """Bad-output controls prove the independent assertions reject key regressions.

    No production code is mutated.  Literal small examples also protect the
    expected-value calculations from sign, orientation and rounding mistakes.
    """
    require(polynomial([Fraction(1, 4), Fraction(5, 4)]) ==
            [Fraction(1), Fraction(-3, 2), Fraction(5, 16)], "literal polynomial oracle")
    for value, digits, expected in ((Fraction(1, 4), 1, "0.3"),
                                    (Fraction(-1, 4), 1, "-0.3"),
                                    (Fraction(2, 3), 3, "0.667"),
                                    (Fraction(-1, 1000), 1, "0.0")):
        require(rounded_decimal(value, digits) == expected, "literal rounding oracle")

    def coefficient_text(values, multiplier=2):
        return f"# L={multiplier};\n" + "".join(
            f"{2 - i}\t{value}\n" for i, value in enumerate(values))

    # det(C)=3, L=2, roots 1 and 2: 3^2 * 2^2 * (z^2-3z+2).
    expected = [Fraction(36), Fraction(-108), Fraction(72)]
    require([value * 3**2 * 2**2 for value in polynomial([Fraction(1), Fraction(2)])]
            == expected, "literal scaled determinant oracle")
    check_coefficients(coefficient_text(expected), expected, 2)
    controls = {
        "missing det(C)^2": coefficient_text([4, -12, 8]),
        "missing L^n": coefficient_text([9, -27, 18]),
        "wrong reported L": coefficient_text(expected, 3),
        "missing coefficient": coefficient_text(expected[:-1]),
        "wrong power": coefficient_text(expected).replace("1\t", "0\t", 1),
    }
    for label, text in controls.items():
        expect_rejected(label, lambda text=text: check_coefficients(text, expected, 2))

    roots = [Fraction(51, 50), Fraction(101, 100)]
    rows = [RootRow(0, "1.0", roots[0], roots[0], Fraction(1, 50)),
            RootRow(1, "1.0", roots[1], roots[1], Fraction(1, 100))]
    check_roots(rows, roots, 1)
    mutations = {
        "wrong exact root": [replace(rows[0], lower=Fraction(103, 100),
                                      upper=Fraction(103, 100)), rows[1]],
        "wrong index": [replace(rows[0], index=1), rows[1]],
        "wrong rounding": [replace(rows[0], decimal="1.1"), rows[1]],
        "noncanonical decimal": [replace(rows[0], decimal="1.00"), rows[1]],
        "small radius": [replace(rows[0], radius=Fraction(0)), rows[1]],
        "overlapping bounds": [rows[0], replace(rows[1], upper=roots[0], radius=Fraction(1, 50))],
        "missing root": rows[:-1],
    }
    for label, changed in mutations.items():
        expect_rejected(label, lambda changed=changed: check_roots(changed, roots, 1))
    return list(controls) + list(mutations)


class Runner:
    def __init__(self, executable, timeout):
        self.executable = executable
        self.timeout = timeout
        self.commands = []

    def run(self, label, arguments):
        command = [str(self.executable), *map(str, arguments)]
        started = time.perf_counter()
        entry = {"label": label, "command": command}
        self.commands.append(entry)
        try:
            result = subprocess.run(command, text=True, capture_output=True, timeout=self.timeout)
        except subprocess.TimeoutExpired as error:
            entry["wall_seconds"] = time.perf_counter() - started
            entry["timed_out"] = True
            raise AssertionError(f"{label}: timeout after {self.timeout}s\n"
                                 f"stdout:\n{error.stdout}\nstderr:\n{error.stderr}") from error
        entry.update(wall_seconds=time.perf_counter() - started, returncode=result.returncode)
        require(result.returncode == 0, f"{label}: exit {result.returncode}\n"
                f"command: {command}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
        return result


def run_case(runner, directory, family, order, variant):
    label = f"{family}-n{order}-c{variant}"
    started = time.perf_counter()
    directory.mkdir()
    roots, digits = spectrum(family, order)
    change, determinant = congruence(order, variant)
    stiffness, mass, multiplier = pencil(roots, change)
    expected = [coefficient * determinant**2 * multiplier**order
                for coefficient in polynomial(roots)]
    require(all(coefficient.denominator == 1 for coefficient in expected),
            f"{label}: expected cleared polynomial must be integral")
    write_matrix(directory / "K.mtx", stiffness)
    write_matrix(directory / "M.mtx", mass)
    inputs = ["--matrix", directory / "K.mtx", "--mass", directory / "M.mtx"]
    outputs = {}
    command_start = len(runner.commands)
    for method, proof in METHODS.items():
        for workers in (1, 4):
            out = directory / f"solve-{method}-w{workers}"
            command_label = f"{label}/solve/{method}/workers={workers}"
            execution = runner.run(command_label, ["solve", *inputs, "--digits", digits,
                "--workers", workers, "--certification", method, "--out", out])
            try:
                result = json.loads((out / "result.json").read_text(encoding="utf-8"))
                require(result["certified"] and result["coefficient_identity_verified"],
                        "solve must certify and independently check coefficient identity")
                require(result["certification_method"] == proof, "requested proof mode")
                require(result["coefficient_oracle"] == "modular" and result["oracle_moduli"] >= 1,
                        "modular coefficient oracle")
                # Each exact phase clamps the requested workers to its own jobs.
                require(result["oracle_workers"] == min(workers, result["oracle_moduli"]) and
                        result["recurrence_workers"] == min(workers, order),
                        "requested worker configuration")
                require(result["certified_eigenvalues"] == order, "complete result root count")
                require(int(result["denominator_clearing_multiplier"]) == multiplier,
                        "result denominator-clearing multiplier")
                check_coefficients((out / "coefficients.tsv").read_text(encoding="utf-8"),
                                   expected, multiplier)
                rows = read_certificate((out / "eigenvalues.tsv").read_text(encoding="utf-8"), digits)
                check_roots(rows, roots, digits)
                if family == "cluster":
                    require(rows[-2].decimal == rows[-1].decimal,
                            "two distinct clustered roots must share printed decimals")
                runner.commands[-1]["solver_times_seconds"] = result["times_seconds"]
            except (AssertionError, KeyError, ValueError, OSError) as error:
                raise AssertionError(f"{command_label}: {error}\nstdout:\n{execution.stdout}"
                                     f"\nstderr:\n{execution.stderr}") from error
            outputs[method, workers] = out
        for artifact in ("coefficients.tsv", "eigenvalues.tsv"):
            require((outputs[method, 1] / artifact).read_bytes() ==
                    (outputs[method, 4] / artifact).read_bytes(),
                    f"{label}/{method}: serial/parallel {artifact} differs")

    # Each proof mode replays the OTHER proof mode's certificate, with the
    # opposite worker count.  Thus all four produced certificates are replayed.
    for method, proof in METHODS.items():
        other = "inertia" if method == "auto" else "auto"
        for workers in (1, 4):
            source = outputs[other, 4 if workers == 1 else 1] / "eigenvalues.tsv"
            out = directory / f"verify-{method}-w{workers}"
            command_label = f"{label}/verify/{method}/workers={workers}"
            execution = runner.run(command_label, ["verify", *inputs, "--cert", source,
                "--workers", workers, "--certification", method, "--out", out])
            try:
                result = json.loads((out / "result.json").read_text(encoding="utf-8"))
                require(execution.stdout.split()[0] == "VERIFIED", "successful replay status")
                require(result["status"] == "VERIFIED" and result["certified"], "certified replay")
                require(result["certification_method"] == proof, "requested replay proof mode")
                # Inertia replay has no coefficient-oracle jobs to parallelize.
                require(result["oracle_workers"] ==
                        (min(workers, result["oracle_moduli"]) if method == "auto" else 0),
                        "replay coefficient-oracle worker count")
                require(result["coefficient_oracle"] == ("modular" if method == "auto" else "none"),
                        "replay coefficient oracle follows proof mode")
                require(result["coefficient_identity_verified"] == (method == "auto"),
                        "replay coefficient identity follows proof mode")
                require(result["certified_eigenvalues"] == order, "complete replay root count")
                require(int(result["denominator_clearing_multiplier"]) == multiplier,
                        "replay denominator-clearing multiplier")
                runner.commands[-1]["solver_times_seconds"] = result["times_seconds"]
            except (AssertionError, IndexError, KeyError, ValueError, OSError) as error:
                raise AssertionError(f"{command_label}: {error}\nstdout:\n{execution.stdout}"
                                     f"\nstderr:\n{execution.stderr}") from error
    elapsed = time.perf_counter() - started
    print(f"PASS {label}: {order} exact roots, 4 solves, 4 replays, {elapsed:.3f}s", flush=True)
    return {"case": label, "family": family, "order": order, "congruence": variant,
            "decimal_places": digits, "det_C": str(determinant), "L": str(multiplier),
            "known_roots_descending": list(map(str, roots)), "wall_seconds": elapsed,
            "command_labels": [item["label"] for item in runner.commands[command_start:]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--orders", type=int, nargs="+", choices=(8, 16, 32), default=[8])
    parser.add_argument("--families", nargs="+", choices=FAMILIES, default=list(FAMILIES))
    parser.add_argument("--report", type=Path, help="new JSON timing report; never overwritten")
    parser.add_argument("--artifacts-dir", type=Path, help="new directory to retain all generated runs")
    parser.add_argument("--timeout", type=float, default=180, help="per-command timeout in seconds")
    options = parser.parse_args()
    require(options.executable.is_file(), f"missing executable: {options.executable}")
    require(options.timeout > 0, "timeout must be positive")
    for path in (options.report, options.artifacts_dir):
        require(path is None or not path.exists(), f"refusing to overwrite {path}")
    require(len(set(options.orders)) == len(options.orders), "duplicate orders")
    require(len(set(options.families)) == len(options.families), "duplicate families")
    runner = Runner(options.executable.resolve(), options.timeout)
    started = time.perf_counter()
    report = {"status": "RUNNING", "orders": options.orders, "families": options.families,
              "cases": [], "commands": runner.commands}
    try:
        report["rejected_oracle_mutations"] = check_oracle_sensitivity()
        with ExitStack() as stack:
            if options.artifacts_dir:
                directory = options.artifacts_dir.resolve()
                directory.mkdir(parents=True)
                report["artifacts_directory"] = str(directory)
            else:
                directory = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix="known-spectrum-")))
            for order in options.orders:
                for family in options.families:
                    for variant in ("a", "b"):
                        case_directory = directory / f"{family}-n{order}-c{variant}"
                        report["cases"].append(run_case(runner, case_directory, family, order, variant))
        report["status"] = "PASS"
    except Exception as error:
        report["status"] = "FAIL"
        report["error"] = str(error)
        raise
    finally:
        report["total_wall_seconds"] = time.perf_counter() - started
        if options.report:
            options.report.parent.mkdir(parents=True, exist_ok=True)
            with options.report.open("x", encoding="utf-8") as output:
                json.dump(report, output, indent=2)
                output.write("\n")
    print(f"PASS {len(report['cases'])} known-spectrum pencils, {len(runner.commands)} commands, "
          f"{len(report['rejected_oracle_mutations'])} rejected oracle mutations, "
          f"{report['total_wall_seconds']:.3f}s")


if __name__ == "__main__":
    main()
