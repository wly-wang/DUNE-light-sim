#!/usr/bin/env python3
import argparse
import math


def read_unique_vals(path: str, column: int, tol: float = 0.5):
    """
    Reads pd_coords.txt lines: idx x y z
    column: 1=x, 2=y, 3=z
    Clusters values within tol (cm).
    """
    vals = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                vals.append(float(parts[column]))
            except ValueError:
                continue

    vals.sort()
    reps = []
    for v in vals:
        if not reps or abs(v - reps[-1]) > tol:
            reps.append(v)
    return reps


def parse_list(s: str):
    return [float(x) for x in s.replace(",", " ").split() if x.strip()]


def fill_gaps_min_points(y_sorted, dy_max: float):
    """
    Given sorted Y list, insert the minimum number of points so that
    consecutive points differ by <= dy_max.
    """
    if dy_max <= 0:
        raise ValueError("dy_max must be > 0")

    y_sorted = sorted(y_sorted)
    filled = [y_sorted[0]]

    for a, b in zip(y_sorted[:-1], y_sorted[1:]):
        gap = b - a
        if gap <= dy_max + 1e-12:
            filled.append(b)
            continue

        n_segments = int(math.ceil(gap / dy_max))  # at least 2
        step = gap / n_segments

        # add intermediate points (minimum count = n_segments - 1)
        for k in range(1, n_segments):
            filled.append(a + k * step)
        filled.append(b)

    # de-duplicate tiny numerical artifacts
    out = []
    for v in filled:
        if not out or abs(v - out[-1]) > 1e-6:
            out.append(v)
    return out


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--pd", default="pd_coords.txt", help="Input pd_coords.txt (idx x y z).")
    ap.add_argument("--tol", type=float, default=0.5, help="Clustering tolerance for PD rows (cm).")

    ap.add_argument("--z-step", type=int, default=3, help="Use every Nth PD Z slice.")
    ap.add_argument("--y-step", type=int, default=1, help="Use every Nth PD Y row BEFORE filling gaps.")

    ap.add_argument("--dy-max", type=float, default=20.0,
                    help="Max allowed spacing between adjacent Y points after filling (cm). "
                         "Smaller => more Y points; larger => fewer.")

    ap.add_argument("--y-max", type=float, default=600.0, help="Top detector boundary (cm). Keep only 0<y<=y_max.")

    ap.add_argument("--d-centers",
                    default="25,75,125,175,225,275,325,375,425,475,525,575,625,675,725,775,825,875,925,975")
    ap.add_argument("--cos-centers",
                    default="0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95")

    ap.add_argument("--x-min", type=float, default=0.0)
    ap.add_argument("--x-max", type=float, default=350.0)

    ap.add_argument("--out", default="origin_grid_lists.txt", help="Output file with X,Y,Z lists.")
    args = ap.parse_args()

    d_centers = parse_list(args.d_centers)
    cos_centers = parse_list(args.cos_centers)

    # ---- Read PD Y and Z from file ----
    y_all = read_unique_vals(args.pd, column=2, tol=args.tol)
    z_all = read_unique_vals(args.pd, column=3, tol=args.tol)

    # pre-filter to top quarter (y>0) and boundary (<=y_max)
    y_all = [y for y in y_all if 0.0 < y <= args.y_max]

    # optional subsample BEFORE filling
    y_used = y_all[::max(1, args.y_step)]
    z_used = z_all[::max(1, args.z_step)]

    # ---- Fill gaps in Y with minimum points ----
    if len(y_used) < 2:
        raise RuntimeError("Not enough Y rows after filtering to fill gaps.")
    y_filled = fill_gaps_min_points(sorted(y_used), dy_max=args.dy_max)

    # ---- Build X ----
    X = set()
    for d in d_centers:
        if d <= 0:
            continue
        for c in cos_centers:
            if not (0.0 <= c <= 1.0):
                continue
            x = d * c
            if args.x_min <= x <= args.x_max:
                X.add(round(x, 3))
    X = sorted(X)

    Y = [round(y, 3) for y in y_filled]
    Z = [round(z, 3) for z in sorted(z_used)]

    # ---- Save lists ----
    with open(args.out, "w") as f:
        f.write("X = " + str(X) + "\n\n")
        f.write("Y = " + str(Y) + "\n\n")
        f.write("Z = " + str(Z) + "\n\n")

    print("================================")
    print(f"PD Y rows (unique, filtered): {len(y_all)}")
    print(f"PD Y rows used pre-fill (step={args.y_step}): {len(y_used)}")
    print(f"PD Y rows after fill (dy_max={args.dy_max} cm): {len(Y)}")

    print(f"PD Z slices total: {len(z_all)}")
    print(f"PD Z slices used (step={args.z_step}): {len(Z)}")

    print(f"|X| = {len(X)}")
    print(f"|Y| = {len(Y)}")
    print(f"|Z| = {len(Z)}")
    print(f"TOTAL ORIGINS = {len(X) * len(Y) * len(Z)}")
    print(f"Wrote lists to: {args.out}")
    print("================================")


if __name__ == "__main__":
    main()