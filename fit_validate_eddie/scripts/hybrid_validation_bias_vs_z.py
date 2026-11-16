#!/usr/bin/env python3

"""Plot source-level hybrid validation bias and RMS versus source Z."""

import argparse
import csv
import math
from array import array
from pathlib import Path

import ROOT


LOW_Z_JOIN_CM = 300.0
HIGH_Z_JOIN_CM = 5500.0
X_REGION_EDGES_CM = (0.0, 350.0 / 3.0, 700.0 / 3.0, 350.0)
Y_REGION_EDGES_CM = (0.0, 200.0, 400.0, 600.0)
STITCH_WINDOWS = (
    ("lowZ", 0.0, 600.0, LOW_Z_JOIN_CM),
    ("highZ", 5200.0, 5800.0, HIGH_Z_JOIN_CM),
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Plot truth-weighted hybrid prediction bias and RMS versus source Z "
            "for the full X/Y sample and each 3x3 X/|Y| region."
        )
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("plots_hybrid/source_population_debug.csv"),
        help="Validation CSV written by FastOpticalFitValidate",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("plots_hybrid/stitch_validation_bias_vs_z"),
        help="Directory for PDF, PNG, and summary CSV outputs",
    )
    parser.add_argument("--y-min", type=float, default=-0.5)
    parser.add_argument("--y-max", type=float, default=0.5)
    return parser.parse_args()


def load_validation_rows(path):
    required = {"X", "Y", "Z", "total_truth", "total_pred"}
    rows = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path} is missing required columns: {', '.join(sorted(missing))}"
            )
        for raw in reader:
            row = {name: float(raw[name]) for name in required}
            if row["total_truth"] <= 0.0:
                continue
            row["residual"] = (
                row["total_pred"] - row["total_truth"]
            ) / row["total_truth"]
            rows.append(row)
    if not rows:
        raise RuntimeError(f"No valid source rows found in {path}")
    return rows


def prediction_region(z):
    if z <= LOW_Z_JOIN_CM or z >= HIGH_Z_JOIN_CM:
        return "radial"
    return "vertical"


def in_region(value, low, high, include_high):
    return low <= value <= high if include_high else low <= value < high


def select_rows(rows, z_low, z_high, x_bounds=None, y_bounds=None):
    selected = []
    for row in rows:
        if not z_low <= row["Z"] <= z_high:
            continue
        if x_bounds is not None:
            x_low, x_high, include_x_high = x_bounds
            if not in_region(row["X"], x_low, x_high, include_x_high):
                continue
        if y_bounds is not None:
            y_low, y_high, include_y_high = y_bounds
            if not in_region(abs(row["Y"]), y_low, y_high, include_y_high):
                continue
        selected.append(row)
    return selected


def summarize_by_z(rows):
    by_z = {}
    for row in rows:
        by_z.setdefault(row["Z"], []).append(row)

    points = []
    for z in sorted(by_z):
        z_rows = by_z[z]
        total_weight = sum(row["total_truth"] for row in z_rows)
        bias = sum(
            row["total_truth"] * row["residual"] for row in z_rows
        ) / total_weight
        variance = sum(
            row["total_truth"] * (row["residual"] - bias) ** 2
            for row in z_rows
        ) / total_weight
        points.append(
            {
                "z": z,
                "n_sources": len(z_rows),
                "sum_truth": total_weight,
                "bias": bias,
                "rms": math.sqrt(max(variance, 0.0)),
                "prediction_region": prediction_region(z),
            }
        )
    return points


def z_errors(points, z_low, z_high, stitch_z):
    z_values = [point["z"] for point in points]
    boundaries = [z_low]
    for left, right in zip(z_values[:-1], z_values[1:]):
        boundary = 0.5 * (left + right)
        if left <= stitch_z < right:
            boundary = stitch_z
        boundaries.append(boundary)
    boundaries.append(z_high)
    low_errors = [z - boundaries[i] for i, z in enumerate(z_values)]
    high_errors = [boundaries[i + 1] - z for i, z in enumerate(z_values)]
    return low_errors, high_errors


def format_edge(value):
    return str(int(round(value)))


def make_plot(points, output_stem, title, z_low, z_high, stitch_z, y_min, y_max):
    if not points:
        print(f"Skipping {output_stem.name}: no validation sources in this selection")
        return

    x_values = [point["z"] for point in points]
    bias_values = [point["bias"] for point in points]
    rms_values = [point["rms"] for point in points]
    x_low_errors, x_high_errors = z_errors(points, z_low, z_high, stitch_z)
    zeros = [0.0] * len(points)

    x = array("d", x_values)
    ex_low = array("d", x_low_errors)
    ex_high = array("d", x_high_errors)
    zero = array("d", zeros)

    graph_rms = ROOT.TGraphAsymmErrors(
        len(points), x, array("d", rms_values), ex_low, ex_high, zero, zero
    )
    graph_bias = ROOT.TGraphAsymmErrors(
        len(points), x, array("d", bias_values), ex_low, ex_high, zero, zero
    )

    canvas = ROOT.TCanvas(f"c_{output_stem.name}", "", 960, 760)
    canvas.SetGrid()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.13)

    graph_rms.SetTitle(title)
    graph_rms.SetMarkerStyle(ROOT.kFullSquare)
    graph_rms.SetMarkerSize(1.4)
    graph_rms.SetMarkerColor(ROOT.kBlack)
    graph_rms.SetLineColor(ROOT.kBlack)
    graph_rms.GetXaxis().SetTitle("source Z [cm]")
    graph_rms.GetYaxis().SetTitle("(N_{pred} - N_{Geant4}) / N_{Geant4}")
    graph_rms.GetXaxis().SetLimits(z_low, z_high)
    graph_rms.GetYaxis().SetRangeUser(y_min, y_max)
    graph_rms.GetXaxis().SetTitleOffset(1.1)
    graph_rms.GetYaxis().SetTitleOffset(1.25)
    graph_rms.Draw("AP")

    graph_bias.SetMarkerStyle(ROOT.kFullCircle)
    graph_bias.SetMarkerSize(1.4)
    graph_bias.SetMarkerColor(ROOT.kRed)
    graph_bias.SetLineColor(ROOT.kRed)
    graph_bias.Draw("P SAME")

    zero_line = ROOT.TLine(z_low, 0.0, z_high, 0.0)
    zero_line.SetLineColor(ROOT.kGray + 2)
    zero_line.SetLineStyle(2)
    zero_line.Draw("SAME")

    stitch_line = ROOT.TLine(stitch_z, y_min, stitch_z, y_max)
    stitch_line.SetLineColor(ROOT.kMagenta + 1)
    stitch_line.SetLineStyle(2)
    stitch_line.SetLineWidth(2)
    stitch_line.Draw("SAME")

    legend = ROOT.TLegend(0.63, 0.72, 0.92, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(graph_rms, "RMS", "P")
    legend.AddEntry(graph_bias, "Bias", "P")
    legend.AddEntry(stitch_line, "stitch boundary", "L")
    legend.Draw()

    canvas.RedrawAxis()
    canvas.SaveAs(f"{output_stem}.pdf")
    canvas.SaveAs(f"{output_stem}.png")


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFontSize(0.04)

    rows = load_validation_rows(args.input_csv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    selections = [("allXY", "all X and |Y|", None, None)]
    for x_index in range(3):
        x_low = X_REGION_EDGES_CM[x_index]
        x_high = X_REGION_EDGES_CM[x_index + 1]
        x_bounds = (x_low, x_high, x_index == 2)
        for y_index in range(3):
            y_low = Y_REGION_EDGES_CM[y_index]
            y_high = Y_REGION_EDGES_CM[y_index + 1]
            y_bounds = (y_low, y_high, y_index == 2)
            key = (
                f"x{format_edge(x_low)}_{format_edge(x_high)}_"
                f"absy{format_edge(y_low)}_{format_edge(y_high)}"
            )
            label = (
                f"X = {format_edge(x_low)}-{format_edge(x_high)} cm, "
                f"|Y| = {format_edge(y_low)}-{format_edge(y_high)} cm"
            )
            selections.append((key, label, x_bounds, y_bounds))

    summary_rows = []
    for side, z_low, z_high, stitch_z in STITCH_WINDOWS:
        for key, label, x_bounds, y_bounds in selections:
            selected = select_rows(rows, z_low, z_high, x_bounds, y_bounds)
            points = summarize_by_z(selected)
            title = f"{side[:-1].capitalize()}-Z hybrid validation: {label}"
            output_stem = args.output_dir / f"validation_bias_vs_z_{side}_{key}"
            make_plot(
                points, output_stem, title, z_low, z_high, stitch_z,
                args.y_min, args.y_max,
            )
            for point in points:
                summary_rows.append(
                    {
                        "plot": output_stem.name,
                        "selection": label,
                        "side": side,
                        **point,
                    }
                )

    summary_path = args.output_dir / "validation_bias_rms_vs_z.csv"
    with summary_path.open("w", newline="") as stream:
        fieldnames = (
            "plot", "selection", "side", "z", "prediction_region",
            "n_sources", "sum_truth", "bias", "rms",
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Read {len(rows)} validation sources from {args.input_csv}")
    print(f"Wrote {len(selections) * len(STITCH_WINDOWS)} plots to {args.output_dir}")
    print(f"Numerical summary: {summary_path}")


if __name__ == "__main__":
    main()
