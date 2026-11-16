#!/usr/bin/env python3

"""Generate a dense fast-simulation scan across the two hybrid Z joins."""

import argparse
import csv
import re
from array import array
from pathlib import Path

import numpy as np
import ROOT


N_ANGLE_BINS = 9
ANGLE_BIN_WIDTH_DEG = 10.0
ABSORPTION_LENGTH_CM = 2000.0
ARAPUCA_WIDTH_CM = 10.0
ARAPUCA_HEIGHT_CM = 47.75
LOW_Z_JOIN_CM = 300.0
HIGH_Z_JOIN_CM = 5500.0
X_REGION_EDGES_CM = np.asarray((0.0, 350.0 / 3.0, 700.0 / 3.0, 350.0))
Y_REGION_EDGES_CM = np.asarray((0.0, 200.0, 400.0, 600.0))


def parse_number_list(value):
    return np.asarray([float(item) for item in value.split(",") if item.strip()], dtype=float)


def read_parameter_array(text, name):
    match = re.search(
        rf"double\s+{re.escape(name)}\s*\[9\]\s*=\s*\{{([^}}]+)\}}\s*;",
        text,
    )
    if not match:
        raise RuntimeError(f"Could not find parameter array '{name}'")
    values = np.asarray(
        [float(value.strip()) for value in match.group(1).split(",")], dtype=float
    )
    if values.size != N_ANGLE_BINS:
        raise RuntimeError(f"Parameter array '{name}' has {values.size} entries, expected 9")
    return values


def load_parameter_set(text, suffix):
    return {
        "p1": read_parameter_array(text, f"p1_{suffix}"),
        "p2": read_parameter_array(text, f"p2_{suffix}"),
        "p3": read_parameter_array(text, f"p3_{suffix}"),
        "p4": read_parameter_array(text, f"p4_{suffix}"),
        "s1": read_parameter_array(text, f"slopes1_{suffix}"),
        "s2": read_parameter_array(text, f"slopes2_{suffix}"),
        "s3": read_parameter_array(text, f"slopes3_{suffix}"),
    }


def load_detector_map(path, x_target=0.0, x_window=1.0):
    detector_map = np.loadtxt(path, dtype=float)
    selected = detector_map[np.abs(detector_map[:, 1] - x_target) <= x_window]
    if selected.size == 0:
        raise RuntimeError(f"No optical detectors selected from {path}")
    return selected[:, 1], selected[:, 2], selected[:, 3]


def rectangle_omega(a, b, d):
    aa = a / (2.0 * d)
    bb = b / (2.0 * d)
    argument = (1.0 + aa * aa + bb * bb) / (
        (1.0 + aa * aa) * (1.0 + bb * bb)
    )
    return 4.0 * np.arccos(np.sqrt(np.clip(argument, 0.0, 1.0)))


def rectangle_solid_angle(dx, dy, dz):
    d1 = np.abs(dy)
    d2 = np.abs(dx)
    d3 = np.abs(dz)
    half_width = 0.5 * ARAPUCA_WIDTH_CM
    half_height = 0.5 * ARAPUCA_HEIGHT_CM
    result = np.zeros_like(d2)

    centered = np.isclose(d1, 0.0) & np.isclose(d3, 0.0)
    result[centered] = rectangle_omega(
        ARAPUCA_WIDTH_CM, ARAPUCA_HEIGHT_CM, d2[centered]
    )

    outside_y = d1 > half_width
    outside_z = d3 > half_height

    mask = (~centered) & outside_y & outside_z
    a = d1[mask] - half_width
    b = d3[mask] - half_height
    result[mask] = 0.25 * (
        rectangle_omega(2.0 * (a + ARAPUCA_WIDTH_CM),
                        2.0 * (b + ARAPUCA_HEIGHT_CM), d2[mask])
        - rectangle_omega(2.0 * a,
                          2.0 * (b + ARAPUCA_HEIGHT_CM), d2[mask])
        - rectangle_omega(2.0 * (a + ARAPUCA_WIDTH_CM),
                          2.0 * b, d2[mask])
        + rectangle_omega(2.0 * a, 2.0 * b, d2[mask])
    )

    mask = (~centered) & (~outside_y) & (~outside_z)
    a = -d1[mask] + half_width
    b = -d3[mask] + half_height
    result[mask] = 0.25 * (
        rectangle_omega(2.0 * (ARAPUCA_WIDTH_CM - a),
                        2.0 * (ARAPUCA_HEIGHT_CM - b), d2[mask])
        + rectangle_omega(2.0 * a,
                          2.0 * (ARAPUCA_HEIGHT_CM - b), d2[mask])
        + rectangle_omega(2.0 * (ARAPUCA_WIDTH_CM - a),
                          2.0 * b, d2[mask])
        + rectangle_omega(2.0 * a, 2.0 * b, d2[mask])
    )

    mask = (~centered) & outside_y & (~outside_z)
    a = d1[mask] - half_width
    b = -d3[mask] + half_height
    result[mask] = 0.25 * (
        rectangle_omega(2.0 * (a + ARAPUCA_WIDTH_CM),
                        2.0 * (ARAPUCA_HEIGHT_CM - b), d2[mask])
        - rectangle_omega(2.0 * a,
                          2.0 * (ARAPUCA_HEIGHT_CM - b), d2[mask])
        + rectangle_omega(2.0 * (a + ARAPUCA_WIDTH_CM),
                          2.0 * b, d2[mask])
        - rectangle_omega(2.0 * a, 2.0 * b, d2[mask])
    )

    mask = (~centered) & (~outside_y) & outside_z
    a = -d1[mask] + half_width
    b = d3[mask] - half_height
    result[mask] = 0.25 * (
        rectangle_omega(2.0 * (ARAPUCA_WIDTH_CM - a),
                        2.0 * (b + ARAPUCA_HEIGHT_CM), d2[mask])
        - rectangle_omega(2.0 * (ARAPUCA_WIDTH_CM - a),
                          2.0 * b, d2[mask])
        + rectangle_omega(2.0 * a,
                          2.0 * (b + ARAPUCA_HEIGHT_CM), d2[mask])
        - rectangle_omega(2.0 * a, 2.0 * b, d2[mask])
    )

    return result


def gaisser_hillas(distance, p1, p2, p3, p4):
    difference = p2 - p4
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        term = np.power((distance - p4) / difference, difference / p3)
        exponential = np.exp((p2 - distance) / p3)
        correction = p1 * term * exponential
    return np.where(np.isfinite(correction) & (correction >= 0.0), correction, 0.0)


def predict_for_parameter_set(
    source_x,
    source_y,
    source_z,
    detector_x,
    detector_y,
    detector_z,
    border_distance,
    parameters,
    generated_photons,
):
    dx = source_x[:, None] - detector_x[None, :]
    dy = source_y[:, None] - detector_y[None, :]
    dz = np.broadcast_to(source_z - detector_z[None, :], dx.shape)
    distance = np.sqrt(dx * dx + dy * dy + dz * dz)
    cosine = np.abs(dx) / distance
    theta = np.degrees(np.arccos(np.clip(cosine, 0.0, 1.0)))
    angle_bin = np.clip(
        np.floor(theta / ANGLE_BIN_WIDTH_DEG).astype(int), 0, N_ANGLE_BINS - 1
    )

    solid_angle = rectangle_solid_angle(dx, dy, dz)
    geometric_hits = (
        np.exp(-distance / ABSORPTION_LENGTH_CM)
        * solid_angle
        * generated_photons
        / (4.0 * np.pi)
    )

    p1 = parameters["p1"][angle_bin] + parameters["s1"][angle_bin] * border_distance[:, None]
    p2 = parameters["p2"][angle_bin] + parameters["s2"][angle_bin] * border_distance[:, None]
    p3 = parameters["p3"][angle_bin] + parameters["s3"][angle_bin] * border_distance[:, None]
    p4 = parameters["p4"][angle_bin]
    correction = gaisser_hillas(distance, p1, p2, p3, p4)

    with np.errstate(divide="ignore", invalid="ignore"):
        predicted_hits = correction * geometric_hits / cosine
    predicted_hits = np.where(np.isfinite(predicted_hits) & (predicted_hits >= 0.0), predicted_hits, 0.0)

    # The C++ validator rounds each optical-detector prediction before summing.
    return np.rint(predicted_hits).sum(axis=1)


def run_scan(
    transition_z,
    z_values,
    source_x,
    source_y,
    detector_x,
    detector_y,
    detector_z,
    radial_parameters,
    vertical_parameters,
    generated_photons,
):
    radial_mean = []
    radial_std = []
    vertical_mean = []
    vertical_std = []
    hybrid_mean = []
    hybrid_std = []
    radial_by_source = []
    vertical_by_source = []
    hybrid_by_source = []
    stitch_radial_by_source = None
    stitch_vertical_by_source = None

    vertical_distance = np.abs(source_y)
    for z in z_values:
        radial_distance = np.hypot(source_y, z - transition_z)
        radial = predict_for_parameter_set(
            source_x, source_y, z,
            detector_x, detector_y, detector_z,
            radial_distance, radial_parameters, generated_photons,
        )
        vertical = predict_for_parameter_set(
            source_x, source_y, z,
            detector_x, detector_y, detector_z,
            vertical_distance, vertical_parameters, generated_photons,
        )

        use_radial = z <= LOW_Z_JOIN_CM if transition_z == LOW_Z_JOIN_CM else z >= HIGH_Z_JOIN_CM
        hybrid = radial if use_radial else vertical

        if np.isclose(z, transition_z):
            stitch_radial_by_source = radial.copy()
            stitch_vertical_by_source = vertical.copy()

        radial_mean.append(radial.mean())
        radial_std.append(radial.std())
        vertical_mean.append(vertical.mean())
        vertical_std.append(vertical.std())
        hybrid_mean.append(hybrid.mean())
        hybrid_std.append(hybrid.std())
        radial_by_source.append(radial)
        vertical_by_source.append(vertical)
        hybrid_by_source.append(hybrid)

    return {
        "z": np.asarray(z_values),
        "radial_mean": np.asarray(radial_mean),
        "radial_std": np.asarray(radial_std),
        "vertical_mean": np.asarray(vertical_mean),
        "vertical_std": np.asarray(vertical_std),
        "hybrid_mean": np.asarray(hybrid_mean),
        "hybrid_std": np.asarray(hybrid_std),
        "radial_by_source": np.asarray(radial_by_source),
        "vertical_by_source": np.asarray(vertical_by_source),
        "hybrid_by_source": np.asarray(hybrid_by_source),
        "stitch_radial_by_source": stitch_radial_by_source,
        "stitch_vertical_by_source": stitch_vertical_by_source,
    }


def make_graph(name, x_values, y_values, color, style, marker=0):
    graph = ROOT.TGraph(
        len(x_values),
        array("d", np.asarray(x_values, dtype=float)),
        array("d", np.asarray(y_values, dtype=float)),
    )
    graph.SetName(name)
    graph.SetLineColor(color)
    graph.SetLineWidth(3)
    graph.SetLineStyle(style)
    if marker:
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(0.45)
    return graph


def draw_scan(results_by_transition, output_dir, generated_photons):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleFontSize(0.045)

    root_file = ROOT.TFile(
        str(output_dir / "hybrid_stitch_simple.root"), "RECREATE"
    )

    for transition_z, region_name in (
        (LOW_Z_JOIN_CM, "low_z"),
        (HIGH_Z_JOIN_CM, "high_z"),
    ):
        result = results_by_transition[transition_z]
        boundary_index = int(np.argmin(np.abs(result["z"] - transition_z)))
        relative_photons = result["radial_mean"] / result["vertical_mean"]

        if transition_z == LOW_Z_JOIN_CM:
            radial_mask = result["z"] <= transition_z
        else:
            radial_mask = result["z"] >= transition_z
        vertical_mask = ~radial_mask

        ratio_span = max(0.01, float(np.ptp(relative_photons)))
        ratio_min = min(1.0, float(relative_photons.min())) - 0.15 * ratio_span
        ratio_max = max(1.0, float(relative_photons.max())) + 0.15 * ratio_span

        ratio_canvas = ROOT.TCanvas(f"c_ratio_{region_name}", "", 900, 650)
        ratio_canvas.SetGrid()
        ratio_canvas.SetLeftMargin(0.14)
        ratio_canvas.SetRightMargin(0.04)
        ratio_canvas.SetBottomMargin(0.13)
        ratio_canvas.SetTopMargin(0.10)
        ratio_frame = ratio_canvas.DrawFrame(
            result["z"].min(), ratio_min, result["z"].max(), ratio_max
        )
        ratio_frame.SetTitle(
            f"Radial/vertical prediction across the {'low-Z' if transition_z == LOW_Z_JOIN_CM else 'high-Z'} stitch"
        )
        ratio_frame.GetXaxis().SetTitle("source Z [cm]")
        ratio_frame.GetYaxis().SetTitle(
            "Relative Number of Photons Predicted (Radial/Vertical)"
        )
        ratio_frame.GetYaxis().SetTitleOffset(1.45)

        ratio_graph = make_graph(
            f"g_radial_over_vertical_{int(transition_z)}", result["z"],
            relative_photons, ROOT.kRed + 1, 1, 20,
        )

        unity = ROOT.TLine(result["z"].min(), 1.0, result["z"].max(), 1.0)
        unity.SetLineColor(ROOT.kGray + 1)
        unity.SetLineStyle(2)
        unity.SetLineWidth(2)
        unity.Draw()
        ratio_graph.Draw("LP SAME")

        ratio_boundary = ROOT.TLine(
            transition_z, ratio_min, transition_z, ratio_max
        )
        ratio_boundary.SetLineColor(ROOT.kMagenta + 1)
        ratio_boundary.SetLineStyle(9)
        ratio_boundary.SetLineWidth(3)
        ratio_boundary.Draw()

        ratio_legend = ROOT.TLegend(0.16, 0.17, 0.47, 0.31)
        ratio_legend.SetBorderSize(0)
        ratio_legend.SetFillStyle(0)
        ratio_legend.AddEntry(ratio_graph, "radial / vertical", "lp")
        ratio_legend.AddEntry(ratio_boundary, "stitch boundary", "l")
        ratio_legend.Draw()

        jump_percent = 100.0 * (relative_photons[boundary_index] - 1.0)
        ratio_annotation = ROOT.TLatex()
        ratio_annotation.SetNDC()
        ratio_annotation.SetTextSize(0.036)
        ratio_annotation.DrawLatex(0.67, 0.22, f"at boundary: {jump_percent:+.2f}%")

        ratio_stem = f"hybrid_stitch_{region_name}_relative"
        ratio_canvas.SaveAs(str(output_dir / f"{ratio_stem}.pdf"))
        ratio_canvas.SaveAs(str(output_dir / f"{ratio_stem}.png"))

        photon_span = max(1.0, float(np.ptp(result["hybrid_mean"])))
        photon_min = max(0.0, float(result["hybrid_mean"].min()) - 0.15 * photon_span)
        photon_max = float(result["hybrid_mean"].max()) + 0.15 * photon_span

        photon_canvas = ROOT.TCanvas(f"c_photons_{region_name}", "", 900, 650)
        photon_canvas.SetGrid()
        photon_canvas.SetLeftMargin(0.14)
        photon_canvas.SetRightMargin(0.04)
        photon_canvas.SetBottomMargin(0.13)
        photon_canvas.SetTopMargin(0.10)
        photon_frame = photon_canvas.DrawFrame(
            result["z"].min(), photon_min, result["z"].max(), photon_max
        )
        photon_frame.SetTitle(
            f"Hybrid photon prediction across the {'low-Z' if transition_z == LOW_Z_JOIN_CM else 'high-Z'} stitch"
        )
        photon_frame.GetXaxis().SetTitle("source Z [cm]")
        photon_frame.GetYaxis().SetTitle("Number of Photons Predicted")
        photon_frame.GetYaxis().SetTitleOffset(1.45)
        photon_frame.GetYaxis().SetMaxDigits(4)

        radial = make_graph(
            f"g_hybrid_radial_photons_{int(transition_z)}", result["z"][radial_mask],
            result["hybrid_mean"][radial_mask], ROOT.kRed + 1, 1, 20,
        )
        vertical = make_graph(
            f"g_hybrid_vertical_photons_{int(transition_z)}", result["z"][vertical_mask],
            result["hybrid_mean"][vertical_mask], ROOT.kBlue + 1, 1, 20,
        )
        radial.Draw("LP SAME")
        vertical.Draw("LP SAME")

        photon_boundary = ROOT.TLine(
            transition_z, photon_min, transition_z, photon_max
        )
        photon_boundary.SetLineColor(ROOT.kMagenta + 1)
        photon_boundary.SetLineStyle(9)
        photon_boundary.SetLineWidth(3)
        photon_boundary.Draw()

        photon_legend = ROOT.TLegend(0.16, 0.70, 0.48, 0.88)
        photon_legend.SetBorderSize(0)
        photon_legend.SetFillStyle(0)
        photon_legend.AddEntry(radial, "radial end-cap", "lp")
        photon_legend.AddEntry(vertical, "vertical middle", "lp")
        photon_legend.AddEntry(photon_boundary, "stitch boundary", "l")
        photon_legend.Draw()

        photon_stem = f"hybrid_stitch_{region_name}_photons"
        photon_canvas.SaveAs(str(output_dir / f"{photon_stem}.pdf"))
        photon_canvas.SaveAs(str(output_dir / f"{photon_stem}.png"))

        both_canvas = ROOT.TCanvas(f"c_both_models_{region_name}", "", 900, 650)
        both_canvas.SetGrid()
        both_canvas.SetLeftMargin(0.14)
        both_canvas.SetRightMargin(0.04)
        both_canvas.SetBottomMargin(0.13)
        both_canvas.SetTopMargin(0.10)
        both_frame = both_canvas.DrawFrame(
            result["z"].min(), photon_min, result["z"].max(), photon_max
        )
        both_frame.SetTitle(
            f"Radial and vertical predictions across the {'low-Z' if transition_z == LOW_Z_JOIN_CM else 'high-Z'} stitch"
        )
        both_frame.GetXaxis().SetTitle("source Z [cm]")
        both_frame.GetYaxis().SetTitle("Number of Photons Predicted")
        both_frame.GetYaxis().SetTitleOffset(1.45)
        both_frame.GetYaxis().SetMaxDigits(4)

        radial_selected = make_graph(
            f"g_radial_selected_{int(transition_z)}", result["z"][radial_mask],
            result["radial_mean"][radial_mask], ROOT.kRed + 1, 1, 20,
        )
        vertical_selected = make_graph(
            f"g_vertical_selected_{int(transition_z)}", result["z"][vertical_mask],
            result["vertical_mean"][vertical_mask], ROOT.kBlue + 1, 1, 20,
        )
        radial_extrapolated = make_graph(
            f"g_radial_extrapolated_{int(transition_z)}", result["z"][vertical_mask],
            result["radial_mean"][vertical_mask], ROOT.kRed + 1, 2,
        )
        vertical_extrapolated = make_graph(
            f"g_vertical_extrapolated_{int(transition_z)}", result["z"][radial_mask],
            result["vertical_mean"][radial_mask], ROOT.kBlue + 1, 2,
        )
        radial_extrapolated.Draw("L SAME")
        vertical_extrapolated.Draw("L SAME")
        radial_selected.Draw("LP SAME")
        vertical_selected.Draw("LP SAME")

        both_boundary = ROOT.TLine(
            transition_z, photon_min, transition_z, photon_max
        )
        both_boundary.SetLineColor(ROOT.kMagenta + 1)
        both_boundary.SetLineStyle(9)
        both_boundary.SetLineWidth(3)
        both_boundary.Draw()

        if transition_z == LOW_Z_JOIN_CM:
            both_legend = ROOT.TLegend(0.16, 0.70, 0.50, 0.88)
        else:
            both_legend = ROOT.TLegend(0.16, 0.16, 0.50, 0.34)
        both_legend.SetBorderSize(0)
        both_legend.SetFillStyle(0)
        both_legend.AddEntry(radial_selected, "radial parameterization", "lp")
        both_legend.AddEntry(vertical_selected, "vertical parameterization", "lp")
        both_legend.AddEntry(both_boundary, "stitch boundary", "l")
        both_legend.Draw()

        both_stem = f"hybrid_stitch_{region_name}_both_models"
        both_canvas.SaveAs(str(output_dir / f"{both_stem}.pdf"))
        both_canvas.SaveAs(str(output_dir / f"{both_stem}.png"))

        root_file.cd()
        ratio_graph.Write()
        ratio_canvas.Write()
        radial.Write()
        vertical.Write()
        photon_canvas.Write()
        radial_selected.Write()
        vertical_selected.Write()
        radial_extrapolated.Write()
        vertical_extrapolated.Write()
        both_canvas.Write()

    root_file.Close()


def region_mask(values, edges, index):
    if index == len(edges) - 2:
        return (values >= edges[index]) & (values <= edges[index + 1])
    return (values >= edges[index]) & (values < edges[index + 1])


def draw_xy_region_scans(results_by_transition, output_dir, source_x, source_y):
    region_dir = output_dir / "xy_regions"
    region_dir.mkdir(parents=True, exist_ok=True)
    root_file = ROOT.TFile(str(region_dir / "xy_region_stitch_scans.root"), "RECREATE")
    csv_path = region_dir / "xy_region_stitch_scans.csv"

    with csv_path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow([
            "x_bin", "abs_y_bin", "x_min_cm", "x_max_cm",
            "abs_y_min_cm", "abs_y_max_cm", "source_positions",
            "transition_z_cm", "source_z_cm", "hybrid_mean_photons",
            "radial_mean_photons", "vertical_mean_photons",
        ])

        for x_index in range(3):
            x_mask = region_mask(source_x, X_REGION_EDGES_CM, x_index)
            for y_index in range(3):
                y_mask = region_mask(np.abs(source_y), Y_REGION_EDGES_CM, y_index)
                source_mask = x_mask & y_mask
                source_count = int(np.count_nonzero(source_mask))
                if source_count == 0:
                    continue

                x_low = X_REGION_EDGES_CM[x_index]
                x_high = X_REGION_EDGES_CM[x_index + 1]
                y_low = Y_REGION_EDGES_CM[y_index]
                y_high = Y_REGION_EDGES_CM[y_index + 1]

                for transition_z, transition_name in (
                    (LOW_Z_JOIN_CM, "low_z"),
                    (HIGH_Z_JOIN_CM, "high_z"),
                ):
                    result = results_by_transition[transition_z]
                    z_values = result["z"]
                    hybrid_mean = result["hybrid_by_source"][:, source_mask].mean(axis=1)
                    radial_mean = result["radial_by_source"][:, source_mask].mean(axis=1)
                    vertical_mean = result["vertical_by_source"][:, source_mask].mean(axis=1)

                    if transition_z == LOW_Z_JOIN_CM:
                        radial_z_mask = z_values <= transition_z
                    else:
                        radial_z_mask = z_values >= transition_z
                    vertical_z_mask = ~radial_z_mask

                    y_span = max(1.0, float(np.ptp(hybrid_mean)))
                    plot_min = max(0.0, float(hybrid_mean.min()) - 0.15 * y_span)
                    plot_max = float(hybrid_mean.max()) + 0.15 * y_span

                    tag = f"x{x_index + 1}_y{y_index + 1}_{transition_name}"
                    canvas = ROOT.TCanvas(f"c_xy_{tag}", "", 900, 650)
                    canvas.SetGrid()
                    canvas.SetLeftMargin(0.14)
                    canvas.SetRightMargin(0.04)
                    canvas.SetBottomMargin(0.13)
                    canvas.SetTopMargin(0.10)
                    frame = canvas.DrawFrame(
                        z_values.min(), plot_min, z_values.max(), plot_max
                    )
                    frame.SetTitle(
                        f"{'Low-Z' if transition_z == LOW_Z_JOIN_CM else 'High-Z'} stitch: "
                        f"X = {x_low:.0f}-{x_high:.0f} cm, |Y| = {y_low:.0f}-{y_high:.0f} cm"
                    )
                    frame.GetXaxis().SetTitle("source Z [cm]")
                    frame.GetYaxis().SetTitle("Number of Photons Predicted")
                    frame.GetYaxis().SetTitleOffset(1.45)
                    frame.GetYaxis().SetMaxDigits(4)

                    radial = make_graph(
                        f"g_xy_radial_{tag}", z_values[radial_z_mask],
                        hybrid_mean[radial_z_mask], ROOT.kRed + 1, 1, 20,
                    )
                    vertical = make_graph(
                        f"g_xy_vertical_{tag}", z_values[vertical_z_mask],
                        hybrid_mean[vertical_z_mask], ROOT.kBlue + 1, 1, 20,
                    )
                    radial.Draw("LP SAME")
                    vertical.Draw("LP SAME")

                    boundary = ROOT.TLine(
                        transition_z, plot_min, transition_z, plot_max
                    )
                    boundary.SetLineColor(ROOT.kMagenta + 1)
                    boundary.SetLineStyle(9)
                    boundary.SetLineWidth(3)
                    boundary.Draw()

                    if transition_z == LOW_Z_JOIN_CM:
                        legend = ROOT.TLegend(0.16, 0.70, 0.48, 0.88)
                    else:
                        legend = ROOT.TLegend(0.16, 0.16, 0.48, 0.34)
                    legend.SetBorderSize(0)
                    legend.SetFillStyle(0)
                    legend.AddEntry(radial, "radial end-cap", "lp")
                    legend.AddEntry(vertical, "vertical middle", "lp")
                    legend.AddEntry(boundary, "stitch boundary", "l")
                    legend.Draw()

                    count_label = ROOT.TLatex()
                    count_label.SetNDC()
                    count_label.SetTextSize(0.034)
                    count_label.DrawLatex(
                        0.72, 0.84, f"{source_count} X/Y positions"
                    )

                    output_stem = f"hybrid_stitch_{tag}_photons"
                    canvas.SaveAs(str(region_dir / f"{output_stem}.pdf"))
                    canvas.SaveAs(str(region_dir / f"{output_stem}.png"))

                    root_file.cd()
                    radial.Write()
                    vertical.Write()
                    canvas.Write()

                    for z_index, z_value in enumerate(z_values):
                        writer.writerow([
                            x_index + 1, y_index + 1, x_low, x_high,
                            y_low, y_high, source_count, transition_z, z_value,
                            hybrid_mean[z_index], radial_mean[z_index],
                            vertical_mean[z_index],
                        ])

    root_file.Close()
    return region_dir, csv_path


def write_outputs(results_by_transition, output_dir, step_cm, source_x, source_y):
    csv_path = output_dir / "hybrid_stitch_scan.csv"
    with csv_path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow([
            "transition_z_cm", "source_z_cm", "hybrid_mean_pe", "hybrid_std_pe",
            "radial_mean_pe", "radial_std_pe", "vertical_mean_pe", "vertical_std_pe",
            "radial_over_vertical_minus_one_percent",
        ])
        for transition_z, result in results_by_transition.items():
            for index, z in enumerate(result["z"]):
                writer.writerow([
                    transition_z, z,
                    result["hybrid_mean"][index], result["hybrid_std"][index],
                    result["radial_mean"][index], result["radial_std"][index],
                    result["vertical_mean"][index], result["vertical_std"][index],
                    100.0 * (result["radial_mean"][index] /
                             result["vertical_mean"][index] - 1.0),
                ])

    source_csv_path = output_dir / "hybrid_stitch_jump_by_source.csv"
    with source_csv_path.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow([
            "source_x_cm", "source_y_cm",
            "low_z_radial_pe", "low_z_vertical_pe", "low_z_difference_percent",
            "high_z_radial_pe", "high_z_vertical_pe", "high_z_difference_percent",
        ])
        low = results_by_transition[LOW_Z_JOIN_CM]
        high = results_by_transition[HIGH_Z_JOIN_CM]
        low_difference = 100.0 * (
            low["stitch_radial_by_source"] / low["stitch_vertical_by_source"] - 1.0
        )
        high_difference = 100.0 * (
            high["stitch_radial_by_source"] / high["stitch_vertical_by_source"] - 1.0
        )
        for index in range(source_x.size):
            writer.writerow([
                source_x[index], source_y[index],
                low["stitch_radial_by_source"][index],
                low["stitch_vertical_by_source"][index],
                low_difference[index],
                high["stitch_radial_by_source"][index],
                high["stitch_vertical_by_source"][index],
                high_difference[index],
            ])

    summary_path = output_dir / "hybrid_stitch_summary.txt"
    with summary_path.open("w") as output:
        for transition_z, result in results_by_transition.items():
            stitch_index = int(np.argmin(np.abs(result["z"] - transition_z)))
            radial = result["radial_mean"][stitch_index]
            vertical = result["vertical_mean"][stitch_index]
            parameter_jump = 100.0 * (radial / vertical - 1.0)
            source_differences = 100.0 * (
                result["stitch_radial_by_source"] /
                result["stitch_vertical_by_source"] - 1.0
            )

            if transition_z == LOW_Z_JOIN_CM:
                radial_side_index = stitch_index
                vertical_side_index = min(stitch_index + 1, len(result["z"]) - 1)
            else:
                vertical_side_index = max(stitch_index - 1, 0)
                radial_side_index = stitch_index

            radial_side = result["hybrid_mean"][radial_side_index]
            vertical_side = result["hybrid_mean"][vertical_side_index]
            sampled_jump = 100.0 * (radial_side / vertical_side - 1.0)

            output.write(f"Z stitch: {transition_z:.1f} cm\n")
            output.write(f"  radial mean at stitch: {radial:.6g} PE\n")
            output.write(f"  vertical mean at stitch: {vertical:.6g} PE\n")
            output.write(f"  radial/vertical difference at same Z: {parameter_jump:+.4f}%\n")
            output.write(
                "  per-source difference min/median/max: "
                f"{source_differences.min():+.4f}% / "
                f"{np.median(source_differences):+.4f}% / "
                f"{source_differences.max():+.4f}%\n"
            )
            output.write(
                f"  sampled hybrid radial-side/vertical-side difference "
                f"({step_cm:g} cm step): {sampled_jump:+.4f}%\n\n"
            )

    return csv_path, source_csv_path, summary_path


def main():
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--parameters",
        type=Path,
        default=base_dir / "GH_parameters_hybrid_2set.txt",
    )
    parser.add_argument(
        "--detector-map",
        type=Path,
        default=base_dir / "protodune_optical_mapping.txt",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=base_dir / "plots_hybrid" / "stitch_scan",
    )
    parser.add_argument("--generated-photons", type=float, default=1.0e7)
    parser.add_argument("--window-cm", type=float, default=300.0)
    parser.add_argument("--step-cm", type=float, default=2.0)
    parser.add_argument("--source-x", default="20,100,200,300,350")
    parser.add_argument("--source-y", default="31.135,171.243,311.352,451.461,591.569")
    parser.add_argument(
        "--full-source-grid",
        action="store_true",
        help="Use every unique X and Y in plots_hybrid/source_population_debug.csv",
    )
    parser.add_argument(
        "--source-grid-csv",
        type=Path,
        default=base_dir / "plots_hybrid" / "source_population_debug.csv",
    )
    parser.add_argument(
        "--xy-region-plots",
        action="store_true",
        help="Write low-Z and high-Z hybrid plots for a 3x3 X/|Y| partition",
    )
    args = parser.parse_args()

    parameter_text = args.parameters.read_text()
    radial_parameters = load_parameter_set(parameter_text, "radial_z_edges")
    vertical_parameters = load_parameter_set(parameter_text, "vertical_z300_5500")
    detector_x, detector_y, detector_z = load_detector_map(args.detector_map)

    if args.full_source_grid:
        with args.source_grid_csv.open() as source_file:
            source_rows = list(csv.DictReader(source_file))
        x_values = np.asarray(
            sorted({float(row["X"]) for row in source_rows}), dtype=float
        )
        y_values = np.asarray(
            sorted({float(row["Y"]) for row in source_rows}), dtype=float
        )
    else:
        x_values = parse_number_list(args.source_x)
        y_values = parse_number_list(args.source_y)
    source_x, source_y = np.meshgrid(x_values, y_values, indexing="ij")
    source_x = source_x.ravel()
    source_y = source_y.ravel()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    results_by_transition = {}
    for transition_z in (LOW_Z_JOIN_CM, HIGH_Z_JOIN_CM):
        z_values = np.arange(
            transition_z - args.window_cm,
            transition_z + args.window_cm + 0.5 * args.step_cm,
            args.step_cm,
        )
        results_by_transition[transition_z] = run_scan(
            transition_z, z_values, source_x, source_y,
            detector_x, detector_y, detector_z,
            radial_parameters, vertical_parameters, args.generated_photons,
        )

    csv_path, source_csv_path, summary_path = write_outputs(
        results_by_transition, args.output_dir, args.step_cm, source_x, source_y
    )
    draw_scan(results_by_transition, args.output_dir, args.generated_photons)
    region_dir = None
    region_csv_path = None
    if args.xy_region_plots:
        region_dir, region_csv_path = draw_xy_region_scans(
            results_by_transition, args.output_dir, source_x, source_y
        )

    print(f"Selected optical detectors: {detector_x.size}")
    print(f"Source X/Y combinations averaged at each Z: {source_x.size}")
    print(f"CSV: {csv_path}")
    print(f"Per-source CSV: {source_csv_path}")
    print(f"Summary: {summary_path}")
    if region_dir is not None:
        print(f"3x3 X/Y region plots: {region_dir}")
        print(f"3x3 X/Y region CSV: {region_csv_path}")
    print(summary_path.read_text())


if __name__ == "__main__":
    main()
