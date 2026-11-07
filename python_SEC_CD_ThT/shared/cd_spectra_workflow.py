from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
import pandas as pd


@dataclass
class SpectrumRecord:
    variant: str
    fraction: int
    run: int
    path: Path
    data: pd.DataFrame


class CDSpectraWorkflow:
    """
    Importable workflow for loading and comparing CD spectra in this project.

    The class is designed for notebook use:

        from cd_spectra_workflow import CDSpectraWorkflow
        workflow = CDSpectraWorkflow(".")
        workflow.summary()
        workflow.plot_fraction_grid([19, 30], variants=("WT", "N27D"))

    Scaling modes
    -------------
    raw:
        Plot the instrument output directly in mdeg.
    normalized:
        Divide each spectrum by a scale factor computed from a chosen wavelength
        window. The default ``trough`` normalization sets the most negative
        point in that window to ``-1`` so shape differences are easy to compare.
    mre:
        Convert to mean residue ellipticity once concentration, path length,
        and residue count are known.
    """

    WT_FOLDER = "WT_022326"
    DEAM_FOLDER = "Deam_021926"

    COLORS = {
        "WT": "#1f4e79",
        "N27D": "#c0392b",
        "Q15E": "#1e8449",
    }

    WT_PATTERN = re.compile(
        r"WT\s+f(?P<fraction>\d{2}).*?(?:run(?P<run>\d+))?\s+avg\.csv$",
        re.IGNORECASE,
    )
    DEAM_PATTERN = re.compile(
        r"(?P<variant>N27D|Q15E)\s+(?P<fraction>\d+)fraction\s+avg(?:\s+cleandata)?\.csv$",
        re.IGNORECASE,
    )

    def __init__(
        self,
        root_dir: str | Path = ".",
        wt_folder: str = WT_FOLDER,
        deam_folder: str = DEAM_FOLDER,
    ) -> None:
        self.root_dir = Path(root_dir).resolve()
        self.wt_dir = self.root_dir / wt_folder
        self.deam_dir = self.root_dir / deam_folder
        self.records = self._discover_records()

    def _read_cd_csv(self, path: Path) -> pd.DataFrame:
        df = pd.read_csv(path, header=None, names=["wavelength_nm", "cd_mdeg"])
        df = df.apply(pd.to_numeric, errors="coerce").dropna()
        return df.sort_values("wavelength_nm").reset_index(drop=True)

    def _discover_records(self) -> dict[str, dict[int, list[SpectrumRecord]]]:
        records: dict[str, dict[int, list[SpectrumRecord]]] = {}
        seen_keys: set[tuple[str, int, int]] = set()

        for path in sorted(self.wt_dir.glob("*.csv")):
            match = self.WT_PATTERN.match(path.name)
            if not match:
                continue

            fraction = int(match.group("fraction"))
            run = int(match.group("run") or 1)
            key = ("WT", fraction, run)
            if key in seen_keys:
                continue

            seen_keys.add(key)
            record = SpectrumRecord(
                variant="WT",
                fraction=fraction,
                run=run,
                path=path,
                data=self._read_cd_csv(path),
            )
            records.setdefault("WT", {}).setdefault(fraction, []).append(record)

        for path in sorted(self.deam_dir.glob("*.csv")):
            match = self.DEAM_PATTERN.match(path.name)
            if not match:
                continue

            variant = match.group("variant").upper()
            fraction = int(match.group("fraction"))
            key = (variant, fraction, 1)
            if key in seen_keys:
                continue

            seen_keys.add(key)
            record = SpectrumRecord(
                variant=variant,
                fraction=fraction,
                run=1,
                path=path,
                data=self._read_cd_csv(path),
            )
            records.setdefault(variant, {}).setdefault(fraction, []).append(record)

        for fractions in records.values():
            for fraction_records in fractions.values():
                fraction_records.sort(key=lambda record: record.run)

        return records

    def summary(self, variants: Sequence[str] | None = None) -> pd.DataFrame:
        variants = tuple(variant.upper() for variant in (variants or ("WT", "N27D", "Q15E")))
        rows: list[dict[str, object]] = []

        for variant in variants:
            for fraction, fraction_records in sorted(self.records.get(variant, {}).items()):
                rows.append(
                    {
                        "variant": variant,
                        "fraction": fraction,
                        "runs": len(fraction_records),
                        "files": ", ".join(record.path.name for record in fraction_records),
                    }
                )

        return pd.DataFrame(rows)

    def available_fractions(
        self,
        variants: Sequence[str] | None = None,
        require_all: bool = False,
    ) -> list[int]:
        variants = tuple(variant.upper() for variant in (variants or tuple(self.records)))
        fraction_sets = [set(self.records.get(variant, {})) for variant in variants]
        if not fraction_sets:
            return []

        fractions = set.intersection(*fraction_sets) if require_all else set.union(*fraction_sets)
        return sorted(fractions)

    def get_records(
        self,
        fraction: int,
        variants: Sequence[str] = ("WT", "N27D"),
    ) -> dict[str, list[SpectrumRecord]]:
        selected: dict[str, list[SpectrumRecord]] = {}

        for variant in variants:
            variant_name = variant.upper()
            fraction_records = self.records.get(variant_name, {}).get(fraction, [])
            if fraction_records:
                selected[variant_name] = fraction_records

        return selected

    def _average_records(self, records: Sequence[SpectrumRecord]) -> pd.DataFrame:
        merged: pd.DataFrame | None = None

        for record in records:
            renamed = record.data.rename(columns={"cd_mdeg": f"run_{record.run}"})
            merged = renamed if merged is None else merged.merge(renamed, on="wavelength_nm", how="inner")

        if merged is None:
            raise ValueError("No spectra were provided for averaging.")

        signal_columns = [column for column in merged.columns if column != "wavelength_nm"]
        averaged = merged.copy()
        averaged["cd_mdeg"] = averaged[signal_columns].mean(axis=1)
        return averaged[["wavelength_nm", "cd_mdeg"]]

    def _select_variant_data(
        self,
        records: Sequence[SpectrumRecord],
        average_replicates: bool = True,
    ) -> tuple[pd.DataFrame, str]:
        if not records:
            raise ValueError("No records were provided.")

        if average_replicates and len(records) > 1:
            label_suffix = f" mean ({len(records)} runs)"
            return self._average_records(records), label_suffix

        return records[0].data.copy(), ""

    def _filter_range(
        self,
        data: pd.DataFrame,
        wavelength_range: tuple[float, float] | None,
    ) -> pd.DataFrame:
        if wavelength_range is None:
            return data.copy()

        low, high = wavelength_range
        filtered = data[(data["wavelength_nm"] >= low) & (data["wavelength_nm"] <= high)].copy()
        if filtered.empty:
            raise ValueError(f"No data remain after filtering to range {wavelength_range}.")
        return filtered

    def _prepare_signal(
        self,
        data: pd.DataFrame,
        mode: str,
        wavelength_range: tuple[float, float] | None,
        normalize_range: tuple[float, float],
        normalization: str,
        concentration_uM: float | None,
        path_length_cm: float,
        num_residues: int | None,
    ) -> tuple[pd.DataFrame, str]:
        mode = mode.lower()

        if mode == "raw":
            filtered = self._filter_range(data, wavelength_range)
            filtered["signal"] = filtered["cd_mdeg"]
            return filtered, "CD signal (mdeg)"

        if mode == "normalized":
            filtered = self._filter_range(data, wavelength_range)
            norm_window = self._filter_range(data, normalize_range)

            if normalization == "trough":
                scale = abs(norm_window["cd_mdeg"].min())
                if np.isclose(scale, 0.0):
                    scale = float(norm_window["cd_mdeg"].abs().max())
            elif normalization == "maxabs":
                scale = float(norm_window["cd_mdeg"].abs().max())
            else:
                raise ValueError("normalization must be 'trough' or 'maxabs'.")

            if np.isclose(scale, 0.0):
                raise ValueError("Normalization scale factor is zero.")

            filtered["signal"] = filtered["cd_mdeg"] / scale
            return filtered, "Normalized CD"

        if mode == "mre":
            if concentration_uM is None or num_residues is None:
                raise ValueError(
                    "MRE plotting requires concentration_uM and num_residues."
                )

            filtered = self._filter_range(data, wavelength_range)
            concentration_molar = concentration_uM * 1e-6
            scale = 10.0 * concentration_molar * path_length_cm * num_residues
            if np.isclose(scale, 0.0):
                raise ValueError("MRE scale factor is zero.")

            filtered["signal"] = filtered["cd_mdeg"] / scale
            return filtered, "[theta] (deg cm^2 dmol^-1)"

        raise ValueError("mode must be 'raw', 'normalized', or 'mre'.")

    def _lighten_color(self, color: str, amount: float) -> tuple[float, float, float]:
        amount = float(np.clip(amount, 0.0, 1.0))
        rgb = np.array(to_rgb(color))
        return tuple(rgb + (1.0 - rgb) * amount)

    def plot_fraction(
        self,
        fraction: int,
        variants: Sequence[str] = ("WT", "N27D"),
        mode: str = "raw",
        average_replicates: bool = True,
        wavelength_range: tuple[float, float] | None = None,
        normalize_range: tuple[float, float] = (200, 250),
        normalization: str = "trough",
        concentration_uM: float | None = None,
        path_length_cm: float = 0.1,
        num_residues: int | None = None,
        ax: plt.Axes | None = None,
        title: str | None = None,
        fraction_label: str | None = None,
        linewidth: float = 2.5,
    ) -> plt.Axes:
        selected = self.get_records(fraction, variants=variants)
        if not selected:
            raise FileNotFoundError(f"No spectra were found for fraction {fraction}.")

        created_figure = False
        if ax is None:
            _, ax = plt.subplots(figsize=(6, 4.5))
            created_figure = True

        ylabel = "CD signal (mdeg)"
        for variant in variants:
            variant_name = variant.upper()
            records = selected.get(variant_name)
            if not records:
                continue

            data, label_suffix = self._select_variant_data(
                records,
                average_replicates=average_replicates,
            )
            prepared, ylabel = self._prepare_signal(
                data=data,
                mode=mode,
                wavelength_range=wavelength_range,
                normalize_range=normalize_range,
                normalization=normalization,
                concentration_uM=concentration_uM,
                path_length_cm=path_length_cm,
                num_residues=num_residues,
            )

            ax.plot(
                prepared["wavelength_nm"],
                prepared["signal"],
                color=self.COLORS.get(variant_name, None),
                linewidth=linewidth,
                label=f"{variant_name}{label_suffix}",
            )

        label_text = fraction_label.strip() if fraction_label else ""
        if title is None:
            title = f"Fraction {fraction}"
            if label_text:
                title = f"{title} - {label_text}"

        ax.set_title(title, fontsize=13, weight="bold")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25, linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(frameon=False)

        if wavelength_range is not None:
            ax.set_xlim(*wavelength_range)

        if created_figure:
            plt.tight_layout()

        return ax

    def plot_state_overlay(
        self,
        state_fractions: dict[str, int],
        variants: Sequence[str] = ("WT", "N27D"),
        mode: str = "raw",
        average_replicates: bool = True,
        wavelength_range: tuple[float, float] | None = (195, 250),
        normalize_range: tuple[float, float] = (195, 250),
        normalization: str = "trough",
        concentration_uM: float | None = None,
        path_length_cm: float = 0.1,
        num_residues: int | None = None,
        ax: plt.Axes | None = None,
        title: str | None = None,
        linewidth: float = 2.6,
        figsize: tuple[float, float] = (8.2, 5.6),
        variant_colors: dict[str, str] | None = None,
        trace_labels: dict[str, str] | None = None,
        monomer_lighten: float = 0.35,
        title_fontsize: float = 13,
        label_fontsize: float = 12,
        tick_fontsize: float = 11,
        legend_fontsize: float = 11,
        legend_loc: str = "best",
        title_weight: str = "bold",
    ) -> plt.Axes:
        if not state_fractions:
            raise ValueError("state_fractions must contain at least one state-to-fraction mapping.")

        created_figure = False
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
            created_figure = True

        ylabel = "CD signal (mdeg)"
        variant_colors = {key.upper(): value for key, value in (variant_colors or {}).items()}
        trace_labels = trace_labels or {}

        for variant in variants:
            variant_name = variant.upper()
            base_color = variant_colors.get(variant_name, self.COLORS.get(variant_name, "#333333"))

            for state_name, fraction in state_fractions.items():
                selected = self.get_records(fraction, variants=(variant_name,))
                records = selected.get(variant_name)
                if not records:
                    continue

                data, _ = self._select_variant_data(
                    records,
                    average_replicates=average_replicates,
                )
                prepared, ylabel = self._prepare_signal(
                    data=data,
                    mode=mode,
                    wavelength_range=wavelength_range,
                    normalize_range=normalize_range,
                    normalization=normalization,
                    concentration_uM=concentration_uM,
                    path_length_cm=path_length_cm,
                    num_residues=num_residues,
                )

                state_key = state_name.strip().lower()
                is_monomer = "monomer" in state_key
                line_style = "--" if is_monomer else "-"
                color = self._lighten_color(base_color, monomer_lighten if is_monomer else 0.0)
                default_label = f"{variant_name} {state_name.lower()}"

                ax.plot(
                    prepared["wavelength_nm"],
                    prepared["signal"],
                    color=color,
                    linestyle=line_style,
                    linewidth=linewidth,
                    label=trace_labels.get(default_label, default_label),
                )

        if title is None:
            title = "CD State Overlay"

        ax.set_title(title, fontsize=title_fontsize, weight=title_weight)
        ax.set_xlabel("Wavelength (nm)", fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.grid(True, alpha=0.25, linewidth=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=tick_fontsize)
        ax.legend(frameon=False, fontsize=legend_fontsize, loc=legend_loc)

        if wavelength_range is not None:
            ax.set_xlim(*wavelength_range)

        if created_figure:
            plt.tight_layout()

        return ax

    def plot_fraction_grid(
        self,
        fractions: Sequence[int],
        variants: Sequence[str] = ("WT", "N27D"),
        mode: str = "raw",
        average_replicates: bool = True,
        wavelength_range: tuple[float, float] | None = None,
        normalize_range: tuple[float, float] = (200, 250),
        normalization: str = "trough",
        concentration_uM: float | None = None,
        path_length_cm: float = 0.1,
        num_residues: int | None = None,
        fraction_labels: dict[int, str] | None = None,
        figsize: tuple[float, float] | None = None,
        sharey: bool = False,
    ) -> tuple[plt.Figure, np.ndarray]:
        if not fractions:
            raise ValueError("fractions must contain at least one fraction number.")

        width = 5.5 * len(fractions)
        figsize = figsize or (width, 4.8)
        fig, axes = plt.subplots(
            1,
            len(fractions),
            figsize=figsize,
            sharey=sharey,
            constrained_layout=True,
        )

        if len(fractions) == 1:
            axes = np.array([axes])

        for ax, fraction in zip(axes, fractions):
            self.plot_fraction(
                fraction=fraction,
                variants=variants,
                mode=mode,
                average_replicates=average_replicates,
                wavelength_range=wavelength_range,
                normalize_range=normalize_range,
                normalization=normalization,
                concentration_uM=concentration_uM,
                path_length_cm=path_length_cm,
                num_residues=num_residues,
                ax=ax,
                fraction_label=(fraction_labels or {}).get(fraction),
            )

        return fig, axes

    def plot_variant_overview(
        self,
        variants: Sequence[str] = ("WT", "N27D"),
        wavelength_range: tuple[float, float] | None = (195, 250),
        figsize: tuple[float, float] | None = None,
    ) -> tuple[plt.Figure, np.ndarray]:
        if not variants:
            raise ValueError("variants must contain at least one variant.")

        figsize = figsize or (7.0 * len(variants), 5.0)
        fig, axes = plt.subplots(
            1,
            len(variants),
            figsize=figsize,
            constrained_layout=True,
            sharey=True,
        )

        if len(variants) == 1:
            axes = np.array([axes])

        for ax, variant in zip(axes, variants):
            variant_name = variant.upper()
            variant_records = self.records.get(variant_name, {})
            if not variant_records:
                raise FileNotFoundError(f"No spectra were found for variant {variant_name}.")

            color_steps = np.linspace(0.45, 0.95, len(variant_records))
            cmap_name = "Blues" if variant_name == "WT" else "Reds" if variant_name == "N27D" else "Greens"
            colors = plt.get_cmap(cmap_name)(color_steps)

            for color, (fraction, records) in zip(colors, sorted(variant_records.items())):
                data, label_suffix = self._select_variant_data(records, average_replicates=True)
                filtered = self._filter_range(data, wavelength_range)
                ax.plot(
                    filtered["wavelength_nm"],
                    filtered["cd_mdeg"],
                    linewidth=2.2,
                    color=color,
                    label=f"F{fraction}{label_suffix}",
                )

            ax.set_title(f"{variant_name} Averaged Spectra", fontsize=13, weight="bold")
            ax.set_xlabel("Wavelength (nm)")
            ax.grid(True, alpha=0.25, linewidth=0.8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if wavelength_range is not None:
                ax.set_xlim(*wavelength_range)
            ax.legend(frameon=False, fontsize=9)

        axes[0].set_ylabel("CD signal (mdeg)")
        return fig, axes

    @staticmethod
    def apply_sec1_style(
        ax: plt.Axes,
        *,
        x_major_ticks: float | None = None,
        y_major_ticks: float | None = None,
        x_minor_ticks: float | None = None,
        y_minor_ticks: float | None = None,
        major_tick_length: float = 10,
        minor_tick_length: float = 6,
        major_tick_width: float = 2,
        minor_tick_width: float = 1,
        tick_labelsize: float = 20,
        bold_ticks: bool = True,
        spine_linewidth: float = 2.0,
        show_box: bool = True,
        font_family: str = "Arial",
    ) -> plt.Axes:
        """
        Apply the sec1.py publication style to a matplotlib Axes.

        Mirrors the SECData.list_plot conventions: bold Arial tick labels,
        ticks on all four sides pointing outward, thick spines, optional
        custom major/minor tick spacing (AutoMinorLocator if not set), and
        no grid.
        """
        plt.rcParams["font.family"] = font_family

        if x_major_ticks is not None:
            ax.xaxis.set_major_locator(MultipleLocator(x_major_ticks))
        if x_minor_ticks is not None:
            ax.xaxis.set_minor_locator(MultipleLocator(x_minor_ticks))
        else:
            ax.xaxis.set_minor_locator(AutoMinorLocator())

        if y_major_ticks is not None:
            ax.yaxis.set_major_locator(MultipleLocator(y_major_ticks))
        if y_minor_ticks is not None:
            ax.yaxis.set_minor_locator(MultipleLocator(y_minor_ticks))
        else:
            ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.tick_params(
            axis="both", which="major",
            labelsize=tick_labelsize, length=major_tick_length,
            direction="out", top=True, bottom=True, left=True, right=True,
            width=major_tick_width,
        )
        ax.tick_params(
            axis="both", which="minor",
            length=minor_tick_length, direction="out",
            top=True, bottom=True, left=True, right=True,
            width=minor_tick_width,
        )

        if bold_ticks:
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_weight("bold")

        for spine in ax.spines.values():
            spine.set_visible(show_box)
            if show_box:
                spine.set_linewidth(spine_linewidth)

        ax.grid(False)
        ax.margins(x=0.02, y=0.02)
        return ax

    def plot_publication_overlay(
        self,
        spectra: dict[str, pd.DataFrame],
        *,
        colors: dict[str, str] | None = None,
        wavelength_range: tuple[float, float] = (185, 250),
        normalize: bool = False,
        normalize_range: tuple[float, float] = (195, 250),
        signal_column: str = "cd_mdeg",
        wavelength_column: str = "wavelength_nm",
        title: str = "",
        xlabel: str = "Wavelength (nm)",
        ylabel: str | None = None,
        figsize: tuple[float, float] = (10, 7),
        dpi: int = 300,
        linewidth: float = 2.8,
        linestyles: dict[str, str] | None = None,
        x_major_ticks: float | None = 10,
        y_major_ticks: float | None = None,
        x_minor_ticks: float | None = 2,
        y_minor_ticks: float | None = None,
        tick_labelsize: float = 20,
        label_fontsize: float = 22,
        title_fontsize: float = 22,
        legend_fontsize: float = 20,
        legend_loc: str = "best",
        legend_box: bool = False,
        zero_line: bool = False,
    ) -> tuple[plt.Figure, plt.Axes]:
        """
        Plot a {variant: dataframe} dict in sec1.py publication style.

        Each DataFrame must contain ``wavelength_column`` and
        ``signal_column``. When ``normalize=True`` each spectrum is divided
        by ``abs(min)`` over ``normalize_range`` so the deepest trough maps
        to -1 (shape comparison).
        """
        if not spectra:
            raise ValueError("spectra must contain at least one variant.")

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        linestyles = linestyles or {}
        default_y = "Normalized CD" if normalize else "CD Signal (mdeg)"

        for variant, df in spectra.items():
            sub = df[
                (df[wavelength_column] >= wavelength_range[0])
                & (df[wavelength_column] <= wavelength_range[1])
            ].sort_values(wavelength_column)

            y = sub[signal_column].to_numpy(dtype=float)
            if normalize:
                norm_mask = (
                    (df[wavelength_column] >= normalize_range[0])
                    & (df[wavelength_column] <= normalize_range[1])
                )
                window = df.loc[norm_mask, signal_column]
                scale = abs(window.min())
                if np.isclose(scale, 0.0):
                    scale = float(window.abs().max())
                if np.isclose(scale, 0.0):
                    raise ValueError(f"Normalization scale is zero for {variant}.")
                y = y / scale

            color = (colors or {}).get(variant) or self.COLORS.get(variant.upper(), None)
            ax.plot(
                sub[wavelength_column].to_numpy(),
                y,
                color=color,
                linewidth=linewidth,
                linestyle=linestyles.get(variant, "-"),
                label=variant,
            )

        if zero_line:
            ax.axhline(0, color="black", linewidth=1.0, linestyle=":")

        ax.set_xlim(*wavelength_range)
        ax.set_xlabel(xlabel, fontsize=label_fontsize, fontweight="normal",
                      fontfamily="Arial", labelpad=10)
        ax.set_ylabel(ylabel or default_y, fontsize=label_fontsize,
                      fontweight="normal", fontfamily="Arial", labelpad=10)
        if title:
            ax.set_title(title, fontsize=title_fontsize, fontweight="bold",
                         fontfamily="Arial", pad=20)

        self.apply_sec1_style(
            ax,
            x_major_ticks=x_major_ticks, y_major_ticks=y_major_ticks,
            x_minor_ticks=x_minor_ticks, y_minor_ticks=y_minor_ticks,
            tick_labelsize=tick_labelsize,
        )

        legend = ax.legend(fontsize=legend_fontsize, frameon=legend_box,
                           loc=legend_loc)
        if legend:
            for text in legend.get_texts():
                text.set_weight("bold")

        fig.tight_layout()
        return fig, ax

    def save_figure(
        self,
        fig: plt.Figure,
        output_stem: str,
        output_dir: str | Path = "publication_plots",
        formats: Sequence[str] = ("png", "svg"),
        dpi: int = 600,
        transparent: bool = False,
    ) -> list[Path]:
        """
        Save a matplotlib figure in one or more publication-friendly formats.

        Args:
            fig: Matplotlib figure to save.
            output_stem: File stem without extension.
            output_dir: Folder relative to the workflow root, or an absolute path.
            formats: Output extensions, for example ``('png', 'svg')``.
            dpi: Export dpi. SVG is vector, but matplotlib still accepts dpi
                for any rasterized elements.
            transparent: Whether to export with a transparent background.

        Returns:
            list[Path]: The saved file paths.
        """
        destination = Path(output_dir)
        if not destination.is_absolute():
            destination = self.root_dir / destination
        destination.mkdir(parents=True, exist_ok=True)

        saved_paths: list[Path] = []
        for fmt in formats:
            extension = fmt.lower().lstrip(".")
            path = destination / f"{output_stem}.{extension}"
            fig.savefig(
                path,
                format=extension,
                dpi=dpi,
                bbox_inches="tight",
                transparent=transparent,
            )
            saved_paths.append(path)

        return saved_paths
