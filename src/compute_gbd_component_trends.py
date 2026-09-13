"""Create descriptive GBD component trend summaries for the US CKM analysis."""
from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parent
DATA_DIR = Path(os.environ.get("CKM_GBD_DATA_DIR", str(REPOSITORY_ROOT / "data" / "gbd"))).expanduser().resolve()
OUTPUT_DIR = Path(os.environ.get("CKM_OUTPUT_DIR", str(REPOSITORY_ROOT / "outputs"))).expanduser().resolve()
OUTPUT = OUTPUT_DIR / "gbd_component_trends_1990_2023.csv"

CAUSES = [
    ("Ischemic_heart_disease", "IHD"),
    ("Stroke", "Stroke"),
    ("Atrial_fibrillation_and_flutter", "AF/flutter"),
    ("Lower_extremity_peripheral_arterial_disease", "PAD"),
    ("Diabetes_Mellitus", "Diabetes"),
    ("Chronic_Kidney_Disease", "CKD"),
]
MEASURES = [
    ("Prevalence", "Age-standardized prevalence per 100,000"),
    ("Incidence", "Age-standardized incidence per 100,000"),
    ("Deaths", "Age-standardized deaths per 100,000"),
    ("DALY", "Age-standardized DALYs per 100,000"),
]


def _percent(coefficient: float) -> float:
    return 100.0 * (np.exp(coefficient) - 1.0)


def _interval(coefficient: float, standard_error: float, degrees_of_freedom: int) -> tuple[float, float]:
    critical = stats.t.ppf(0.975, degrees_of_freedom)
    return _percent(coefficient - critical * standard_error), _percent(coefficient + critical * standard_error)


def _fit_hac(design: np.ndarray, outcome: np.ndarray):
    return sm.OLS(outcome, design).fit(cov_type="HAC", cov_kwds={"maxlags": 1})


def trend_summary(values: pd.DataFrame) -> dict[str, float | int]:
    """Return AAPC, segment slopes, global sensitivity, and diagnostics."""
    year = values["Year"].to_numpy(dtype=float)
    log_rate = np.log(values["Value"].to_numpy(dtype=float))
    centered_year = year - year.mean()
    n = len(values)

    linear = _fit_hac(sm.add_constant(centered_year, has_constant="add"), log_rate)
    linear_params = np.asarray(linear.params, dtype=float)
    linear_cov = np.asarray(linear.cov_params(), dtype=float)
    linear_slope = float(linear_params[1])
    linear_se = float(np.sqrt(linear_cov[1, 1]))
    linear_lower, linear_upper = _interval(linear_slope, linear_se, n - 2)

    quadratic = _fit_hac(np.column_stack([np.ones(n), centered_year, centered_year**2]), log_rate)
    quadratic_p = float(np.asarray(quadratic.pvalues, dtype=float)[2])
    ljung_box_p = float(acorr_ljungbox(linear.resid, lags=[1], return_df=True)["lb_pvalue"].iloc[0])

    candidates = range(int(year.min()) + 8, int(year.max()) - 7)
    candidate_fits = []
    for breakpoint in candidates:
        hinge = np.maximum(year - breakpoint, 0.0)
        design = sm.add_constant(np.column_stack([centered_year, hinge]), has_constant="add")
        candidate_fits.append((sm.OLS(log_rate, design).fit().bic, breakpoint))
    best_bic, best_breakpoint = min(candidate_fits, key=lambda item: item[0])

    hinge = np.maximum(year - best_breakpoint, 0.0)
    segmented = _fit_hac(sm.add_constant(np.column_stack([centered_year, hinge]), has_constant="add"), log_rate)
    parameters = np.asarray(segmented.params, dtype=float)
    covariance = np.asarray(segmented.cov_params(), dtype=float)
    early_slope = float(parameters[1])
    late_slope = float(parameters[1] + parameters[2])
    early_se = float(np.sqrt(covariance[1, 1]))
    late_variance = covariance[1, 1] + covariance[2, 2] + 2.0 * covariance[1, 2]
    late_se = float(np.sqrt(max(late_variance, 0.0)))

    late_weight = float((year.max() - best_breakpoint) / (year.max() - year.min()))
    aapc_slope = float(early_slope + late_weight * parameters[2])
    aapc_variance = covariance[1, 1] + late_weight**2 * covariance[2, 2] + 2.0 * late_weight * covariance[1, 2]
    aapc_se = float(np.sqrt(max(aapc_variance, 0.0)))
    aapc_lower, aapc_upper = _interval(aapc_slope, aapc_se, n - 3)
    early_lower, early_upper = _interval(early_slope, early_se, n - 3)
    late_lower, late_upper = _interval(late_slope, late_se, n - 3)

    return {
        "global_log_linear_percent": _percent(linear_slope),
        "global_log_linear_lower_percent": linear_lower,
        "global_log_linear_upper_percent": linear_upper,
        "aapc_percent": _percent(aapc_slope),
        "aapc_lower_percent": aapc_lower,
        "aapc_upper_percent": aapc_upper,
        "aapc_break_year": int(best_breakpoint),
        "aapc_early_segment_percent": _percent(early_slope),
        "aapc_early_segment_lower_percent": early_lower,
        "aapc_early_segment_upper_percent": early_upper,
        "aapc_late_segment_percent": _percent(late_slope),
        "aapc_late_segment_lower_percent": late_lower,
        "aapc_late_segment_upper_percent": late_upper,
        "aapc_bic": float(best_bic),
        "aapc_delta_bic_vs_global_linear": float(best_bic - linear.bic),
        "durbin_watson": float(durbin_watson(linear.resid)),
        "ljung_box_p_lag1": ljung_box_p,
        "quadratic_term_p": quadratic_p,
    }


def _read_source(measure: str, cause_stem: str) -> pd.DataFrame:
    path = DATA_DIR / f"SourceData_GBD_{measure}_{cause_stem}.csv"
    data = pd.read_csv(path)
    required = {"Year", "Value", "Lower bound", "Upper bound"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
    data = data.loc[data["Year"].between(1990, 2023)].sort_values("Year").copy()
    if len(data) != 34 or data["Value"].isna().any() or (data["Value"] <= 0).any():
        raise ValueError(f"{path.name} must contain 34 positive annual values from 1990 through 2023")
    return data


def main() -> None:
    rows = []
    for measure, metric_label in MEASURES:
        for cause_stem, label in CAUSES:
            data = _read_source(measure, cause_stem)
            first, last = data.iloc[0], data.iloc[-1]
            rows.append({
                "measure": metric_label,
                "component": label,
                "rate_1990": first["Value"],
                "rate_2023": last["Value"],
                "rate_2023_lower": last["Lower bound"],
                "rate_2023_upper": last["Upper bound"],
                "absolute_change_1990_2023": last["Value"] - first["Value"],
                "percent_change_1990_2023": 100.0 * (last["Value"] / first["Value"] - 1.0),
                **trend_summary(data),
            })
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT, index=False)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
