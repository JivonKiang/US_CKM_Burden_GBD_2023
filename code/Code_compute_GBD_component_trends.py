"""Create reproducible GBD component trend tables for the US CKM manuscript."""

from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.stats.stattools import durbin_watson


PROJECT = Path(os.environ.get("CKM_PROJECT_ROOT", Path.cwd()))
OUTPUT = PROJECT / "Figures_Tables" / "Table_GBD_US_CKM_component_trends_1990_2023.csv"

CAUSES = [
    ("Ischemic heart disease", "IHD"),
    ("Stroke", "Stroke"),
    ("Atrial fibrillation and flutter", "AF/flutter"),
    ("Lower extremity peripheral arterial disease", "PAD"),
    ("Diabetes Mellitus", "Diabetes"),
    ("Chronic Kidney Disease", "CKD"),
]


def eapc(values: pd.DataFrame) -> dict[str, float]:
    """Return HAC-robust EAPC plus diagnostics and a two-segment sensitivity fit."""
    year = values["Year"].to_numpy(dtype=float)
    centered_year = year - year.mean()
    log_rate = np.log(values["Value"].to_numpy(dtype=float))
    design = sm.add_constant(centered_year)
    ordinary = sm.OLS(log_rate, design).fit()
    robust = ordinary.get_robustcov_results(cov_type="HAC", maxlags=1)
    n = len(values)
    critical = stats.t.ppf(0.975, n - 2)
    slope = float(robust.params[1])
    stderr = float(robust.bse[1])
    lower = slope - critical * stderr
    upper = slope + critical * stderr
    quadratic = sm.OLS(log_rate, np.column_stack([np.ones(n), centered_year, centered_year**2])).fit(
        cov_type="HAC", cov_kwds={"maxlags": 1}
    )
    ljung_box_p = float(acorr_ljungbox(ordinary.resid, lags=[1], return_df=True)["lb_pvalue"].iloc[0])
    convert = lambda coefficient: 100 * (np.exp(coefficient) - 1)

    # Sensitivity analysis for nonlinearity: fit continuous two-segment
    # log-linear models over prespecified interior break years and retain the
    # BIC-minimising breakpoint. This is not presented as a formal Joinpoint
    # software result; it is a reproducible check on whether one EAPC masks a
    # materially different early-versus-late slope.
    candidates = range(int(year.min()) + 8, int(year.max()) - 7)
    piecewise = []
    for breakpoint in candidates:
        hinge = np.maximum(year - breakpoint, 0.0)
        design_pw = sm.add_constant(np.column_stack([year - year.mean(), hinge]))
        fit_pw = sm.OLS(log_rate, design_pw).fit()
        piecewise.append((fit_pw.bic, breakpoint, fit_pw))
    best_bic, best_breakpoint, best_fit = min(piecewise, key=lambda item: item[0])
    early_slope = float(best_fit.params[1])
    late_slope = float(best_fit.params[1] + best_fit.params[2])
    return {
        "eapc_percent": convert(slope),
        "eapc_lower_percent": convert(lower),
        "eapc_upper_percent": convert(upper),
        "durbin_watson": float(durbin_watson(ordinary.resid)),
        "ljung_box_p_lag1": ljung_box_p,
        "quadratic_term_p": float(quadratic.pvalues[2]),
        "piecewise_break_year": int(best_breakpoint),
        "piecewise_eapc_early_percent": convert(early_slope),
        "piecewise_eapc_late_percent": convert(late_slope),
        "piecewise_bic": float(best_bic),
        "piecewise_delta_bic_vs_linear": float(best_bic - ordinary.bic),
    }


def summarize(measure: str, metric_label: str) -> list[dict]:
    rows = []
    for filename, label in CAUSES:
        data = pd.read_csv(PROJECT / "data" / measure / f"{filename}.csv")
        data = data.loc[data["Year"].between(1990, 2023)].sort_values("Year")
        first, last = data.iloc[0], data.iloc[-1]
        diagnostics = eapc(data)
        rows.append(
            {
                "measure": metric_label,
                "component": label,
                "rate_1990": first["Value"],
                "rate_2023": last["Value"],
                "rate_2023_lower": last["Lower bound"],
                "rate_2023_upper": last["Upper bound"],
                "absolute_change_1990_2023": last["Value"] - first["Value"],
                "percent_change_1990_2023": 100 * (last["Value"] / first["Value"] - 1),
                **diagnostics,
            }
        )
    return rows


def main() -> None:
    rows = summarize("Prevalence", "Age-standardized prevalence per 100,000")
    rows += summarize("Incidence", "Age-standardized incidence per 100,000")
    rows += summarize("Deaths", "Age-standardized deaths per 100,000")
    rows += summarize("DALY", "Age-standardized DALYs per 100,000")
    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT, index=False)
    print(f"Wrote {OUTPUT}")
    print(result.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
