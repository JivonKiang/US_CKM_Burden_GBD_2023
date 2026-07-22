"""Register the NHANES source files required for the US CKM analysis.

The historical XPT files supplied by the investigator remain read-only in
their original directory.  Only files absent from that archive are retrieved
from the official CDC public-use static-data endpoint.  This script does not
derive CKM stages; it creates an auditable input manifest for later analysis.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
from urllib.request import urlopen


PROJECT = Path(os.environ.get("CKM_PROJECT_ROOT", Path.cwd()))
LOCAL_ROOT_ENV = os.environ.get("CKM_NHANES_SOURCE_ROOT")
if not LOCAL_ROOT_ENV:
    raise RuntimeError("Set CKM_NHANES_SOURCE_ROOT to a local, authorized NHANES source directory before running this input-registration script.")
LOCAL_ROOT = Path(LOCAL_ROOT_ENV)
DOWNLOAD_ROOT = PROJECT / "data" / "NHANES" / "official_missing_modules"
MANIFEST = PROJECT / "data" / "NHANES" / "nhanes_ckm_input_manifest_20260716.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_xpt(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(url, timeout=120) as response, destination.open("wb") as output:
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
    # SAS transport files begin with a textual header; reject HTML/error pages.
    prefix = destination.read_bytes()[:200].decode("ascii", errors="ignore")
    if "HEADER RECORD" not in prefix:
        destination.unlink(missing_ok=True)
        raise RuntimeError(f"CDC response was not a SAS transport file: {url}")


def local_file(folder: str, filename: str) -> tuple[str, Path, str]:
    path = LOCAL_ROOT / folder / filename
    return "investigator_archive", path, ""


def official_file(cycle: str, filename: str) -> tuple[str, Path, str]:
    url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{cycle}/DataFiles/{filename}"
    return "cdc_public_use", DOWNLOAD_ROOT / filename, url


def cycle_filename(module: str, suffix: str, extension: str = "XPT") -> str:
    """NHANES 2017-March 2020 files use P_ prefix without a suffix."""
    return f"P_{module}.{extension}" if suffix == "P" else f"{module}_{suffix}.{extension}"


def resolve_preferred_file(
    folder: str, filename: str, cdc_year: str
) -> tuple[str, Path, str]:
    source, path, url = local_file(folder, filename)
    if path.exists():
        return source, path, url
    return official_file(cdc_year, filename.lower())


def main() -> None:
    records = []
    cycles = {
        "2011-2012": {"suffix": "G", "cdc_year": "2011", "prefix": ""},
        "2013-2014": {"suffix": "H", "cdc_year": "2013", "prefix": ""},
        "2015-2016": {"suffix": "I", "cdc_year": "2015", "prefix": ""},
        "2017-March 2020": {"suffix": "P", "cdc_year": "2017", "prefix": "P_"},
        "2021-2023 candidate": {"suffix": "L", "cdc_year": "2021", "prefix": ""},
    }
    local_modules = {
        "DEMO": "DEMO",
        "BMX": "BMX",
        "BPQ": "BPQ",
        "DIQ": "DIQ",
        "MCQ": "MCQ",
        "GHB": "Glycohemoglobin",
        "GLU": "Glucose",
        "HDL": "Cholesterol",
        "TCHOL": "Cholesterol",
        "TRIGLY": "Cholesterol",
        "SMQ": "SMQ",
    }

    for cycle, info in cycles.items():
        suffix = info["suffix"]
        for module, folder in local_modules.items():
            filename = cycle_filename(module, suffix)
            source, path, url = resolve_preferred_file(folder, filename, info["cdc_year"])
            records.append((cycle, module, source, path, url, "historical local source"))

        bp_filename = (
            "P_BPXO.XPT" if suffix == "P" else "BPXO_L.XPT" if suffix == "L" else f"BPX_{suffix}.XPT"
        )
        source, path, url = resolve_preferred_file("BPX", bp_filename, info["cdc_year"])
        records.append((cycle, "BPX", source, path, url, "blood pressure device differs for P/L; harmonisation required"))

        if suffix != "L":
            source, path, url = resolve_preferred_file("UA", cycle_filename("BIOPRO", suffix), info["cdc_year"])
            records.append((cycle, "BIOPRO", source, path, url, "serum creatinine available locally"))
        else:
            source, path, url = official_file(info["cdc_year"], "BIOPRO_L.xpt")
            records.append((cycle, "BIOPRO", source, path, url, "missing from investigator archive"))

        alb_filename = cycle_filename("ALB_CR", suffix, extension="xpt")
        source, path, url = official_file(info["cdc_year"], alb_filename)
        records.append((cycle, "ALB_CR", source, path, url, "missing from investigator archive"))

    manifest_rows = []
    for cycle, module, source, path, url, note in records:
        if source == "cdc_public_use" and not path.exists():
            print(f"Downloading {path.name}")
            download_xpt(url, path)
        status = "available" if path.exists() else "missing"
        manifest_rows.append(
            {
                "cycle": cycle,
                "module": module,
                "source": source,
                "path": str(path),
                "official_url": url,
                "status": status,
                "bytes": path.stat().st_size if path.exists() else "",
                "sha256": sha256(path) if path.exists() else "",
                "note": note,
            }
        )

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)
    unavailable = [row for row in manifest_rows if row["status"] != "available"]
    print(f"Wrote {MANIFEST}")
    print(f"Registered {len(manifest_rows)} inputs; unavailable: {len(unavailable)}")


if __name__ == "__main__":
    main()
